# Content described in ADIPEC paper SPE-229932-MS, 2025 titled:
# "A Universal, EOS-Based Correlation for Z-Factor, Viscosity and Enthalpy For Hydrocarbon and H₂, N₂, CO₂, H₂S Gas Mixtures"
# M. W. Burgoyne, Santos; M. H. Nielsen, Whitson AS; M. Stanko, Whitson AS

# BNS_10.py - Peng-Robinson Z-Factor, Thermodynamic, and Viscosity Properties

import numpy as np
from typing import Dict, Any
from dataclasses import dataclass

@dataclass(frozen=True)
class Constants:
    # Universal and empirical constants (field units unless otherwise noted)
    R: float = 10.731577089016         # ft³·psia/(lb-mol·R) (universal gas constant for PVT units)
    R_THERMO: float = 1.98588          # Btu/(lb-mol·R) (universal gas constant for Thermal units)
    MW_AIR: float = 28.97              # Molecular weight of air
    DEG_F_TO_R: float = 459.67         # °F to °R offset
    FT3_PSIA_TO_BTU: float = 5.403     # Conversion for (ft³·psia) to Btu
    mwCH4: float = 16.0425             # Methane MW
    TcCH4: float = 343.008             # Methane critical temperature (deg R)
    PcCH4: float = 667.029             # Methane critical pressure (psia)
    VcZcCH4: float = R * TcCH4 / PcCH4 # Methane Vc/Zc (ft3/lb-mol)

CONSTS = Constants()

# Cp polynomial coefficients for components (order: a, b, c, d, e)
# Each row is a component: [CO2, H2S, N2, H2, Gas (C1+)]
# Fitted to NIST data at zero psia (Burgoyne, 2025)
Cp_poly_coeffs = np.array([
    [2.725473196,  0.004103751,  1.5602E-05, -4.19321E-08,  3.10542E-11],  # CO2
    [4.446031265, -0.005296052,  2.0533E-05, -2.58993E-08,  1.25555E-11],  # H2S
    [3.423811591,  0.001007461, -4.58491E-06,  8.4252E-09, -4.38083E-12],  # N2
    [1.421468418,  0.018192108, -6.04285E-05,  9.08033E-08, -5.18972E-11], # H2
    [5.369051342, -0.014851371,  4.86358E-05, -3.70187E-08,  1.80641E-12], # Methane
])

# All component properties, by index (CO2, H2S, N2, H2, Gas)
PROPS = {
    'mws'   : np.array([44.01, 34.082, 28.014, 2.016, 0.0]),
    'tcs'   : np.array([547.416, 672.120, 227.160, 47.430, 1.0]),                 # Custom Tc for H2 (Neislen, 2023)
    'pcs'   : np.array([1069.51, 1299.97, 492.84, 187.53, 1.0]),
    'ACF'   : np.array([0.12253, 0.04909, 0.037, -0.217, -0.03899]),              # Retuned Burgoyne, 2025
    'VSHIFT': np.array([-0.27607, -0.22901, -0.21066, -0.36270, -0.19076]),       # Retuned Burgoyne, 2025
    'OmegaA': np.array([0.427671, 0.436725, 0.457236, 0.457236, 0.457236]),       # Retuned CO2 and H2S: Burgoyne, 2025
    'OmegaB': np.array([0.0696397, 0.0724345, 0.0777961, 0.0777961, 0.0777961]),  # Retuned CO2 and H2S: Burgoyne, 2025
    'VCVIS' : np.array([1.46352, 1.46808, 1.35526, 0.68473,  0.0]),               # Retuned Burgoyne, 2025
    'Cp_poly': Cp_poly_coeffs,
}

def tc_ag(x):
    a, b, c = 2695.14765, 274.341701, CONSTS.TcCH4
    return a * x / (b + x) + c

def tc_gc(x):
    a, b, c = 1098.10948, 101.529237, CONSTS.TcCH4
    return a * x / (b + x) + c

def pc_fn(x, vc_slope, tc_p):
    vc_on_zc = vc_slope * x + CONSTS.VcZcCH4
    return CONSTS.R * tc_p / vc_on_zc

def pseudo_critical(sg_hc, AG = False):
    """
    Calculates pseudo-critical temperature and pressure for the pseudo-hydrocarbon component.
    - Uses custom linear fits for gas condensates (AG = False) and associated gas (AG = True) (Burgoyne, 2025).
    - These relations are derived from property fitting per paper reference.
    """
    x = max(0, CONSTS.MW_AIR * sg_hc - CONSTS.mwCH4)
    if AG:
        tpc_hc = tc_ag(x)
        vc_slope = 0.177497835
    else:
        tpc_hc = tc_gc(x)
        vc_slope = 0.170931432
    ppc_hc = pc_fn(x, vc_slope, tpc_hc)
    return tpc_hc, ppc_hc   


def update_hydrocarbon(props, sg_hc, AG):
    """
    Updates the hydrocarbon pseudo-component's critical props, MW, and viscosity parameters
    based on current pseudo-component MW (from overall gas specific gravity).
    Ensures no inplace mutation of the shared props dict.
    """
    tcs, pcs, mws, VCVIS = [arr.copy() for arr in (props['tcs'], props['pcs'], props['mws'], props['VCVIS'])]
    hc_mw = sg_hc * CONSTS.MW_AIR
    tcs[-1], pcs[-1] = pseudo_critical(sg_hc, AG)
    mws[-1] = hc_mw
    # VCVIS: Empirical fit for viscosity pseudo-critical volume (Burgoyne, 2025)
    VCVIS[-1] = 0.0576710 * (hc_mw - CONSTS.mwCH4) + 1.44383
    newprops = dict(props)
    newprops.update({'tcs': tcs, 'pcs': pcs, 'mws': mws, 'VCVIS': VCVIS})
    return newprops

def methane_adjust(cp_coeffs: np.ndarray, hc_mw: float) -> np.ndarray:
    """
    Adjusts the Cp polynomial coefficients for the 'Gas' (C1+) pseudo-component
    based on the MW of the hydrocarbon fraction.
    Each Cp term gets a different quadratic scaling.
    Only the last row (C1+) is changed.
    Custom function (Burgoyne, 2025)
    """
    a0 = np.array([7.8570E-04, 1.3123E-03, 9.8133E-04, 1.6463E-03, 1.7306E-02])
    a1 = np.array([-8.1649E-03, 5.5485E-03, 8.3258E-02, 2.0635E-01, 2.5551E+00])
    x = hc_mw - CONSTS.mwCH4
    scale = a0 * x ** 2 + a1 * x + 1  # Apply to each polynomial coefficient
    cp_adj = cp_coeffs.copy()
    cp_adj[-1, :] *= scale  # Only for 'Gas'
    return cp_adj

def get_z_fractions(co2, h2s, n2, h2) -> np.ndarray:
    """
    Returns array of mole fractions in fixed order [CO2, H2S, N2, H2, Gas].
    Fraction for 'Gas' is 1 minus all others.
    """
    frac_hc = 1.0 - (co2 + h2s + n2 + h2)
    return np.array([co2, h2s, n2, h2, frac_hc])

def hydrocarbon_sg(sg, zf, mws) -> float:
    """
    Given overall gas specific gravity, non-HC fractions, and MWs,
    solves for the specific gravity of the hydrocarbon pseudo-component.
    """
    frac_hc = zf[-1]
    sum_nonhc = np.dot(zf[:-1], mws[:-1])
    if frac_hc > 0:
        return (sg - sum_nonhc / CONSTS.MW_AIR) / frac_hc
    else:
        return 0.75  # fallback for near-zero HC content

def fugacity_coefficient(Z, A, B):
    """
    Returns the fugacity coefficient for given Z, A, B (PR EOS, pure component form).
    """
    sqrt2 = np.sqrt(2)
    arg = (Z + (1 + sqrt2) * B) / (Z + (1 - sqrt2) * B)
    # Safe log: handle division by zero or negative input gracefully
    if arg <= 0 or (Z - B) <= 0:
        return np.inf
    ln_phi = (Z - 1) - np.log(Z - B) - (A / (2 * sqrt2 * B)) * np.log(arg)
    return np.exp(ln_phi)


def cubic_root(coeffs, A=None, B=None):
    """
    Finds all real, positive roots of the cubic.
    - If more than one, selects the root with minimum fugacity coefficient (lowest free energy).
    - If only one, returns it.
    - If none, returns Z=1.0 and sets ideal_gas_flag = True.
    Returns (Z, ideal_gas_flag).
    """
    roots = np.roots(coeffs)
    # Only consider real, positive roots
    real_positive_roots = [r.real for r in roots if abs(r.imag) < 1e-8 and r.real > 0]
    if not real_positive_roots:
        return 1.0, True  # Fallback: ideal gas
    if A is None or B is None or len(real_positive_roots) == 1:
        return max(real_positive_roots), False  # Default to vapor root
    # Select root with minimum fugacity coefficient
    fugacities = [fugacity_coefficient(z, A, B) for z in real_positive_roots]
    min_idx = int(np.argmin(fugacities))
    return real_positive_roots[min_idx], False


def pr_m_params(ACF: np.ndarray) -> np.ndarray:
    """
    Computes 'm' parameter for PR EOS alpha(T) function for each component.
    """
    return 0.37464 + 1.54226 * ACF - 0.26992 * (ACF ** 2)

def calc_bips(degR, tpc_hc):
    """
    Temperature-dependent Binary Interaction Parameters (BIPs) for all pairs,
    with analytic first and second T-derivatives for use in PR analytic derivatives.
    - Fitted forms: constant + slope/Tr, Tr based on component crit temp.
    - Returns: kij, dkij_dT, d2kij_dT2 (all NxN)
    - Tuned to experimental VLE data Burgoyne, 2025
    """
    components = ['CO2', 'H2S', 'N2', 'H2', 'Gas']
    bip_parameters = {
        ("Gas", "CO2"): {"constant": -0.145561 ,  "Tr_slope": 0.276572 ,  "tc": tpc_hc  },
        ("Gas", "H2S"): {"constant": 0.16852   ,  "Tr_slope": -0.122378,  "tc": tpc_hc  },
        ("Gas", "N2"):  {"constant": -0.108    ,  "Tr_slope": 0.0605506,  "tc": tpc_hc  },
        ("Gas", "H2"):  {"constant": -0.0620119,  "Tr_slope": 0.0427873,  "tc": tpc_hc  },
        ("CO2", "H2S"): {"constant": 0.248638  ,  "Tr_slope": -0.138185,  "tc": 547.416 },
        ("CO2", "N2"):  {"constant": -0.25     ,  "Tr_slope": 0.11602  ,  "tc": 547.416 },
        ("CO2", "H2"):  {"constant": -0.247153 ,  "Tr_slope": 0.16377  ,  "tc": 547.416 },
        ("H2S", "N2"):  {"constant": -0.204414 ,  "Tr_slope": 0.234417 ,  "tc": 672.12  },
        ("H2S", "H2"):  {"constant": 0         ,  "Tr_slope": 0        ,  "tc": 672.12  },
        ("N2",  "H2"):  {"constant": -0.166253 ,  "Tr_slope": 0.0788129,  "tc": 227.16  },
    }
    def lookup_key(i, j):
        return (i, j) if (i, j) in bip_parameters else (j, i)
    n = len(components)
    kij = np.zeros((n, n))
    dkij_dT = np.zeros((n, n))
    d2kij_dT2 = np.zeros((n, n))
    for i, ci in enumerate(components):
        for j, cj in enumerate(components):
            if ci == cj:
                kij[i, j] = 0.0
                dkij_dT[i, j] = 0.0
                d2kij_dT2[i, j] = 0.0
            else:
                params = bip_parameters[lookup_key(ci, cj)]
                tc = params["tc"]
                const = params["constant"]
                slope = params["Tr_slope"]
                Tr = degR / tc
                kij[i, j] = const + slope / Tr
                dkij_dT[i, j] = -slope * tc / (degR**2)
                d2kij_dT2[i, j] = 2*slope*tc / (degR**3)
    return kij, dkij_dT, d2kij_dT2

def stiel_thodos_viscosity(t_rankine, molecular_weights, tcs_, pcs_):
    """
    Stiel-Thodos dilute gas viscosity (cP) for each component.
    - Used as the basis for compositional viscosity mixing.
    """
    Tr = t_rankine / tcs_
    Tc_k = tcs_ * (5.0 / 9.0)
    Pc_atm = pcs_ / 14.696
    eta_factor = Tc_k ** (1.0 / 6.0) / (np.sqrt(molecular_weights) * Pc_atm ** (2.0 / 3.0))
    mu_vals = np.where(
        Tr <= 1.5,
        34e-5 * (Tr ** 0.94) / eta_factor,
        17.78e-5 * ((4.58 * Tr) - 1.67) ** (5.0 / 8.0) / eta_factor
    )
    return mu_vals

def lbc_viscosity(Z, degf, psia, sg, co2=0.0, h2s=0.0, n2=0.0, h2=0.0, AG=False):
    """
    Computes LBC viscosity (cP) using dilute gas mixing + dense phase polynomial in reduced density.
    - Used by PR property package as a field engineering default.
    """
    zf = get_z_fractions(co2, h2s, n2, h2)
    props = dict(PROPS)
    sg_hc = hydrocarbon_sg(sg, zf, props['mws'])
    sg_hc = max(sg_hc, CONSTS.mwCH4 / CONSTS.MW_AIR)
    hc_mw = sg_hc * CONSTS.MW_AIR
    props = update_hydrocarbon(props, sg_hc, AG)
    mws, tcs, pcs, VCVIS = props['mws'], props['tcs'], props['pcs'], props['VCVIS']
    degR = degf + CONSTS.DEG_F_TO_R
    mu_components = stiel_thodos_viscosity(degR, mws, tcs, pcs)
    sqrt_mw = np.sqrt(mws)
    numerator = np.dot(zf, mu_components * sqrt_mw)
    denominator = np.dot(zf, sqrt_mw)
    mu_dilute = numerator / denominator
    # Reduced density mixing
    rhoc = 1.0 / np.dot(VCVIS, zf)
    gas_density = psia / (Z * CONSTS.R * degR)
    rhor = gas_density / rhoc
    # Polynomial correction for dense phase
    a = [0.1023, 0.023364, 0.058533, -0.0392852, 0.00926279]            # Retuned Burgoyne, 2025
    lhs = a[0] + a[1]*rhor + a[2]*rhor**2 + a[3]*rhor**3 + a[4]*rhor**4
    Tc_mix = np.dot(zf, tcs) * (5.0 / 9.0)
    Pc_mix = np.dot(zf, pcs) / 14.696
    Mw_mix = np.dot(zf, mws)
    eta_mix = (Tc_mix ** (1.0 / 6.0)) / (np.sqrt(Mw_mix) * (Pc_mix ** (2.0 / 3.0)))
    viscosity = (lhs ** 4 - 1e-4)/eta_mix + mu_dilute
    return viscosity
             

def pr_properties(
    temp: float, pres: float, sg: float,
    co2: float=0.0, h2s: float=0.0, n2: float=0.0, h2: float=0.0, AG: bool=False,
    viscosity: bool=False, density: bool=False, thermo: bool=False, Metric: bool=False, verbose: bool=False
    ) -> Dict[str, Any]:
    """
    Main Peng-Robinson property function.
    - Handles composition, EOS property calculation, Z-factor, viscosity, and all closed-form
      thermodynamic departures and analytic derivatives.
    - Only calculates analytic departures if thermo==True for speed.
    
    With Metric = True, temperature is in deg C, pressure in MPa and results are also in Metric units per below
    With Metric = False (or omitted), temperature in deg F, pressure in psia, and results return in Field units
    
    Inputs:
    -------
    temp: degF or degC
    pres: psia or MPa
    sg: gas mixture specific gravity relative to air. Used to calculate the hydrocarbon MW
        - If 100% inert mixture, then sg input is ignored
        - If implied hydrocarbon gas sg < methane, then sg input is ignored and hydrocarbon MW set to methane.
    co2, h2s, n2, h2: Mole fractions of respective inert components of the mixture. Defaults to zero if undefined
    AG: Boolean flag that controls the Tc, Pc relationship to use for hydrocarbon pseudocomponent. Defaults to False (for gas condensate). True used for Associated gas
    viscosity: Boolean flag that controls whether viscosity is calculated (default False)
    density: Boolean flag that controls whether density is calculated (default False)
    thermo: Boolean flag that controls whether thermal properties are calculated (default False)
    verbose: Boolean flag that controls whether additional intermediate values and derivatives - always in Field units - are also returned (default False)

    Outputs:
    --------
    Returns a dictionary of results - depending on requested values - that always includes Z-factor and optionally;
    {'Z':        Z-Factor                                     (Dimensionless)          
    'Density':   Gas Density                                  (lbm/cuft or kg/m3)
    'H':         Enthalpy relative to 60 degF and 14.606 psia (Btu/(lb-mol) or kJ/(kmol))
    'Cp':        Isobaric heat capacity                       (Btu/(lb-mol·R) or kJ/(kmol K))
    'JT':        Joule-Thomson Coefficient                    (degF/psi or degC/MPa)
    'Viscosity': Gas viscosity                                (cP or mPa·s)} 

    """
    
    # Note that Cv has known innacuracies as it relies on second derivatives poorly characterized with cubic EOS
    # I have left the Cv calculation in the code, but do not return it
    # 'Cv':        Isochoric heat capacity                      (Btu/(lb-mol·R) or kJ/(kmol K)) - 

    if (co2 + h2s + n2 + h2) > 1.0:
        raise ValueError(
            f"Invalid composition: sum of CO₂({co2}) + H₂S({h2s}) + N₂({n2}) + H₂({h2}) = "
            f"{co2 + h2s + n2 + h2:.6f} > 1.0"
        )    
    
    if Metric:
        degF = temp * 9 / 5 + 32
        psia = pres * 145.0377
    else:
        degF = temp
        psia = pres
        
        
    zf = get_z_fractions(co2, h2s, n2, h2)
    props = dict(PROPS)
    sg_hc = hydrocarbon_sg(sg, zf, props['mws'])
    sg_hc = max(sg_hc, CONSTS.mwCH4 / CONSTS.MW_AIR)
    hc_mw = sg_hc * CONSTS.MW_AIR
    props = update_hydrocarbon(props, sg_hc, AG)
    mws, tcs, pcs, ACF, VSHIFT, OmegaA, OmegaB, VCVIS = (props[k] for k in ('mws','tcs','pcs','ACF','VSHIFT','OmegaA','OmegaB','VCVIS'))
    mws[-1] = hc_mw
    degR = degF + CONSTS.DEG_F_TO_R
    tpc_hc, _ = pseudo_critical(sg_hc)
    trs = degR / tcs
    prs = psia / pcs
    kij, dkij_dT, d2kij_dT2 = calc_bips(degR, tpc_hc)
    m_i = pr_m_params(ACF)
    # PR alpha(T): temperature correction for attraction parameter
    alpha = (1.0 + m_i * (1.0 - np.sqrt(trs))) ** 2
    a_c_i = OmegaA * CONSTS.R ** 2 * tcs ** 2 / pcs
    a_i = a_c_i * alpha
    b_i = OmegaB * CONSTS.R * tcs / pcs
    aij = np.sqrt(np.outer(a_i, a_i)) * (1.0 - kij)
    a_mix = np.sum(np.outer(zf, zf) * aij)
    b_mix = np.dot(zf, b_i)
    RT = CONSTS.R * degR
    A = a_mix * psia / (RT ** 2)
    B = b_mix * psia / RT
    coeffs = [1.0, -(1.0 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3)]
    z_eos, ideal_gas_flag = cubic_root(coeffs, A, B)
    Bi = OmegaB * (prs / trs)
    z_vshift = z_eos - np.dot(zf, VSHIFT * Bi)
    result = {"Z": float(z_vshift)}
    
    if density:
        mwt = float(np.dot(zf, mws))
        # density in lbm/ft  (Field units)
        result['Density'] = float(mwt * psia / (z_vshift * CONSTS.R * degR))
        

    if thermo:
        # ----- Analytic derivatives -----
        d_alpha_dTr = -m_i * (1.0 + m_i * (1.0 - np.sqrt(trs))) / np.sqrt(trs)
        d_alpha_dT = d_alpha_dTr / tcs
        da_i_dT = a_c_i * d_alpha_dT  # d(ai)/dT for each component

        # Mixing rules: analytic da_mix/dT (includes BIP T-derivative terms)
        daij_dT = np.zeros_like(aij)
        for i in range(5):
            for j in range(5):
                sqrt_ai_aj = np.sqrt(a_i[i] * a_i[j])
                N = da_i_dT[i] * a_i[j] + a_i[i] * da_i_dT[j]
                daij_dT[i,j] = -dkij_dT[i,j] * sqrt_ai_aj + (1.0 - kij[i,j]) * 0.5 * (N / sqrt_ai_aj)
        da_mix_dT = np.sum(np.outer(zf, zf) * daij_dT)

        # d2a_mix/dT²: needed for analytic Cp (in departure)
        d2alpha_dTr2 = 2*(-m_i/(2*np.sqrt(trs)))**2 + 2*(1.0 + m_i*(1.0 - np.sqrt(trs))) * (m_i/(4*trs**1.5))
        d2alpha_dT2 = d2alpha_dTr2 / (tcs**2)
        d2a_i_dT2 = a_c_i * d2alpha_dT2
        d2aij_dT2 = np.zeros_like(aij)
        for i in range(5):
            for j in range(5):
                sqrt_ai_aj = np.sqrt(a_i[i] * a_i[j])
                N = da_i_dT[i] * a_i[j] + a_i[i] * da_i_dT[j]
                cross_2_da = 2 * da_i_dT[i] * da_i_dT[j]
                d2a_term = d2a_i_dT2[i]*a_i[j] + cross_2_da + a_i[i]*d2a_i_dT2[j]
                D = sqrt_ai_aj
                term1 = - (N / D) * dkij_dT[i,j]
                term2 = (1.0 - kij[i,j]) * 0.5 * ((d2a_term / D) - (N**2) / (2 * D**3))
                term3 = - D * d2kij_dT2[i,j]
                d2aij_dT2[i,j] = term1 + term2 + term3
        d2a_mix_dT2 = np.sum(np.outer(zf, zf) * d2aij_dT2)

        def compute_dz_dT_constP(T, P, a_mix_val, da_mix_dT_val, b_mix_val, z_val):
            """
            Closed-form analytic dZ/dT at constant pressure for PR EOS.
            Used in Cp and JT calculations. 
            """
            A_ = a_mix_val * P / (CONSTS.R * T) ** 2
            B_ = b_mix_val * P / (CONSTS.R * T)
            dB_dT_ = -B_ / T
            dA_dT_ = (P / (CONSTS.R * T) ** 2) * da_mix_dT_val - 2 * a_mix_val * P / (CONSTS.R ** 2 * T ** 3)
            dF_dz = 3 * z_val ** 2 - 2 * (1 - B_) * z_val + (A_ - 3 * B_ ** 2 - 2 * B_)
            dF_dA = z_val - B_
            dF_dB = z_val ** 2 - (6 * B_ + 2) * z_val - A_ + 2 * B_ + 3 * B_ ** 2
            dF_dT_ = dF_dA * dA_dT_ + dF_dB * dB_dT_
            dz_dT_ = - dF_dT_ / dF_dz
            return dz_dT_
        dz_dT = compute_dz_dT_constP(degR, psia, a_mix, da_mix_dT, b_mix, z_eos)

        # -- Ideal gas Cp (polynomial, componentwise, with MW-adjusted C1+) --
        Cp_poly = methane_adjust(props['Cp_poly'], hc_mw)
        T_K = degR * 5.0 / 9.0
        Cp_vals = np.array([np.polyval(Cp_poly[i, ::-1], T_K) for i in range(5)])  # Each component's Cp(T)
        Cp_IG = np.dot(zf, Cp_vals) * CONSTS.R_THERMO

        # --- PR departure and analytic property calculations ---
        sqrt2 = np.sqrt(2)
        Bdim = (b_mix * psia) / (CONSTS.R * degR)
        def log_term(z):
            return np.log((z + (sqrt2 + 1) * Bdim) / (z - (sqrt2 - 1) * Bdim))
        H_dep_ft3psia = (CONSTS.R * degR * (z_eos - 1) + (degR * da_mix_dT - a_mix) / (2 * sqrt2 * b_mix) * log_term(z_eos))
        H_dep_Btu = H_dep_ft3psia / CONSTS.FT3_PSIA_TO_BTU
        # Empirical reference enthalpy at 60°F for each component - From the results of this function at 60 degF and 14.696 psia (Burgoyne, 2025)
        H0_ = [-16.6022, -21.5512, -3.57757, 0.008054, 0]
        H0_[-1] = -0.015774 * (hc_mw - CONSTS.mwCH4)**2 - 0.646645 * (hc_mw - CONSTS.mwCH4) - 8.2551915 # Fitted to The results of this function at different HC MW's (Burgoyne, 2025)
        H0 = np.dot(zf, H0_)
        
        # Ideal-gas enthalpy change (from reference T to T)
        T_ref_F = 60.0
        T_ref_K = (T_ref_F + CONSTS.DEG_F_TO_R) * 5.0 / 9.0
        cp = Cp_poly
        H_IG_Btu = CONSTS.R_THERMO * 9 / 5 * np.dot(zf, (
            cp[:,0] * (T_K - T_ref_K)
          + cp[:,1]/2 * (T_K**2 - T_ref_K**2)
          + cp[:,2]/3 * (T_K**3 - T_ref_K**3)
          + cp[:,3]/4 * (T_K**4 - T_ref_K**4)
          + cp[:,4]/5 * (T_K**5 - T_ref_K**5)
        ))
        H_total_Btu = H_IG_Btu + H_dep_Btu - H0

        def Hdep_deriv(T, z, dzdT):
            """
            Analytic derivative of PR departure enthalpy wrt temperature (constant P).
            - Used for Cp departure.
            """
            dF1dT = CONSTS.R * (z - 1) + CONSTS.R * T * dzdT
            N = z + (sqrt2 + 1) * Bdim
            D = z - (sqrt2 - 1) * Bdim
            dBdim_dT = - Bdim / T
            dN_dT = dzdT + (sqrt2 + 1) * dBdim_dT
            dD_dT = dzdT - (sqrt2 - 1) * dBdim_dT
            dln_dT = (1.0 / N) * dN_dT - (1.0 / D) * dD_dT
            X = (T * da_mix_dT - a_mix) / (2 * sqrt2 * b_mix)
            dX_dT = (da_mix_dT + T * d2a_mix_dT2 - da_mix_dT) / (2 * sqrt2 * b_mix)
            dF2dT = dX_dT * np.log(N / D) + X * dln_dT
            return (dF1dT + dF2dT) / CONSTS.FT3_PSIA_TO_BTU

        Cp_total_Btu = Cp_IG + Hdep_deriv(degR, z_eos, dz_dT)

        # Analytic PR volume
        V = z_eos * CONSTS.R * degR / psia
        dP_dT_constV = CONSTS.R / (V - b_mix) - da_mix_dT / (V**2 + 2*b_mix*V - b_mix**2)
        dV_dT = (CONSTS.R / psia) * (z_eos + degR * dz_dT)
        den = (V ** 2 + 2 * b_mix * V - b_mix ** 2) ** 2
        dP_dV_constT = (- CONSTS.R * degR / (V - b_mix) ** 2 + 2 * a_mix * (V + b_mix) / den)
        Cv_total_Btu = Cp_total_Btu + (degR * dV_dT ** 2 * dP_dV_constT) / CONSTS.FT3_PSIA_TO_BTU

        # JT coefficient (R/psia)
        mu_JT = (degR * dV_dT - V) / (Cp_total_Btu * CONSTS.FT3_PSIA_TO_BTU)

        # ----- Verbose Diagnostics -----
        if verbose:
            print("---- VERBOSE PR PROPERTIES ----")
            print(f"Inputs: degF={degF:.3f}, psia={psia:.3f}, sg={sg:.3f}, co2={co2:.3f}, h2s={h2s:.3f}, n2={n2:.3f}, h2={h2:.3f}")
            print(f"degR: {degR:.5f}")
            print("z_fractions:", zf)
            print("tcs_local:", tcs)
            print("pcs_local:", pcs)
            print("trs:", trs)
            print("prs:", prs)
            print("kij:\n", kij)
            print("dkij_dT:\n", dkij_dT)
            print("d2kij_dT2:\n", d2kij_dT2)
            print("a_c_i:", a_c_i)
            print("m_i:", m_i)
            print("alpha:", alpha)
            print("a_i:", a_i)
            print("b_i:", b_i)
            print("a_mix:", a_mix)
            print("b_mix:", b_mix)
            print("A:", A)
            print("B:", B)
            print("z_eos:", z_eos)
            print("z_vshift:", z_vshift)
            print("da_mix_dT:", da_mix_dT)
            print("d2a_mix_dT2:", d2a_mix_dT2)
            print("dz_dT:", dz_dT)
            print("Bdim:", Bdim)
            print("Cp_IG:", Cp_IG)
            print("H_dep_Btu:", H_dep_Btu)
            print("Cp_total_Btu:", Cp_total_Btu)
            print("V (ft³/lb-mol):", V)
            print("dP_dT_constV:", dP_dT_constV)
            print("dV_dT:", dV_dT)
            print("Cv_total_Btu:", Cv_total_Btu)
            print("JT:", mu_JT)
        result.update({
            "H": float(H_total_Btu),
            "Cp": float(Cp_total_Btu),
            #"Cv": float(Cv_total_Btu),
            "JT": float(mu_JT)
        })
    if viscosity:
        result['Viscosity'] = float(lbc_viscosity(z_vshift, degF, psia, sg, co2, h2s, n2, h2, AG))
        if verbose:
            print("Viscosity (cP):", result["Viscosity"])
    
    if Metric:
        if density:
            result['Density'] *= 16.01846337396 # kg/m 
        if thermo:
            result["H"] *= 2.326 # kJ/(kmol)
            result["Cp"] *= 4.186800585 # kJ/(kmol·K)
            result["Cv"] *= 4.186800585 # kJ/(kmol·K)
            result["JT"] *= 80.576521 # degC/MPa
    
    return result