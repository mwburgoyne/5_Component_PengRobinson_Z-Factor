import numpy as np

###############################################################################
# Global Constants and Arrays
###############################################################################
R      = 10.731577089016
mwAir  = 28.97
degF2R = 459.67

names   =          ['CO2',     'H2S',      'N2',     'H2',     'Hydrocarbon']
mws     = np.array([44.01,     34.082,    28.014,    2.016,     0.0         ])
tcs     = np.array([547.416,   672.120,   227.160,   47.430,    1.0         ])  # H2 Tc modified
pcs     = np.array([1069.51,   1299.97,   492.84,    187.53,    1.0         ])
ACF     = np.array([0.12256,   0.04916,   0.037,    -0.21700,  -0.03899     ])
VSHIFT  = np.array([-0.27593, -0.22896,  -0.21066,  -0.32400,  -0.19076     ])
OmegaA  = np.array([0.427705,  0.436743,  0.457236,  0.457236,  0.457236    ])
OmegaB  = np.array([0.0696460, 0.0724373, 0.0777961, 0.0777961, 0.0777961   ])
VCVIS   = np.array([1.46020,   1.46460,   1.35422,   0.67967,   0.0         ])  # cuft/lbmol

###############################################################################
# 1) Analytic Cubic Solver 
###############################################################################
def cubic_root(a, flag=0):
    """
    Returns real root(s) of the cubic polynomial a[0]*Z^3 + a[1]*Z^2 + a[2]*Z + a[3] = 0.
      - If flag=1, returns the maximum real root
      - If flag=-1, returns the minimum real root
      - If flag=0, returns all real roots (as a NumPy array)
    """
    if a[0] != 1:
        a = np.asarray(a, dtype=float) / a[0]

    a0, a1, a2, a3 = a  # after normalization, a0 = 1
    p = (3*a2 - a1**2) / 3.0
    q = (2*a1**3 - 9*a1*a2 + 27*a3) / 27.0
    disc = (q**2 / 4.0) + (p**3 / 27.0)

    if disc < 0:
        m = 2.0 * np.sqrt(-p/3.0)
        qpm = 3.0*q / (p*m) if abs(p*m) > 1e-14 else 0.0
        qpm = np.clip(qpm, -1, 1)   # clamp for arccos domain
        theta = np.arccos(qpm) / 3.0

        roots = m * np.array([
            np.cos(theta),
            np.cos(theta + 4.0*np.pi/3.0),
            np.cos(theta + 2.0*np.pi/3.0),
        ])
        roots -= a1 / 3.0
    else:
        def real_cubic(x):
            return np.sign(x) * abs(x)**(1.0/3.0)

        temp1 = -q/2.0 + np.sqrt(disc)
        temp2 = -q/2.0 - np.sqrt(disc)
        P = real_cubic(temp1)
        Q = real_cubic(temp2)
        roots = np.array([P + Q]) - a1 / 3.0

    if flag == -1:
        return np.min(roots)
    if flag == 1:
        return np.max(roots)
    return roots


###############################################################################
# 2) Pseudo-Critical Temperature/Pressure 
###############################################################################
def tc_pc(sg_hc):
    """
    Computes (Tpc, Ppc) from a custom correlation for a hydrocarbon SG.
    """
    v1, v2, v3 = 0.000229975, 0.186415901, 2.470903632
    offset, pl, vl = -0.032901049, 42.6061669, 3007.108548

    mw_gas   = mwAir * sg_hc
    vc_on_zc = v1 * mw_gas**2 + v2*mw_gas + v3
    tpc_hc   = (offset + vc_on_zc) * vl / (offset + vc_on_zc + pl)
    ppc_hc   = tpc_hc * R / vc_on_zc
    return tpc_hc, ppc_hc


###############################################################################
# 3) calc_bips 
###############################################################################
def calc_bips(hc_mw, degf):
    """
    Returns a 5x5 BIP matrix for [CO2, H2S, N2, H2, HC] at specified hc_mw, degF.
    """
    degR = degf + degF2R

    intcpts     = np.array([0.386557, 0.267007, 0.486589, 0.776917])
    mw_slopes   = np.array([-0.00219806, -0.00396541, -0.00316789, 0.0106061])
    degR_slopes = np.array([-158.333, -58.611, -226.239, -474.283])

    hc_bips   = intcpts + degR_slopes/degR + mw_slopes*hc_mw
    inert_bips = np.array([0.0600319, -0.229807, -0.18346, 0.646796, 0.65, 0.369087])
    bips      = np.concatenate([hc_bips, inert_bips])

    bip_pairs = [
        (0, 4), (1, 4), (2, 4), (3, 4),
        (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)
    ]

    bip_matrix = np.zeros((5, 5))
    for idx, (i, j) in enumerate(bip_pairs):
        bip_matrix[i, j] = bips[idx]
        bip_matrix[j, i] = bips[idx]

    return bip_matrix


###############################################################################
# 4) Peng-Robinson Z-Factor
###############################################################################
def peng_robinson_z(degf, psia, sg, co2=0.0, h2s=0.0, n2=0.0, h2=0.0):
    """
    Computes the Peng-Robinson Z-factor for a gas mixture:
      y = [co2, h2s, n2, h2, hydrocarbon].
    Returns the maximum real root minus volume shift.
    """
    # Basic checks
    if any(x < 0 for x in (co2, h2s, n2, h2)) or (co2 + h2s + n2 + h2 > 1):
        return None

    # Mole fraction array
    y = np.array([co2, h2s, n2, h2, 1 - (co2 + h2s + n2 + h2)])
    degR = degf + degF2R

    # inert_sg = sum(inert fraction * MW_i / MW_air)
    # sum_inerts = co2 + h2s + n2 + h2 = np.sum(y[:-1])
    sum_inerts = np.sum(y[:-1])
    inert_sg   = np.sum(y[:-1] * mws[:-1]) / mwAir  # partial SG from inerts

    # Effective HC SG
    if sum_inerts < 1.0:
        sg_hc = (sg - inert_sg) / (1.0 - sum_inerts)
    else:
        sg_hc = 0.75                            # Dummy value for when inert fraction == 1
    sg_hc = min(max(sg_hc, 0.553779772), 1.75)  # Methane lower limit. 1.75 upper limit for hydrocarbon gas
    hc_mw = sg_hc * mwAir

    # Update hydrocarbon Tcs, Pcs
    tcs[-1], pcs[-1] = tc_pc(sg_hc)

    # Reduced T/P
    trs = degR / tcs
    prs = psia / pcs

    # BIPs
    kij = calc_bips(hc_mw, degf)

    # PR alpha factors
    m_factors = 0.37464 + 1.54226*ACF - 0.26992*(ACF**2)
    alpha     = (1.0 + m_factors*(1.0 - np.sqrt(trs)))**2

    # Ai, Bi
    Ai = OmegaA * alpha * prs / (trs**2)      # Eq 4.12c of Phase Behavior Monograph (Whitson et al)
    Bi = OmegaB * prs / trs                   # Eq 4.12d of Phase Behavior Monograph (Whitson et al)

    # A, B per Eq 4.16 of Phase Behavior Monograph (Whitson et al)
    A = np.sum(y[:, None] * y * np.sqrt(np.outer(Ai, Ai)) * (1.0 - kij))
    B = np.sum(y * Bi)

    # PR cubic:  [1, -(1-B), (A - 3B^2 - 2B), -(A*B - B^2 - B^3)]
    coeffs = [
        1.0,
        -(1.0 - B),
        A - 3.0*B**2 - 2.0*B,
        -(A*B - B**2 - B**3)
    ]

    # Solve and apply volume shift
    z_factor = cubic_root(coeffs, flag=1)
    vol_shift = np.sum(y * VSHIFT * Bi)
    return z_factor - vol_shift


###############################################################################
# 5) LBC (Lorenz-Bray-Clark) Viscosity
###############################################################################
def lbc(Z, degf, psia, sg, co2=0.0, h2s=0.0, n2=0.0, h2=0.0):
    """
    Computes gas viscosity via the Lorenz-Bray-Clark (LBC) correlation, using
    Stiel-Thodos for pure-component viscosities + mixing rules (Herning & Zippener).
    """
    # Basic checks
    if any(x < 0 for x in (co2, h2s, n2, h2)) or (co2 + h2s + n2 + h2 > 1):
        return None

    degR = degf + degF2R
    
    y = np.array([co2, h2s, n2, h2, 1 - (co2 + h2s + n2 + h2)])

    # inert_sg = sum(y[:-1] * mws[:-1]) / mwAir
    sum_inerts = np.sum(y[:-1])
    inert_sg   = np.sum(y[:-1] * mws[:-1]) / mwAir

    if sum_inerts < 1.0:
        sg_hc = (sg - inert_sg) / (1.0 - sum_inerts)
    else:
        sg_hc = 0.75
    sg_hc = min(max(sg_hc, 0.553779772), 1.75)  # Methane lower limit. 1.75 upper limit for hydrocarbon gas
    hc_gas_mw = sg_hc * mwAir

    # Update arrays for the hydrocarbon slot
    mws[-1]  = hc_gas_mw
    tcs[-1], pcs[-1] = tc_pc(hc_gas_mw / mwAir)

    # Update VCVIS for Hydrocarbon
    def vcvis_hc(mw):
        return 0.057511062*mw + 0.478400158
    VCVIS[-1] = vcvis_hc(hc_gas_mw)

    # Stiel-Thodos for pure-component viscosity
    def stiel_thodos(tempR, mol_wts):
        Tr_array = tempR / tcs
        Tc_k     = tcs * 5.0 / 9.0   # Degrees K
        Pc_atm   = pcs / 14.696    # Atmospheres

        base = (Tc_k**(1/6.0)) / (np.sqrt(mol_wts) * (Pc_atm**(2/3)))
        ui_list = []
        for i, Tr in enumerate(Tr_array):
            if Tr <= 1.5:
                ui_list.append(34e-5 * Tr**0.94 / base[i])
            else:
                ui_list.append(17.78e-5 * (4.58*Tr - 1.67)**(5/8) / base[i])
        return np.array(ui_list)

    ui = stiel_thodos(degR, mws)

    # Herning & Zippener "dilute gas" mixture viscosity
    def dilute_visc(y, pure_visc, mol_wts):
        sqrt_mw = np.sqrt(mol_wts)
        return np.sum(y * pure_visc * sqrt_mw) / np.sum(y * sqrt_mw)

    # LBC polynomial (4th and 5th terms have been adjusted)
    a = [0.1023, 0.023364, 0.058533, -3.92835e-02, 9.28591e-03]
    rhoc = 1.0 / np.sum(VCVIS * y)
    rhor = (psia / (Z * R * degR)) / rhoc

    lhs = sum(a[i] * (rhor ** i) for i in range(len(a)))

    Tc_k   = tcs * 5.0/9.0     # Degrees Kelvin
    Pc_atm = pcs / 14.696      # Atmospheres
    sum_Tc = np.sum(y * Tc_k)
    sum_Pc = np.sum(y * Pc_atm)
    sum_mw = np.sum(y * mws)

    mixture_eta = (sum_Tc**(1/6.0)) / (np.sqrt(sum_mw) * (sum_Pc**(2/3)))

    # Combine the LBC correlation
    diluted = dilute_visc(y, ui, mws)
    return (lhs**4 - 1e-4)/mixture_eta + diluted
