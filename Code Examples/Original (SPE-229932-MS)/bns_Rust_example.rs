#![allow(non_snake_case, dead_code, unused_variables, unused_mut)]

use std::f64::consts::SQRT_2;

// ===== Utilities =====
fn dot5(a: &[f64; 5], b: &[f64; 5]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
fn poly5_eval(coeffs_low_to_high: [f64; 5], x: f64) -> f64 {
    let mut p = 0.0;
    let mut xn = 1.0;
    for c in coeffs_low_to_high {
        p += c * xn;
        xn *= x;
    }
    p
}
fn clamp_min(x: f64, lo: f64) -> f64 { if x < lo { lo } else { x } }

// ===== Constants (Field units) =====
struct Constants {
    R: f64,            // ft^3·psia/(lb-mol·R)
    R_THERMO: f64,     // Btu/(lb-mol·R)
    MW_AIR: f64,
    DEG_F_TO_R: f64,
    FT3_PSIA_TO_BTU: f64,
    mwCH4: f64,
    TcCH4: f64,        // deg R
    PcCH4: f64,        // psia
    VcZcCH4: f64,      // ft3/lb-mol
}
const CONSTS: Constants = Constants {
    R: 10.731_577_089_016,
    R_THERMO: 1.98588,
    MW_AIR: 28.97,
    DEG_F_TO_R: 459.67,
    FT3_PSIA_TO_BTU: 5.403,
    mwCH4: 16.0425,
    TcCH4: 343.008,
    PcCH4: 667.029,
    VcZcCH4: (10.731_577_089_016 * 343.008) / 667.029,
};

// ===== Fixed coefficient tables (CO2, H2S, N2, H2, Gas(C1+)) =====
struct Props {
    mws:  [f64; 5],
    tcs:  [f64; 5],
    pcs:  [f64; 5],
    acf:  [f64; 5],
    vshift: [f64; 5],
    omega_a: [f64; 5],
    omega_b: [f64; 5],
    vcvis: [f64; 5],
    cp_poly: [[f64; 5]; 5], // Cp poly coeffs (K powers), low?high
}
fn base_props() -> Props {
    Props {
        mws:   [44.01,   34.082, 28.014, 2.016, 0.0],
        tcs:   [547.416, 672.120, 227.160, 47.430, 1.0],
        pcs:   [1069.51, 1299.97, 492.84, 187.53, 1.0],
        acf:   [0.12253, 0.04909, 0.037,  -0.217, -0.03899],
        vshift:[-0.27607,-0.22901,-0.21066,-0.36270,-0.19076],
        omega_a: [0.427671, 0.436725, 0.457236, 0.457236, 0.457236],
        omega_b: [0.0696397,0.0724345,0.0777961,0.0777961,0.0777961],
        vcvis: [1.46352, 1.46808, 1.35526, 0.68473, 1.44383], // last updated later for HC
        cp_poly: [
            [2.725473196,  0.004103751,  1.5602E-05, -4.19321E-08,  3.10542E-11], // CO2
            [4.446031265, -0.005296052,  2.0533E-05, -2.58993E-08,  1.25555E-11], // H2S
            [3.423811591,  0.001007461, -4.58491E-06,  8.4252E-09, -4.38083E-12], // N2
            [1.421468418,  0.018192108, -6.04285E-05,  9.08033E-08, -5.18972E-11],// H2
            [5.369051342, -0.014851371,  4.86358E-05, -3.70187E-08,  1.80641E-12],// Gas (C1+)
        ],
    }
}

// ===== Pseudocritical relations =====
fn tc_ag(x: f64) -> f64 { 2695.14765 * x / (274.341701 + x) + CONSTS.TcCH4 }
fn tc_gc(x: f64) -> f64 { 1098.10948 * x / (101.529237 + x) + CONSTS.TcCH4 }
fn pc_from_vc_zc(x: f64, vc_slope: f64, tc: f64) -> f64 {
    let vc_over_zc = vc_slope * x + CONSTS.VcZcCH4;
    CONSTS.R * tc / vc_over_zc
}
fn pseudo_critical(sg_hc: f64, ag: bool) -> (f64, f64) {
    let x = (CONSTS.MW_AIR * sg_hc - CONSTS.mwCH4).max(0.0);
    let (tpc, slope) = if ag { (tc_ag(x), 0.177497835) } else { (tc_gc(x), 0.170931432) };
    (tpc, pc_from_vc_zc(x, slope, tpc))
}
fn update_hydrocarbon(mut p: Props, sg_hc: f64, ag: bool) -> Props {
    let (tpc, ppc) = pseudo_critical(sg_hc, ag);
    let hc_mw = sg_hc * CONSTS.MW_AIR;
    p.tcs[4] = tpc;
    p.pcs[4] = ppc;
    p.mws[4] = hc_mw;
    p.vcvis[4] = 0.0576710 * (hc_mw - CONSTS.mwCH4) + 1.44383; // empirical
    p
}

fn get_z_fractions(co2: f64, h2s: f64, n2: f64, h2: f64) -> [f64; 5] {
    let hc = 1.0 - (co2 + h2s + n2 + h2);
    [co2, h2s, n2, h2, hc]
}
fn hydrocarbon_sg(sg_bulk: f64, zf: &[f64; 5], mws: &[f64; 5]) -> f64 {
    let frac_hc = zf[4];
    let sum_nonhc = zf[0]*mws[0] + zf[1]*mws[1] + zf[2]*mws[2] + zf[3]*mws[3];
    if frac_hc > 0.0 {
        (sg_bulk - sum_nonhc / CONSTS.MW_AIR) / frac_hc
    } else {
        0.75
    }
}

// ===== PR cubic helpers =====
fn fugacity_coeff(z: f64, a: f64, b: f64) -> f64 {
    let arg = (z + (1.0 + SQRT_2) * b) / (z + (1.0 - SQRT_2) * b);
    if arg <= 0.0 || (z - b) <= 0.0 { return f64::INFINITY; }
    let ln_phi = (z - 1.0) - (z - b).ln() - (a / (2.0 * SQRT_2 * b)) * arg.ln();
    ln_phi.exp()
}
fn solve_pr_cubic_select(coeffs: [f64; 4], a: f64, b: f64) -> (f64, bool) {
    let c2 = coeffs[1];
    let c1 = coeffs[2];
    let c0 = coeffs[3];

    // Depressed cubic t^3 + pt + q with z = t - c2/3
    let a_div3 = c2 / 3.0;
    let p = c1 - c2*c2/3.0;
    let q = (2.0*c2*c2*c2)/27.0 - (c2*c1)/3.0 + c0;
    let disc = (q*q)/4.0 + (p*p*p)/27.0;

    let mut roots: Vec<f64> = Vec::new();
    if disc > 0.0 {
        let sqrt_disc = disc.sqrt();
        let u = (-q/2.0 + sqrt_disc).cbrt();
        let v = (-q/2.0 - sqrt_disc).cbrt();
        let z = (u + v) - a_div3;
        if z.is_finite() && z > 0.0 { roots.push(z); }
    } else if disc.abs() <= 1e-14 {
        let u = (-q/2.0).cbrt();
        let z1 = (2.0*u) - a_div3;
        let z2 = (-u) - a_div3;
        if z1.is_finite() && z1 > 0.0 { roots.push(z1); }
        if z2.is_finite() && z2 > 0.0 { roots.push(z2); }
    } else {
        let r = (-p*p*p/27.0).sqrt();
        let phi = (-q/2.0) / r;
        let phi = phi.acos();
        let m = 2.0 * (-p/3.0).sqrt();
        for t in [
            m * (phi/3.0).cos(),
            m * ((phi + 2.0*std::f64::consts::PI)/3.0).cos(),
            m * ((phi + 4.0*std::f64::consts::PI)/3.0).cos(),
        ] {
            let z = t - a_div3;
            if z.is_finite() && z > 0.0 { roots.push(z); }
        }
    }

    if roots.is_empty() { return (1.0, true); }
    if roots.len() == 1 || !a.is_finite() || !b.is_finite() { return (roots[0], false); }
    // choose minimum fugacity coefficient (vapor-like root)
    let mut best = roots[0];
    let mut best_phi = fugacity_coeff(best, a, b);
    for &r in roots.iter().skip(1) {
        let phi = fugacity_coeff(r, a, b);
        if phi < best_phi { best = r; best_phi = phi; }
    }
    (best, false)
}
fn pr_m_param(acf: f64) -> f64 {
    0.37464 + 1.54226 * acf - 0.26992 * acf * acf
}

// ===== Temperature-dependent BIPs =====
struct Bips {
    kij: [[f64; 5]; 5],
    dkij_dT: [[f64; 5]; 5],
    d2kij_dT2: [[f64; 5]; 5],
}
fn calc_bips(degR: f64, tpc_hc: f64) -> Bips {
    #[derive(Clone, Copy)]
    struct Pair { c: f64, s: f64, tc: f64 }
    let mut out = Bips {
        kij: [[0.0;5];5],
        dkij_dT: [[0.0;5];5],
        d2kij_dT2: [[0.0;5];5],
    };
    let get = |i: usize, j: usize| -> Pair {
        match (i, j) {
            (4, 0)|(0, 4) => Pair{ c:-0.145561,  s: 0.276572,  tc: tpc_hc },
            (4, 1)|(1, 4) => Pair{ c: 0.16852,  s:-0.122378,  tc: tpc_hc },
            (4, 2)|(2, 4) => Pair{ c:-0.108,    s: 0.0605506, tc: tpc_hc },
            (4, 3)|(3, 4) => Pair{ c:-0.0620119,s: 0.0427873, tc: tpc_hc },
            (0, 1)|(1, 0) => Pair{ c: 0.248638, s:-0.138185,  tc: 547.416 },
            (0, 2)|(2, 0) => Pair{ c:-0.25,     s: 0.11602,   tc: 547.416 },
            (0, 3)|(3, 0) => Pair{ c:-0.247153, s: 0.16377,   tc: 547.416 },
            (1, 2)|(2, 1) => Pair{ c:-0.204414, s: 0.234417,  tc: 672.12 },
            (1, 3)|(3, 1) => Pair{ c: 0.0,      s: 0.0,       tc: 672.12 },
            (2, 3)|(3, 2) => Pair{ c:-0.166253, s: 0.0788129, tc: 227.16 },
            _ => Pair{ c: 0.0, s: 0.0, tc: 1.0 },
        }
    };
    for i in 0..5 {
        for j in 0..5 {
            if i == j { continue; }
            let p = get(i, j);
            let tr = degR / p.tc;
            out.kij[i][j] = p.c + p.s / tr;
            out.dkij_dT[i][j] = -p.s * p.tc / (degR*degR);
            out.d2kij_dT2[i][j] = 2.0 * p.s * p.tc / (degR*degR*degR);
        }
    }
    out
}

// ===== Stiel–Thodos dilute viscosity =====
fn stiel_thodos_viscosity(
    t_rankine: f64,
    mws: &[f64;5],
    tcs: &[f64;5],
    pcs: &[f64;5]
) -> [f64;5] {
    let mut mu = [0.0;5];
    for i in 0..5 {
        let tr = t_rankine / tcs[i];
        let tc_k = tcs[i] * (5.0/9.0);
        let pc_atm = pcs[i] / 14.696;
        let eta = tc_k.powf(1.0/6.0) / (mws[i].sqrt() * pc_atm.powf(2.0/3.0));
        mu[i] = if tr <= 1.5 {
            34e-5 * tr.powf(0.94) / eta
        } else {
            17.78e-5 * ((4.58*tr) - 1.67).powf(5.0/8.0) / eta
        };
    }
    mu
}

// ===== Methane-adjust C1+ Cp polynomial vs HC MW =====
fn methane_adjust(cp: &[[f64;5];5], hc_mw: f64) -> [[f64;5];5] {
    let a0 = [7.8570E-04, 1.3123E-03, 9.8133E-04, 1.6463E-03, 1.7306E-02];
    let a1 = [-8.1649E-03, 5.5485E-03, 8.3258E-02, 2.0635E-01, 2.5551E+00];
    let x = hc_mw - CONSTS.mwCH4;
    let mut out = *cp;
    for k in 0..5 {
        let scale = a0[k]*x*x + a1[k]*x + 1.0;
        out[4][k] = cp[4][k] * scale;
    }
    out
}

// ===== LBC viscosity (mixture) =====
fn lbc_viscosity(
    z: f64, degF: f64, psia: f64, sg: f64,
    co2: f64, h2s: f64, n2: f64, h2: f64, ag: bool
) -> f64 {
    let zf = get_z_fractions(co2, h2s, n2, h2);
    let mut props = base_props();
    let sg_hc = {
        let sg_hc0 = hydrocarbon_sg(sg, &zf, &props.mws);
        clamp_min(sg_hc0, CONSTS.mwCH4 / CONSTS.MW_AIR)
    };
    let hc_mw = sg_hc * CONSTS.MW_AIR;
    props = update_hydrocarbon(props, sg_hc, ag);
    let mws = props.mws;
    let tcs = props.tcs;
    let pcs = props.pcs;
    let vcvis = props.vcvis;

    let degR = degF + CONSTS.DEG_F_TO_R;

    // dilute-gas mixture viscosity
    let mu_i = stiel_thodos_viscosity(degR, &mws, &tcs, &pcs);
    let sqrt_mw: Vec<f64> = mws.iter().map(|x| x.sqrt()).collect();
    let numerator: f64 = zf.iter().enumerate().map(|(i, zi)| zi * mu_i[i] * sqrt_mw[i]).sum();
    let denom: f64 = zf.iter().enumerate().map(|(i, zi)| zi * sqrt_mw[i]).sum();
    let mu_dilute = numerator / denom;

    // reduced density mixing
    let rhoc = 1.0 / dot5(&vcvis, &zf);
    let gas_density = psia / (z * CONSTS.R * degR);
    let rhor = gas_density / rhoc;

    // dense-phase polynomial correction
    let a = [0.1023, 0.023364, 0.058533, -0.0392852, 0.00926279];
    let lhs = a[0] + a[1]*rhor + a[2]*rhor.powi(2) + a[3]*rhor.powi(3) + a[4]*rhor.powi(4);

    let tc_mix_k = dot5(&tcs, &zf) * (5.0/9.0);
    let pc_mix_atm = dot5(&pcs, &zf) / 14.696;
    let mw_mix = dot5(&mws, &zf);
    let eta_mix = tc_mix_k.powf(1.0/6.0) / (mw_mix.sqrt() * pc_mix_atm.powf(2.0/3.0));

    (lhs.powi(4) - 1e-4)/eta_mix + mu_dilute // cP
}

// ===== Result struct =====
#[derive(Debug, Clone, Default)]
struct PrResult {
    z: f64,
    density: Option<f64>,
    h: Option<f64>,
    cp: Option<f64>,
    cv: Option<f64>,
    jt: Option<f64>,
    viscosity: Option<f64>,
}

// ===== Main API =====
#[allow(clippy::too_many_arguments)]
fn pr_properties(
    temp: f64,       // degF (Metric=false) or degC (Metric=true)
    pres: f64,       // psia or MPa
    sg: f64,
    co2: f64, h2s: f64, n2: f64, h2: f64,
    ag: bool,
    viscosity: bool,
    density: bool,
    thermo: bool,
    metric: bool,
    verbose: bool
) -> Result<PrResult, String> {
    if co2 + h2s + n2 + h2 > 1.0 + 1e-12 {
        return Err(format!(
            "Invalid composition: CO2({})+H2S({})+N2({})+H2({}) > 1.0",
            co2, h2s, n2, h2
        ));
    }
    // Normalize to Field
    let (degF, psia) = if metric {
        (temp * 9.0/5.0 + 32.0, pres * 145.0377)
    } else { (temp, pres) };

    let zf = get_z_fractions(co2, h2s, n2, h2);
    let mut props = base_props();

    // Hydrocarbon sg and MW
    let mut sg_hc = hydrocarbon_sg(sg, &zf, &props.mws);
    sg_hc = clamp_min(sg_hc, CONSTS.mwCH4 / CONSTS.MW_AIR);
    let hc_mw = sg_hc * CONSTS.MW_AIR;
    props = update_hydrocarbon(props, sg_hc, ag);

    // Unpack
    let mut mws = props.mws;
    let tcs = props.tcs;
    let pcs = props.pcs;
    let acf = props.acf;
    let vshift = props.vshift;
    let omega_a = props.omega_a;
    let omega_b = props.omega_b;

    mws[4] = hc_mw;

    let degR = degF + CONSTS.DEG_F_TO_R;

    let (tpc_hc, _ppc_hc) = pseudo_critical(sg_hc, ag);
    let mut trs = [0.0; 5];
    let mut prs_r = [0.0; 5];
    for i in 0..5 {
        trs[i] = degR / tcs[i];
        prs_r[i] = psia / pcs[i];
    }

    let bips = calc_bips(degR, tpc_hc);

    // PR alpha and parameters
    let mut m_i = [0.0;5];
    for i in 0..5 { m_i[i] = pr_m_param(acf[i]); }
    let mut alpha = [0.0;5];
    for i in 0..5 { alpha[i] = (1.0 + m_i[i] * (1.0 - trs[i].sqrt())).powi(2); }
    let mut a_c_i = [0.0;5];
    let mut a_i = [0.0;5];
    let mut b_i = [0.0;5];
    for i in 0..5 {
        a_c_i[i] = omega_a[i] * CONSTS.R * CONSTS.R * tcs[i] * tcs[i] / pcs[i];
        a_i[i]   = a_c_i[i] * alpha[i];
        b_i[i]   = omega_b[i] * CONSTS.R * tcs[i] / pcs[i];
    }

    // a_mix, b_mix with kij
    let mut a_mix = 0.0;
    for i in 0..5 {
        for j in 0..5 {
            let aij = (a_i[i]*a_i[j]).sqrt() * (1.0 - bips.kij[i][j]);
            a_mix += zf[i] * zf[j] * aij;
        }
    }
    let b_mix = dot5(&b_i, &zf);

    // A, B and cubic
    let rt = CONSTS.R * degR;
    let a_dim = a_mix * psia / (rt*rt);
    let b_dim = b_mix * psia / rt;

    // PR cubic (monic): z^3 - (1-B) z^2 + (A - 3B^2 - 2B) z - (A B - B^2 - B^3) = 0
    let c2 = -(1.0 - b_dim);
    let c1 = a_dim - 3.0*b_dim*b_dim - 2.0*b_dim;
    let c0 = -(a_dim*b_dim - b_dim*b_dim - b_dim*b_dim*b_dim);
    let (z_eos, _ideal_fallback) = solve_pr_cubic_select([1.0, c2, c1, c0], a_dim, b_dim);

    // volume-shifted Z
    let mut bi = [0.0;5];
    for i in 0..5 { bi[i] = omega_b[i] * (prs_r[i] / trs[i]); }
    let mut shift = 0.0;
    for i in 0..5 { shift += zf[i] * vshift[i] * bi[i]; }
    let z_vshift = z_eos - shift;

    let mut res = PrResult { z: z_vshift, ..Default::default() };

    // Density (lbm/ft^3)
    if density {
        let mwt = dot5(&mws, &zf);
        let rho = mwt * psia / (z_vshift * CONSTS.R * degR);
        res.density = Some(rho);
    }

    if thermo {
        // Derivatives for departures
        let mut d_alpha_dT = [0.0;5];
        let mut d2_alpha_dT2 = [0.0;5];
        for i in 0..5 {
            let sqrt_tr = trs[i].sqrt();
            let d_alpha_dTr = -m_i[i] * (1.0 + m_i[i]*(1.0 - sqrt_tr)) / sqrt_tr;
            d_alpha_dT[i] = d_alpha_dTr / tcs[i];

            let term1 = 2.0 * (-m_i[i]/(2.0*sqrt_tr)).powi(2);
            let term2 = 2.0 * (1.0 + m_i[i]*(1.0 - sqrt_tr)) * (m_i[i] / (4.0*trs[i].powf(1.5)));
            let d2alpha_dTr2 = term1 + term2;
            d2_alpha_dT2[i] = d2alpha_dTr2 / (tcs[i]*tcs[i]);
        }
        let mut da_i_dT = [0.0;5];
        let mut d2a_i_dT2 = [0.0;5];
        for i in 0..5 {
            da_i_dT[i]  = a_c_i[i] * d_alpha_dT[i];
            d2a_i_dT2[i]= a_c_i[i] * d2_alpha_dT2[i];
        }

        let mut da_mix_dT = 0.0;
        let mut d2a_mix_dT2 = 0.0;
        for i in 0..5 {
            for j in 0..5 {
                let sqrt_ai_aj = (a_i[i]*a_i[j]).sqrt();
                let n = da_i_dT[i]*a_i[j] + a_i[i]*da_i_dT[j];

                let daij_dT = -bips.dkij_dT[i][j]*sqrt_ai_aj
                    + (1.0 - bips.kij[i][j]) * 0.5 * (n / sqrt_ai_aj);
                da_mix_dT += zf[i]*zf[j]*daij_dT;

                let cross_2_da = 2.0 * da_i_dT[i] * da_i_dT[j];
                let d2a_term = d2a_i_dT2[i]*a_i[j] + cross_2_da + a_i[i]*d2a_i_dT2[j];
                let d = sqrt_ai_aj;
                let term1 = - (n / d) * bips.dkij_dT[i][j];
                let term2 = (1.0 - bips.kij[i][j]) * 0.5 * ((d2a_term / d) - (n*n)/(2.0*d*d*d));
                let term3 = - d * bips.d2kij_dT2[i][j];
                let d2aij_dT2 = term1 + term2 + term3;
                d2a_mix_dT2 += zf[i]*zf[j]*d2aij_dT2;
            }
        }

        // dZ/dT|P via implicit differentiation (non-shifted z_eos)
        let a_ = a_mix;
        let b_ = b_mix;
        let aA = a_ * psia / (CONSTS.R * degR).powi(2);
        let bB = b_ * psia / (CONSTS.R * degR);

        let dB_dT = -bB / degR;
        let dA_dT = (psia / (CONSTS.R * degR).powi(2)) * da_mix_dT
            - 2.0 * a_ * psia / (CONSTS.R*CONSTS.R * degR.powi(3));

        let dF_dz = 3.0*z_eos*z_eos - 2.0*(1.0 - bB)*z_eos + (aA - 3.0*bB*bB - 2.0*bB);
        let dF_dA = z_eos - bB;
        let dF_dB = z_eos*z_eos - (6.0*bB + 2.0)*z_eos - aA + 2.0*bB + 3.0*bB*bB;
        let dF_dT = dF_dA * dA_dT + dF_dB * dB_dT;
        let dz_dT = - dF_dT / dF_dz;

        // Ideal-gas Cp from polynomials (with methane-adjusted C1+)
        let cp_poly = methane_adjust(&props.cp_poly, hc_mw);
        let t_k = degR * 5.0/9.0;

        let mut cp_vals = [0.0;5];
        for i in 0..5 { cp_vals[i] = poly5_eval(cp_poly[i], t_k); }
        let cp_ig = dot5(&cp_vals, &zf) * CONSTS.R_THERMO; // Btu/(lb-mol·R)

        // Enthalpy departure
        let bdim = (b_mix * psia) / (CONSTS.R * degR);
        let log_term = |z: f64| ((z + (SQRT_2 + 1.0)*bdim) / (z - (SQRT_2 - 1.0)*bdim)).ln();
        let h_dep_ft3psia =
            CONSTS.R*degR*(z_eos - 1.0)
          + (degR*da_mix_dT - a_mix) / (2.0*SQRT_2*b_mix) * log_term(z_eos);
        let h_dep_btu = h_dep_ft3psia / CONSTS.FT3_PSIA_TO_BTU;

        // Reference enthalpy @ 60°F
        let mut h0 = [-16.6022, -21.5512, -3.57757, 0.008054, 0.0];
        h0[4] = -0.015774*(hc_mw - CONSTS.mwCH4).powi(2)
              - 0.646645*(hc_mw - CONSTS.mwCH4)
              - 8.2551915;
        let h0_mix = dot5(&h0, &zf);

        // Ideal-gas enthalpy change from 60°F (integrate Cp)
        let t_ref_f = 60.0;
        let t_ref_k = (t_ref_f + CONSTS.DEG_F_TO_R) * 5.0/9.0;
        let mut h_ig_btu = 0.0;
        for i in 0..5 {
            let a = cp_poly[i];
            let term =
                a[0]*(t_k - t_ref_k)
              + a[1]/2.0*(t_k*t_k - t_ref_k*t_ref_k)
              + a[2]/3.0*(t_k.powi(3) - t_ref_k.powi(3))
              + a[3]/4.0*(t_k.powi(4) - t_ref_k.powi(4))
              + a[4]/5.0*(t_k.powi(5) - t_ref_k.powi(5));
            h_ig_btu += zf[i] * CONSTS.R_THERMO * 9.0/5.0 * term;
        }
        let h_total_btu = h_ig_btu + h_dep_btu - h0_mix;

        // Cp departure = d(H_dep)/dT
        let dBdim_dT = - bdim / degR;
        let n = z_eos + (SQRT_2 + 1.0)*bdim;
        let d = z_eos - (SQRT_2 - 1.0)*bdim;
        let dln_dT = (dz_dT + (SQRT_2 + 1.0)*dBdim_dT)/n - (dz_dT - (SQRT_2 - 1.0)*dBdim_dT)/d;
        let x = (degR*da_mix_dT - a_mix) / (2.0*SQRT_2*b_mix);
        let dx_dT = (degR * d2a_mix_dT2) / (2.0*SQRT_2*b_mix);
        let dHdep_dT_btu = ( CONSTS.R*(z_eos - 1.0) + CONSTS.R*degR*dz_dT
                           + dx_dT * (n/d).ln() + x * dln_dT ) / CONSTS.FT3_PSIA_TO_BTU;

        let cp_total_btu = cp_ig + dHdep_dT_btu;

        // Cv and JT
        let v = z_eos * CONSTS.R * degR / psia;
        let dV_dT = (CONSTS.R / psia) * (z_eos + degR*dz_dT);
        let den = (v*v + 2.0*b_mix*v - b_mix*b_mix).powi(2);
        let dP_dV_constT = - CONSTS.R*degR / (v - b_mix).powi(2) + 2.0*a_mix*(v + b_mix) / den;
        let cv_total_btu = cp_total_btu + (degR * dV_dT*dV_dT * dP_dV_constT) / CONSTS.FT3_PSIA_TO_BTU;

        let mu_jt = (degR * dV_dT - v) / (cp_total_btu * CONSTS.FT3_PSIA_TO_BTU);

        res.h  = Some(h_total_btu);
        res.cp = Some(cp_total_btu);
        res.cv = Some(cv_total_btu);
        res.jt = Some(mu_jt);
        if verbose {
            eprintln!("DEBUG: z_eos={:.8}, z_vshift={:.8}, Cp_IG={:.6}, H_dep={:.6}", z_eos, z_vshift, cp_ig, h_dep_btu);
        }
    }

    if viscosity {
        let mu = lbc_viscosity(z_vshift, degF, psia, sg, co2, h2s, n2, h2, ag);
        res.viscosity = Some(mu);
    }

    // Metric output conversion (keep inputs normalized already)
    if metric {
        if let Some(rho) = res.density.as_mut() {
            *rho *= 16.018_463_373_96; // kg/m^3
        }
        if let Some(h) = res.h.as_mut()   { *h  *= 2.326; }              // kJ/kmol
        if let Some(cp)= res.cp.as_mut()  { *cp *= 4.186_800_585; }      // kJ/(kmol·K)
        if let Some(cv)= res.cv.as_mut()  { *cv *= 4.186_800_585; }      // kJ/(kmol·K)
        if let Some(jt)= res.jt.as_mut()  { *jt *= 80.576_521; }         // °C/MPa
        // viscosity already cP ~ mPa·s
    }

    Ok(res)
}

// ===== Demo run =====
fn main() {
    let out = pr_properties(
        /*temp*/ 88.0,     /*pres*/ 1080.0,     /*sg*/ 0.65,
        /*co2*/ 0.99, /*h2s*/ 0.0, /*n2*/ 0.00, /*h2*/ 0.0,
        /*AG*/ false,
        /*viscosity*/ true, /*density*/ true, /*thermo*/ true,
        /*Metric*/ false,   /*verbose*/ false
    ).expect("pr_properties failed");

    println!("Z   = {}", out.z);
    println!("rho = {:?} (lbm/ft^3)", out.density.unwrap());
    println!("mu  = {:?} (cP)", out.viscosity.unwrap());
    println!("H  = {:?} (Btu/(lb-mol))", out.h.unwrap());
    println!("Cp  = {:?} (Btu/(lb-mol·R))", out.cp.unwrap());
    println!("Cv  = {:?} (Btu/(lb-mol·R))", out.cv.unwrap());
    println!("JT  = {:?} (°F/psia)", out.jt.unwrap());
}