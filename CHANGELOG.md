# Changes since SPE-229932-MS

The paper (Burgoyne, Nielsen and Stanko, ADIPEC 2025, [SPE-229932-MS](https://doi.org/10.2118/229932-MS)) describes the five-component model. `Code Examples/Original (SPE-229932-MS)/` holds the code exactly as published with it, the repository as it stood on 7 November 2025 (commit `9454963`), defects included. `Code Examples/Latest (with Helium)/` carries every change below. Newest first; each entry says what changed, why, and what it does to the numbers.

## October 2026: implementations aligned with the regression basis

Porting helium into all four languages meant cross-checking them case by case for the first time. They did not fully agree with each other, nor with the regression that defines the model (`Data/04 - Hydrocarbon Tc & Pc fitting.zip`). The defects below were fixed in `Latest` only; `Original` keeps them, as published.

| Implementation | Defect in `Original` | Fix in `Latest` |
|---|---|---|
| Python | With `AG=True`, the hydrocarbon-inert BIPs used the gas-condensate hydrocarbon Tc while the EOS used the associated-gas Tc | BIPs use the same AG-dependent Tc as the EOS |
| Fortran | The hydrocarbon-inert BIPs used methane's Tc (343 °R) for every gas | BIPs use the hydrocarbon pseudo-critical Tc |
| VBA | CO2 acentric factor 0.12256 (regression value 0.12253); the hydrocarbon Tc/Pc correlation used an air MW of 28.9625 (28.97 elsewhere) and lacked the floor at methane; the LBC hydrocarbon SG floor was 16.043/28.97 rather than 16.0425/28.97 | All three set to the regression values |
| Rust | None found | |

The regression decks settle the BIP convention: in all 115 regressed samples the BIPs were built from the AG-dependent hydrocarbon Tc. The Rust and VBA already did this.

**Effect on the numbers.** Python results change only for `AG=True` with inerts present. For gas SG 0.65-1.0 (60-300 °F, 500-15,000 psia), Z changes by up to 1.7% (median 0.3-0.7% at SG 0.8-1.0), viscosity by up to 5.1%, and H by up to 50 Btu/lb-mol. At SG 1.4 near the phase boundary, the corrected BIPs can switch the selected cubic root, giving step changes. The Fortran changes for every gas heavier than methane that contains inerts: up to 0.22% in Z and 2.3% in H on the test grid. The VBA corrections are small: up to 0.023% in Z and density, 0.043% in viscosity and 0.46% in H. A fourth VBA item, the methane Vc/Zc constant now computed as R·Tc/Pc instead of rounded to 5.51852587149056, changes results by about 1e-10.

**Result.** On a 252-case grid covering compositions with and without helium, SG 0.65-0.8, AG on and off, 60-300 °F and 14.7-15,000 psia, Python, Rust and Fortran agree to 1e-12 relative in Z, density, Cp, Cv, JT and viscosity, and to 1e-12 Btu/lb-mol in H. The VBA cannot be run outside Excel. A line-by-line Python translation of the `Latest` module agrees with the same grid to the same tolerance, but the module itself has not yet been compiled and run in Excel.

## October 2026: helium as a sixth component

**What.** Helium is added as a sixth component, giving the order [CO2, H2S, N2, H2, He, Gas]. Every implementation (Python, Rust, Fortran, VBA module) takes a helium mole fraction as a new last argument (`he`; an `Optional` argument in VBA), so existing calls are unchanged. All helium binary interaction parameters are zero; only pure-helium properties were regressed.

**Helium constants.** MW 4.003, Tc 6.35 °R, Pc 32.9236 psia, acentric factor −0.17984, volume shift −0.078082, PR default Ω_A and Ω_B, VcVis 0.76778 ft³/lb-mol, ideal-gas Cp/R = 5/2 exactly, reference enthalpy offset 0.425547 Btu/lb-mol.

> **The helium critical temperature is deliberately non-standard: 6.35 °R, not the NIST 9.35 °R (5.195 K). It is not a typo; do not "correct" it.** In this model Tc also sets the Stiel-Thodos dilute-gas viscosity inside LBC. At the true Tc that term is about 25% low for helium and VcVis cannot repair it, while the density fit is insensitive to Tc because the volume shift and acentric factor absorb it. The helium Tc, acentric factor, volume shift and VcVis are one regressed set. H2's 47.43 °R (against NIST 59.7 °R) is the same kind of effective constant. The acentric factor is also bounded so that the PR alpha term stays well behaved to 500 °F.

**Basis.** Regressed in PhazeComp 1.81 to 29,400 NIST WebBook points for pure helium (49 isotherms, 60-300 °F; 14.7-14,990 psia), volume shift and acentric factor to molar density first, then VcVis to viscosity with the LBC coefficients held at their published values. Decks and inputs: `Data/06 - Helium EOS and VcVis Regression.zip`.

| Pure helium vs NIST (mean / max absolute error) | Tc 6.35 °R (used) | Tc 9.34 °R (standard) |
|---|---|---|
| Density | 0.31% / 0.84% | 0.32% / 0.86% |
| Viscosity | 0.70% / 3.1% | 18.6% / 57% |
| Cp | 0.23% / 0.60% | |

Thermal outputs were not fitted: against NIST the pressure dependence of helium enthalpy is about 13% low (475 against 543 Btu/lb-mol at 60 °F and 14,990 psia, relative to 60 °F and 14.7 psia), and the Joule-Thomson coefficient is 7-19% low in magnitude.

**Effect on existing results: none.** Adding helium did not change any result when helium is zero: 320 reference cases across composition, 50-300 °F, 14.7-15,000 psia and both hydrocarbon correlations reproduce to 5e-15 relative. The implementation corrections above are separate.

**Implementations.** Python, Rust, Fortran and the VBA module agree to 1e-12 (see the entry above). The VBA module is exported as `bns_VBA_Module1.bas`; the workbook in `Latest` is not provided because VBA cannot be rewritten safely from outside Excel. Import the module by hand (see the README). The same helium model ships in pyResToolbox 3.8.3 (`he=` on the gas functions).

## August 2026: volume shift carried into enthalpy and the Joule-Thomson coefficient

**What was wrong.** The published implementations applied the Peneloux volume shift to Z-factor and density, but computed enthalpy, Cp, Cv and the Joule-Thomson coefficient from the untranslated EOS root.

**The correction.** The shift reduces to a constant molar volume offset, `c = SUM(z_i * VSHIFT_i * b_i)`. A constant translation leaves Cp, Cv and entropy unchanged, but shifts enthalpy by `−c·p` and the JT coefficient by `+c/Cp`. Enthalpy is reported relative to 60 °F and 14.696 psia, so the shift enters as `−c·(p − 14.696)` and the reference-state offsets stay valid.

**Effect on the numbers.** Z, density, Cp, Cv and viscosity are unchanged to machine precision. **Enthalpy and the JT coefficient change.** In the README's field-units worked example (120 °F, 2000 psia, sg 0.8 with CO2, H2S, N2 and H2), H moved by 5.8% and JT by 3.4%. Against reference EOS (CoolProp) over 60-300 °F and 100-10,000 psia, the mean JT bias went from +4.2% to +0.7% for pure CO2 and from +12.0% to +5.9% for pure methane, and the mean absolute enthalpy-departure error for methane fell from 65 to 31 Btu/lb-mol. Enthalpy or JT values computed with the `Original` code, or with any copy from before August 2026, will differ accordingly.

**Implementations.** Applied in Python, Rust, Fortran and the VBA module, and in the workbook (commits `8dcda71`, `06bf561`, `c8a2f08`, worked examples regenerated in `e5d02e5`). A derivation of the volume shift as applied to Z-factor is in `Code Examples/Derivation of Volume Shift Applied to Z-Factor.pdf`.
