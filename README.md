# A Universal, EOS-Based Correlation for Z-Factor, Viscosity and Thermal Properties For Hydrocarbon and H<sub>2</sub>, N<sub>2</sub>, CO<sub>2</sub>, H<sub>2</sub>S Gas Mixtures

**Author**: Mark Burgoyne  
**First released**: 05-04-2024  
**Most recent update**: 13-08-2025

Example implementations in Python and Excel VBA of method outlined in ADIPEC 2025 paper, SPE-229932-MS

**Paper Abstract:** Accurate estimation of gas compressibility (Z-Factor), viscosity, and enthalpy is vital across conventional hydrocarbon production, greenhouse gas storage, and hydrogen energy applications. Legacy Z-factor correlations like Dranchuk-Abou-Kassem (DAK) and Hall-Yarborough (HY) perform well for typical natural gases but lose accuracy - and sometimes fail - at low temperatures, high pressures, or with inert-rich mixtures. Advanced frameworks such as GERG2008 address a broader range of compositions but are often impractical for spreadsheets and standard engineering workflows due to software, complexity, and compositional constraints.

We present a unified thermodynamic model combining a volume-translated Peng-Robinson Equation of State (EOS) for Z-factor, density, enthalpy, Cp, and Joule-Thomson coefficient with the Lohrenz-Bray-Clark (LBC) method for viscosity. The model covers one hydrocarbon and four inert components (H₂, N₂, CO₂, H₂S), with parameters tuned to Standing-Katz, NIST, binary vapor liquid equilibrium (VLE) and publicly available experimental Z-Factor data. Validation shows accuracy equal to or better than traditional correlations from 50-300 °F and 14.7-15,000 psia, including gas condensates and associated hydrocarbon gases as well as inert and hydrogen-rich mixtures. Viscosity results closely match reference data.

This approach streamlines gas property calculations with a single, analytically solved, spreadsheet-ready framework - facilitating integration into reservoir simulators and meeting the needs of both conventional and emerging carbon dioxide and hydrogen storage workflows.

Additionally, in its simplest form, this analytically solvable cubic EOS-based method can be a compact, direct and robust replacement for iterative legacy methods (e.g., DAK, HY), eliminating numerical convergence issues and simplifying integration into standard engineering workflows.

Calculations, regression data, and example implementations in Python and Excel VBA are available here.


## Example usage
```python
import bns as bns

# With Metric = True, temperature is in deg C, pressure in MPa and results are also in Metric units per below comments
bns.pr_properties(temp=48.88889, pres=13.78948965, sg=0.8,  co2=0.2, h2s=0.1, n2=0.02, h2=0.1, viscosity=True, density=True, thermo=True, Metric = True)

{'Z': 0.7941023149604872,
 'Density': 150.3029810272768,       # kg/m3    
 'H': -1240.8929074628306,           # kJ/(kmol)
 'Cp': 53.92887147310112,            # kJ/(kmol K)
 'Cv': 29.640351748197087,           # kJ/(kmol K)
 'JT': 2.9028170987704853,           # degC/MPa
 'Viscosity': 0.01790543522353429}   # mPa·s
 
# With Metric = False (or omitted), temperature in deg F, pressure in psia, and results return in Field units
bns.pr_properties(temp=120, pres=2000, sg=0.8,  co2=0.2, h2s=0.1, n2=0.02, h2=0.1, viscosity=True, density=True, thermo=True, Metric = False)

{'Z': 0.7941021413708897,          
 'Density': 9.383130066514621,       # lbm/cuft
 'H': -533.4898678489947,            # Btu/(lb-mol)
 'Cp': 12.880695477802188,           # Btu/(lb-mol·R)
 'Cv': 7.079475893321194,            # Btu/(lb-mol·R)
 'JT': 0.03602552422338762,          # degF/psia
 'Viscosity': 0.017905453206830867}  # cP
 
# With other flags to False, or omitted, just Z-Factor is calculated & returned
# If 100% inert mixture, then gas sg input is ignored
# If implied hydrocarbon gas sg < methane, then sg input is ignored and hydrocarbon MW set to methane.
# Any inerts with undefined mole fractions default to zero mole fraction
bns.pr_properties(temp=60, pres=2000, sg=0.75,  co2=1.0)   # 100% CO2
{'Z': 0.2778261815694524}


```

## Helium extension (October 2026)

> **⚠ Helium uses a deliberately non-standard critical temperature: Tc = 6.35 °R, not the NIST 9.35 °R (5.195 K). This is not a typo. Do not "correct" it.** The helium constants (Tc, AF, VSHIFT, VcVis) were regressed together and only work as a set.

`bns.py` now accepts `he=` (mole fraction) alongside `co2`, `h2s`, `n2` and `h2`. Helium is a sixth component with all helium BIPs set to zero; only pure-helium properties were regressed, and no previously fitted parameter was changed. With `he=0` every output is identical to the five-component version (320 reference cases across composition, 50-300 °F and 14.7-15,000 psia agree to 5e-15 relative). Only the Python implementation carries helium so far; the VBA, workbook, Rust and Fortran versions do not yet.

Why the non-standard Tc: in this model Tc also sets the Stiel-Thodos dilute-gas viscosity inside LBC. With the true Tc that term is about 25% low for helium (N2 −1.9%, CO2 −0.1%), and VcVis cannot repair it because LBC collapses to the dilute term at low density. At helium's reduced temperatures the density fit is insensitive to Tc (VSHIFT and AF absorb it), so an effective Tc was regressed with them, as the custom H2 Tc (47.43 °R against NIST 59.7 °R) was earlier.

| Helium constants | MW | Tc (°R) | Pc (psia) | AF | VSHIFT | Ω_A, Ω_B | VcVis (ft³/lb-mol) |
|---|---|---|---|---|---|---|---|
| | 4.003 | **6.35 (non-standard)** | 32.9236 | −0.17984 | −0.078082 | PR defaults | 0.76778 |

Accuracy against 29,400 NIST WebBook points for pure helium (49 isotherms 60-300 °F, 14.7-14,990 psia), mean / max absolute error:

| | Tc 6.35 °R (used) | Tc 9.34 °R (standard) |
|---|---|---|
| Density | 0.31% / 0.84% | 0.32% / 0.86% |
| Viscosity | 0.70% / 3.1% | 18.6% / 57% |
| Cp | 0.23% / 0.60% | |

Thermal outputs were not fitted (as for the other components). Against NIST the pressure dependence of enthalpy is about 13% low (H at 60 °F, 14,990 psia: 475 against 543 Btu/lb-mol relative to 60 °F and 14.7 psia), and the JT coefficient is 7-19% low in magnitude. The AF value also keeps the PR alpha term well-behaved (1 + m(1 − √Tr) ≥ 0) up to 500 °F; above that, helium results are an extrapolation.

```python
bns.pr_properties(temp=150, pres=3000, sg=0.70, he=0.05, co2=0.03, n2=0.10, viscosity=True, density=True)
```

## Evolution of work
- Original Single component PR EOS model for hydrocarbon gas in reduced temperature and pressure space (per [Linkedin post 5th April 2024](https://www.linkedin.com/pulse/z-factors-natural-gas-simple-eos-based-approach-mark-burgoyne-aazrc))
- Fist update for inerts (per [Linkedin post 27th April 2024](https://www.linkedin.com/pulse/improving-single-component-peng-robinson-z-factor-inerts-burgoyne-zfxcc)) fitted constant BIP's between inert and hydrocarbon pairs.  
- Numerous updates since until - after discussion with Markus Hayes Neilsen and Milan Stanko - we decided to publish this approach as a paper. Additional functionality was added (support for Hydrogen storage, and inclusion of thermal calculations), and methodology writeup is contained in the paper [SPE-229932-MS](https://doi.org/10.2118/229932-MS).

## Data Sources used

- **Standing & Katz Z-Factors for pure hydrocarbon gas**
  - 5,940 Z-Factor data points from Appendix A, Table A-2 of Handbook of Natural Gas Engineering, Katz et al.  
- **Pure inert critical parameters**
  - NIST data for pure CO2, H2S, N2, and H2 across temperatures of 50-300°F and pressures of 14.7-15,000 psia. Densities were converted into Z-Factors, while viscosities and thermal properties tabulated.
- **Mixture properties for fitting temperature dependent inert:inert BIP pairs**
  - Publically available experimental VLE data was used to tune BIP pairs.
- **Mixture properties for fitting hydrocarbon critical property relationships**
  - 1,062 Z-factor measurements were digitized from 90 samples detailed in Wichert's 1970 thesis, covering mixtures containing 0 - 54.5% CO2, 0 - 73.9% H2S, and 0 - 25.2% N2
  - Publicly available data from Geoscience Australia’s NOPIMS website with 61 samples and >1,500 Z-Factor data points covering 0.58 - 2.58 SG gas
- **Pure Viscosities for LBC Regression**
  - Pure NIST viscosities from each of the inerts and methane were used, along with synthetic natural gas viscosities using the Lee Gonzalez and Eakin correlation over the range of 0.6 - 2.0 SG, 14.7 - 5,000 psia and 60 - 300 degF
- **Ideal gas heat capacity**
  - Pure NIST heat capacities at zero pressure were used to (re)fit Cp_ig with Riazi's polynomial form. There was no additional tuning to match thermal results, simply applying accepted departure cubic EOS model calculations leveraging our EOS model.




## Volume-shift correction to enthalpy and Joule-Thomson (August 2026)

Earlier versions of every implementation in this repository applied the Peneloux volume
shift to Z-factor and density but computed enthalpy, Cp, Cv and the Joule-Thomson
coefficient from the untranslated EOS root. Because the shift reduces to a constant molar
volume offset `c = SUM(z_i * VSHIFT_i * b_i)`, a constant translation leaves Cp, Cv and
entropy unchanged but shifts enthalpy by `-c*p` and the JT coefficient by `+c/Cp`. The
correction is now applied consistently in the Python, Fortran, Rust and VBA implementations
(all cross-checked against each other).

Measured against reference EOS (CoolProp) over 60-300 degF and 100-10,000 psia, the
correction takes the mean JT bias from +4.2% to +0.7% for pure CO2 and from +12.0% to
+5.9% for pure methane, and reduces mean absolute enthalpy-departure error for methane
from 65 to 31 Btu/lb-mol. Z, density, Cp, Cv and viscosity are unchanged to machine
precision. If you have H or JT values computed with an earlier copy of this code, expect
them to differ accordingly; the worked examples above reflect the corrected outputs.

## Note on the Excel workbook

`Code Examples/bns_VBA_Module1.bas` is Module1 exported from
`5-Component Universal EOS.xlsm`. The workbook and the .bas are in step, and both carry the
volume-shift correction to enthalpy and the Joule-Thomson coefficient, matching the Python,
Fortran and Rust implementations.

Keep them in step by hand. Excel owns the copy inside the workbook, so edit the VBA in the VBE,
then export Module1 over the .bas; patching the workbook from outside Excel is not safe, because
Excel can keep running the cached p-code and silently ignore patched source.

## Additional Resources

Datasets used for these regressions have been uploaded for public access. To replicate these findings, you will need to download and install Aaron Zicks' [PhazeComp](https://www.zicktech.com/phazecomp.html) software, which will run these models with the free functionality. I strongly encourage those interested in deepening their understanding of EOS modelling to invest time in mastering this software.
