# A Universal, EOS-Based Correlation for Z-Factor, Viscosity and Thermal Properties For Hydrocarbon and H<sub>2</sub>, N<sub>2</sub>, CO<sub>2</sub>, H<sub>2</sub>S Gas Mixtures

**Author**: Mark Burgoyne  
**First released**: 05-04-2024  
**Most recent update**: 03-10-2026 (helium added; see [CHANGELOG.md](CHANGELOG.md))

Example implementations in Python, Excel VBA, Rust and Fortran of the method outlined in ADIPEC 2025 paper, SPE-229932-MS, and of its later extensions.

## Original and Latest

`Code Examples/` holds two versions of every implementation:

- **`Original (SPE-229932-MS)/`**: the five-component model exactly as described in the paper (the repository on 7 November 2025). Use it to reproduce the paper.
- **`Latest (with Helium)/`**: the current model. It adds helium as a sixth component, applies the volume shift to enthalpy and the Joule-Thomson coefficient (August 2026), and corrects implementation defects that made the four languages disagree with each other and with the regression (October 2026). **Results from `Latest` therefore differ from `Original` even without helium**: enthalpy and JT in every language (by 5.8% and 3.4% in the field-units worked example below), and Z, density and viscosity for associated gas with inerts in Python, for gases with inerts in Fortran, and slightly in VBA. Python, Rust and Fortran in `Latest` agree to 1e-12.

Every change since the paper, with its effect on the numbers, is listed in [CHANGELOG.md](CHANGELOG.md).

**Paper Abstract:** Accurate estimation of gas compressibility (Z-Factor), viscosity, and enthalpy is vital across conventional hydrocarbon production, greenhouse gas storage, and hydrogen energy applications. Legacy Z-factor correlations like Dranchuk-Abou-Kassem (DAK) and Hall-Yarborough (HY) perform well for typical natural gases but lose accuracy - and sometimes fail - at low temperatures, high pressures, or with inert-rich mixtures. Advanced frameworks such as GERG2008 address a broader range of compositions but are often impractical for spreadsheets and standard engineering workflows due to software, complexity, and compositional constraints.

We present a unified thermodynamic model combining a volume-translated Peng-Robinson Equation of State (EOS) for Z-factor, density, enthalpy, Cp, and Joule-Thomson coefficient with the Lohrenz-Bray-Clark (LBC) method for viscosity. The model covers one hydrocarbon and four inert components (H₂, N₂, CO₂, H₂S), with parameters tuned to Standing-Katz, NIST, binary vapor liquid equilibrium (VLE) and publicly available experimental Z-Factor data. Validation shows accuracy equal to or better than traditional correlations from 50-300 °F and 14.7-15,000 psia, including gas condensates and associated hydrocarbon gases as well as inert and hydrogen-rich mixtures. Viscosity results closely match reference data.

This approach streamlines gas property calculations with a single, analytically solved, spreadsheet-ready framework - facilitating integration into reservoir simulators and meeting the needs of both conventional and emerging carbon dioxide and hydrogen storage workflows.

Additionally, in its simplest form, this analytically solvable cubic EOS-based method can be a compact, direct and robust replacement for iterative legacy methods (e.g., DAK, HY), eliminating numerical convergence issues and simplifying integration into standard engineering workflows.

Calculations, regression data, and example implementations in Python and Excel VBA are available here.


## Example usage

The examples use `Code Examples/Latest (with Helium)/bns.py`.

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

## Helium (October 2026)

> **⚠ Helium uses a deliberately non-standard critical temperature: Tc = 6.35 °R, not the NIST 9.35 °R (5.195 K). This is not a typo. Do not "correct" it.** In this model Tc also sets the Stiel-Thodos dilute-gas viscosity inside LBC, which is about 25% low for helium at the true Tc; the density fit is insensitive to Tc. The helium Tc, acentric factor, volume shift and VcVis were regressed together and only work as a set.

Every implementation in `Latest` takes a helium mole fraction as a new last argument (`he`), with all helium BIPs zero. Against 29,400 NIST points for pure helium (60-300 °F, 14.7-14,990 psia) the model reproduces density to 0.31% mean / 0.84% max and viscosity to 0.70% / 3.1%. Constants, basis, accuracy and limitations are in [CHANGELOG.md](CHANGELOG.md).

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
- **Helium (Latest only)**
  - NIST WebBook densities and viscosities for pure helium, 49 isotherms at 60-300 °F and 14.7-14,990 psia (29,400 points); PhazeComp decks and inputs in `Data/06 - Helium EOS and VcVis Regression.zip`.




## Volume-shift correction to enthalpy and Joule-Thomson (August 2026)

The paper's implementations applied the volume shift to Z-factor and density but not to enthalpy and the Joule-Thomson coefficient. `Latest` corrects this: enthalpy shifts by `−c·(p − 14.696)` and JT by `+c/Cp`, with `c = SUM(z_i * VSHIFT_i * b_i)`, while Z, density, Cp, Cv and viscosity are unchanged. If you have H or JT values from the `Original` code or any earlier copy, expect them to differ. Details and measured effect: [CHANGELOG.md](CHANGELOG.md).

## Note on the Excel workbook

`Original (SPE-229932-MS)/` contains the workbook `5-Component Universal EOS.xlsm` as published, with its Module1 also exported as `bns_VBA_Module1.bas` for reading and diffing.

`Latest (with Helium)/` contains the updated module, `bns_VBA_Module1.bas`, but no workbook: VBA cannot be rewritten safely from outside Excel, because Excel can keep running the cached p-code and silently ignore patched source. To build the Latest workbook, open a copy of the Original workbook, press Alt+F11, right-click Module1 and choose Remove Module (no export needed), then File > Import File and select `Latest (with Helium)/bns_VBA_Module1.bas`, and save. The helium argument is optional, so existing worksheet formulas keep working.

## Additional Resources

Datasets used for these regressions have been uploaded for public access. To replicate these findings, you will need to download and install Aaron Zicks' [PhazeComp](https://www.zicktech.com/phazecomp.html) software, which will run these models with the free functionality. I strongly encourage those interested in deepening their understanding of EOS modelling to invest time in mastering this software.
