# A Universal, EOS-Based Correlation for Z-Factor, Viscosity and Enthalpy For Hydrocarbon and H2, N2, CO2, H2S Gas Mixtures

**Author**: Mark Burgoyne  
**First released**: 05-04-2024  
**Most recent release**: 15-06/2025
 

Accurate estimation of gas compressibility (Z-Factor), viscosity and enthalpy are essential in applications spanning conventional hydrocarbon production, greenhouse gas (GHG) storage, and hydrogen-based energy systems. Industry-standard legacy Z-Factor correlations, such as Dranchuk-Abou-Kassem (DAK) and Hall-Yarborough (HY), are accurate for predominantly hydrocarbon gases but can struggle to converge at low temperatures or high pressures, with accuracy degrading significantly when using standard critical property correlations for mixtures containing high inert fractions. Conversely, advanced frameworks such as GERG2008 handle diverse mixtures comprehensively but face practical challenges including limited availability in software, complexity for spreadsheet-based implementation, and predefined compositional constraints.

The thermodynamic model developed consists of a volume-translated Peng-Robinson EOS for prediction of gas compressibility factor, density, enthalpy and specific heat (using departure function calculations), and the Lohrenz-Bray-Clark (LBC) method for viscosity prediction, delivering an all-in-one solution applicable for both conventional petroleum engineering as well as new storage applications. Five components were considered; one representing hydrocarbons and four inerts (H2, N2, CO2, H2S). Component properties (including temperature dependent binary interaction parameters) and model constants were regressed to reproduce the Standing-Katz chart and inert component data from the National Institute of Standards and Technology (NIST) and experimental binary VLE data. Validation confirms that the proposed model matches or surpasses the accuracy of traditional Z-factor correlations across temperatures from 50-300 F, and pressures from 14.7 - 15,000 psia, significantly outperforming them for inert mixtures, and introduces support for hydrogen-rich mixtures not previously available. Viscosity predictions similarly align closely with established correlations and reference data.

This unified method simplifies scoping level gas property calculations by using one coherent, physics-based framework, analytically solved, enhancing accuracy and usability. Its versatility enables straightforward spreadsheet implementation for most applications and easy integration into reservoir simulation software, providing a practical alternative suitable for traditional petroleum engineering workflows and emerging applications in carbon capture, utilization, and storage (CCUS) and underground hydrogen storage.

## Where are the files?
I will update the final version, complete with (a) data files used to replicate this work, (b) Working python code example, (c) Working Excel example, (d) link to Techncial paper describing the methodology once completed. In the mean time, contact me directly if you woud like an interim version: mark dot w dot burgoyne at gmail dot com

## Example usage
```python
import bns as bns

# With Metric = True, temperature is in deg C, pressure in MPa and results are also in Metric units per below comments
bns.pr_properties(temp=20, pres=1, sg=16.0425/28.97,  co2=0.2, h2s=0.1, n2=0.02, h2=0, viscosity=True, density=True, thermo=True, Metric = True, verbose=False)

{'Z': 0.9774927909395308,
 'Density': 6.733398040501448,      # kg/m3
 'H': -9.291398955385844,           # kJ/(kmol K)
 'Cp': 20.29653028742051,           # kJ/(kmol K)
 'Cv': 15.180692602056988,          # kJ/(kmol K)
 'JT': 5.1202038900231805,          # degC/MPa
 'Viscosity': 0.01109866095887059}  # cP
 
# With Metric = False (or omitted), temperature in deg F and pressure in psia
bns.pr_properties(temp=120, pres=2000, sg=16.0425/28.97,  co2=0.2, h2s=0.1, n2=0.02, h2=0, viscosity=True, density=True, thermo=True, Metric = False, verbose=False)

{'Z': 0.7659286647434862,
 'Density': 9.939570654063326,      # lbm/cuft
 'H': -701.2889499297443,           # Btu/(lb-mol·R)
 'Cp': 13.57348973097375,           # Btu/(lb-mol·R)
 'Cv': 7.053849310183024,           # Btu/(lb-mol·R)
 'JT': 0.03997877139841768,         # degF/psia
 'Viscosity': 0.01858423360291733}  # cP
 
# With other flags to False, or omitted, just Z-Factor is calculated & returned
# If 100% inert mixture, then gas sg input is ignored
# If implied hydrocarbon gas sg < methane, then sg input is ignored and hydrocarbon MW set to methane.
# Any inerts with undefined mole fractions default to zero mole fraction
bns.pr_properties(temp=60, pres=2000, sg=0.75,  co2=1.0)   # 100% CO2
{'Z': 0.2778261815694524}


```

## Evolution of work
- Original Single component PR EOS model for hydrocarbon gas in reduced temperature and pressure space (per [Linkedin post 5th April 2024](https://www.linkedin.com/pulse/z-factors-natural-gas-simple-eos-based-approach-mark-burgoyne-aazrc))
- Fist update for inerts (per [Linkedin post 27th April 2024](https://www.linkedin.com/pulse/improving-single-component-peng-robinson-z-factor-inerts-burgoyne-zfxcc)) fitted constant BIP's between inert and hydrocarbon pairs.  
- Numerous updates since until - after discussion with Markus Hayes Neilsen and Milan Stanko - we decided to publish this approach as a paper. Additional functionality was added (support for Hydrogen storage, and inclusion of thermal calculations), and methodology writeup will be much clearer in the paper (hence it's all been removed from this Github site).

## Data Sources used

- **Standing & Katz Z-Factors for pure hydrocarbon gas**
  - 5,940 Z-Factor data points from Appendix A, Table A-2 of Handbook of Natural Gas Engineering, Katz et al.  
- **Pure inert critical parameters**
  - NIST data for pure CO2, H2S, N2, and H2 across temperatures of 50-300°F and pressures of 14.7-15,000 psia. Densities were converted into Z-Factors, while viscosities and thermal properties tabulated.
- **Mixture properties for fitting temperature dependent inert:inert BIP pairs**
  - Publically available experimental VLE data was used to tune BIP pairs.
- **Mixture properties for fitting hydrocarbon critical property relationships**
  - 1,061 Z-factor measurements were digitized from 89 samples detailed in Wichert's 1970 thesis, which covered mixtures containing 0-54.5% CO2, 0-73.9% H2S, and 0-25.2% N2
  - Additional data points were also generated using the GERG2008 model for synthetic hydrocarbon gas mixtures with higher MW's to better constrain / describe fits with higher MW hydrocarbon gas.
- **Pure Viscosities for LBC Regression**
  - Pure NIST viscosities from each of the inerts and methane were used, along with synthetic natural gas viscosities using the Lee Gonzalez and Eakin correlation over the range of 0.6 - 2.0 SG, 14.7 - 5,000 psia and 60 - 300 degF
- **Ideal gas heat capacity**
  - Pure NIST heat capacities at zero pressure were used to (re)fit Cp_ig with Riazi's polynomial form. There was no additional tuning to match thermal results, simply applying accepted departure model calculations leveraging our EOS model.




## Additional Resources

All the datasets used for these regressions have been uploaded for public access. To replicate these findings, you will need to download and install Aaron Zicks' [PhazeComp](https://www.zicktech.com/phazecomp.html) software, which will run these models with the free functionality. I strongly encourage those interested in deepening their understanding of EOS modelling to invest time in mastering this software.
