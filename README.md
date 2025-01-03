# Improving the Single Component Peng-Robinson Z-Factor Approach for Inerts

**Author**: Mark Burgoyne  
**First released**: 05-04-2024  
**First Update**: 27-04-2024  
**Second Update**: 05-05-2024  
**Third Update**: 19-05-2024 -> 05-06-2024  
**Fourth Update**: 15-06-2024  
**Fifth Update**: 14-07-2024  
**Sixth Update**: 29-09-2024  
**Seventh Update**: 28-12-2024
 

Following feedback from Curtis Whitson and Simon Tortike, I explored the potential of extending the single-component Peng-Robinson Z-Factor method to explicitly incorporate inerts. This was driven by two main considerations: 
1. The accuracy of my original single-component model was inherently limited by the choice of critical pressure and temperature correlation.
2. As we increasingly encounter scenarios such as CCUS with high inert concentrations, a simplified yet accurate approach is needed to handle up to 100% inerts - beyond the range tested with approaches such as Wichert & Aziz.

Additional discussions with Markus Hayes Neilsen prompted me to extend the approach to incorporate Hydrogen as well

A ‘simple’ tuned Peng Robinson model scheme, complete with binary interaction coefficients as a function of temperature and hydrocarbon gas MW, and hydrocarbon Tc and Pc as a function of hydrocarbon gas MW.
By assuming single phase gas, we do not need to do flash calculations, solving only for maximum real root of the cubic equation.
Instead of using conventional correlations such as Sutton & Wichert Aziz or PMC that estimate Tc and Pc from MW of the gas mixture and applying corrections to Tc & Pc correlations for inerts, this approach treats hydrocarbon gas, CO2, H2S, N2 and H2 as distinct species, permitting more reliable application at inert fractions up to and including 1.0.

## Key differentiators include
  1.	Non iterative, intrinsically stable solution scheme well suited to analytical applications
  2.	Broad applicability to up to 100% Nitrogen, CO2, H2S or H2, making it well suited to not only standard oilfield gas, but also scoping storage injection studies.
  3.	Based upon industry standard cubic EOS, enabling greater confidence to use over wide P, T range.
  4.	Results can be replicated with any industry standard applications that support PR EOS with volume translation.
  5.	Complete with a paired tuned LBC viscosity model that is similarly applicable over wide ranges of P, T and composition.

## Evolution of work
- Original Single component PR EOS model for hydrocarbon gas in reduced temperature and pressure space (per [Linkedin post 5th April 2024](https://www.linkedin.com/pulse/z-factors-natural-gas-simple-eos-based-approach-mark-burgoyne-aazrc))
- Fist update for inerts (per [Linkedin post 27th April 2024](https://www.linkedin.com/pulse/improving-single-component-peng-robinson-z-factor-inerts-burgoyne-zfxcc)) fitted constant BIP’s between inert and hydrocarbon pairs.  
- Second update of this work (5th May 2024), delivered additional functionality/accuracy including:
  1. All BIP pairs temperature dependent
  2. Hydrocarbon MW dependent adjustments for Hydrocarbon : Inert BIP pairs
  3. Implemented LBC viscosity calculations, tuned for up to 100% mole fraction of natural gas or inerts
- Third update of this work (19th May - 5th June 2024, with current GitHub content), delivered additional functionality/accuracy including:
  1. Reverted to a single component gas EOS fit that avoided changing Omega values
  2. Refitted all inert critical properties to GERG2008, including efforts to minimize OmegaA/B deviations from default for CO2
  3. Refitted inert temperature dependant pairs to GERG2008 data (previously used noisy Wichert data)
  4. Updated all knock-on coefficients and dependancies. Augmented data to fit against with four synthetic gas samples at various richnesses and inert fractions using the GERG2008 EOS.
- Fourth update of this work (15th June):
  1. Retuned inerts and hydrocarbon single component critical parameters to minimize deviance from standard critical properties
  2. Added support for Hydrogen (critical paremeters and BIP pairs)
- Fifth update of this work (23rd June 2024):
  1. Re-digitized Standing & Katz data to address slight errors, and re-regressed the 'Gas' component with default Omegas
  2. Simplified the Inert:Inert BIP pairs to constants.
  3. Moved to H2 critical properties outlined in SPE214437
  4. Moved to alternate Tc and Pc form that extrapolates more inline with Twu at higher MW's
- Sixth update of this work (2nd September 2024):
  1. Refit 'Gas' model to Standing & Katz data from Table A2 of Natural Gas Engineering Handbook from Katz et al. Delivers very minor change, however now data tied to published information and with higher granularity.
  2. Pure inert component fits left unchanged. Refit subsequent BIP and Viscosity parameters
  3. Fit Tc as a function of Vc/Zc rather than Gas MW
- Seventh update of this work (28th December 2024, with current GitHub content):
  1. Refit viscosities: (a) Introduce pure methane viscosity from NIST, (b) Limit maximum pressure to 10,000 psia from NIST, (c) Limit Lee Gonzalez & Eakin to 5,000 psia for pure hydrocarbons heavier than methane due to accuracy concerns above this

## Data Sources used

- **Standing & Katz Z-Factors for pure hydrocarbon gas**
  - 5,940 Z-Factor data points from Appendix A, Table A-2 of Handbook of Natural Gas Engineering, Katz et al.  
- **Pure inert critical parameters**
  - Generated GERG2008 EOS data for pure CO2, H2S, N2, and H2 across temperatures of 50-300°F and pressures of 14.7-15,000 psia. Densities were converted into Z-Factors, and along with viscosities tabulated.
- **Mixture properties for fitting temperature dependent inert:inert BIP pairs**
  - A GERG2008 multicomponent model was used to create 4,990 mixture Z-Factor data points over the temperature range 90 – 300 degF and 14.7 – 15,000 psia, and various CO2, H2S, N2, H2 mole fractions
- **Mixture properties for fitting MW and temperature dependent HC:inert BIP pairs**
  - 1,061 Z-factor measurements were digitized from 89 samples detailed in Wichert’s 1970 thesis, which covered mixtures containing 0-54.5% CO2, 0-73.9% H2S, and 0-25.2% N2
  - An additional ~45,000 data points were also generated using the GERG2008 model for synthetic hydrocarbon gas mixtures with MW's of 16, 26, 34 and 45 lb/lb-mol, paired with 0.25 - 0.75 mole fraction inerts in order to better constrain / describe fits to (a) higher hydrocarbon gas MW's and (b) higher inert fractions.
- **Pure Viscosities for LBC Regression**
  - Pure NIST viscosities from each of the inerts and methane were used, along with synthetic natural gas viscosities using the Lee Gonzalez and Eakin correlation over the range of 0.6 - 2.0 SG, 14.7 – 5,000 psia and 60 – 300 degF

## Steps

1.	From Standing & Katz data, using PR & Peneloupe EOS, fit a single component ‘Gas’ model to reduced pressure and temperature data by regressing on ACF and VTRANS. OmegA and OmegaB
2.	With data generated from GERG2008 EOS, regress on critical parameters for pure CO2, H2S, N2 and H2. 
  a.	Changing ACF, VTRAN, OmegaA and OmegaB for CO2
  b.	VTRAN, OmegaA and OmegaB for H2S
  c.	VTRAN only for N2 and H2
3.	Using (a) real data published by Wichert, augmented with (b) synthetic GERG2008 data to anchor scenarios with elevated hydrocarbon MW, regress on (i) hydrocarbon gas : inert BIP’s as a function of temperature and hydrocarbon MW, and inert:inert BIP's and (ii) Tc and Pc behaviour as a function of hydrocarbon gas MW to match Z-Factors
4.	Tune VCVIS and LBC coefficients to match single component inert (& methane) viscosity (from NIST to 10,000 psia) and Lee Gonzalez & Eakin viscosity (to 5,000 psia) for pure hydrocarbon gas.


## Results

**Critical parameters**  

| Comp | MW    | Tc (R)  | Pc (psia) | ACF     | VTRAN   | OmegaA  | OmegaB   | VcVis (ft³/lbmol)|
|------|-------|---------|-----------|---------|---------|---------|----------|------------------|
| CO2  | 44.01 | 547.416 | 1069.51   | 0.12256 | -0.27593| 0.427705| 0.0696460| 1.46020          |
| H2S  | 34.082| 672.12  | 1299.97   | 0.04916 | -0.22896| 0.436743| 0.0724373| 1.46460          |
| N2   | 28.014| 227.16  | 492.84    | 0.03700 | -0.21067| 0.457236| 0.0777961| 1.35422          |
| H2   |  2.016|  47.430 | 187.53    | -0.2170 | -0.32400| 0.457236| 0.0777961| 0.67967          |
| Gas  | *     | *       | *         | -0.03899| -0.19076| 0.457236| 0.0777961| *                |

**Properties are MW dependent*  

| Variable   | Tc              | Vc/Zc        |
|------------|-----------------|--------------|
| A          | 3007.108547898  | 2.2997507E-4 |
| B          | 42.6061669      | 0.1864159    |
| C          | -0.032901049    | 2.470903632  |
    
mw_hc = Inert free hydrocarbon gas MW

`Hydrocarbon Gas Vc/Zc (ft³/lbmol) = A * mw_hc**2 + B * mw_hc + C`  
`Hydrocarbon Gas Tc (R) = A * ((Vc/Zc) + C) /(B + (Vc/Zc) + C)`  
`Hydrocarbon Gas Pc (psia) = R * Tc / (Vc/Zc),  where R = 10.7316`  
`Hydrocarbon Gas VcVis (ft³/lbmol) = 0.057511062 *  mw_hc + 0.478400158`  
`LBC P3, P4 = -3.92835e-02,  9.28591e-03`  


| Gas:Inert BIP Parameters | A             |  B           |   C        |
|--------------------------|---------------|--------------|------------|
| CO2                      |  0.386557     |-0.00219806	  |  -158.333  |
| H2S                      |  0.267007     |-0.00396541	  |  -58.611   |
| N2                       |  0.486589     |-0.00316789	  |  -226.239  |
| H2                       |  0.776917     | 0.0106061	   |  -474.283  |

`Hydrocarbon:Inert BIPs = A + B * mw_hc + C / DegR`  
    
| BIP Pair Values  |  CO2     |  H2S    |  N2    | 
|------------------|----------|---------|--------|
| H2S              | 0.0600319|         |        | 
| N2               |-0.229807 |-0.18346 |        | 
| H2               | 0.646796 | 0.65000 | 0.36917| 

# Pure Hydrocarbon Gas Residual Error Plots

| ![Relative Error Map](images/S&K_Rel_Error_vs_PR.png)                 |


# Pure CO2 Residual Error Plots

| ![Cross Plot Calculated vs GERG2008 Z-Factors](images/CO2_1.png) | ![Relative Z-Factor Error](images/CO2_2.png)     |
|------------------------------------------------------------------|--------------------------------------------------|
| ![Relative Error Map](images/CO2_3.png)                          | ![Relative Molar Volume Error](images/CO2_4.png) |

# Pure H2S Residual Error Plots

| ![Cross Plot Calculated vs GERG2008 Z-Factors](images/H2S_1.png) | ![Relative Z-Factor Error](images/H2S_2.png)     |
|------------------------------------------------------------------|--------------------------------------------------|
| ![Relative Error Map](images/H2S_3.png)                          | ![Relative Molar Volume Error](images/H2S_4.png) |

# Pure N2 Residual Error Plots

| ![Cross Plot Calculated vs GERG2008 Z-Factors](images/N2_1.png) | ![Relative Z-Factor Error](images/N2_2.png)     |
|-----------------------------------------------------------------|-------------------------------------------------|
| ![Relative Error Map](images/N2_3.png)                          | ![Relative Molar Volume Error](images/N2_4.png) |


# Wichert Z-Factor Cross Plots

| ![Cross Plot Calculated vs Reported Z-Factors & comparison with DAK](images/Crossplot_Z-PR_DAK.png) | ![Cross Plot Calculated vs Wichert & GERG2008 Z-Factors to higher richness](images/Crossplot_Z-PR_GERG&Wichert.png)    |
|-----------------------------------------------------------------------------------------------------|--------------------------------------------------|
| ![Cross Plot Calculated vs Wichert & Extended GERG Z-Factors & comparison with DAK](images/Crossplot_Z-PR_DAK_extended.png)                         |   |



| Method (vs Wichert data only) | Avg Rel. Error | Max Rel. Error | 95% of rel. errors < |
|-------------------------------|----------------|----------------|----------------------|
| Peng Robinson                 | -0.003         | 0.052          | 0.0243               |
| DAK + Sutton & Wichert        | 0.007          | 0.097          | 0.0336               |
| DAK + PMC                     | 0.004          | 0.157          | 0.0349               |

| Method (vs Wichert & GERG data) | Avg Rel. Error | Max Rel. Error | 95% of rel. errors < |
|---------------------------------|----------------|----------------|----------------------|
| Peng Robinson                   | -0.001         | 0.082          | 0.0277               |
| DAK + Sutton & Wichert          | 0.006          | 0.431          | 0.0999               |
| DAK + PMC                       | 0.032          | 0.654          | 0.1427               |


| ![Relative Z-Factor Error](images/rel_wichert.png) | ![Wichert Data Correlation Matrix for inputs and residual relative error](images/corel_wichert.png) |


## Additional Resources

All the datasets used for these regressions have been uploaded for public access. To replicate these findings, you will need to download and install Aaron Zicks’ [PhazeComp](https://www.zicktech.com/phazecomp.html) software, which will run these models with the free functionality. I strongly encourage those interested in deepening their understanding of EOS modelling to invest time in mastering this software.
