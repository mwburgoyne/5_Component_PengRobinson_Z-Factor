# Improving the Single Component Peng-Robinson Z-Factor Approach for Inerts

**Author**: Mark Burgoyne  
**First released**: 05-04-2024  
**First Update**: 27-04-2024  
**Second Update**: 05-05-2024  
**Third Update**: 19-05-2024 -> 05-06-2024  
**Fourth Update**: 15-06-2024

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
- Fourth update of this work (15th June, with current GitHub content):
  1. Retuned inerts and hydrocarbon single component critical parameters to minimize deviance from standard critical properties
  2. Added support for Hydrogen (critical paremeters and BIP pairs)

## Data Sources used

- **Standing & Katz Z-Factors for pure hydrocarbon gas**
  - Digitized 515 Z-Factor data points from original Standing & Katz plot found in SPE-942140 - "Density of Natural Gases"  
- **Pure inert critical parameters**
  - Generated GERG2008 EOS data for pure CO2, H2S, N2, and H2 across temperatures of 50-300°F and pressures of 14.7-15,000 psia. Densities were converted into Z-Factors, and along with viscosities tabulated.
- **Mixture properties for fitting temperature dependent inert:inert BIP pairs**
  - A GERG2008 multicomponent model was used to create 4,990 mixture Z-Factor data points over the temperature range 90 – 300 degF and 14.7 – 15,000 psia, and various CO2, H2S, N2, H2 mole fractions
- **Mixture properties for fitting MW and temperature dependent HC:inert BIP pairs**
  - 1,061 Z-factor measurements were digitized from 89 samples detailed in Wichert’s 1970 thesis, which covered mixtures containing 0-54.5% CO2, 0-73.9% H2S, and 0-25.2% N2
  - An additional ~45,000 data points were also generated using the GERG2008 model for synthetic hydrocarbon gas mixtures with MW's of 16, 26, 34 and 45 lb/lb-mol, paired with 0.25 - 0.75 mole fraction inerts in order to better constrain / describe fits to (a) higher hydrocarbon gas MW's and (b) higher inert fractions.
- **Pure Viscosities for LBC Regression**
  - Pure NIST viscosities from each of the inerts were used, along with 24,000 synthetic natural gas viscosities using the Lee Gonzalez and Eakin correlation over the range of 0.5 - 2.0 SG, 14.7 – 15,000 psia and 60 – 300 degF

## Steps

1.	From digitized Standing & Katz plot, using PR & Peneloupe EOS, fit a single component ‘Gas’ model to reduced pressure and temperature data by regressing on ACF, VTRANS OmegA and OmegaB
2.	With data generated from GERG2008 EOS, regress on critical parameters for pure CO2, H2S, N2 and H2. 
  a.	Changing ACF, VTRAN, OmegaA and OmegaB for CO2
  b.	VTRAN, OmegaA and OmegaB for H2S
  c.	VTRAN only for N2 and H2
3.	Create data with GERG2008 of synthetic 25%, 50% and 75% binary mixtures of inerts over a range of pressures and temperatures, fitting BIP’s as a function of temperature for each pair.
4.	Using (a) real data published by Wichert, augmented with (b) synthetic GERG2008 data to anchor scenarios with elevated hydrocarbon MW, regress on (i) hydrocarbon gas : inert BIP’s as a function of temperature and hydrocarbon MW and (ii) Tc and Pc behaviour as a function of hydrocarbon gas MW to match Z-Factors
5.	Tune VCVIS and LBC coefficients to match single component inert viscosity (from NIST) and Lee Gonzalez & Eakin viscosity for pure hydrocarbon gas.


## Results

**Critical parameters**

| Comp | MW    | Tc (R)  | Pc (psia) | ACF     | VTRAN   | OmegaA  | OmegaB   | VcVis (ft³/lbmol)|
|------|-------|---------|-----------|---------|---------|---------|----------|------------------|
| CO2  | 44.01 | 547.416 | 1069.51   | 0.10560 | -0.19240| 0.440853| 0.0730166| 1.55342          |
| H2S  | 34.082| 672.12  | 1299.97   | 0.01839 | -0.18906| 0.441796| 0.0739200| 1.57142          |
| N2   | 28.014| 227.16  | 492.84    | 0.03700 | -0.21067| 0.457236| 0.0777961| 1.41785          |
| H2   |  2.016|  59.364 | 187.53    | -0.2170 | -0.12899| 0.457236| 0.0777961| 1.18168          |
| Gas  | *     | *       | *         | -0.04051| -0.19185| 0.457236| 0.0777961| *                |

**Properties are MW dependent*
mw_hc = Inert free hydrocarbon gas MW  

| Variable   | Tc             | Pc              |
|------------|----------------|-----------------|
| A          | 4.76601        |   -1.18331      |
| B          | 45.4709        | 1104.13         |
| C          |  0.00244823    | 0.00702049875   |
| D          |  0.3219414605  |    1.50107866   |

`Gas Tc (R) = (A * mw_hc + B)/(C * mw_hc + D)`  
`Gas Pc (psia) = (A * mw_hc + B)/(C * mw_hc + D)`  
`Gas VcVis (ft³/lbmol) = 0.070294775 * mw_hc + 0.333582429`  
`LBC P3, P4 = -4.01411e-02, 9.00617e-03`  


| BIP Pair Parameters: A | H2S                             | N2                     | H2                   |Hydrocarbon Gas      |
|------------------------|---------------------------------|------------------------|----------------------|---------------------|
| CO2                    | -0.000811630909090909           | -0.00112389            | -3.949E-05           | 0.356872            |
| H2S                    |                                 | 0.00113879071428571    | 0.00377893           | 0.13333             |
| N2                     |                                 |                        | 0.000227832857142857 | -0.0429868          |
| H2                     |                                 |                        |                      | 0.712742            |

| BIP Pair Parameters: B | H2S                             | N2                     | H2                   |Hydrocarbon Gas      |
|------------------------|---------------------------------|------------------------|----------------------|---------------------|
| CO2                    | 0.244400581818182               | -0.1389305             | 0.661337             | -179.465            |
| H2S                    |                                 | 0.251688807142857      | 0.5633018            | -24.7713            |
| N2                     |                                 |                        | 0.199822661904762    | -3.39882            |
| H2                     |                                 |                        |                      | -349.25             |

| BIP Pair Parameters: C | H2S                             | N2                     | H2                   |Hydrocarbon Gas      |
|------------------------|---------------------------------|------------------------|----------------------|---------------------|
| CO2                    | 0.000199935714285714            | 0.000648015753811406   | 0.00116141530208922  | -0.00128858         |
| H2S                    |                                 | 0.00193186545454545    | 0.00119598402597403  | -0.00141921         |
| N2                     |                                 |                        | 0.000141143636363636 | 0.000728009         |
| H2                     |                                 |                        |                      | 0.00110919          |

| BIP Pair Parameters: D | H2S                             | N2                     | H2                   | 
|------------------------|---------------------------------|------------------------|----------------------| 
| CO2                    | 0.0843990523809524              | -0.294767026538679     | 0.530504334274421    | 
| H2S                    |                                 | 0.0794874545454544     | 0.831735004329004    | 
| N2                     |                                 |                        | 0.216575909090909    | 

`H2S:CO2, CO2:N2, H2S:N2, CO2:H2 BIPs = max((A * degF + B),(C * degF + D))`  ,
`H2S:H2, N2:H2 BIPs = min((A * degF + B),(C * degF + D))`  ,
`Hydrocarbon:Inert BIPs = A + B / DegR + C * mw_hc`  

# Pure Hydrocarbon Gas Residual Error Plots

| ![Calculated vs Standing & Katz Z-Factors](images/z-factors_S&K.png) |
|----------------------------------------------------------------------|
| ![Relative Error Map](images/Relative_Error_Map.png)                 |
[Uploading Inerts_PR_Match_VSHIFT_Only_AB_CO2_BigLoop_BIPS_Temp_reduced_MW_Slope4_Linear_v3.phz…]()



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
| Peng Robinson                 | 0.002          | 0.074          | 0.0321               |
| DAK + Sutton & Wichert        | 0.007          | 0.097          | 0.0336               |
| DAK + PMC                     | 0.004          | 0.157          | 0.0349               |

| Method (vs Wichert & GERG data) | Avg Rel. Error | Max Rel. Error | 95% of rel. errors < |
|---------------------------------|----------------|----------------|----------------------|
| Peng Robinson                   | 0.001          | 0.074          | 0.0277               |
| DAK + Sutton & Wichert          | 0.003          | 0.359          | 0.0618               |
| DAK + PMC                       | 0.024          | 0.651          | 0.122                |


| ![Relative Z-Factor Error](images/rel_wichert.png) | ![Wichert Data Correlation Matrix for inputs and residual relative error](images/corel_wichert.png) |


## Additional Resources

All the datasets used for these regressions have been uploaded for public access. To replicate these findings, you will need to download and install Aaron Zicks’ PhazeComp software, which will run these models with the free functionality. I strongly encourage those interested in deepening their understanding of EOS modelling to invest time in mastering this software.
