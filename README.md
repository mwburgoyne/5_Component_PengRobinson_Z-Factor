# Improving the Single Component Peng-Robinson Z-Factor Approach for Inerts

**Author**: Mark Burgoyne  
**Version 1**: 27-04-2024  
**Version 2**: 5/5/2024

Following feedback from Curtis Whitson and Simon Tortike, I explored the potential of extending the single-component Peng-Robinson Z-Factor method to explicitly incorporate inerts. This was driven by two main considerations: 
1. The accuracy of my original single-component model was inherently limited by the choice of critical pressure and temperature correlation.
2. As we increasingly encounter scenarios such as CCUS with high inert concentrations, a simplified yet accurate approach is needed to handle up to 100% inerts - beyond the range tested with approaches such as Wichert & Aziz.

Original Single component PR EOS model for hydrocarbon gas in reduced temperature and pressure space (per [LinkedIn post](https://www.linkedin.com/pulse/z-factors-natural-gas-simple-eos-based-approach-mark-burgoyne-aazrc) 5th April 2024)
Version 1 of this work (per [LinkedIn post](https://www.linkedin.com/pulse/improving-single-component-peng-robinson-z-factor-inerts-burgoyne-zfxcc) 27th April 2024) fitted constant BIP’s between inert and hydrocarbon pairs.  
Version 2 of this work (5th May 2024), delivered additional functionality/accuracy including:
1. All BIP pairs temperature dependent
2. Hydrocarbon MW dependent adjustments for Hydrocarbon : Inert BIP pairs
3. Implemented LBC viscosity calculations, tuned for up to 100% mole fraction of natural gas or inerts

## Data Sources used

- **Standing & Katz Z-Factors for pure hydrocarbon gas**
  - Digitized 515 Z-Factor data points from original Standing & Katz plot found in SPE-942140 - "Density of Natural Gases"  
- **Pure inert critical parameters and viscosity**
  - Gathered 68,668 density data points from the NIST database for pure vapor and supercritical states of CO2, H2S, and N2, across temperatures of 50-300°F and pressures of 14.7-15,000 psia. Densities were converted into Z-Factors, and along with viscosities tabulated.
- **Mixture properties for fitting temperature dependent inert:inert BIP pairs**
  - A GERG2008 multicomponent model was used to create 4,990 mixture Z-Factor data points over the temperature range 90 – 300 degF and 14.7 – 15,000 psia, and various CO2, H2S, N2 mole fractions
- **Mixture properties for fitting MW and temperature dependent HC:inert BIP pairs**
  - 1,061 Z-factor measurements were digitized from 89 samples detailed in Wichert’s 1970 thesis, which covered mixtures containing 0-54.5% CO2, 0-73.9% H2S, and 0-25.2% N2
- **Pure Viscosities for LBC Regression**
  - Paired NIST viscosities from each of the inert Z-Factor data points were used, along with 24,000 synthetic natural gas viscosities using the Lee Gonzalez and Eakin correlation over the range of 0.5 - 2.0 SG, 14.7 – 15,000 psia and 60 – 300 degF

## Steps

1. With Peng-Robinson equation of state (EOS), regress single component hydrocarbon gas critical properties to digitized Standing & Katz Z-Factor data in reduced temperature and pressure space.
2. Regress Volume Shift for CO2, N2, and H2S, and adjusting the OmegaA and OmegaB parameters for CO2 to match NIST inerts Z-Factors. CO2 was the most challenging to model accurately, yet the approach achieved better than 1% average error and maintained less than 5% error for 99% of the data points, except near the critical point.
3. Regress Inert:Inert temperature dependent BIP pairs to the GERG2008 mixture data, using a form as follows: `BIP = A + B / DegR`.
4. Simultaneously regress (a) HC:Inert temperature and MW dependent BIP pairs and (b) Sutton critical property coefficients to the Wichert Z-Factor data, using a form as follows: `BIP = A + B / DegR + C * HC_MW`.
5. Regress VCVIS of inerts and hydrocarbon gas at differing hydrocarbon MW’s to both the NIST data (inerts) as well as synthetic LGE viscosity data for pure hydrocarbon, while also investigating custom 3rd and 4th LBC coefficients.

## Results

**Critical parameters**

| Comp | MW    | Tc (R)  | Pc (psia) | ACF    | VTRAN   | OmegaA  | OmegaB  | VcVis (ft³/lbmol) |
|------|-------|---------|-----------|--------|---------|---------|---------|------------------|
| CO2  | 44.01 | 547.416 | 1069.51   | 0.225  | -0.27025| 0.441273| 0.070128| 1.3548           |
| H2S  | 34.082| 672.12  | 1299.97   | 0.09   | -0.11575| 0.457236| 0.077796| 1.37424          |
| N2   | 28.014| 227.16  | 492.84    | 0.037  | -0.20976| 0.457236| 0.077796| 1.25466          |
| Gas  | *     | *       | *         | -0.04896| -0.3949 | 0.429188| 0.069255| *                |

**Properties are MW dependent*

mw_hc = Inert free hydrocarbon gas MW  
`Gas Tc (R) = coefic_tc[0] * mw_hc ** 2 + coefic_tc[1] * mw_hc + coefic_tc[2]`  
`Gas Pc (psia) = coefic_pc[0] * mw_hc ** 2 + coefic_pc[1] * mw_hc + coefic_pc[2]`  
`Gas VcVis (ft³/lbmol) = 0.053212 * mw_hc + 0.487479`  
`coefic_tc = [-9.91619530e-02, 1.18144133e+01, 1.78887597e+02]`  
`coefic_pc = [-3.66048577e-03, -2.73719309e+00, 7.01192731e+02]`  


| BIP Pair Parameters: A | H2S          | N2            | Gas          |
|------------------------|--------------|---------------|--------------|
| CO2                    | -4.427770E-01| -2.305710E-01 | 9.215660E-01 |
| H2S                    |              | 2.657560E+00  | 8.314670E-03 |
| N2                     |              |               | -5.405930E+00|

| BIP Pair Parameters: B | H2S          | N2            | Gas          |
|------------------------|--------------|---------------|--------------|
| CO2                    | 3.741680E+02 | 1.330530E+02  | 1.892310E+01 |
| H2S                    |              | -1.601900E+03 | 2.093930E+01 |
| N2                     |              |               | 2.136390E+03 |

| BIP Pair Parameters: C | H2S          | N2            | Gas          |
|------------------------|--------------|---------------|--------------|
| CO2                    | 0            | 0             | -5.247030E-02|
| H2S                    |              | 0             | 5.168520E-03 |
| N2                     |              |               | 7.717860E-02 |

`BIP = A + B / DegR + C * mw_hc`

## Additional Resources

All the datasets used for these regressions have been uploaded for public access. To replicate these findings, you will need to download and install Aaron Zicks’ PhazeComp software, which will run these models with the free functionality. I strongly encourage those interested in deepening their understanding of EOS modelling to invest time in mastering this software.
