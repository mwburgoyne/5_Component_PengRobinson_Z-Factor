**Improving the single component Peng-Robinson Z-Factor Approach for Inerts**

Mark Burgoyne, 27-04-2024

Following feedback from Curtis Whitson and Simon Tortike, I explored the potential of extending the single-component Peng-Robinson Z-Factor method to explicitly incorporate inerts. This was driven by two main considerations: First, the accuracy of my original single-component model was inherently limited by the choice of critical pressure and temperature correlation. Second, as we increasingly encounter scenarios such as CCUS with high inert concentrations, a simplified yet accurate approach is needed to handle up to 100% inerts - beyond the range tested with approaches such as Wichert & Aziz.

**Methodology**

**Step 1:** I gathered 68,668 density data points from the NIST database for pure vapor and supercritical states of CO2, H2S, and N2, across temperatures of 50-300°F and pressures of 14.7-15,000 psia. These densities were then converted into Z-Factors.

**Step 2:** I applied regression analysis to these data points using a Peng-Robinson equation of state (EOS), modifying the Volume Shift for CO2, N2, and H2S, and adjusting the OmegaA and OmegaB parameters for CO2. CO2 was the most challenging to model accurately, yet the approach achieved better than 1% average error and maintained less than 5% error for 99% of the data points, except near the critical point.

**Step 3:** I digitized 1,052 Z-factor measurements from 88 samples detailed in Wichert’s 1970 thesis, which covered mixtures containing 0-54.5% CO2, 0-73.9% H2S, and 0-25.2% N2. I then developed a four-component Peng-Robinson model by altering the coefficients of the Sutton critical property correlation for the hydrocarbon gas and adjusting the binary interaction coefficients among all four components. This optimization was performed through a Python-driven workflow, minimizing the overall root mean square error by systematically tweaking the Sutton coefficients with Python and the binary interaction parameters with PhazeComp.

**Additional Resources**

All the datasets used for these regressions have been uploaded for public access. To replicate these findings, you will need to download and install Aaron Zicks’ PhazeComp software, which will run these models with the free functionality. I strongly encourage those interested in deepening their understanding of EOS modelling to invest time in mastering this software.
