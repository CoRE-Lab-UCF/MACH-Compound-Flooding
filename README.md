# A multivariate statistical framework for mixed populations in compound flood analysis

This directory contains MATLAB scripts used for the analysis of compound flooding events in coastal regions, focusing on the distinct contributions from tropical cyclones (TCs) and non-tropical cyclones (non-TCs). The methodology and results are detailed in the paper "[A multivariate statistical framework for mixed populations in compound flood analysis](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-1122)."

<br>

>Maduwantha, P., Wahl, T., Santamaria-Aguilar, S., Jane, R. A., Booth, J. F., Kim, H., and Villarini, G.: A multivariate statistical framework for mixed populations in compound flood analysis, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2024-1122, 2024.

<br>

## Overview of Scripts

### 1. Bias Correction
- **`Bias_Correction.m`**
  - **Description**: This script applies bias correction to rainfall data using quantile mapping. The method adjusts the hourly rainfall gauge data to match the basin-average values derived from the Analysis of Record for Calibration (AORC) data. 

### 2. Data Checking
- **`Checking_TC_and_ET_cylcone_tracks.m`**
  - **Description**: This script validates the tracks of tropical cyclones (TCs) and extratropical cyclones (ETCs) used in the study. 

### 3. Data Accumulation
- **`Creating_Accumulated_RF.m`**
  - **Description**: This script generates accumulated rainfall data over specified time intervals. It processes rainfall data to create cumulative totals, which are then used to analyze rainfall patterns during and after storm events. 

### 4. Isoline Creation
- **`Creating_Isolines_From_Two_Populations.m`**
  - **Description**: This script produces isolines based on statistical data from two distinct populations. Isolines are lines of equal value that help in comparing flood risk profiles between different data sets, such as those derived from different storm types or different time periods.

### 5. Peak Over Threshold (POT) Analysis
- **`Creating_POT_Extremes_NTR.m`**
  - **Description**: This script identifies peak-over-threshold (POT) extremes for non-tidal residuals (NTR). It is used to determine extreme NTR events that are crucial for understanding the upper tails of the NTR distribution, which significantly impact flood risk assessments.
- **`Creating_POT_Extremes_RF.m`**
  - **Description**: This script identifies POT extremes for rainfall data. It focuses on extreme rainfall events that are critical for flood risk analysis, helping to identify periods of intense rainfall that could lead to flooding.

### 6. Distribution Fitting
- **`Fitting_Distributions.m`**
  - **Description**: This script fits various statistical distributions to the observed data, which is essential for probabilistic forecasting and simulation. It evaluates different distributions to find the best fit, aiding in the accurate modeling of extreme events.

### 7. Sampling and Stratification
- **`Sampling_Con_NTR.m`**
  - **Description**: This script conducts conditional sampling on non-tidal residuals (NTR) data. It selects NTR data based on specific conditions for use in simulations, which is important for detailed flood risk analysis.
- **`Sampling_Con_RF.m`**
  - **Description**: This script conducts conditional sampling on rainfall data. It selects rainfall data under specified conditions, enhancing the analysis of flood risks related to extreme rainfall events.
- **`Stratification_Con_NTR.m`**
  - **Description**: This script stratifies non-tidal residuals (NTR) data. It organizes NTR data into different strata based on certain criteria, allowing for a more nuanced analysis under varying conditions.
- **`Stratification_Con_RF.m`**
  - **Description**: This script stratifies rainfall data. It organizes rainfall data into different strata, which is crucial for detailed analysis of rainfall patterns and their impacts on flooding.

### 8. Utility Scripts
- **`DKW_conf_int.m`**
  - **Description**: This script calculates Dvoretzky-Kiefer-Wolfowitz (DKW) confidence intervals for empirical distribution functions. It provides statistical bounds for the distribution of data, which is important for validating the results of extreme value analyses.
- **`Dummy`**
  - **Description**: This script serves as a placeholder or template for new scripts. It can be used for testing purposes or as a starting point for developing new analysis scripts.
- **`GPD_EST.m`**
  - **Description**: This script estimates the parameters of the Generalized Pareto Distribution (GPD). It is used to fit extreme value data to the GPD, which is commonly used in the analysis of POT extremes.
- **`Gamma_EST.m`**
  - **Description**: This script estimates the parameters for the Gamma distribution. It fits data to a Gamma distribution, which is often used to model the distribution of rainfall and other hydrological variables.
- **`Logistic_EST.m`**
  - **Description**: This script estimates the parameters for the Logistic distribution. It fits data to a Logistic distribution, which can be useful in modeling data with certain characteristics.
- **`ut_reconstr.m`**
  - **Description**: This script provides utilities for reconstructing datasets. It includes functions for data reconstruction, which are useful in preparing datasets for analysis.
- **`ut_solv.m`**
  - **Description**: This script provides utility functions for solving mathematical problems involved in data analysis. It includes functions for solving equations and performing calculations that are required in other scripts.

## Getting Started

1. **MATLAB Installation**: Ensure MATLAB is installed and set up on your system.
2. **Data Preparation**: Download the required datasets as specified in the associated paper.
3. **Script Execution**: Open MATLAB and run each script as needed. Modify input parameters based on your specific data and analysis requirements.


## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References

Refer to the associated paper for detailed context and methodologies: [A multivariate statistical framework for mixed populations in compound flood analysis](https://doi.org/10.5194/egusphere-2024-1122).






