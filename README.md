# A multivariate statistical framework for mixed populations in compound flood analysis

This directory contains MATLAB scripts used for the analysis of compound flooding events in coastal regions, focusing on the distinct contributions from tropical cyclones (TCs) and non-tropical cyclones (non-TCs). The methodology and results are detailed in the paper "[A multivariate statistical framework for mixed populations in compound flood analysis](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-1122)."


## Overview of Scripts

# MACH-Compound Flooding - Scripts Directory

This directory contains MATLAB scripts used in the analysis of compound flooding events as detailed in the paper "A multivariate statistical framework for mixed populations in compound flood analysis" by Pravin Maduwantha et al. These scripts are essential for processing, analyzing, and simulating data related to compound flooding events driven by multiple factors like storm surges, rainfall, and different storm types.

## Overview of Scripts

### 1. Bias Correction
- **`Bias_Correction.m`**
  - **Purpose**: Applies bias correction to rainfall data.
  - **Description**: Uses quantile mapping to adjust hourly rainfall gauge data to match basin-average values derived from AORC data. Ensures that local measurements are representative of the entire catchment area.

### 2. Data Checking
- **`Checking_TC_and_ET_cylcone_tracks.m`**
  - **Purpose**: Validates tropical cyclone (TC) and extratropical cyclone (ETC) tracks.
  - **Description**: Ensures that the storm tracks used in the analysis are accurate and reliable.

### 3. Data Accumulation
- **`Creating_Accumulated_RF.m`**
  - **Purpose**: Generates accumulated rainfall data.
  - **Description**: Accumulates rainfall data over different time intervals to analyze rainfall patterns during and after storm events.

### 4. Isoline Creation
- **`Creating_Isolines_From_Two_Populations.m`**
  - **Purpose**: Generates isolines from two data populations.
  - **Description**: Compares flood risk profiles by creating isolines from two distinct data sets.

### 5. Peak Over Threshold (POT) Analysis
- **`Creating_POT_Extremes_NTR.m`**
  - **Purpose**: Identifies POT extremes for non-tidal residuals (NTR).
  - **Description**: Identifies extreme NTR events crucial for understanding the upper tails of the distributions.
- **`Creating_POT_Extremes_RF.m`**
  - **Purpose**: Identifies POT extremes for rainfall.
  - **Description**: Focuses on extreme rainfall events for flood risk assessment.

### 6. Distribution Fitting
- **`Fitting_Distributions.m`**
  - **Purpose**: Fits statistical distributions to data.
  - **Description**: Evaluates various statistical distributions to find the best fit for observed data, aiding in probabilistic forecasting.

### 7. Sampling and Stratification
- **`Sampling_Con_NTR.m`**
  - **Purpose**: Conducts conditional sampling on NTR data.
  - **Description**: Samples NTR data based on specific conditions for use in simulations.
- **`Sampling_Con_RF.m`**
  - **Purpose**: Conducts conditional sampling on rainfall data.
  - **Description**: Samples rainfall data under specified conditions to enhance flood risk analysis.
- **`Stratification_Con_NTR.m`**
  - **Purpose**: Stratifies NTR data.
  - **Description**: Organizes NTR data for detailed analysis under varying conditions.
- **`Stratification_Con_RF.m`**
  - **Purpose**: Stratifies rainfall data.
  - **Description**: Organizes rainfall data for nuanced analysis.

### 8. Utility Scripts
- **`DKW_conf_int.m`**
  - **Purpose**: Calculates Dvoretzky-Kiefer-Wolfowitz confidence intervals.
  - **Description**: Computes confidence intervals for empirical distribution functions.
- **`Dummy`**
  - **Purpose**: Placeholder script.
  - **Description**: Used for testing or as a template for new scripts.
- **`GPD_EST.m`**
  - **Purpose**: Estimates Generalized Pareto Distribution parameters.
  - **Description**: Fits data to a GPD for extreme value analysis.
- **`Gamma_EST.m`**
  - **Purpose**: Estimates parameters for the Gamma distribution.
  - **Description**: Fits data to a Gamma distribution.
- **`Logistic_EST.m`**
  - **Purpose**: Estimates parameters for the Logistic distribution.
  - **Description**: Fits data to a Logistic distribution.
- **`ut_reconstr.m`**
  - **Purpose**: Reconstructs datasets.
  - **Description**: Provides utilities for data reconstruction.
- **`ut_solv.m`**
  - **Purpose**: Solves mathematical problems.
  - **Description**: Provides utility functions for solving equations involved in data analysis.

## Getting Started

1. **MATLAB Installation**: Ensure MATLAB is installed and set up on your system.
2. **Data Preparation**: Download the required datasets as specified in the main project README.
3. **Script Execution**: Open MATLAB and run each script as needed. Modify input parameters based on your specific data and analysis requirements.

## Detailed Usage

1. **Bias Correction**:
   - Open `Bias_Correction.m` in MATLAB.
   - Ensure that input data paths are correctly specified.
   - Run the script to apply bias correction to your rainfall data.

2. **Data Validation**:
   - Use `Checking_TC_and_ET_cylcone_tracks.m` to validate storm track data.
   - Ensure accurate TC and ETC tracks are used in subsequent analysis.

3. **Data Accumulation**:
   - Execute `Creating_Accumulated_RF.m` to generate accumulated rainfall datasets.
   - Adjust time intervals as required.

4. **POT Analysis**:
   - Run `Creating_POT_Extremes_NTR.m` and `Creating_POT_Extremes_RF.m` to identify extreme events.
   - Use the identified extremes for further statistical analysis.

5. **Distribution Fitting**:
   - Use `Fitting_Distributions.m` to fit various distributions to your data.
   - Evaluate the best fit using the provided statistical criteria.

6. **Sampling and Stratification**:
   - Execute sampling scripts (`Sampling_Con_NTR.m` and `Sampling_Con_RF.m`) to conditionally sample your data.
   - Stratify the data using `Stratification_Con_NTR.m` and `Stratification_Con_RF.m` for detailed analysis.

7. **Utility Scripts**:
   - Use utility scripts (`DKW_conf_int.m`, `GPD_EST.m`, etc.) as needed for specific statistical calculations and data manipulations.


## License

This project is licensed under the MIT License. See the [LICENSE](../LICENSE) file for details.

## References

For more details, refer to the associated paper: [A multivariate statistical framework for mixed populations in compound flood analysis](https://doi.org/10.5194/egusphere-2024-1122).




