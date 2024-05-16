# A multivariate statistical framework for mixed populations in compound flood analysis

This directory contains MATLAB scripts used for the analysis of compound flooding events in coastal regions, focusing on the distinct contributions from tropical cyclones (TCs) and non-tropical cyclones (non-TCs). The methodology and results are detailed in the paper "[A multivariate statistical framework for mixed populations in compound flood analysis](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-1122)."


## Overview of Scripts

### 1. Bias_Correction.m
- **Description**: Adjusts rainfall data to reflect larger catchment characteristics using quantile mapping methods, aligning point measurements with broader basin averages.

### 2. Checking_TC_and_ET_cylcone_tracks.m
- **Description**: Ensures the accuracy of cyclone track data utilized in simulations, which is critical for reliable analysis outcomes.

### 3. Creating_Accumulated_RF.m
- **Description**: Compiles rainfall data into accumulated datasets over set intervals, aiding in the analysis of precipitation patterns throughout and following storm events.

### 4. Creating_Isolines_From_Two_Populations.m
- **Description**: Produces isolines for comparative analysis between two distinct data sets, facilitating risk profile comparison across different flood scenarios.

### 5. Creating_POT_Extremes_NTR.m
- **Description**: Identifies extreme non-tidal residual events using peak over threshold methods to analyze upper distribution tails.

### 6. Creating_POT_Extremes_RF.m
- **Description**: Focuses on extracting extreme rainfall events that surpass predefined thresholds, crucial for detailed flood risk assessments.

### 7. Fitting_Distributions.m
- **Description**: Employs various statistical distributions to model observed data, supporting probabilistic forecasting efforts by identifying the best fit based on statistical criteria.

### 8. Sampling_Con_NTR.m
- **Description**: Implements conditional sampling on non-tidal residuals, selecting data under specified conditions for model input.

### 9. Sampling_Con_RF.m
- **Description**: Similarly, this script conditionally samples rainfall data, preparing it for detailed flood risk analysis.

### 10. Stratification_Con_NTR.m
- **Description**: Stratifies non-tidal residual data, organizing it for granular analysis under diverse environmental conditions.

### 11. Stratification_Con_RF.m
- **Description**: Stratifies rainfall data to enable nuanced analysis, essential for understanding varied flood risks.

### 12. Utility Scripts
- **DKW_conf_int.m**: Computes confidence intervals for empirical distribution functions, enhancing data reliability.
- **Dummy**: Serves as a placeholder or template, useful for testing new procedures or scripts.
- **GPD_EST.m**: Estimates parameters for the Generalized Pareto Distribution, ideal for modeling extreme values.
- **Gamma_EST.m**: Fits data to a Gamma distribution, commonly used in hydrological data analysis.
- **Logistic_EST.m**: Fits data to a Logistic distribution, useful in various statistical analyses for flood modeling.
- **ut_reconstr.m**: Assists in data reconstruction, aiding in data set preparation and correction.
- **ut_solv.m**: Solves complex mathematical problems that arise during data handling.

## Getting Started

Ensure MATLAB is installed and data paths are configured as outlined in the main project documentation. Each script is executable within MATLAB by adjusting input parameters to suit specific data sets and analysis needs.

## License

Licensed under the MIT License. See the [LICENSE](../LICENSE) file for details.

## References

Refer to the associated paper for detailed context and methodologies: [A multivariate statistical framework for mixed populations in compound flood analysis](https://doi.org/10.5194/egusphere-2024-1122).






