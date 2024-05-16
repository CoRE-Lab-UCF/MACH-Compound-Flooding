# A multivariate statistical framework for mixed populations in compound flood analysis

This directory contains MATLAB scripts used for the analysis of compound flooding events in coastal regions, focusing on the distinct contributions from tropical cyclones (TCs) and non-tropical cyclones (non-TCs). The methodology and results are detailed in the paper "[A multivariate statistical framework for mixed populations in compound flood analysis](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-1122)."


## Overview of Scripts

### Data Preparation and Correction

- `Bias_Correction.m`: Corrects biases in meteorological data, making it more accurate for model simulations.
- `Checking_TC_and_ET_cylcone_tracks.m`: Validates the tracks of tropical and extratropical cyclones used in the study to ensure accuracy in simulation inputs.
- `Creating_Accumulated_RF.m`: Generates accumulated rainfall datasets over specified time intervals to analyze during and after storm events.

### Statistical Analysis

- `Creating_Isolines_From_Two_Populations.m`: Produces isolines based on statistical data from two distinct populations to compare their flood risk profiles.
- `Creating_POT_Extremes_NTR.m` and `Creating_POT_Extremes_RF.m`: Identify peak over threshold (POT) extremes for non-tidal residuals and rainfall, respectively, crucial for understanding the upper tails of the distributions.
- `Fitting_Distributions.m`: Fits various statistical distributions to the observed data, which is essential for probabilistic forecasting and simulation.

### Simulation and Sampling

- `Sampling_Con_NTR.m` and `Sampling_Con_RF.m`: Conducts conditional sampling on non-tidal residuals and rainfall data for use in simulations.
- `Stratification_Con_NTR.m` and `Stratification_Con_RF.m`: Stratifies data to allow for detailed analysis under varying conditions, facilitating more accurate flood risk assessments.

### Utility Scripts

- `ut_reconstr.m` and `ut_solv.m`: Provide utility functions such as reconstructing data sets and solving mathematical problems involved in other scripts.

## Getting Started

Ensure MATLAB is installed and set up appropriately. Data files and dependencies must be configured as per the project specifications.

## Usage

Each script can be executed from the MATLAB command window or integrated into a larger MATLAB project. Parameters may need to be adjusted based on specific project requirements.



