# A multivariate statistical framework for mixed populations in compound flood analysis

This directory contains MATLAB and R scripts used for estimating the joint probabilities in coastal regions, focusing on the mixed populations (tropical cyclones (TCs) and non-tropical cyclones (non-TCs)). The methodology and results are detailed in the paper "[A multivariate statistical framework for mixed populations in compound flood analysis](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-1122)."

<br>

>Maduwantha, P., Wahl, T., Santamaria-Aguilar, S., Jane, R. A., Booth, J. F., Kim, H., and Villarini, G.: A multivariate statistical framework for mixed populations in compound flood analysis, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2024-1122, 2024.

<br>

## Overview of Scripts
**The input and output files are generated by the MATLAB or R scripts (except the raw data inputs). Under the "Input" files, the file name is given first following the format of the file**

### 1. Bias Correction
- **`Bias_Correction.m`**
  - **Description**: This script applies bias correction to rainfall data using quantile mapping. The method adjusts the hourly rainfall gauge data to match the basin-average values derived from the Analysis of Record for Calibration (AORC) data.-

**Input**:
- `Pr_Measured.mat`:- [date-time, hourly precipitation]
- `Pr_AORC.mat`:- [date-time, hourly precipitation]

**output**: 
- `Bias Corrected MS at {location}_With_AORC_of_GC.mat`:- [date-time, hourly precepitation]

 
### 2. Calculating accumulated rainfall 
- **`Creating_Accumulated_RF.m`**
  - **Description**: The script calculates the accumulated rainfall from 1 to 48 hours suing the hourly rainfall data.
    
**Input**: 
- `Bias Corrected MS at “location”_With_AORC_of_GC.mat`:- [date-time, hourly precipitation]

**output**:  
- `Hourly_accumulation_Bias_Corrected_RF_data_”location”.mat`:- [date-time, hourly accumulated Rf (nx48)]



   

### 3. 	Defining peak over threshold extreme events
- **`Scripts/Creating_POT_Extremes_NTR.m`**
  - **Description**: The script uses the Peak over threshold (POT) approach to define events over thresholds for non-tidal residuals (NTR).
Functions: ut_solv, ut_reconstr are needed and should be downloaded through U-tide package (https://www.mathworks.com/matlabcentral/fileexchange/46523-utide-unified-tidal-analysis-and-prediction-functions). The years with substantial data gaps are removed manually. 

**Input**:
- `Water level data` :- [date-time, hourly water levels]

**output**: 
- `POT_NTR_and_NTR_timeseries.mat` :- structure[NTR.POT.WL_raw.Tide.MovA_WL.Threshold]


- **`Scripts/Creating_POT_Extremes_RF.m`**
  - **Description**: The script uses the Peak over threshold (POT) approach to define events over thresholds for basin average rainfall. The code calculates the threshold exceedances for rainfall accumulation times from 1 to 48.

**Input**:
- `Hourly_accumulation_Bias_Corrected_RF_data_”location”.mat`:- [date-time, hourly accumulated Rf (nx48)]

**output**: 
- `POT_Events_for_each_RF_Accumulation_time.mat`: -structure[events[POT.Threshold]




### 4. Two way sampling
4. 1 Conditioned on NTR
  - **Description**: The script finds the maximum accumulated rainfall (for all the hourly accumulations from 1 to 48) within a given time window around the selected POT event when conditioning NTR.
- **`Scripts/Sampling_Con_NTR.m`**

**Input**:
- `Hourly_accumulation_Bias_Corrected_RF_data_”location”.mat`:- [date-time, hourly accumulated Rf (nx48)]
- `POT_NTR_and_NTR_timeseries.mat`:- structure[NTR.POT.WL_raw.Tide.MovA_WL.Threshold]

**output**:  
- `Maximum_RF_events_for_each_POT_NTR.ma`:-structure[Acumulation[Time_NTR.POT_NTR.Time_RF.Max_RF]

4. 2 Conditioned on RF
  - **Description**: The script finds the maximum NTR (for all the selected accumulation time) within a given time window around the selected POT event when conditioning RF.
- **`Scripts/Sampling_Con_RF.m`**

**Input**:
- `POT_NTR_and_NTR_timeseries.mat`:- structure[NTR.POT.WL_raw.Tide.MovA_WL.Threshold]
- `POT_Events_for_each_RF_Accumulation_time.mat`: -structure[events[POT.Threshold]


**output**: 
- `Maximum_NTR_events_for_each_POT_RF_for_”Selected_accumulation_time”_RF_acc.mat`:-structure [Time_NTR.Max_NTR.Time_RF.POT_RF]

## 5. Stratification
  - **Description**: The two conditional samples are stratified in to two sets as: 1. The events caused by Tropical cyclones 2. The events that were not caused by tropical cyclones. The Hurdat 2  data set is used to identify the events induced by tropical cylcones.
    
5.1 Conditioned on NTR
- **`Scripts/Stratification_Con_NTR.m`**

**Input**: 
- `Maximum_RF_events_for_each_POT_NTR.mat`:-structure[Acumulation[Time_NTR.POT_NTR.Time_RF.Max_RF]
- `Cyclone_Track_data_from_1850.mat`:- [date-time, lat, lon, speed, distance from the city]

**output**: 
- `TC_events_conditioning_POT_NTR.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]
- `ETC_events_conditioning_POT_NTR.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]

5.2 Conditioned on RF
- **`Scripts/Stratification_Con_RF.m`**

**Input**:
- `Maximum_NTR_events_for_each_POT_RF_for_”Selected_accumulation_time”_RF_acc.mat`:-structure [Time_NTR.Max_NTR.Time_RF.POT_RF]
- `Cyclone_Track_data_from_1850.mat`:- [date-time, lat, lon, speed, distance from the city]

**output**: 
- `TC_events_conditioning_POT_RF.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]
- `ETC_events_conditioning_POT_RF.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]

## 6. Fitting distributions
  - **Description**: The two stratified conditional samples are fitted in to different parametric distributions to check the most appropriate distribution. The selection of the distribution is done based on the AIC and the “Multihazard R package” was used (https://rdrr.io/github/rjaneUCF/MultiHazard/f/README.md). Following codes and functions will calculate the confidence interval and will fit the selected distributions to the respective samples.
- **`Scripts/Fitting_Ditributions.m`**
- **`Scripts/DKW_conf_int.m`**


## 7. Fitting Copulas
  - **Description**: The two stratified conditional samples are fitted in to different Copula families to check the most appropriate copula family. This is done based on the AIC and the “Multihazard R package” was used (https://rdrr.io/github/rjaneUCF/MultiHazard/f/README.md). 

## 8. Bivariate analysis and estimating Joint probability distributions
  - **Description**: The two stratified conditional samples are used to:
1.	Calculate the annual exceedance probabilities for all the selected combinations of NTR and RF in the parametric space.
2.	Generate N number of realizations from the fitted copulas.
Here the following functions (modified from “Multihazard R package” (https://rdrr.io/github/rjaneUCF/MultiHazard/f/README.md) are used for the calculations separately for tropical cyclones and non-tropical cyclones. 

For Tropical Cyclone (TC) events:
- **`Scripts/ Copula_TC_BC. R`**

**Input**:
- `Data`:Hourly time series of NTR and RF [NTR, RF] , 
- `Data_Con1`:TC POT events conduiting NTR for last 30 years, 
- `Data_Con2`: TC POT events conduiting NTR for last 30 years, 
- `Thres1` :Threshold for NTR, 
- `Thres2` :Threshold for RF, 
- `Copula_Family1`: Character array of selected copula _Family_NTR,
- `Copula_Family2`: Character array of selected copula _Family_RF, 
- `Marginal_Dist1`:Selected Marginal_Dist. For NTR, 
- `Marginal_Dist2`: Selected Marginal_Dist. For RF,
- `RP`: an array of return periods interested (c(5,10,20,50,100)),
- `N`:Number of realizations expected,
- `N_Ensemble`:Number of realizations expected in each isoline,
- `mu` : Annual rate of data availability, 
- `Con1`:Character array for the conditioning sample 1 ("NTR"),
- `Con2`:Character array for the conditioning sample 2 ("RF"),
- `GPD_con1`:Selected GPD_NTR,
- `GPD_con2`:Selected GPD_RF,
- `Data_Con1_M`:TC POT sample Conditioning_NTR [POT_NTR max_RF],
- `Data_Con2_M:` TC POT sample Conditioning_RF[Max_NTR POT RF],
- `EL_con1`:Average inter arrival time of TC threshold exceedances conditioning  NTR,
- `EL_con2`: Average inter arrival time of TC threshold exceedances conditioning  RF,

**output**: 
- `RP_TC.csv` : CSV file containing calculated return periods for all the selected combinations of NTR and RF in the parametric space based on TC samples. [NTR, RF, return period cond. NTR, return period cond. RF]
- `TC_Cop_Sample.csv`: CSV file containing N number of realizations based on the probability distribution of TC samples. [NTR, RF]

For non-Tropical Cyclone events:
- **`Scripts/ Copula_TC_BC. R`**

**Input**:
- `Data`:Hourly time series of NTR and RF [NTR, RF] , 
- `Data_Con1`:non-TC POT events conduiting NTR for last 30 years[POT_NTR max_RF],
- `Data_Con2`: non-TC POT events conduiting NTR for last 30 years[Max_NTR POT RF],
- `Thres1` :Threshold for NTR, 
- `Thres2` :Threshold for RF, 
- `Copula_Family1`: Character array of selected copula _Family_NTR,
- `Copula_Family2`: Character array of selected copula _Family_RF, 
- `Marginal_Dist1`:Selected Marginal_Dist. For NTR, 
- `Marginal_Dist2`: Selected Marginal_Dist. For RF,
- `RP`: an array of return periods interested (c(5,10,20,50,100)),
- `N`:Number of realizations expected,
- `N_Ensemble`:Number of realizations expected in each isoline,
- `mu` : Annual rate of data availability, 
- `Con1`:Character array for the conditioning sample 1 ("NTR"),
- `Con2`: Character array for the conditioning sample 2 ("RF"),
- `GPD_con1`:Selected GPD_NTR,
- `GPD_con2`:Selected GPD_RF,
- `Data_Con1_M`:non-TC POT sample Conditioning_NTR [POT_NTR max_RF],
- `Data_Con2_M`: non-TC POT sample Conditioning_RF [Max_NTR POT RF],
- `EL_con1`:Average inter arrival time of non-TC threshold exceedances conditioning  NTR,
- `EL_con2`: Average inter arrival time of non-TC threshold exceedances conditioning  RF,

**output**: 
- `RP_ETC.csv` : CSV file containing calculated return periods for all the selected combinations of NTR and RF in the parametric space based on non-TC samples. [NTR, RF, return period cond. NTR, return period cond. RF]
- `ETC_Cop_Sample.csv`: CSV file containing N number of realizations based on the probability distribution. Of non-TC samples [NTR, RF]

## 9. Combining popuations
  - **Description**: The separately estimated return periods (annual exceedance probabilities) of TC samples and non-TC samples are combined:
1.	Selecting the maximum annual exceedance probability from two conditioned samples of each population
2.	Calculating total annual exceedance probability from both populations.
- **`Scripts/Creating_Isolines_From_Two_Populations.m`**

**Input**:
- `RP` : character array of the return periods of interest {'5','10','20','50','100'}
- `Rp`: numeric array of return periods of interest [5 10 20 50 100]
- `TC` :RP_TC.csv
- `ETC` : RP_ETC.csv
- `TC_ext_CON_ntr_yrs` : TC_events_conditioning_POT_NTR.mat
- `TC_ext_CON_RF_yrs` : TC_events_conditioning_POT_RF.mat
- `Non_TC_ext_CON_ntr_yrs` :ETC_events_conditioning_POT_NTR.mat
- `Non_TC_ext_CON_RF_yrs` : ETC_events_conditioning_POT_RF.mat
- `COP_sample_TC` : TC_Cop_Sample.csv
- `COP_sample_ETC` : ETC_Cop_Sample.csv

Following variables are only for plotting.
- `c_map` : [ linspace(0.8,1,512)', linspace(0,0.9,512)', linspace(0,0.2,512)']; % Colour Map
- `F_Size` :15; % font size of the figures
- `Con_NTR_LW` : 1; % Lower limit of the NTR axis
- `Con_NTR_M_Size` : 6; % maximum of the NTR axis
- `Con_RF_LW` : 1;% Lowe limit of the RF axis
- `Con_RF_M_Size` : 4;% maximum of the RF axis
- `width` : 20; % width of figrues
- `height` : 20; % height of figures
- `resolution` : 1000; % Number of grid cells in each direction of discretizing
- `x0`:0.2;
- `y0`:0.2;
- `x_U_lim` : 8; % limits of discretizing
- `x_L_lim` : -1.2; % limits of discretizing
- `Y_U_ilm` : 300;% limits of discretizing
- `Y_L_lim`: 0;% limits of discretizing



## Getting Started

1. **MATLAB Installation**: Ensure MATLAB is installed and set up on your system.
2. **Data Preparation**: Download the required datasets as specified in the associated paper.
3. **Script Execution**: Open MATLAB and run each script as needed. Modify input parameters based on your specific data and analysis requirements.


## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References

Refer to the associated paper for detailed context and methodologies: [A multivariate statistical framework for mixed populations in compound flood analysis](https://doi.org/10.5194/egusphere-2024-1122).






