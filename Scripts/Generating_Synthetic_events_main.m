%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Pravin
%
% This script generates synthetic storm events using peak NTR–RF
% combinations sampled from fitted copulas.
%
% NOTE: File paths are currently specific to the author’s directory. Update
% them as needed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% Loading historical extreme event metadata
% These include [Peak_Time Peak Duration Intensity Rising_Duration Falling_Duration Start_Ind Ending_Ind Tidal_Lag_time Zero_St_Ind Zero_end_Ind RF_Peak_Lag_Time]
% NTR events
load("non_TC_NTR_event_data.mat");
load("TC_NTR_event_data.mat");

% RF events
load("non_TC_RF_event_data.mat");
load("TC_RF_event_data.mat");

Model_Runs = [];

% Supporting data: Monthly tidal signal segments, MSL distributions, monthly storm occurrence
Monthly_Tidal_curves = load('Monthly_Tidal_Curves.mat');
MSL = load("Monthly_MSL_Distributions.mat");
TC_monthly_dis = table2array(struct2table(load('PDF_TC_Monthly_disrtibution.mat')));
non_TC_monthly_dis = table2array(struct2table(load('PDF_non_TC_Monthly_disrtibution.mat')));

% rainfall field data
AORC = load("AORC_Gloucester_city_cliped_data.mat");
RF_fields = AORC.AORC_GC_rectang_clip.PRCP;

%% Load simulated peak NTR and RF values from fitted copulas

Ob_TC = 43;
Ob_Non_TC = 517;

% Simulated copula samples for TC and non-TC events
TC_Sim = table2array(readtable("TC_Cop_Sample.csv"));
scatter(TC_Sim(:,1), TC_Sim(:,2), [], 'green', 'filled'); hold on;

Non_TC_Sim = table2array(readtable("ETC_Cop_Sample.csv"));
scatter(Non_TC_Sim(:,1), Non_TC_Sim(:,2), [], 'filled');

Full_Sample = [TC_Sim; Non_TC_Sim];

% Load computed return periods of the grid points in the parametric space
load("Return_period_data.mat");
X = RP_data(:,1);
Y = RP_data(:,2);

% Interpolating return periods of each grid point to the simulated peaks
RP_Combined = griddata(X, Y, RP_data(:,5), Full_Sample(:,1), Full_Sample(:,2));
RP_TC = griddata(X, Y, RP_data(:,3), Full_Sample(:,1), Full_Sample(:,2));
RP_non_TC = griddata(X, Y, RP_data(:,4), Full_Sample(:,1), Full_Sample(:,2));

% Assign simulated NTR and RF peaks
DES_NTR = Full_Sample(:,1);
DES_RF = Full_Sample(:,2);

%% Combine historical event metadata

% Convert structs to tables and merge TC + non-TC events
NTR_Data = [struct2table(TC_NTR_event_data); struct2table(ETC_NTR_event_data)];
RF_Data = [struct2table(TC_RF_event_data); struct2table(ETC_RF_event_data)];

% Round lag times to nearest hour (not requred)
NTR_Data.RF_Peak_Lag_Time = round(NTR_Data.RF_Peak_Lag_Time);
RF_Data.RF_Peak_Lag_Time = round(RF_Data.RF_Peak_Lag_Time);

%% Set up NTR time series
NTR_thres = 0.63; % Threshold to define POT events
POT_NTR_data = load("POT_NTR_and_NTR_timeseries_for_Philli_airport.mat");
NTR = [POT_NTR_data.POT_NTR_and_NTR_timeseries.Time, ...
       POT_NTR_data.POT_NTR_and_NTR_timeseries.NTR];

%% Main loop — generate one synthetic compound event per sample
for i = 1:length(DES_NTR)

    % Step 1: Select an observed NTR event close to simulated peak
    [NTR_select] = Samp_NTR_curves(DES_NTR(i), NTR_Data.Peak, ...
                    NTR_Data.RF_Peak_Lag_Time, 1, 1, 35);

    % Step 2: Scale the selected NTR curve to match simulated peak
    [NTR_scaled] = Scaling_NTR_curves_M_2(DES_NTR(i), NTR, ...
                    NTR_Data, NTR_select, NTR_thres);

    % Step 3: Add monthly MSL and tidal curve to get full storm tide hydrograph
    [Storm_Tide_hydrograph] = Comb_NTR_MSL_TIDE(RP_Combined(i), i, ...
                    length(TC_Sim), length(Non_TC_Sim), ...
                    NTR_scaled, Monthly_Tidal_curves, ...
                    MSL, TC_monthly_dis, non_TC_monthly_dis);

    % Save storm tide hydrograph to event structure
    Event(i).Storm_Tide_hydrograph = Storm_Tide_hydrograph;

    % Step 4: Select RF event and scale it
    [RF_select] = Samp_RF_Events(DES_RF(i), RF_Data, ...
                     RF_Data.RF_Peak_Lag_Time, 1, 1, 35);

    [RF_Scaled] = Scaling_RF_events(Acc, DES_RF(i), RF_Data, ...
                    RF_fields, RF_select);

    % Step 5: Combine water level and rainfall into one compound forcing
    [Run] = Combining_RF_WL(DES_NTR(i), DES_RF(i), RF_Scaled, Storm_Tide_hydrograph);

    % Save rainfall data and update WL hydrograph
    RF_Scaled.Compound_RF_field = Run.RF;
    RF_Scaled.Selected_RF_NTR_lag = Run.lag_time;
    Event(i).Design_RF_Field = RF_Scaled;
    Event(i).Storm_Tide_hydrograph.Storm_Tide_HG = Run.WL;

end

% Save the full event structure
save("Generated_Events.mat", "Event");

%% Save results for plotting and analysis

% Convert structure into simplified vector format
Sim_events_data = struct;
for i = 1:length(Event)
    Sim_events_data(i,1).Design_NTR = Event(i).Storm_Tide_hydrograph.Design_NTR;
    Sim_events_data(i,1).Design_RF = Event(i).Design_RF_Field.Design_RF;
    Sim_events_data(i,1).NTR_NTR_Event_No = Event(i).Storm_Tide_hydrograph.Event_No;
    Sim_events_data(i,1).NTR_Month = Event(i).Storm_Tide_hydrograph.Month;
    Sim_events_data(i,1).NTR_MSL = Event(i).Storm_Tide_hydrograph.MSL;
    Sim_events_data(i,1).NTR_Scaled_Falling_dur = Event(i).Storm_Tide_hydrograph.Scaled_Falling_dur;
    Sim_events_data(i,1).NTR_Scaled_Rising_dur = Event(i).Storm_Tide_hydrograph.Scaled_Rising_dur;
    Sim_events_data(i,1).NTR_Scaling_factor = Event(i).Storm_Tide_hydrograph.Scaling_factor;
    Sim_events_data(i,1).RF_scaled_hourly_peak = Event(i).Design_RF_Field.scaled_hourlu_peak;
    Sim_events_data(i,1).RF_Duration = Event(i).Design_RF_Field.Duration;
    Sim_events_data(i,1).RF_Scaling_fac = Event(i).Design_RF_Field.Scaling_fac;
    Sim_events_data(i,1).RF_Original_hourly_peak = Event(i).Design_RF_Field.Original_Peak;
    Sim_events_data(i,1).Selected_lag_time = Event(i).Design_RF_Field.Selected_RF_NTR_lag;
    Sim_events_data(i,1).RF_RF_Event_No = Event(i).Design_RF_Field.Event_No;
    Sim_events_data(i,1).NTR_Tidal_lag = Event(i).Storm_Tide_hydrograph.Original_Tidal_Lag;
end

% Save both as structure and array
save("Sim_events_data.mat", "Sim_events_data");
Sim_events_data_INT = table2array(struct2table(Sim_events_data));
save("Sim_events_data_INT.mat", "Sim_events_data_INT");
