

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This conde will generate synthetic compound events using the peak NTR-RF combinations simulated from fitted copulas  
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
clc

% Loadinfg the selecetd extreme events extracted from the time series from 1979 to 2021
load("ETC_NTR_event_data.mat");% 
load("TC_NTR_event_data.mat");% 
load("ETC_RF_event_data.mat");% 
load("TC_RF_event_data.mat");% 

Model_Runs=[];

Monthly_Tidal_curves=load('Monthly_Tidal_Curves.mat'); % loading monthly tidal curves
MSL = load("Monthly_MSL_Distributions.mat");
TC_monthly_dis=load('PDF_TC_Monthly_disrtibution.mat');% Monthly distribution of occurence
TC_monthly_dis=table2array(struct2table(TC_monthly_dis));

ETC_monthly_dis=load('PDF_non_TC_Monthly_disrtibution.mat');% Monthly distribution of occurence
ETC_monthly_dis=table2array(struct2table(ETC_monthly_dis));

Iso_data = load('Isoline_data.mat');% Isoline data
AORC = load("AORC_Gloucester_city_cliped_data.mat");
numb_TC = 42;
numb_ETC = 518;

RF_fields=AORC.AORC_GC_rectang_clip.PRCP;

%%  Loading the simulated data

%  The below section is for genrating events. If the events are alredy
%%%%%  generated the below secton is not needed to be run %%%%%

Ob_TC = 43; % Number of TC events
Ob_Non_TC = 517; % Number of non-TC events

TC_Sim = table2array(readtable("TC_Cop_Sample.csv"));
scatter(TC_Sim(:,1),TC_Sim(:,2),[],'green','filled','.','MarkerEdgeColor','flat');hold on;
Non_TC_Sim = table2array(readtable("ETC_Cop_Sample.csv"));
scatter(Non_TC_Sim(:,1),Non_TC_Sim(:,2),[],'filled','.','MarkerEdgeColor','flat'); 
Full_Sample = [TC_Sim;Non_TC_Sim];

% calculating the return period of each event
load("Return_period_data.mat");
X=RP_data(:,1);
Y=RP_data(:,2);

RP_Combined = griddata(X, Y, RP_data(:,5), Full_Sample(:,1), Full_Sample(:,2)); % Combined population generated return period for each event in the sample
RP_TC = griddata(X, Y, RP_data(:,3), Full_Sample(:,1), Full_Sample(:,2)); % TC generated return period for each event in the sample
RP_non_TC = griddata(X, Y, RP_data(:,4), Full_Sample(:,1), Full_Sample(:,2)); % Non-Tc generated return period for each event in the sample

scatter(Full_Sample(:,1),Full_Sample(:,2),[],'red','filled','.','MarkerEdgeColor','flat'); 

DES_NTR=Full_Sample(:,1);
DES_RF=Full_Sample(:,2);



%% Creating files merging two populations
NTR_Data = [struct2table(TC_NTR_event_data);struct2table(ETC_NTR_event_data)];
RF_Data = [struct2table(TC_RF_event_data);struct2table(ETC_RF_event_data)];

% Rounding the lag time nearest integer
NTR_Data.RF_Peak_Lag_Time = round(NTR_Data.RF_Peak_Lag_Time);
RF_Data.RF_Peak_Lag_Time = round(RF_Data.RF_Peak_Lag_Time);


%% Selecting a NTR curve
% loadnig the NTR time series
NTR_thres=0.63; % Threshold used for defining NTR POT events
POT_NTR_data = load("POT_NTR_and_NTR_timeseries_for_Philli_airport.mat");

for i=1:length(DES_NTR)

% DES_NTR = simulated target NTR peak from fitted copulas
% N = number of events needed
% Power = The power of euclidian distance 
[NTR_select]=Samp_NTR_curves(DES_NTR(i),NTR_Data.Peak,NTR_Data.RF_Peak_Lag_Time,1,1,35);



% Scaling of NTR Curve


% The scaling is done only along the y axis
% NTR = [time NTR] the time series of NTR
% NTR_select = Selected NTR [value Index]
% NTR_thres = the threshold used to calculate the duration
% NTR_Data=the data file contain NTR obervation information should be
% columns in order

NTR = [POT_NTR_data.POT_NTR_and_NTR_timeseries.Time, POT_NTR_data.POT_NTR_and_NTR_timeseries.NTR];
[NTR_scaled]=Scaling_NTR_curves_M_2(DES_NTR(i),NTR,NTR_Data,NTR_select,NTR_thres);




% Combinign Scaled NTR with MSL and Tide ( The observed NTR pak lag time with next high tide
% will be used)
% RP_Comb = Combined retunr period of the event
% len_TC_sim =  Number of TC simulatins
% len_non_TC_sim =  Number of non-TC simulations
% NTR_scaled = Scaled NTR data (vector should be come from previous
% function)
% Monthly_Tidal_curves = Tidal segemnt for each month
% MSL = MSL distribution
% numb_TC =  Number of TCs recorded
% numb_ETC =  Number of ETCs reocorded
% TC_monthly_dis =  Frequency of occurence of TCs in monthly basis
% ETC_monthly_dis =  frequency of occurenece of non-TCs on monthly basis


[Storm_Tide_hydrograph]=Comb_NTR_MSL_TIDE(RP_Combined(i),i,length(TC_Sim),length(Non_TC_Sim),NTR_scaled,Monthly_Tidal_curves,MSL,TC_monthly_dis,ETC_monthly_dis);


% Assigin the output to a vector
Event(i).Storm_Tide_hydrograph=Storm_Tide_hydrograph;


% Selecting a RF event
% N= number of events slected for each target value
% Power = the power for inverse distance metrix to pick nearby events
% lag_thres = maximum allowable lag time between RF peak and NTR peak
% RF_Lag_time = array with lag times for each observed Rf events in the RF_
% data vector
% [RF_select]=Samp_RF_Events(Des_RF,RF_Data,RF_Lag_time,N,Power,lag_thres)
[RF_select]=Samp_RF_Events(DES_RF(i),RF_Data,RF_Data.RF_Peak_Lag_Time,1,1,35);


% Scaling the RF Event
% RF_fields = the evctor of observed Rf fields (3D metrix, 2d array for
% each hour)
%[RF_Scaled]=Scaling_RF_events(RF_Acc,DES_RF,RF_Data,RF_fields,RF_select)
[RF_Scaled]=Scaling_RF_events(18,DES_RF(i),RF_Data,RF_fields,RF_select);


% Combinnig Storm tide hydrograph with RFevent
[Run]=Combining_RF_WL(DES_NTR(i),DES_RF(i),RF_Scaled,Storm_Tide_hydrograph);

% assigning the data to create data file for numerical model runs
RF_Scaled.Compound_RF_field = Run.RF;
RF_Scaled.Selected_RF_NTR_lag = Run.lag_time;

% Assigin the output to a vector
Event(i).Design_RF_Field=RF_Scaled;

% % Assigning the NTR curve agin from RUN vector since some of Storm-tide hydrographs might
% be changed for independedent events
Event(i).Storm_Tide_hydrograph.Storm_Tide_HG=Run.WL;
end

save("Generated_Events.mat","Event");


%% Saving the events data in a structure file for plotting
% Note: when plotting lag times, the events without lag times (independednt
% events were removed)

Sim_events_data=struct;
for i=1:length(Event)
    Sim_events_data(i,1).Design_NTR=Event(i).Storm_Tide_hydrograph.Design_NTR;
    Sim_events_data(i,1).Design_RF=Event(i).Design_RF_Field.Design_RF;
    Sim_events_data(i,1).NTR_NTR_Event_No=Event(i).Storm_Tide_hydrograph.Event_No;
    Sim_events_data(i,1).NTR_Month=Event(i).Storm_Tide_hydrograph.Month;
    Sim_events_data(i,1).NTR_MSL=Event(i).Storm_Tide_hydrograph.MSL;
    Sim_events_data(i,1).NTR_Scaled_Falling_dur=Event(i).Storm_Tide_hydrograph.Scaled_Falling_dur;
    Sim_events_data(i,1).NTR_Scaled_Rising_dur=Event(i).Storm_Tide_hydrograph.Scaled_Rising_dur;
    Sim_events_data(i,1).NTR_Scaling_factor=Event(i).Storm_Tide_hydrograph.Scaling_factor;
    Sim_events_data(i,1).RF_scaled_hourly_peak=Event(i).Design_RF_Field.scaled_hourlu_peak;
    Sim_events_data(i,1).RF_Duration=Event(i).Design_RF_Field.Duration;
    Sim_events_data(i,1).RF_Scaling_fac=Event(i).Design_RF_Field.Scaling_fac;
    Sim_events_data(i,1).RF_Original_hourly_peak=Event(i).Design_RF_Field.Original_Peak;
    Sim_events_data(i,1).Selected_lag_time=Event(i).Design_RF_Field.Selected_RF_NTR_lag;
    Sim_events_data(i,1).RF_RF_Event_No=Event(i).Design_RF_Field.Event_No;
    Sim_events_data(i,1).NTR_Tidal_lag=Event(i).Storm_Tide_hydrograph.Original_Tidal_Lag;
    
end
save("Sim_events_data.mat","Sim_events_data");
Sim_events_data_INT = table2array(struct2table(Sim_events_data));
save("Sim_events_data_INT.mat","Sim_events_data_INT");
