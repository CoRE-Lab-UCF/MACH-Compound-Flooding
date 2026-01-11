%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Pravin
%
% IMPORTANT: The paths included in the script are according to the
% author's directory. Please change them accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This function creates storm-tide hydrographs by combining the scaled NTR
% curve with MSL and a tidal segment.
% The observed NTR peak lag time relative to the next high tide is used.

% RP_Comb            = Combined return period of the event
% event_no           = Event number in the simulated sample
% len_TC_sim         = Number of TC simulations from copulas
% len_non_TC_sim     = Number of non-TC simulations from copulas
% NTR_scaled         = Scaled NTR data structure from the previous function
% Monthly_Tidal_curves = Tidal signal segments around each high tide, separated by month
% MSL                = MSL distribution separated by calendar month
% TC_monthly_dis     = Monthly occurrence distribution for TCs
% ETC_monthly_dis    = Monthly occurrence distribution for non-TCs
% Storm_Tide_hydrograph = Structure containing the storm-tide hydrograph and related parameters

function[Storm_Tide_hydrograph]=Comb_NTR_MSL_TIDE(RP_Comb,event_no,len_TC_sim,len_non_TC_sim,NTR_scaled,Monthly_Tidal_curves,MSL,TC_monthly_dis,ETC_monthly_dis)

NTR_time_series_Scaled=NTR_scaled.Scaled_NTR;

%%%% Selecting the distribution probabilities for selecting a month
month=[];
Sel_dis=0;
if event_no <= len_TC_sim % Use TC monthly distribution for TC simulations
    Sel_dis=1;
    month=randsample(TC_monthly_dis(:,1),1,true,TC_monthly_dis(:,2)); % Sample month using the probability distribution
else
    month=randsample(ETC_monthly_dis(:,1),1,true,ETC_monthly_dis(:,2)); % Sample month using the probability distribution
    Sel_dis=2;
end


% Load monthly tidal curves and select one curve randomly
Tidal_curves = Monthly_Tidal_curves.Monthly_Tidal_Curves(month).Month;

% Taking a random curve
Curve_no = randi([1 size(Tidal_curves,1)],1,1);

% Selecting a curve
Tide = Tidal_curves(Curve_no,:);

% Select lag time (NTR peak lag relative to the next high tide)
lag_T = NTR_scaled.Original_Tidal_Lag;

% Load monthly MSL distribution and sample one value
monthly_MSL = MSL.Monthly_MSL(month).Month;
Selected_MSL = randsample(monthly_MSL,1);

%     % We don't check this at this point
%     % Check if the NTR signal has periodic variability
%     % Finding the local peaks
%     [pks,locs]=findpeaks(NTR_time_series);
%     A=find(diff(locs)==12);
%     if length(A)>4
%         k
%         lag_TIME_of_NTR=4;
%     end   


% Combine scaled NTR, tide, and MSL

% Interpolate tidal signal to hourly resolution (original is 15-min resolution)
hourly_TIDE = interp1(1:0.25:(length(Tide)/4+0.75),Tide,(1:1:(length(Tide)/4+0.75)));

% Add MSL and apply the lag time for the peak surge
WL_predicted = [NaN(lag_T,1) ; hourly_TIDE']+Selected_MSL; % Add lag time to align the tide

WL_predicted_5day = WL_predicted(1:length(hourly_TIDE),1); % Select the WL prediction for the 5-day window

% Find the peak location of the scaled NTR curve
[~,max_ind_NTR]=max(NTR_time_series_Scaled);

% Align the scaled NTR time series to the 5-day WL window (peak near the center)
% Peak NTR is aligned close to the mid-point of WL_predicted_5day

% If both rising and falling durations are <= 2.5 days
%figure
if max_ind_NTR <= (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)<=(length(WL_predicted_5day)-1)/2 
    NTR_Scaled_mod = [zeros(((length(WL_predicted_5day)-1)/2)+1-max_ind_NTR,1);NTR_time_series_Scaled];
    NTR_Scaled_mod = [NTR_Scaled_mod;zeros(length(WL_predicted_5day)-length(NTR_Scaled_mod),1)];
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);

% If rising duration is <= 2.5 days and falling duration is > 2.5 days
elseif  max_ind_NTR <= (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)>(length(WL_predicted_5day)-1)/2 
    NTR_Scaled_mod = [zeros(((length(WL_predicted_5day)-1)/2)+1-max_ind_NTR,1);NTR_time_series_Scaled];
    NTR_Scaled_mod(length(WL_predicted_5day)+1:end)=[];
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);

% If rising duration is > 2.5 days and falling duration is <= 2.5 days
elseif  max_ind_NTR > (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)<=(length(WL_predicted_5day)-1)/2 
    
    NTR_time_series_Scaled(1:(max_ind_NTR-((length(WL_predicted_5day)-1)/2)-1))=[];
    NTR_Scaled_mod= [NTR_time_series_Scaled;zeros(length(WL_predicted_5day)-length(NTR_time_series_Scaled),1)];
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);

% If both rising and falling durations are > 2.5 days
elseif max_ind_NTR > (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)>((length(WL_predicted_5day)-1)/2) 
    NTR_time_series_Scaled(1:(max_ind_NTR-((length(WL_predicted_5day)-1)/2)-1))=[];
    NTR_time_series_Scaled(length(WL_predicted_5day)+1:end)=[];
    NTR_Scaled_mod= NTR_time_series_Scaled;
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);
else
    length(NTR_time_series_Scaled)-max_ind_NTR; % checking for any errors
    (length(WL_predicted_5day)-1)/2 ; % checking for any errors
    x=5; % checking for any errors

end

% Checking the maximum WL occurs in the middle area (model requirement)
%%%%%%%%%% Don't check this here %%%%%%%%%%
% if max(Storm_Tide_hydrograph(:,k))== max(Storm_Tide_hydrograph(60:108,k))
%     k=k+1;
% elseif Design_NTR < 0
%     k=k+1;
% end
% 

% Save outputs to the structure
Storm_Tide_hydrograph=NTR_scaled;
Storm_Tide_hydrograph.TC_or_Non_TC_Dis =Sel_dis;
Storm_Tide_hydrograph.Month = month;
Storm_Tide_hydrograph.MSL = Selected_MSL;
Storm_Tide_hydrograph.Five_day_Tide = WL_predicted_5day-Selected_MSL;
Storm_Tide_hydrograph.Storm_Tide_HG = storm_tide;
Storm_Tide_hydrograph.RP = RP_Comb;

end
