%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% This function Creates storm-Tide-hydrograps combining NTR_curves, Lag
% time, Tidal curves
% Combinign Scaled NTR with MSL and ide ( The observed NTR pak lag time with next high tide
% will be used)
% RP_Comb = Combined retunr period of the event
% RP_TC =  The Return period from the TC
% RP_non_TC =  The return period from the ETc
% NTR_scaled = Scaled NTR data (vector should be come from previous
% function)
% Monthly_Tidal_curves = Tidal segemnt for each month
% MSL = MSL distribution
% numb_TC =  Number of TCs recorded
% numb_ETC =  Number of ETCs reocorded
% TC_monthly_dis =  Frequency of occurence of TCs in monthly basis
% ETC_monthly_dis =  frequency of occurenece of non-TCs on monthly basis




function[Storm_Tide_hydrograph]=Comb_NTR_MSL_TIDE(RP_Comb,event_no,len_TC_sim,len_non_TC_sim,NTR_scaled,Monthly_Tidal_curves,MSL,TC_monthly_dis,ETC_monthly_dis)

NTR_time_series_Scaled=NTR_scaled.Scaled_NTR;

%%%% Selecting the distribution probabilities for selectig a month
month=[];
Sel_dis=0;
if event_no <= len_TC_sim % Using TC probability as a basis to get the month
    Sel_dis=1;
    month=randsample(TC_monthly_dis(:,1),1,true,TC_monthly_dis(:,2)); % sampling given the probability 
else
    month=randsample(ETC_monthly_dis(:,1),1,true,ETC_monthly_dis(:,2));
    Sel_dis=2;
end

%%% old method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if Rnd <= Clc_prob_TC % Using TC probability as a basis to get the month
%     Sel_dis = 1;
%     C = rand(1,1); % this random number is to select the month
%     for i=1:length(cdf_TC)-1
%         if C<cdf_TC(1,1)
%             month = 1;
%         elseif C>cdf_TC(i,1) && C < cdf_TC(i+1,1)
%             month = i;
%         end
%     end
%     
% else
%     Sel_dis = 2;
%     C = rand(1,1);
%     for i=1:length(cdf_ETC)-1
%         if C<cdf_ETC(1,1)
%             month = 1;
%         elseif C>cdf_ETC(i,1) && C < cdf_ETC(i+1,1)
%             month = i;
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Loading the Monthly Tidal curves distribution and findnig a curve
Tidal_curves = Monthly_Tidal_curves.Monthly_Tidal_Curves(month).Month;

% taking a random curve
Curve_no = randi([1 size(Tidal_curves,1)],1,1);

% Selecting a curve
Tide = Tidal_curves(Curve_no,:);

% select lag time
lag_T = NTR_scaled.Original_Tidal_Lag;

% Loading MSL distribution and select one
monthly_MSL = MSL.Monthly_MSL(month).Month;
Selected_MSL = randsample(monthly_MSL,1);

%     % We dont check this at this point
%     % Check if the NTR signal has a periodic variability
%     % finding the local peaks
%     [pks,locs]=findpeaks(NTR_time_series);
%     A=find(diff(locs)==12);
%     if length(A)>4
%         k
%         lag_TIME_of_NTR=4;
%     end   


% Combining all together

% intrepolate tidal signal to hourly resolution (its in 15 min
% resolution)
hourly_TIDE = interp1(1:0.25:(length(Tide)/4+0.75),Tide,(1:1:(length(Tide)/4+0.75)));

% Adding WL to MSL
WL_predicted = [NaN(lag_T,1) ; hourly_TIDE']+Selected_MSL; % add the lag time of peak surge


WL_predicted_5day = WL_predicted(1:length(hourly_TIDE),1); % Selectnig the WL_prediction data for 5 days

% Find the heighest point of Scaled NTR curve
[~,max_ind_NTR]=max(NTR_time_series_Scaled);


% making the NTR array with same length as WL_predicted_5day
% Peak NTR should come to the point 85th location

% If the rising duration and falling duration both are less than 2.5
% days
%figure
if max_ind_NTR <= (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)<=(length(WL_predicted_5day)-1)/2 
    NTR_Scaled_mod = [zeros(((length(WL_predicted_5day)-1)/2)+1-max_ind_NTR,1);NTR_time_series_Scaled];
    NTR_Scaled_mod = [NTR_Scaled_mod;zeros(length(WL_predicted_5day)-length(NTR_Scaled_mod),1)];
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);

% If the rising duration higher than 2.5 and falling duration both are less than 2.5
% days
elseif  max_ind_NTR <= (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)>(length(WL_predicted_5day)-1)/2 
    NTR_Scaled_mod = [zeros(((length(WL_predicted_5day)-1)/2)+1-max_ind_NTR,1);NTR_time_series_Scaled];
    NTR_Scaled_mod(length(WL_predicted_5day)+1:end)=[];
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);

elseif  max_ind_NTR > (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)<=(length(WL_predicted_5day)-1)/2 
    
    NTR_time_series_Scaled(1:(max_ind_NTR-((length(WL_predicted_5day)-1)/2)-1))=[];
    NTR_Scaled_mod= [NTR_time_series_Scaled;zeros(length(WL_predicted_5day)-length(NTR_time_series_Scaled),1)];
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);

% If the rising duration and falling duration both are higher than 2.5
% days

elseif max_ind_NTR > (length(WL_predicted_5day)-1)/2 && (length(NTR_time_series_Scaled)-max_ind_NTR)>((length(WL_predicted_5day)-1)/2) 
    NTR_time_series_Scaled(1:(max_ind_NTR-((length(WL_predicted_5day)-1)/2)-1))=[];
    NTR_time_series_Scaled(length(WL_predicted_5day)+1:end)=[];
    NTR_Scaled_mod= NTR_time_series_Scaled;
    storm_tide = NTR_Scaled_mod+WL_predicted_5day;
    %plot(NTR_Scaled_mod);hold on;plot(WL_predicted_5day);hold on;plot(storm_tide);
else
    length(NTR_time_series_Scaled)-max_ind_NTR
    (length(WL_predicted_5day)-1)/2 
    x=5

end


% Checking the maximum WL occures in the middle are. This was a
% requirement for the flood model
%%%%%%%%%% Dont check this here %%%%%%%%%%
% if max(Storm_Tide_hydrograph(:,k))== max(Storm_Tide_hydrograph(60:108,k))
%     k=k+1;
% elseif Design_NTR < 0
%     k=k+1;
% end
% 

Storm_Tide_hydrograph=NTR_scaled;
Storm_Tide_hydrograph.TC_or_Non_TC_Dis =Sel_dis;
Storm_Tide_hydrograph.Month = month;
Storm_Tide_hydrograph.MSL = Selected_MSL;
Storm_Tide_hydrograph.Five_day_Tide = WL_predicted_5day-Selected_MSL;
Storm_Tide_hydrograph.Storm_Tide_HG = storm_tide;
Storm_Tide_hydrograph.RP = RP_Comb;


end