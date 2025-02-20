%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Creating a prbability distribution based on the inverse distance of two
% points
% Des_NTR = Design NTR value
% NTR = NTR_time series [Time NTR]
% Strt = Starting of the event
% ending =  Ending of the event
% Rising ratio = Rising Ratio
% Falling ratio = Falling Ratio
% NTR_Data = The structure file with NTR data
% NTR_select = Selected NTR [value Index]

function [RF_Scaled]=Scaling_RF_events(RF_Acc,DES_RF,RF_Data,RF_fields,RF_select)
    RF_Scaled=struct;
    BA_Scaled_RF_Field_sec=[];
    ind = RF_select(1,2);
    st = RF_Data.Start_Ind(ind);
    end_t = RF_Data.Ending_Ind(ind);

    scaling_fac = DES_RF/RF_Data.Peak_18_hr_RF(ind);


    RF_Fields_sec = RF_fields(:,:,st:end_t);
    Scaled_RF_Field_sec = RF_Fields_sec.*scaling_fac;

    BA_Scaled_RF_Field_sec(:,1) = mean(Scaled_RF_Field_sec,[1,2]);

    scaled_peak_18_hr = max(movsum(BA_Scaled_RF_Field_sec,RF_Acc));

    RF_Scaled.Original_Peak_Time=RF_Data.Peak_Time(ind);
    RF_Scaled.Original_Peak=RF_Data.Peak(ind);
    RF_Scaled.Event_No=ind;
    RF_Scaled.Original_Start_ind=RF_Data.Start_Ind(ind);
    RF_Scaled.Original_end_ind=RF_Data.Ending_Ind(ind);
    RF_Scaled.Original_peak_18_hr = RF_Data.Peak_18_hr_RF(ind);
    RF_Scaled.Original_intensity = RF_Data.Intensity(ind);
    RF_Scaled.Original_RF_peak_lag = RF_Data.RF_Peak_Lag_Time(ind);
    RF_Scaled.Original_RF_field = RF_Fields_sec;
    RF_Scaled.Duration = RF_Data.Duartion(ind);
    RF_Scaled.Design_RF=DES_RF;
    RF_Scaled.Scaling_fac= scaling_fac;
    RF_Scaled.Scaled_RF_field = Scaled_RF_Field_sec;
    RF_Scaled.BA_Scaled_RF_Field_sec= BA_Scaled_RF_Field_sec;
    RF_Scaled.scaled_peak_18_hr= scaled_peak_18_hr;
    RF_Scaled.scaled_hourlu_peak= max(BA_Scaled_RF_Field_sec);

  
end