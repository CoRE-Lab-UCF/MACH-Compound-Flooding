%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Pravin
%
% IMPORTANT: The paths included in the script are according to the
% author's directory. Please change them accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function scales the selected rainfall field to match the target
% 18-hour RF peak, and saves the scaled fields and basin-average series.

% RF_Acc    = Accumulation window (hours) used for moving-sum calculations
% DES_RF    = Target RF value derived from fitted copulas
% RF_Data   = Structure (or table) containing RF event metadata
% RF_fields = Rainfall field data (x, y, t)
% RF_select = Selected RF event [value, index]

function [RF_Scaled]=Scaling_RF_events(RF_Acc,DES_RF,RF_Data,RF_fields,RF_select)
    RF_Scaled=struct;
    BA_Scaled_RF_Field_sec=[];
    ind = RF_select(1,2);
    st = RF_Data.Start_Ind(ind);
    end_t = RF_Data.Ending_Ind(ind);

    % Scaling factor based on the 18-hour peak RF
    scaling_fac = DES_RF/RF_Data.Peak_18_hr_RF(ind);

    % Extract and scale the RF field for the selected event window
    RF_Fields_sec = RF_fields(:,:,st:end_t);
    Scaled_RF_Field_sec = RF_Fields_sec.*scaling_fac;

    % Basin-average RF time series (mean over area)
    BA_Scaled_RF_Field_sec(:,1) = mean(Scaled_RF_Field_sec,[1,2]);

    % Scaled 18-hour peak based on basin-average moving accumulation
    scaled_peak_18_hr = max(movsum(BA_Scaled_RF_Field_sec,RF_Acc));

    % Save outputs and metadata
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
