%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Pravin
%
% IMPORTANT: The paths included in the script are according to the
% author's directory. Please change them accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates a probability distribution based on the inverse
% distance, and selects an RF event close to the target peak.
%
% N         = Number of events selected for each target value
% Power     = Power used for the inverse-distance metric to favor nearby events
% lag_thres = Maximum allowable lag time between RF peak and NTR peak [hours]
% RF_Lag_time = Array of lag times for each observed RF event in RF_Data

function [RF_select]=Samp_RF_Events(Des_RF,RF_Data,RF_Lag_time,N,Power,lag_thres)

    % Filter extreme lag times (outside ±lag_thres)
    RF_Lag_time(RF_Lag_time>lag_thres,:)=NaN; % Set extreme positive lag times to NaN
    RF_Lag_time(RF_Lag_time<lag_thres*(-1),:)=NaN; % Set extreme negative lag times to NaN

    % Compute inverse-distance weights (closer peaks get higher probability)
    dis = abs(RF_Data.Peak_18_hr_RF-Des_RF).^Power;
    inv_dis = 1./dis;

    % Build sampling weights (normalized by the maximum weight)
    Probs_RF=[RF_Data.Peak_18_hr_RF (inv_dis./max(inv_dis))];
    
    % Sample an event subject to scaling and lag constraints
    lag_sel = NaN; % Initialize to run the while loop
    if Des_RF > 0 % For compound events, use only compound observations
        while isnan(lag_sel) % Keep sampling until a valid lag is selected
            RF_sel=-0.1; % Dummy value to run the while loop (do not use zero)
            while RF_sel*4 < Des_RF % Allow up to 4× scaling
                RF_sel = randsample(RF_Data.Peak_18_hr_RF,N,true,Probs_RF(:,2));
            end
            ind = find(RF_Data.Peak_18_hr_RF(:,1) == RF_sel);
            RF_select = [RF_sel(1,1) ind(1,1)];
            lag_sel=RF_Lag_time(ind(1,1));
        end
    else % zero rainfall events (yet use the same method, but will be zero at the end )
        RF_sel=-0.1; % Dummy value to run the while loop (do not use zero)
        while RF_sel*4 < Des_RF % Allow up to 4× scaling
            RF_sel = randsample(RF_Data.Peak_18_hr_RF,N,true,Probs_RF(:,2));
        end
        ind = find(RF_Data.Peak_18_hr_RF(:,1) == RF_sel);
        RF_select = [RF_sel(1,1) ind(1,1)];

    end

end
