%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Pravin
%
% IMPORTANT: The paths included in the script are according to the
% author’s directory. Please change them accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates a probability distribution based on the inverse
% distance between the target NTR and observed NTR peaks.

% Des_NTR     = Target NTR value derived from fitted copulas
% Peak_NTR    = Peak NTR values of all extracted NTR events [Date Value]
% RF_Lag_time = RF lag times of all extracted NTR events [hours]
% N           = Number of events to be sampled
% Power       = Power used for the inverse distance (typically 1)
% lag_thres   = Lag-time threshold to identify independent events [hours] (here 36)
% NTR_select  = Selected NTR event [Peak value, Peak index]

function [NTR_select]=Samp_NTR_curves(Des_NTR,Peak_NTR,RF_Lag_time, N,Power,lag_thres)
    
    % Filter extreme lag times (outside ±lag_thres)
    RF_Lag_time(RF_Lag_time>lag_thres,:)=NaN; % Set extreme positive lag times to NaN
    RF_Lag_time(RF_Lag_time<lag_thres*(-1),:)=NaN; % Set extreme negative lag times to NaN

    % Compute inverse-distance weights (closer peaks get higher probability)
    if Des_NTR>=0 % For positive target NTR values
        dis = abs(Des_NTR-Peak_NTR).^Power;
        inv_dis = 1./dis;
    else % For negative target NTR values
        dis = abs(Peak_NTR+Des_NTR).^Power;
        inv_dis = 1./dis;
    end

    % Build sampling weights (normalized by the maximum weight)
    Probs_NTR=[Peak_NTR (inv_dis./max(inv_dis))];
    
    % Sample an event subject to scaling and lag constraints
    lag_sel = NaN; % Initialize to run the while loop
    if Des_NTR > 0 % For compound events, use only compound observations
        while isnan(lag_sel) % Keep sampling until a valid lag is selected
            NTR_sel=-10; % Dummy value to run the while loop (must be lower than min Des_NTR)
            while NTR_sel*4 < Des_NTR % Allow up to 4× scaling
                NTR_sel = randsample(Peak_NTR,N,true,Probs_NTR(:,2));
            end
            ind = find(Peak_NTR == NTR_sel);
            NTR_select = [NTR_sel ind(1,1)];
            lag_sel = RF_Lag_time(ind);
        end
    else    % For independent events, allow any observation
        NTR_sel=-10; % Dummy value to run the while loop (must be lower than min Des_NTR)
        while NTR_sel*4 < Des_NTR % Allow up to 4× scaling
            NTR_sel = randsample(Peak_NTR,N,true,Probs_NTR(:,2));
        end
        ind = find(Peak_NTR == NTR_sel);
        NTR_select = [NTR_sel ind(1,1)];
        lag_sel = RF_Lag_time(ind);
    end
   
end
