%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Creating a prbability distribution based on the inverse distance of two
% points
% Selecting a RF event
% N= number of events slected for each target value
% Power = the power for inverse distance metrix to pick nearby events
% lag_thres = maximum allowable lag time between RF peak and NTR peak
% RF_Lag_time = array with lag times for each observed Rf events in the RF_
% data vector
% [RF_select]=Samp_RF_Events(Des_RF,RF_Data,RF_Lag_time,N,Power,lag_thres)

function [RF_select]=Samp_RF_Events(Des_RF,RF_Data,RF_Lag_time,N,Power,lag_thres)

    RF_Lag_time(RF_Lag_time>lag_thres,:)=NaN; % Making all exrteme lag times NaNs
    RF_Lag_time(RF_Lag_time<lag_thres*(-1),:)=NaN; % Making all exrteme lag times NaNs

    dis = abs(RF_Data.Peak_18_hr_RF-Des_RF).^Power;
    inv_dis = 1./dis;

    Probs_RF=[RF_Data.Peak_18_hr_RF (inv_dis./max(inv_dis))];
    

    lag_sel = NaN; % To run the while loop
    if Des_RF > 0 % For compound events use only the compound observations
        while isnan(lag_sel) % The selected lag should be greater than 36 hours to run the loop
            RF_sel=-0.1; % a dummy value to run the while loop (dont put zero since Des_NTR might == 0)
            while RF_sel*4 < Des_RF % Allow only 5 times of scaling
                RF_sel = randsample(RF_Data.Peak_18_hr_RF,N,true,Probs_RF(:,2));
            end
            ind = find(RF_Data.Peak_18_hr_RF(:,1) == RF_sel);
            RF_select = [RF_sel(1,1) ind(1,1)];
            lag_sel=RF_Lag_time(ind(1,1));
        end
    else
        RF_sel=-0.1; % a dummy value to run the while loop (dont put zero since Des_NTR might == 0)
        while RF_sel*4 < Des_RF % Allow only 5 times of scaling
            RF_sel = randsample(RF_Data.Peak_18_hr_RF,N,true,Probs_RF(:,2));
        end
        ind = find(RF_Data.Peak_18_hr_RF(:,1) == RF_sel);
        RF_select = [RF_sel(1,1) ind(1,1)];

    end

end