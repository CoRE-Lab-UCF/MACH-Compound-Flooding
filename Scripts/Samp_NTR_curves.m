
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Creating a prbability distribution based on the inverse distance of two
% points




function [NTR_select]=Samp_NTR_curves(Des_NTR,Peak_NTR,RF_Lag_time, N,Power,lag_thres)
    
    RF_Lag_time(RF_Lag_time>lag_thres,:)=NaN; % Making all exrteme lag times NaNs
    RF_Lag_time(RF_Lag_time<lag_thres*(-1),:)=NaN; % Making all exrteme lag times NaNs

    if Des_NTR>=0 % for positive des NTR values
        dis = abs(Des_NTR-Peak_NTR).^Power;
        inv_dis = 1./dis;
    else % fro negative des_NTR values
        dis = abs(Peak_NTR+Des_NTR).^Power;
        inv_dis = 1./dis;
    end

    Probs_NTR=[Peak_NTR (inv_dis./max(inv_dis))];
    
    
    lag_sel = NaN; % To run the while loop
    if Des_NTR > 0 % For compound events use only the compound observations
        while isnan(lag_sel) % The selected lag should be greater than 36 hours to run the loop
            NTR_sel=-10; % a dummy value to run the while loop (Should be less than the minimum des_NTR )
            while NTR_sel*4 < Des_NTR % Allow only 5 times of scaling
                NTR_sel = randsample(Peak_NTR,N,true,Probs_NTR(:,2));
            end
            ind = find(Peak_NTR == NTR_sel);
            NTR_select = [NTR_sel ind(1,1)];
            lag_sel = RF_Lag_time(ind);
        end
    else    % For Independend events Use any of the observation
        NTR_sel=-10; % a dummy value to run the while loop (Should be less than the minimum des_NTR )
        while NTR_sel*4 < Des_NTR % Allow only 5 times of scaling
            NTR_sel = randsample(Peak_NTR,N,true,Probs_NTR(:,2));
        end
        ind = find(Peak_NTR == NTR_sel);
        NTR_select = [NTR_sel ind(1,1)];
        lag_sel = RF_Lag_time(ind);
    end
   
end