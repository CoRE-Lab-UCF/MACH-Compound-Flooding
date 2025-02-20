%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Storm_Tide_hydrograph = The vector os Storm Tide hydrographs
% RF_fields = 3D vector with selected RF events
% Dess_NTR = Target NTR values simulated from fitted copulas
% Dess_RF = Target RF values simulated from fitted copulas (for the selected
% Accumulations time)
% Des_RF_Acc = Design_RF_Accumulation
% Run = Structure array with the synthetc rainfall field and still walter
% level time series

function [Run]=Combining_RF_WL(Dess_NTR,Dess_RF,RF_Scaled,Storm_Tide_hydrograph)

NTR_lag = Storm_Tide_hydrograph.Original_RF_Lag;
RF_lag = RF_Scaled.Original_RF_peak_lag;

STH = Storm_Tide_hydrograph.Storm_Tide_HG;
RF = RF_Scaled.Scaled_RF_field;



% finding the peak indices
[~,NTR_pk_ind] = max(Storm_Tide_hydrograph.Scaled_NTR);
[~,RF_pk_ind] = max(RF_Scaled.BA_Scaled_RF_Field_sec);

% Option 1
% Slecting the minimum lag time
% [~,ind] = min([abs(NTR_lag),abs(RF_lag)]);
% if ind ==1
%     lag = NTR_lag;
% else
%     lag = RF_lag;
% end
% 
% if lag >= 0
%     if lag > lag_thres
%         lag=lag_thres;
%     end
% else
%     if lag < lag_thres*(-1)
%         lag=lag_thres*(-1);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Option 2
% Selecting a random lag time out of two
lagss = 1:1:11; % A vector of possible leg times before the high tide
if Dess_NTR > 0 && Dess_RF >0 % Compound events
    LG = [NTR_lag RF_lag];
    lag = randsample(LG,1);
   
else
    lag = randsample(lagss,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% check whether the lag times within the thresholds

% We add the lag time to RF time series    
Z = zeros(size(RF,1),size(RF,2),((length(STH)+1)/2)-lag-RF_pk_ind);

RF=cat(3,Z,RF);

if Dess_NTR >0
    Run.WL=STH;
else
    Run.WL=Storm_Tide_hydrograph.Five_day_Tide+Storm_Tide_hydrograph.MSL;
end

Run.RF=RF;
Run.Des_RF=RF_Scaled.Design_RF;
Run.Des_NTR=Storm_Tide_hydrograph.Design_NTR;
Run.lag_time = lag;

end

