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

function [Design_NTR_time_series]=Scaling_NTR_curves_M_2(Des_NTR,NTR,NTR_Data,NTR_select,NTR_thres)

    Design_NTR_time_series=struct;

    for i=1:size(NTR_select,1)
        ind = NTR_select(i,2);

        Strt = NTR_Data.Zero_St_Ind(ind);
        ending = NTR_Data.Zero_end_Ind(ind);
    

        Design_NTR_time_series(i).Original_Date = NTR_Data.Peak_Time(ind);
        Design_NTR_time_series(i).Original_NTR_Peak = NTR_select(i,1);
        Design_NTR_time_series(i).Original_Rising_dur = NTR_Data.Rising_Duration(ind);
        Design_NTR_time_series(i).Original_Falling_dur = NTR_Data.Falling_Duration(ind);
        Design_NTR_time_series(i).Original_Intensity = NTR_Data.Intensity(ind);
        Design_NTR_time_series(i).Original_Tidal_Lag = NTR_Data.Tidal_Lag_time(ind);
        Design_NTR_time_series(i).Original_RF_Lag = NTR_Data.RF_Peak_Lag_Time(ind);
        Design_NTR_time_series(i).Design_NTR = Des_NTR;
        Design_NTR_time_series(i).Event_No = ind;


        NTR_Sec = NTR(Strt:ending,2);
        Time_Sec = NTR(Strt:ending,1);

        %%%%%%% Finding the duration of the NTR event to match the rising Scaled %
        % NTR event (Only use the rising section Ratio)
       
        Scale_fac = Des_NTR/NTR_select(i,1); % calculating the scaling factor

        % Scaling the time series
        Scaled_NTR = NTR_Sec.*Scale_fac;
        

        % Saving to the data file
      
        Design_NTR_time_series(i).Scaled_NTR = Scaled_NTR;
        Design_NTR_time_series(i).Original_NTR = NTR_Sec;
        Design_NTR_time_series(i).Scaling_factor = Scale_fac;

        [~,time_NTR_ind] = max(Scaled_NTR);
        
        eta_bf=Scaled_NTR(time_NTR_ind);
        k=1;
        while eta_bf > NTR_thres
            eta_bf=Scaled_NTR(time_NTR_ind-k);
            k=k+1;
           if time_NTR_ind-k < 1
               break
           end
        end
        dur_bf = k-1;
    
        eta_aft=Scaled_NTR(time_NTR_ind);
        m=1;
        while eta_aft > NTR_thres
            eta_aft=Scaled_NTR(time_NTR_ind+m);
            m=m+1;
            if time_NTR_ind+m > length(Scaled_NTR)
                break
            end
        end
        dur_aft = m-1;

        Design_NTR_time_series(i).Scaled_Rising_dur = dur_bf;
        Design_NTR_time_series(i).Scaled_Falling_dur = dur_aft;
        
    end

    %plot(Design_NTR_time_series(1).Original_NTR);hold on;
    %plot(Design_NTR_time_series(1).Scaled_NTR);hold on
    

end