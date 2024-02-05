
clear
clc
load('C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa\Two Way Sampling\Con_RF\18_hr_ACC\Maximum_NTR_events_for_each_POT_RF_for_18hr_RF_acc.mat');
load('C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport\Cyclone_Track\Cyclone_Track_data_from_1850_new.mat');

L=18; %Selected accumulation
%% 

distance = km2deg(350); % in degrees 300 km 
Tropical_Cyclones_Con_RF = struct;
ET_Cyclones_Con_RF = struct;
count_ETC =0;
count_TC = 0;
%%
for i = L % only done for the selected accumulation
    m=1;n=1;
    for j =1:580

        Time_RF = Maximum_NTR_events_of_POT_RF(j).Time_RF;
        RF = Maximum_NTR_events_of_POT_RF(j).POT_RF;

        Time_ntr = Maximum_NTR_events_of_POT_RF(j).Time_NTR;
        NTR = Maximum_NTR_events_of_POT_RF(j).Max_NTR;

    
        T_elements = datevec(Time_RF);
        hr = T_elements(1,4);
        hr_new = floor(hr/6)*6; 
        T_elements(1,4)= hr_new;
        Time_NTR_mod = datenum(T_elements);
        strt_hr = Time_NTR_mod-1.5;
        end_hr = Time_NTR_mod+1.75;
        
        index = find(Cyclone_track_data(:,1)>=strt_hr & Cyclone_track_data(:,1)<=end_hr);
    
        if length(index)>0
            if min(Cyclone_track_data(index,5))<= distance
               Tropical_Cyclones(i).Event(m).Time_NTR = Time_ntr;
               Tropical_Cyclones(i).Event(m).Max_NTR = NTR;
               Tropical_Cyclones(i).Event(m).Time_RF = Time_RF;
               Tropical_Cyclones(i).Event(m).POT_RF = RF;
               m=m+1;
               count_TC = count_TC+1;
            else
               ET_Cyclones(i).Event(n).Time_NTR = Time_ntr;
               ET_Cyclones(i).Event(n).Max_NTR = NTR;
               ET_Cyclones(i).Event(n).Time_RF = Time_RF;
               ET_Cyclones(i).Event(n).POT_RF = RF;
               n=n+1;
               count_ETC= count_ETC+1;
            end
        else
           ET_Cyclones(i).Event(n).Time_NTR = Time_ntr;
           ET_Cyclones(i).Event(n).Max_NTR = NTR;
           ET_Cyclones(i).Event(n).Time_RF = Time_RF;
           ET_Cyclones(i).Event(n).POT_RF = RF;
           n=n+1;
           count_ETC = count_ETC+1;

        end
        
    end
end
%%
save('TC_events_conditioning_POT_RF_for_18hr_RF_acc',"Tropical_Cyclones")
save('ET_events_conditioning_POT_RF_for_18hr_RF_acc.mat',"ET_Cyclones")

