%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script checks and plots the TCs and ETCs on a given day
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% input data: Please change accordingly
Time_NTR_mod = datenum(1996,1,20,12,0,0); % change the date and the time when you want to plot cyclones
Cir_1=350; % search radius of the circles
Cir_2=0; % search radius of the circles
Cir_3=1000; % search radius of the circles
T_before_event = 2; % Time window: before the event
T_after_event = 1; % Time window: after the event
Acc = 18; % Selected RF Accumulation time
% Co-ordinates of the center of search radius (basin)
GC_LT=39.888749; % need to enter the corresponding coordinates of the center of the catchment
GC_LN=-75.097364;

%% Loading Cyclone Track data

Data_GC_Con_NTR = load("%%Path%%/ET_events_conditioning_POT_NTR.mat");
Data_AO  = load("ET_cyclone_tracls_from_1983","Data_AO_Complete");

%% Checking for Conditioning NTR 18hr Accumulation %%%%%%%%%%%%%%%
% Taking the time matrix of GC ETC data
Data_GC_2hr_Con_NTR  = Data_GC_Con_NTR.ET_Cyclones(Acc);
for i = 1:544 % Need to change according to the number of POT events in the array
    Time_NTR_Con_NTR (i,1) = Data_GC_2hr_Con_NTR.Event(i).Time_NTR;
end

%% Plotting

Data_AO_Complete=Data_AO.Data_AO_Complete;    

figure
dist = km2deg(350); % The distance threshold in km
count =0;



    strt_hr = Time_NTR_mod-T_before_event; % Assuming a 3-day window please change it accordingly
    end_hr = Time_NTR_mod+T_after_event+0.25; % Assuming a 3-day window and adjustment of 6hr resolution
    
    index_all = find(Data_AO_Complete(:,2)>=strt_hr & Data_AO_Complete(:,2)<=end_hr); % Finding all the ET indices passing during 3 day period
    selected_data = Data_AO_Complete(index_all,:);
    UID = selected_data(:,1);
    UID = unique(UID);

    for n = 1: length(UID)
        UID_one = UID(n,1);
        indinces_single_ET = find(selected_data(:,1) == UID_one);
        cord_of_ET_min = selected_data(indinces_single_ET,4:5);
        cord_of_ET_min(cord_of_ET_min>180)=cord_of_ET_min(cord_of_ET_min>180)-360;
        
        
        geoplot(cord_of_ET_min(:,2),cord_of_ET_min(:,1),'b');
        geolimits([0 90],[-180 180]);
        hold on
    end


hold on

%%  Loading TC track data
load("%%path%%/Cyclone_Track_data_from_1850.mat");

%% finding the corresponding tropical cyclone

Datenum = Cyclone_track_data(:,1);
index_TC = find(Datenum(:,1)>=strt_hr & Datenum(:,1)<=end_hr); % Finding all the TC indices passing during 3 day period

clear M H D A B C

%% Mapping The Hurricanes
dif = diff(index_TC);
breaks = find(dif(:,1)~=1);
breaks = [breaks;length(index_TC)];
st=1;
for bb = 1: length(breaks)
    endd=breaks(bb);
    Index_Tc_single = index_TC(st:endd);

    geoplot(Cylcone_Codinates_Modified(Index_Tc_single,6),Cylcone_Codinates_Modified(Index_Tc_single,7),'r');
    hold on

    st=breaks(bb)+1;
end

%% Creating Circles

radius = km2deg(Cir_1);

startx = GC_LN-radius;
endingx = GC_LN+radius;
k=1;

for x = startx-0.99:0.05:endingx+0.99
    d=1;
    y1=GC_LT;
    while d< radius
        d= distance(GC_LT,GC_LN,y1,x);
        y1=y1+0.05;  
    end

    Y1(k,1)=x;
    Y1(k,2)=y1;
    

    dd=1;
    y2=GC_LT;
    while dd< radius
        dd= distance(GC_LT,GC_LN,y2,x);
        y2=y2-0.05;  
    end

    Y2(k,1)=x;
    Y2(k,2)=y2;
    k=k+1;


end

geoplot(Y1(:,2),Y1(:,1),'r-')
hold on
geoplot(Y2(:,2),Y2(:,1),'r-')
hold on
clearvars y1 y2 radius Y1 Y2
%%%%%%%%%%%%%%%%

%%

radius = km2deg(Cir_2);

startx = GC_LN-radius;
endingx = GC_LN+radius;
k=1;

for x = startx-1.4:0.05:endingx+1.45
    d=1;
    y1=GC_LT;
    while d< radius
        d= distance(GC_LT,GC_LN,y1,x);
        y1=y1+0.05;  
    end

    Y1(k,1)=x;
    Y1(k,2)=y1;
    

    dd=1;
    y2=GC_LT;
    while dd< radius
        dd= distance(GC_LT,GC_LN,y2,x);
        y2=y2-0.05;  
    end

    Y2(k,1)=x;
    Y2(k,2)=y2;
    k=k+1;


end

geoplot(Y1(:,2),Y1(:,1),'g-')
hold on
geoplot(Y2(:,2),Y2(:,1),'g-')
hold on
clearvars Y1 Y2 radius y1 y2
%%%%%%%%%%%%%%%%
%%
radius = km2deg(Cir_3);

startx = GC_LN-radius;
endingx = GC_LN+radius;
k=1;

for x = startx-3:0.05:endingx+3
    d=1;
    y1=GC_LT;
    while d< radius
        d= distance(GC_LT,GC_LN,y1,x);
        y1=y1+0.05;  
    end

    Y1(k,1)=x;
    Y1(k,2)=y1;
    

    dd=1;
    y2=GC_LT;
    while dd< radius
        dd= distance(GC_LT,GC_LN,y2,x);
        y2=y2-0.05;  
    end

    Y2(k,1)=x;
    Y2(k,2)=y2;
    k=k+1;


end

geoplot(Y1(:,2),Y1(:,1),'b-')
hold on
geoplot(Y2(:,2),Y2(:,1),'b-')
hold on
clearvars y1 y2 radius Y1 Y2
%%%%%%%%%%%%%%%%



