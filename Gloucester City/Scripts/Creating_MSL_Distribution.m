%% Fitting distribution to the WL moving average for the last 5 years 
clear 
clc
close all
NTR=load("C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Water Levels\POT_NTR_and_NTR_timeseries_for_Philli_airport.mat");
Time = NTR.POT_NTR_and_NTR_timeseries.Time(end-43800:end);
Mov_avg = NTR.POT_NTR_and_NTR_timeseries.MovA_WL(end-43800:end);
Time(isnan(Mov_avg)==1, :) = [];
Mov_avg(isnan(Mov_avg)==1, :) = [];

x_values = min(Mov_avg):0.01:max(Mov_avg);

histogram(Mov_avg,x_values)

%% Checking monthly MSL variation
monthOrder = {'January','February','March','April','May','June','July', ...
    'August','September','October','November','December'};
tt = datevec(Time);
app_length = 31*24*5;
MSL = nan(app_length,12);
intervalWidth = 0.025;
X=-0.2:intervalWidth:0.55;
t = tiledlayout(6,2);
Monthly_MSL=struct;
for m=1:12
    monht_ind = find(tt(:,2) == m);
    Monthly_data = Mov_avg(monht_ind);
    %ploting an example curve for month
    MSL(1:length(Monthly_data),m)=Monthly_data;
    
    Monthly_data=sort(Monthly_data,'ascend');
    %Storing data to a structure file
    Monthly_MSL(m).Month=Monthly_data;

    ncount=histc(Monthly_data,X); 
    relativefreq = ncount./length(Monthly_data(~isnan(Monthly_data)))*100;
    nexttile
    bar(X-intervalWidth/2, relativefreq,1);title(monthOrder{m});ylim([0 45]);

   
end
t.TileSpacing = 'compact';
title(t,'MSL distribution for 2017 to 2021 ')
xlabel(t,'30 days averaged MSL(m)')
ylabel(t,'Relative frequncy (%)')
saveas(t,"MSL Ditribution.fig")


% ploting a box whishker plot
figure
b1=boxchart(MSL); xlabel('Month');ylabel('30 days averaged MSL (m)');title('MSL variability for 2017 to 2021 ');%ylim([0.65 1.3]);
saveas(b1,"Box_plot_of_high_tide.fig")

save("Monthly_MSL_Distributions",'Monthly_MSL');
