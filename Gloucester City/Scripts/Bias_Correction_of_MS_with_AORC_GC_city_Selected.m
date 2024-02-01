%%%%%% %%%%%%%%
clear
close all
clc


Pr_Measured = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport\Spatial Averaging of RF DATA\Phlli_airport_RF_data_1901_to_2021.mat");
Pr_Measured=Pr_Measured.Philli_airport_RF_data_1901_to_2021;

AORC=load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Spatial Averaging of RF DATA\AORC\AORC_9_grid_loc_Gloucester.mat");

Pr_AORC =mean(AORC.prcpdata(:,2:end),2);
Time_AORC=AORC.prcpdata(:,1);

q=[0:0.0002:1];
q(end)=0.99998;
%% Make the trace precipitation values 0.1 mm
% select the data after 2013
t=datenum(2013,1,1,0,0,0);
I=find(Pr_Measured(:,1)==t);

D1=Pr_Measured(1:I,:);
D2=Pr_Measured(I+1:end,:);

D2(isnan(D2(:,2)),2)=0.1;


Pr_Measured=[D1;D2];
%%
t1 = datetime(1979,02,01,00,00,00,'Format','yyyy MM dd HH mm'); %%insert the starting date manually / NOTE THAT IT STATRTS WITH 02ND OF JAN
t2 = datetime(2021,12,31,23,00,00,'Format','yyyy MM dd HH mm'); %% Enter the end date manually
 
t=t1:hours(1):t2; %create the time series
t=t';
%Dates = datevec(t); %Split the date time 


%% Finding the threshold
%finding the index of begining of AORC data
index=find(Pr_Measured(:,1)==datenum(t1));

MS_Sorted = sort(Pr_Measured(index:end,2),'ascend');
AORC_Sorted = sort(Pr_AORC,'ascend');

index2 = find(MS_Sorted > 0, 1); % find the location of the first element larger than 0
Threshold = AORC_Sorted (index2);

%% Creating Small rainfall days in to dry days

MOD_PR_AORC = Pr_AORC;
MOD_PR_AORC(MOD_PR_AORC < Threshold)= 0; %Make all the lower threshould values zero and this should be taken finally



%%
C=[];
Pr_Measured_1=Pr_Measured(:,2); % We take the entire dta set to fit the distribution

gamma_MS=fitdist(Pr_Measured_1,'gamma');
shp_par_MS = gamma_MS.a;
scl_par_MS = gamma_MS.b;

gamma_AORC=fitdist(MOD_PR_AORC,'gamma');
shp_par_AORC = gamma_AORC.a;
scl_par_AORC = gamma_AORC.b;

y = cdf('gamma',Pr_Measured_1,shp_par_MS,scl_par_MS); % Find the CDF of the data
Bias_corrected_MS(:,2) = gaminv(y,shp_par_AORC,scl_par_AORC); % Find the inverse values of the data

Bias_corrected_MS(:,1)=Pr_Measured(:,1);
%save("Bias Corrected MS at Philli airport_With_AORC_of_GC.mat","Bias_corrected_MS");

%%
% 
% Q_data_M1 =[];
% j=1;
% for q= 0:0.001:1
%  
%     MS_Quantile = quantile (Pr_Measured(:,2),q);
%     AORC_Quantile = quantile (Bias_corrected_MS(:,2),q);
%     Q_data_M1 (j,1)= q;
%     Q_data_M1 (j,2)= MS_Quantile;
%     Q_data_M1(j,3)= AORC_Quantile;
%     j=j+1;
% end
% %% Ploting the Quantiles
% % PLease select the quantile range needed to be plotted
% 
% figure
% scatter(Q_data_M1(:,2),Q_data_M1(:,3),'x','r');
% 
% xlim([0 60]); ylim([0 60]);
% xlabel('Rainfall (mm)- Measured');ylabel('Rainfall (mm) - Chelsa');%legend('Measured','Chelsa') % data 1 and data 2 respectively
% 
% hold on
% x=linspace(0,120,2);
% y=x;
% plot(x,y,'color','black');
% 
% legend('Method 1','x=y line','Location','northeastoutside');
% 
% %annotation('textbox', [0.5, 0.9, 0, 0], 'string');
% 
% hold on
%%

%% The Q_data output files should be created before running this code
fig2=figure;
hold on
Q_data =[];
j=1;

% This is for daily data comparision
% AORC_daily = sum(reshape(Pr_AORC,24,[]));
% BSMS_daily = sum(reshape(Bias_corrected_MS(index:end,2),24,[]));
% MS_daily = sum(reshape(Pr_Measured(index:end,2),24,[]));

AORC_daily = Pr_AORC;
BSMS_daily = Bias_corrected_MS(index:end,2);
MS_daily = Pr_Measured(index:end,2);





MS_Quantile = quantile (MS_daily,q);
AORC_Quantile = quantile (AORC_daily,q);
BC_Quantile = quantile(BSMS_daily,q);

Q_data (:,1)= q;
Q_data (:,3)= MS_Quantile;
Q_data (:,2)= AORC_Quantile;
Q_data (:,4)= BC_Quantile;


a=scatter(Q_data(:,3),Q_data(:,2));hold on;
b=scatter(Q_data(:,4),Q_data(:,2));

%qqplot(MS_daily,AORC_daily)

a.SizeData = 80;b.SizeData = 80;a.LineWidth=1.5; b.LineWidth=1.5;
a.MarkerEdgeColor = [0.83 0.14 0.14]; b.MarkerEdgeColor=[0 0.5 1];


xlim([0 45]); ylim([0 45]);
xlabel('Measured gauge data (mm)');ylabel('Basin average AORC (mm)') % data 1 and data 2 respectively

x=linspace(0,120,2);
y=x;
plot(x,y,'color','black');

legend('Before bias correction','After bias correction','Location','best');
set(gca,'FontSize',20,'FontName','Times New Roman');grid on
set(gca, 'ytick', 0:10:70);set(gca,'xtick',0:10:70)
ax.XTickLabelRotation= 0;
title('')

set(gcf,'units','inches','position',[0,0,8,7])
saveas(fig2,'Bias_conrrection_Gloucester_City')

