%% Load The data
clear
clc
close all

%th = 0.7232; %3 events per year
%th = 0.672; % 4 events per year
th = 0.634; % 5 events per year
%th = 0.575;% 6 events per year
%th = 0.57;% 7 events per year

%th = 0.61; % same threshold

load('C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa\Water Levels\POT_NTR_and_NTR_timeseries_for_Philli_airport.mat');


Time=POT_NTR_and_NTR_timeseries.Time;
%WL_Detrended=POT_NTR_and_NTR_timeseries.WL_Det;
NTR=POT_NTR_and_NTR_timeseries.NTR;

%  Events over threshold
% decluster time

ts=[Time NTR];

%q= 0.995615%[.98 .99 .9975];
%th= quantile(ts(:,2),q);


% events over threshold
EOT= ts(:,2);
EOT(EOT< th)= NaN;
EOT= [ts(:,1) EOT];
EOT(isnan(EOT(:,2)),:)= [];

% Plot the NTR and values over Threshold

hh= figure;
set(hh,'units','centimeters','Position',[0.2 1.5 28 12],...
    'color','w');
hold all
% h1= plot(ts(:,1),ts(:,2),'.-','LineStyle', 'none');
h1= plot(ts(:,1),ts(:,2),'.-');

ylabel('NTR (m)');
set(gca,'FontSize',11,'FontName','Times New Roman');
grid minor
ax = gca;
ax.XTickLabelRotation= 0;

title('Philadelphia Tide Gauge')

axs= axis;
hold all; 
% plot(EOT(:,1),EOT(:,2),',-k','LineStyle', 'none')
plot(EOT(:,1),EOT(:,2),'.k')

% Decluster time
dec_tim = 2.5; % in days
xnan= sum(isnan(EOT(:,2)));
POT= nan(length(EOT),2);

while xnan< size(EOT,1)
    
    for i= 1: size(POT,2)
        
        [~,fmax]= max(EOT(:,2));
        
        POT(fmax,:)= EOT(fmax,:);
        hold all; plot(EOT(fmax,1),EOT(fmax,2),'LineWidth',1);
        
        dec_wind= find(EOT(:,1)>= EOT(fmax,1)-dec_tim & EOT(:,1)<= EOT(fmax,1)+dec_tim);
        
        hold all; hi= plot(EOT(dec_wind,1),EOT(dec_wind,2),...
            '.-','Color','r');
        
        EOT(dec_wind,2)= nan;
        
        xnan= sum(isnan(EOT(:,2)));
    end
    
end

POT(isnan(POT(:,2)),:)= [];
hold all
plot(POT(:,1),POT(:,2),'og')

%% Plot
hh= figure;
set(hh,'units','centimeters','Position',[0.2 1.5 28 12],...
    'color','w');
hold all
plot(POT(:,1),POT(:,2),'.-','LineStyle', 'none')

ylabel('Sea level (m)');
set(gca,'FontSize',11,'FontName','Times New Roman')
grid minor
ax = gca;
ax.XTickLabelRotation= 0;

title('Peaks over threshold (POT) Phladelphia')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Peaks over the threshold (POT)

[x,is]= sort(POT(:,2));
xtime = POT(is,1);

x= x -th;
x(x<=0) = 0.001;

hold all; plot(xtime,x,'.');

dates= datevec(ts(:,1));
uyears= unique(dates(:,1)); 
n= length(uyears);
        
%% Pareto fitting
% To know the POT parameters;
[parm,bounds] = gpfit(x,0.05);

%% check fitting
lowerBnd = min(x);
xmax= max(x)*1.1;
ygrid = linspace(lowerBnd,xmax,100);

% empirical function
[F,yi]= ecdf(x);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],...
    'color','w');

plot(ygrid,gpcdf(ygrid,parm(1),parm(2),th),'-k','LineWidth',2);
hold on;
stairs(yi+th,F,'x','LineWidth',2,'Color',[0.6953 0.1328 0.1328]);
hold off;
xlabel('Sea level (m)');
ylabel('Cumulative Probability');
legend('Fitted Pareto CDF','Empirical CDF','location','Best');
xlim([lowerBnd xmax]);

title('Pareto CDF Philadelphia');
set(gca,'FontSize',12,'FontName','Times New Roman')
grid minor

%% Return level

% return period
% RP= [2 5 10 25 50 100 250 1000];
RP= (2:1:5000);

% vector of non-excedance probs.
mu = n/length(x);
Prob= 1-mu./RP;

% Return levels
RL= gpinv(Prob,parm(1),parm(2),th);
RL_bel= gpinv(Prob,bounds(1,1),bounds(1,2),th);
RL_abo= gpinv(Prob,bounds(2,1),bounds(2,2),th);

% Empirical return period
rank = fliplr(1:1:length(x))';
N = length(x);

GRINGORTEN = (n + 0.12)./(rank - 0.44);
Weibull    = (n+1)./rank;

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],'color','w');

h1= plot(RP,RL,'k','LineWidth',2); hold on;
plot(RP,RL_bel,'--k','LineWidth',1,'Color',h1.Color)
plot(RP,RL_abo,'--k','LineWidth',1,'Color',h1.Color)

hold all
h2= plot(Weibull,x+th,'x','LineWidth',2,'Color',[0.6953 0.1328 0.1328]);
h3= plot(GRINGORTEN,x+th,'o','LineWidth',2,'Color',[0.6953 0.1328 0.1328]);

ylabel('Return level (m)'); xlabel('Return period (years)')

set(gca,'FontSize',11,'FontName','Times New Roman')
grid minor

legend([h1,h2,h3],'Fitted Pareto','Weibull','Gringorten','location','Best');
set(gca, 'XScale', 'log','XTick',[2 5 10 25 50 100 250 1000])


axs= axis;
xlim([1 axs(2)])

%%


POT_NTR_and_NTR_timeseries.NTR=NTR;
POT_NTR_and_NTR_timeseries.POT = POT;
%POT_NTR_and_NTR_timeseries.WL_raw = Data(:,2);
POT_NTR_and_NTR_timeseries.Threshold = th;
save('POT_NTR_and_NTR_timeseries.mat','POT_NTR_and_NTR_timeseries')

%%


load('POT_NTR_and_NTR_timeseries.mat');

load('C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa\Spatial Averaging of RF DATA\Hourly_accumulation_Bias_Corrected_RF_data_Philli_airport.mat');
%load('C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport\Spatial Averaging of RF DATA\Hourly_accumulation_Bias_Corrected_RF_data_Philli_airport.mat');

POT_NTR = POT_NTR_and_NTR_timeseries.POT;

%%
Maximum_RF_events = struct;
for j=1:48
    for i = 1:length(POT_NTR)
        Time_NTR = POT_NTR(i,1);
        [index_1] = find (Data(:,1)== Time_NTR);
        [Max_RF, index_2] = max(Data(index_1-35:index_1+36,j+1)); % taking the 3 day window
        index_2 = index_2+index_1-36; % correcting the index to to find the row in 'Data'
    

        Maximum_RF_events(j).Accumulation(i).Time_NTR = Time_NTR;
        Maximum_RF_events(j).Accumulation(i).POT_NTR = POT_NTR(i,2);
              
        Maximum_RF_events(j).Accumulation(i).Time_RF = Data(index_2,1);
        Maximum_RF_events(j).Accumulation(i).Max_Rainfall = Max_RF;

    end
end
%%
save('Maximum_RF_events_for_each_POT_NTR',"Maximum_RF_events")

%%

load('Maximum_RF_events_for_each_POT_NTR.mat');
load('C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport\Cyclone_Track\Cyclone_Track_data_from_1850_new.mat');


%% 
count_TC= zeros(48,1);
count_ETC =zeros(48,1);
%distance = km2deg(350); % in degrees 350 km 
distance = km2deg(350);
Tropical_Cyclones = struct;
ET_Cyclones = struct;

%% One event will be manually transfered in to the List date: 1995/10/06

X=Maximum_RF_events(1).Accumulation;
X=struct2table(X);
X=X.POT_NTR;
leng=length(X);
%%
for i = 1:48

    m=1;n=1;
    for j =1:leng % This should be changed according the number of Extreme events you have
        Time_ntr = Maximum_RF_events(i).Accumulation(j).Time_NTR;
        NTR = Maximum_RF_events(i).Accumulation(j).POT_NTR;
        Time_RF = Maximum_RF_events(i).Accumulation(j).Time_RF;
        RF = Maximum_RF_events(i).Accumulation(j).Max_Rainfall;
    
    
        T_elements = datevec(Time_ntr); % Resolutions of data are different 6 hr and 1 hr
        hr = T_elements(1,4);
        hr_new = floor(hr/6)*6; 
        T_elements(1,4)= hr_new;
        Time_NTR_mod = datenum(T_elements);
        strt_hr = Time_NTR_mod-2;
        end_hr = Time_NTR_mod+1.75;
        
        index = find(Cyclone_track_data(:,1)>=strt_hr & Cyclone_track_data(:,1)<=end_hr);
    
        if length(index)>0
            if Time_ntr>728938.316666667 && Time_ntr<728939.416666667
                j
    
               Tropical_Cyclones(i).Event(m).Time_NTR = Time_ntr;
               Tropical_Cyclones(i).Event(m).NTR = NTR;
               Tropical_Cyclones(i).Event(m).Time_RF = Time_RF;
               Tropical_Cyclones(i).Event(m).RF = RF;
               m=m+1;
               count_TC(i,1) = count_TC(i,1)+1;
            
            elseif min(Cyclone_track_data(index,5))<= distance
            %if min(Cyclone_track_data(index,5))<= distance
            %if min(Cyclone_track_data(index,5))> distance && min(Cyclone_track_data(index,5))< distance2 
               Tropical_Cyclones(i).Event(m).Time_NTR = Time_ntr;
               Tropical_Cyclones(i).Event(m).NTR = NTR;
               Tropical_Cyclones(i).Event(m).Time_RF = Time_RF;
               Tropical_Cyclones(i).Event(m).RF = RF;
               m=m+1;
               count_TC(i,1) = count_TC(i,1)+1;
            else
               ET_Cyclones(i).Event(n).Time_NTR = Time_ntr;
               ET_Cyclones(i).Event(n).NTR = NTR;
               ET_Cyclones(i).Event(n).Time_RF = Time_RF;
               ET_Cyclones(i).Event(n).RF = RF;
               n=n+1;
               count_ETC(i,1) = count_ETC(i,1)+1;
            end
        else
           ET_Cyclones(i).Event(n).Time_NTR = Time_ntr;
           ET_Cyclones(i).Event(n).NTR = NTR;
           ET_Cyclones(i).Event(n).Time_RF = Time_RF;
           ET_Cyclones(i).Event(n).RF = RF;
           n=n+1;
           count_ETC(i,1) = count_ETC(i,1)+1;

        end
        
    end
end
%%
save('TC_events_conditioning_POT_NTR_350',"Tropical_Cyclones")
save('ET_events_conditioning_POT_NTR_350',"ET_Cyclones")

%%
load('Maximum_RF_events_for_each_POT_NTR.mat');
load('ET_events_conditioning_POT_NTR_350.mat');
load('TC_events_conditioning_POT_NTR_350.mat');

%% Calculation for No_Condition taking all the events
Kendalls=[];

%%
TC=Tropical_Cyclones(1).Event;
TC=struct2table(TC);
TC=table2array(TC);
len_tc=length(TC);

ETC=ET_Cyclones(1).Event;
ETC=struct2table(ETC);
%ETC([1,2],:)=[]; % this is done due to unavailablity of data for 1979 Jan
ETC=table2array(ETC);

len_et=length(ETC);

EV=Maximum_RF_events(1).Accumulation;
EV=struct2table(EV);
EV=table2array(EV);
len_EV=length(EV);



%%
for i = 1:48
    x=[];y=[];
    for j = 1:len_EV
        NTR = Maximum_RF_events(i).Accumulation(j).POT_NTR;
        x(j,1)=NTR;
        RF = Maximum_RF_events(i).Accumulation(j).Max_Rainfall;
        y(j,1)=RF;
    end
    X=isnan(x);
    y(X)=[];
    x(X)=[];
    Y=isnan(y);
    x(Y)=[];
    y(Y)=[];
  
    Kendalls(i,2) = corr(x,y,'type','Kendall');
    Kendalls(i,1)=i;
end

%% Calculation for TC events

for i = 1:48
    x=[];y=[];
    for j = 1:len_tc % Only 37 events% Should be manually entered
        NTR = Tropical_Cyclones(i).Event(j).NTR;
        x(j,1)=NTR;
        RF = Tropical_Cyclones(i).Event(j).RF;
        y(j,1)=RF;
    end
    i
    X=isnan(x);
    y(X)=[];
    x(X)=[];
    Y=isnan(y);
    x(Y)=[];
    y(Y)=[];
    Kendalls(i,3) = corr(x,y,'type','Kendall');
    
end



%% Calculation for ETC events

for i = 1:48
    x=[];y=[];
    for j = 1:len_et % Only 563 events in this scenario
        NTR = ET_Cyclones(i).Event(j).NTR;
        x(j,1)=NTR;
        RF = ET_Cyclones(i).Event(j).RF;
        y(j,1)=RF;
    end
    i
    X=isnan(x);
    y(X)=[];
    x(X)=[];
    Y=isnan(y);
    x(Y)=[];
    y(Y)=[];
    Kendalls(i,4) = corr(x,y,'type','Kendall');
end


%%
save('Kendals_Toe_for_All_TC_ETC_respectively_350',"Kendalls")

%%


hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],'color','w');

h1= plot(Kendalls(:,1),Kendalls(:,2),'LineWidth',1,'Marker','.','MarkerSize',10); hold on;
plot(Kendalls(:,1),Kendalls(:,3),'LineWidth',1,'Color','r','Marker','.','MarkerSize',10)
plot(Kendalls(:,1),Kendalls(:,4),'LineWidth',1,'Color','g','Marker','.','MarkerSize',10)
ylabel("Kendall's T "); xlabel('Rainfall Accumulation (hours)')
set(gca,'FontSize',14,'FontName','Times New Roman','YLim',[-0.1 0.5])
% = [0 1];
grid("minor");
title(['Variation of the Dependancy of Max. RF & POT NTR with the Rainfall duration (Con. NTR)'])


hold all
legend('All','TC','Non TC','Location','northwest');

clearvars x y X Y
%% Calculating the correlation for a moving window 
% 
Acc = 18; % Accumulation hour of interest

TC=Tropical_Cyclones(Acc).Event;
TC=struct2table(TC);
TC=table2array(TC);
len_tc=length(TC);

ETC=ET_Cyclones(Acc).Event;
ETC=struct2table(ETC);
%ETC([1,2],:)=[]; % this is done due to unavailablity of data for 1979 Jan
ETC=table2array(ETC);

len_et=length(ETC);

EV=Maximum_RF_events(Acc).Accumulation;
EV=struct2table(EV);
EV=table2array(EV);
len_EV=length(EV);


%
Acc_TC_data=struct2table(Tropical_Cyclones(Acc).Event);Acc_TC_data=table2array(Acc_TC_data);
Acc_ETC_data=struct2table(ET_Cyclones(Acc).Event);Acc_ETC_data=table2array(Acc_ETC_data);
Acc_All_data=struct2table(Maximum_RF_events(Acc).Accumulation);Acc_All_data=table2array(Acc_All_data);

%
Kendalls=[];
p=[];
i=1;
for yr=1901:1991
    Kendalls(i,1)=yr+15;
    p(i,1)=i;

    t1=datenum(yr,1,1,0,0,0);
    t2=datenum(yr+30,1,1,0,0,0);
    ind_tc=find (Acc_TC_data(:,1)>t1 & Acc_TC_data(:,1)<t2);
    no_tc=length(ind_tc);
    
    ind_etc=find (Acc_ETC_data(:,1)>t1 & Acc_ETC_data(:,1)<t2);
    
    ind_all=find (Acc_All_data(:,1)>t1 & Acc_All_data(:,1)<t2);
    
    x=Acc_TC_data(ind_tc,2);
    y=Acc_TC_data(ind_tc,4);
    [Kendalls(i,2), p(i,2)] = corr(x,y,'type','Kendall');
    %clearvars x y

    x=Acc_ETC_data(ind_etc,2);
    y=Acc_ETC_data(ind_etc,4);
    [Kendalls(i,3), p(i,3)] = corr(x,y,'type','Kendall');
    %clearvars x y

    x=Acc_All_data(ind_all,2);
    y=Acc_All_data(ind_all,4);
    [Kendalls(i,4), p(i,4)] = corr(x,y,'type','Kendall');
    %clearvars x y
    i=i+1;
end

% Sampling with natural variability
A=1901:1:2021;
Years=datevec(Acc_All_data(:,1));


for i=1:10000
    B = randsample(A,30);
    data=[];
    for j=1:length(B)
        ind=find(Years(:,1)==B(1,j));
        data = [data; Acc_All_data(ind,:)];
    end
    Re_sample(i,1)= corr(data(:,2),data(:,4),'type','Kendall');
end

Re_sample=sort(Re_sample,"ascend");
Low_b=Re_sample(floor(length(Re_sample)*0.1),1);
up_b=Re_sample(floor(length(Re_sample)*0.9),1);

%
figure
for i=1:length(p(:,4))
    if p(i,4)<= 0.05
        h1=plot(Kendalls(i,1),Kendalls(i,4),LineStyle="none",Marker="o",Color='black',MarkerFaceColor="r");
        hold on
    else
        h1=plot(Kendalls(i,1),Kendalls(i,4),LineStyle="none",Marker="o",Color='black');
        hold on
    end

end

xlim([1900 2021]);
ylim([0 0.5]);

p1 = patch([1900 1900 2021 2021],[up_b Low_b Low_b up_b ],'black','FaceAlpha',0.2);
p1.EdgeColor = "k";
p1.LineStyle = "--";

ylabel("Kendall's T "); xlabel('30 yr moving window from 1901')
set(gca,'FontSize',14,'FontName','Times New Roman','YLim',[-0.1 0.5])


grid("minor");
title('Changes in the rank correlation between 18hr Rf Acc & NTR extremes');

saveas(h1,"Changes in the rank correlation between 18hr Rf Acc & NTR extremes")

%%
