
clear
clc

path= ['C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1'];
TC_ext_CON_ntr_yrs = load([path,'\Copulas\TC_events_conditioning_POT_NTR_18hr_Acc.txt']);
TC_ext_CON_RF_yrs = load([path,'\Copulas\TC_events_conditioning_POT_RF_for_18hr_RF_acc.txt']);
Non_TC_ext_CON_ntr_yrs = load([path,'\Copulas\ET_events_conditioning_POT_NTR_18hr_Acc.txt']);
Non_TC_ext_CON_RF_yrs = load([path,'\Copulas\ET_events_conditioning_POT_RF_for_18hr_RF_acc.txt']);

LineWidth=2;
Colour=[0.8350 0.0780 0.1840];
width = 5;
height = 5.5;
resolution = 3000;
x0=1;
y0=1;
%% Peaks over the threshold (POT)_ NTR_TC
th = 0.633; % The threshold used 

[x,is]= sort(TC_ext_CON_ntr_yrs(:,2));
xtime = TC_ext_CON_ntr_yrs(is,1);
x= x-th;
x(x<=0) = 0.001;

dates= datevec(TC_ext_CON_ntr_yrs(:,1));
years= unique(dates(:,1)); 
n= length(years);
        
% Pareto fitting
% To know the POT parameters;
[parm,bounds] = gpfit(x,0.05);

% check fitting
lowerBnd = min(x);
xmax= max(x)*1.8;
ygrid = linspace(lowerBnd,xmax,100);


% empirical function
[F,yi]= ecdf(x);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],...
    'color','w');
stairs(yi+th,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;
plot(ygrid,gpcdf(ygrid,parm(1),parm(2),th),'-k','LineWidth',LineWidth);
hold on;
xlabel('NTR (m)');
ylabel('Cumulative Probability');

xlim([th-0.1 3]);

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(TC_ext_CON_ntr_yrs(:,2),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;

legend('Empirical CDF','Fitted Pareto CDF','','','location','southeast');


%title('Pareto CDF of POT NTR ETCs when Conditioning NTR');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on

saveas(hh,'Pareto CDF of POT NTR_TC');
%
exportgraphics(hh,'Pareto CDF of POT NTR_TC.png','Resolution',500)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% Peaks over the threshold (POT)_ RF_TC
th = 35; % The threshold used 

[x,is]= sort(TC_ext_CON_RF_yrs(:,4));
xtime = TC_ext_CON_RF_yrs(is,3);
x= x-th;
x(x<=0) = 0.001;

hold all; 

dates= datevec(TC_ext_CON_RF_yrs(:,3));
years= unique(dates(:,1)); 
n= length(years);
        
% Pareto fitting
% To know the POT parameters;
[parm,bounds] = gpfit(x,0.05);

% check fitting
lowerBnd = min(x);
xmax= max(x)*1.8;
ygrid = linspace(lowerBnd,xmax,100);

% empirical function
[F,yi]= ecdf(x);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],...
    'color','w');


stairs(yi+th,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;
plot(ygrid,gpcdf(ygrid,parm(1),parm(2),th),'-k','LineWidth',LineWidth);
hold on;
xlabel('18 hour rainfall (mm)');
ylabel('Cumulative Probability');

xlim([th-0.1 160]);

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(TC_ext_CON_RF_yrs(:,4),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;


legend('Empirical CDF','Fitted Pareto CDF','','','location','southeast');

%title('Pareto CDF of POT NTR ETCs when Conditioning NTR');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on

saveas(hh,'Pareto CDF of POT RF_TC');
exportgraphics(hh,'Pareto CDF of POT RF_TC.png','Resolution',500)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% Peaks over the threshold (POT)_ NTR_TC
th = 0.633; % The threshold used 

[x,is]= sort(Non_TC_ext_CON_ntr_yrs(:,2));
xtime = Non_TC_ext_CON_ntr_yrs(is,1);
x= x-th;
x(x<=0) = 0.001;

hold all; 

dates= datevec(Non_TC_ext_CON_ntr_yrs(:,1));
years= unique(dates(:,1)); 
n= length(years);
        
% Pareto fitting
% To know the POT parameters;
[parm,bounds] = gpfit(x,0.05);

% check fitting
lowerBnd = min(x);
xmax= max(x)*1.8;
ygrid = linspace(lowerBnd,xmax,100);

% empirical function
[F,yi]= ecdf(x);

% Plot
hh= figure;


stairs(yi+th,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;
plot(ygrid,gpcdf(ygrid,parm(1),parm(2),th),'-k','LineWidth',LineWidth);
hold on;
xlabel('NTR (m)');set(hh,'units','centimeters','Position',[3 1.5 20 12],...
    'color','w');

ylabel('Cumulative Probability');

xlim([th-0.1 2]);

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(Non_TC_ext_CON_ntr_yrs(:,2),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;

legend('Empirical CDF','Fitted Pareto CDF','','','location','southeast');


%title('Pareto CDF of POT NTR ETCs when Conditioning NTR');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on

saveas(hh,'Pareto CDF of POT NTR_ETC');
exportgraphics(hh,'Pareto CDF of POT NTR_ETC.png','Resolution',500)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% Peaks over the threshold (POT)_ RF_TC
th = 35; % The threshold used 

[x,is]= sort(Non_TC_ext_CON_RF_yrs(:,4));
xtime = Non_TC_ext_CON_RF_yrs(is,3);
x= x-th;
x(x<=0) = 0.001;

hold all;

dates= datevec(Non_TC_ext_CON_RF_yrs(:,3));
years= unique(dates(:,1)); 
n= length(years);
        
% Pareto fitting
% To know the POT parameters;
[parm,bounds] = gpfit(x,0.05);

% check fitting
lowerBnd = min(x);
xmax= max(x)*1.8;
ygrid = linspace(lowerBnd,xmax,100);


% empirical function
[F,yi]= ecdf(x);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],...
    'color','w');


stairs(yi+th,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;
plot(ygrid,gpcdf(ygrid,parm(1),parm(2),th),'-k','LineWidth',LineWidth);
hold on;
xlabel('18 hour rainfall (mm)');
ylabel('Cumulative Probability');

xlim([th-0.1 160]);

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(Non_TC_ext_CON_RF_yrs(:,4),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;

legend('Empirical CDF','Fitted Pareto CDF','','','location','southeast');


%title('Pareto CDF of POT NTR ETCs when Conditioning NTR');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on

saveas(hh,'Pareto CDF of POT RF_ETC');
exportgraphics(hh,'Pareto CDF of POT RF_ETC.png','Resolution',500)




%% Fitting the distribution for non conditional samples
% For TC
% for max RF (Gamma dis)


Max_RF = TC_ext_CON_ntr_yrs(:,4);
Max_RF(Max_RF==0)=0.00000001;
Max_NTR_Time = TC_ext_CON_ntr_yrs(:,3);
%
% empirical function
[F,yi]= ecdf(Max_RF);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],'color','w');
stairs(yi,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;

% Gamma distribution
Max_RF = sort(Max_RF,'ascend');
pd = fitdist(Max_RF,'Gamma');
mu = pd.a;
sigma = pd.b;

x=[min(Max_RF):0.01:max(Max_RF)*1.5];
x=x';
CDF_1 = gamcdf(x,mu,sigma);

fig1 = plot(x,CDF_1,'-k','LineWidth',LineWidth);
hold on;

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(TC_ext_CON_ntr_yrs(:,4),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;


xlabel('18 hour rainfall (mm)');
ylabel('Cumulative Probability');
xlim([0 160]);
%title = title('CDF of Max RF of ETCs when Conditioning NTR');
legend('Empirical CDF','Fitted Gamma (2)','','','location','southeast');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on
saveas(hh,'Non_Extreme_CDF of Max_RF_TC');
exportgraphics(hh,'Non_Extreme_CDF of Max_RF_TC.png','Resolution',500)

%% For TC
% for max NTR (logistic)


Max_NTR = TC_ext_CON_RF_yrs(:,2);
%Max_NTR(Max_NTR==0)=0.00000001;
Max_NTR_Time = TC_ext_CON_RF_yrs(:,1);
%

% empirical function
[F,yi]= ecdf(Max_NTR);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],'color','w');
stairs(yi,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;

% logis distribution
Max_NTR = sort(Max_NTR,'ascend');
pd = fitdist(Max_NTR,'Logistic');
mu = pd.mu;
sigma = pd.sigma;

x=[min(Max_NTR):0.01:max(Max_NTR)*1.5];
x=x';
CDF_1 = cdf('Logistic',x,mu, sigma);

fig1 = plot(x,CDF_1,'-k','LineWidth',LineWidth);
hold on;

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(TC_ext_CON_RF_yrs(:,2),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;


xlabel('NTR (m)');
ylabel('Cumulative Probability');
xlim([min(Max_NTR) 2.5]);
%title = title('CDF of Max RF of ETCs when Conditioning NTR');
legend('Empirical CDF','Fitted logistic','','','location','southeast');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on
saveas(hh,'Non_Extreme_CDF of Max_NTR_TC');
exportgraphics(hh,'Non_Extreme_CDF of Max_NTR_TC.png','Resolution',500)


%%
% For ETC
% for max RF (Tweedie dis) % imported from R 0:0.1:300

TWE_cdf = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Fitting Distributios\ETC\ETC_TWEEdie_Sample.csv");

Max_RF = Non_TC_ext_CON_ntr_yrs(:,4);
Max_RF(Max_RF==0)=0.00000001;
Max_NTR_Time = Non_TC_ext_CON_ntr_yrs(:,3);
%
% empirical function
[F,yi]= ecdf(Max_RF);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],'color','w');
stairs(yi,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;

% Tweedie distribution
x=0:0.1:300;
fig1 = plot(x,TWE_cdf.x,'-k','LineWidth',LineWidth);
hold on;

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(Non_TC_ext_CON_ntr_yrs(:,4),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;


xlabel('18 hour rainfall (mm)');
ylabel('Cumulative Probability');
xlim([0 100]);
%title = title('CDF of Max RF of ETCs when Conditioning NTR');
legend('Empirical CDF','Fitted Tweedie','location','southeast');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on
saveas(hh,'Non_Extreme_CDF of Max_RF_ETC');
exportgraphics(hh,'Non_Extreme_CDF of Max_RF_ETC.png','Resolution',500)




%% for ETC
% for max NTR (logistic)


Max_NTR = Non_TC_ext_CON_RF_yrs(:,2);
Max_NTR(Max_NTR==0)=0.00000001;
Max_NTR_Time = Non_TC_ext_CON_RF_yrs(:,1);
%

% empirical function
[F,yi]= ecdf(Max_NTR);

% Plot
hh= figure;
set(hh,'units','centimeters','Position',[3 1.5 20 12],'color','w');
stairs(yi,F,'.','MarkerSize',18,'LineWidth',1,'Color',Colour);
hold on;

% logistic distribution
Max_NTR = sort(Max_NTR,'ascend');
pd = fitdist(Max_NTR,'Logistic');
mu = pd.mu;
sigma = pd.sigma;

x=[min(Max_NTR):0.01:max(Max_NTR)*1.5];
x=x';
CDF_1 = cdf('Logistic',x,mu, sigma);

fig1 = plot(x,CDF_1,'-k','LineWidth',LineWidth);
hold on;

% Here 0.95 is the alpha for the function
[CI_U,CI_L,Xi]= DKW_conf_int(Non_TC_ext_CON_RF_yrs(:,2),0.95);

plot(Xi,CI_U,'-k','LineWidth',1.5,'LineStyle',':');
hold on;
plot(Xi,CI_L,'-k','LineWidth',1.5,'LineStyle',':');
hold on;


xlabel('NTR (m)');
ylabel('Cumulative Probability');
xlim([min(Max_NTR) 2]);
%title = title('CDF of Max RF of ETCs when Conditioning NTR');
legend('Empirical CDF','Fitted logistic','','','location','southeast');
set(gca,'FontSize',18,'FontName','Times New Roman');set(gcf,'units','inches','position',[x0,y0,width,height])
grid on
saveas(hh,'Non_Extreme_CDF of Max_NTR_ETC');
exportgraphics(hh,'Non_Extreme_CDF of Max_NTR_ETC.png','Resolution',500)



% 
% 
% 
% %% Exponential Ditribution
% 
% clear mu sigma pd parm 
% 
% pd = fitdist(Max_RF,'Exponential');
% mu = pd.mu;
% 
% CDF_3 = expcdf(x,mu);
% plot(x,CDF_3,'color',[0.9290 0.6940 0.1250],'LineWidth',1);
% 
% %% Lognormal
% Max_RF(Max_RF==0)=[];
% 
% clear mu sigma pd parm 
% 
% pd = fitdist(Max_RF,'Lognormal');
% mu = pd.mu;
% sigma = pd.sigma;
% 
% CDF_2 = logncdf(x,mu,sigma);
% plot(x,CDF_2,'color',[0.4660 0.6740 0.1880],'LineWidth',1);
% 
% 
% 
% %% Weibull Distribution
% 
% clear mu sigma pd parm 
% 
% pd = fitdist(Max_RF,'Weibull');
% mu = pd.A;
% sigma = pd.B;
% 
% CDF_4 = wblcdf(x,mu);
% plot(x,CDF_4,'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
% 
% 
% %% BirnbaumSaunders Distribution
% 
% clear mu sigma pd parm 
% 
% pd = fitdist(Max_RF,'BirnbaumSaunders');
% mu = pd.beta;
% sigma = pd.gamma;
% 
% CDF_5 = cdf('BirnbaumSaunders',x,mu,sigma);
% plot(x,CDF_5,'Color',[0 0.4470 0.7410],'LineWidth',1);
% legend('Empirical CDF','Gamma','Exponential','Lognormal','Weibull','Birnbaum-Saunder','location','Best');
% 
% 
% saveas(hh,'CDF of Max RF of ETCs when Conditioning NTR')
% saveas(hh,'CDF of Max RF of ETCs when Conditioning NTR.png')
% 
% %% 
% % 
% % Save figures
% 
% % FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% % for iFig = 1:length(FigList)
% %     FigHandle = FigList(iFig);
% %     FigName   = ['Fig ' num2str(get(FigHandle, 'Number')) '.png'];
% %     saveas(FigHandle, fullfile(FSave, FigName));
% % end

