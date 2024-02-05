%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script combines the AEP from two popultions and plot them
%   
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%               The word "ETC" referes to non-TC events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
clc
RP ={'5','10','20','50','100'}; % The return periods of interest
rp=[5 10 20 50 100]; % Array of RPs of interest


% reading the RP array if TC population 
TC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\RP_TC.csv");
%
% reading the RP array if Non-TC populatio
ETC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\RP_ETC.csv");

% Loading the POT extreme events
TC_ext_CON_ntr_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_events_conditioning_POT_NTR_18hr_Acc.txt");
TC_ext_CON_RF_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_events_conditioning_POT_RF_for_18hr_RF_acc.txt");
Non_TC_ext_CON_ntr_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ET_events_conditioning_POT_NTR_18hr_Acc.txt");
Non_TC_ext_CON_RF_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ET_events_conditioning_POT_RF_for_18hr_RF_acc.txt");

% Loadnig the hourly time series of data
Data_1901 = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\NTR_and_18hr_RF_acc_data.mat");
Data_1901.D=Data_1901.C;

% Loading the simulated events from copulas
COP_sample_TC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_Cop_Sample.csv");
COP_sample_ETC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ETC_Cop_Sample.csv");

% List of names of events hit Gloucester City
event = {'Hazel','Unnamed','Sandy','Connie','Isabel','','Irene','','David','Floyd','Ida','','Ernesto','','Agnes','Able','','Groundhog',     'Danielle','Ivan','Bertha','Isaias','Jeanne'};
date ={'10/15/1954','08/23/1933','10/30/2012','08/12/1955','09/19/2003','', '08/29/2011', '','09/06/1979','09/17/1999','09/02/2021','','09/03/2006','','06/23/1972','09/01/1952','','02/02/1952',    '09/27/1992','09/21/2004','07/14/1996','08/15/2020','09/30/2004'};

% Following variables are created
c_map = [ linspace(0.8,1,512)', linspace(0,0.9,512)', linspace(0,0.2,512)']; % Colour Map
F_Size = 15; % font size of the figures
Con_NTR_LW = 1; % Lower limit of the NTR axis
Con_NTR_M_Size = 6; % maximum of the NTR axis
Con_RF_LW = 1;% Lowe limit of the RF axis
Con_RF_M_Size = 4;% maximum of the RF axis
width = 20; % width of figrues
height = 20; % height of figures
resolution = 1000; % Number of grid cells in each direction of discretizing
x0=0.2;
y0=0.2;

x_U_lim = 8; % limits of discretizing
x_L_lim = -1.2; % limits of discretizing
Y_U_ilm = 300;% limits of discretizing
Y_L_lim = 0;% limits of discretizing


%% finding isolines for conditioning parameter one


xv = linspace(x_L_lim, x_U_lim, resolution);
yv = linspace(Y_L_lim, Y_U_ilm, resolution);
[Xm,Ym] = ndgrid(xv, yv);

Zm1 = griddata(TC.Var1, TC.Var2, TC.V3, Xm, Ym);

for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm1,[rp(i) rp(i)]); hold on;
    cord_con1(i).X = M(1,2:end);
    cord_con1(i).Y = M(2,2:end);
end


%% finding isolines for conditioning parameter two

Zm2 = griddata(TC.Var1, TC.Var2, TC.V4, Xm, Ym);


figure
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm2,[rp(i) rp(i)]); hold on;
    cord_con2(i).X = M(1,2:end);
    cord_con2(i).Y = M(2,2:end);
end



%% Plotting together
figure
for i=1:length(rp)
    plot(cord_con1(i).X,cord_con1(i).Y,'r','LineWidth',3); hold on;
    plot(cord_con2(i).X,cord_con2(i).Y,'b','LineWidth',3); 
end
xlabel('NTR (m)'); ylabel('18 hr Rainfall (mm)');


%% 
% Combine the two conditional data sets based on the highest probability

An_ex_con1 = 1./TC.V3;
An_ex_con2 = 1./TC.V4;

% Checking the maximum probability and plot them
for i=1:length(An_ex_con2)
    if An_ex_con1(i)>An_ex_con2(i)
        An_ex_comb(i)=An_ex_con1(i);
    else
        An_ex_comb(i)=An_ex_con2(i);
    end
end

RP_Comb1 = 1./An_ex_comb;
Zm3 = griddata(TC.Var1, TC.Var2, RP_Comb1, Xm, Ym);

cord_con_comb_TC=[];
figure
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm3,[rp(i) rp(i)]); hold on;
    cord_con_comb_TC(i).X = M(1,2:end);
    cord_con_comb_TC(i).Y = M(2,2:end);
end
%% Plotting the isolines for TC sample

figure
for i=1:length(rp)

    prd=[];
    prd = ksdensity([COP_sample_TC.Var1,COP_sample_TC.Var2],[cord_con_comb_TC(i).X',cord_con_comb_TC(i).Y']); % Calculating the prbability density
    
    pr=prd/max(prd); % normalize the relative probability
    fig2=scatter(cord_con_comb_TC(i).X',cord_con_comb_TC(i).Y',[],pr,'filled');hold on;

    % Plotting most lkely event
    [~,ind_MLI] = max(pr);
    MLE = scatter(cord_con_comb_TC(i).X(ind_MLI)',cord_con_comb_TC(i).Y(ind_MLI)','k','filled','^');hold on;

    % Plotting the RP labels
    str = ['RP = ',RP(i)];
    str=strcat(str(1,1),str(1,2));
    m=max(cord_con_comb_TC(i).Y);
    text(0,m+5,str,"FontSize",F_Size,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');

end


%Plotting Extreme event sample: TC
ntr=plot(TC_ext_CON_ntr_yrs(:,2),TC_ext_CON_ntr_yrs(:,4),Color='r',Marker='o',MarkerSize=Con_NTR_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_NTR_LW);hold on;
rf=plot(TC_ext_CON_RF_yrs(:,2),TC_ext_CON_RF_yrs(:,4),Color='b',Marker='^',MarkerSize=Con_RF_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_RF_LW);
hh.Position=[0 0 width height];

%for TC events
loc1 = sortrows(TC_ext_CON_ntr_yrs,2,"descend");


xlabel('NTR (m)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
ylabel('18 Hour Rainfall (mm)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
c = colorbar;c.Label.String = "Relative Probability";c.Label.Rotation = 270;c.Label.VerticalAlignment = "bottom"; c.Label.FontSize = F_Size;
colormap(flipud(c_map))
ax = gca;
ax.XTickLabelRotation= 0;
title('')
hold on
ylim([0 175]);xlim([0 3]);

legend([fig2,MLE,ntr,rf],'JP Isoline','Most likely event','Sample Con NTR','Sample Con RF');
set(gcf,'units','inches','position',[x0,y0,width,height])
saveas(fig2,'Combined_TC_with_Simulated_kde')


%% Calculations for Non_TC populaition %%

% finding isolines for conditioning parameter one

Zm1 = griddata(ETC.Var1, ETC.Var2, ETC.V3, Xm, Ym);
figure
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm1,[rp(i) rp(i)]); hold on;
    cord_con1(i).X = M(1,2:end);
    cord_con1(i).Y = M(2,2:end);
end


%% finding isolines for conditioning parameter two

Zm2 = griddata(ETC.Var1, ETC.Var2, ETC.V4, Xm, Ym);
figure
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm2,[rp(i) rp(i)]); hold on;
    cord_con2(i).X = M(1,2:end);
    cord_con2(i).Y = M(2,2:end);
end



%% Plotting together
figure
for i=1:length(rp)
    plot(cord_con1(i).X,cord_con1(i).Y,'r','LineWidth',3); hold on;
    plot(cord_con2(i).X,cord_con2(i).Y,'b','LineWidth',3); hold on;xlabel('NTR (m)'); ylabel('18 hr Rainfall (mm)');
end

%%
% Combine the two conditional data sets based on the highest probability

An_ex_con1 = 1./ETC.V3;
An_ex_con2 = 1./ETC.V4;

% Checking the maximum probability and plot them
for i=1:length(An_ex_con2)
    if An_ex_con1(i)>An_ex_con2(i)
        An_ex_comb(i)=An_ex_con1(i);
    else
        An_ex_comb(i)=An_ex_con2(i);
    end
end

RP_Comb2 = 1./An_ex_comb;
Zm3 = griddata(ETC.Var1, ETC.Var2, RP_Comb2, Xm, Ym);

% Ploting
figure
cord_con_comb_ETC=[];
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm3,[rp(i) rp(i)]); hold on;
    cord_con_comb_ETC(i).X = M(1,2:end);
    cord_con_comb_ETC(i).Y = M(2,2:end);
end
%%

figure
for i=1:length(rp)
    %plot(cord_con_comb(i).X,cord_con_comb(i).Y,'r','LineWidth',3); hold on;xlabel('NTR (m)'); ylabel('18 hr Rainfall (mm)');
    prd = ksdensity([COP_sample_ETC.Var1,COP_sample_ETC.Var2],[cord_con_comb_ETC(i).X',cord_con_comb_ETC(i).Y']); % Calculating the prbability density
    
    pr=prd/max(prd); % normalize the relative probability
    fig3=scatter(cord_con_comb_ETC(i).X',cord_con_comb_ETC(i).Y',[],pr,'filled');hold on;

    % Plotting most lkely event
    [~,ind_MLI] = max(pr);
    MLE = scatter(cord_con_comb_ETC(i).X(ind_MLI)',cord_con_comb_ETC(i).Y(ind_MLI)','k','filled','^');hold on;


    str = ['RP = ',RP(i)];
    str=strcat(str(1,1),str(1,2));
    m=max(cord_con_comb_ETC(i).Y);
    text(0,m+5,str,"FontSize",F_Size,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');

end
xlabel('NTR (m)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
ylabel('18 Hour Rainfall (mm)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
c = colorbar;c.Label.String = "Relative Probability";c.Label.Rotation = 270;c.Label.VerticalAlignment = "bottom"; c.Label.FontSize = F_Size;
colormap(flipud(c_map))
ax = gca;
ax.XTickLabelRotation= 0;
title('')
hold on
ylim([0 175]);xlim([0 3]);



%Plotting Extreme event sample: TC
ntr=plot(Non_TC_ext_CON_ntr_yrs(:,2),Non_TC_ext_CON_ntr_yrs(:,4),Color='r',Marker='o',MarkerSize=Con_NTR_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_NTR_LW);hold on;
rf=plot(Non_TC_ext_CON_RF_yrs(:,2),Non_TC_ext_CON_RF_yrs(:,4),Color='b',Marker='^',MarkerSize=Con_RF_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_RF_LW);


hh.Position=[0 0 width height];
legend([fig3,MLE,ntr,rf],'JP Isoline','Most likely event','Sample Con NTR','Sample Con RF');
set(gcf,'units','inches','position',[x0,y0,width,height])
saveas(fig3,'Combined_Non_TC_with_Simulated_kde')

%% Combining two populations

An_ex_con1_TC = 1./TC.V3;
An_ex_con2_TC = 1./TC.V4;

An_ex_con1_ETC = 1./ETC.V3;
An_ex_con2_ETC = 1./ETC.V4;

% Checking the maximum probability and add them
for i=1:length(An_ex_con1_TC)
    if An_ex_con1_TC(i)>An_ex_con2_TC(i)
        An_ex_comb_TC(i,1)=An_ex_con1_TC(i);
    else
        An_ex_comb_TC(i,1)=An_ex_con2_TC(i);
    end
    
    if An_ex_con1_ETC(i)>An_ex_con2_ETC(i)
        An_ex_comb_ETC(i,1)=An_ex_con1_ETC(i);
    else
        An_ex_comb_ETC(i,1)=An_ex_con2_ETC(i);
    end

end



% combinnig two populations
An_ex_comb_ppl =1-(1-An_ex_comb_TC).*(1-An_ex_comb_ETC);
RP_An_ex_comb_ppl = 1./An_ex_comb_ppl;

rp=[5 10 20 50 100];

xv = linspace(x_L_lim, x_U_lim, resolution);
yv = linspace(Y_L_lim, Y_U_ilm, resolution);
[Xm,Ym] = ndgrid(xv, yv);

Zm4 = griddata(ETC.Var1, ETC.Var2, RP_An_ex_comb_ppl, Xm, Ym);

figure
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm4,[rp(i) rp(i)]); hold on;
    cord_con_comb(i).X = M(1,2:end);
    cord_con_comb(i).Y = M(2,2:end);
end

RP_data(:,1:2)=[TC.Var1, TC.Var2];
RP_data(:,3)=RP_Comb1';
RP_data(:,4)=RP_Comb2';
RP_data(:,5)=RP_An_ex_comb_ppl;

save("Return_period_data","RP_data");


%% Selecting a simulated sample from copulas

% definnign the copula sample based on the number of events
count_select_TC = round(length(TC_ext_CON_RF_yrs)/(length(TC_ext_CON_RF_yrs)+length(Non_TC_ext_CON_ntr_yrs))*length(table2array(COP_sample_TC))); % How many of events are selected from TC simulations
count_select_ETC = length(table2array(COP_sample_TC))-count_select_TC; % How many of events are selected from ETC simulations

TC_select_ind = randsample(length(table2array(COP_sample_TC)),count_select_TC); % generate random N numbers where N=count_select_TC
ETC_select_ind = randsample(length(table2array(COP_sample_TC)),count_select_ETC); % generate random N numbers where N=count_select_ETC

cop_sample_select_TC = [COP_sample_TC.Var1(TC_select_ind) COP_sample_TC.Var2(TC_select_ind)];
cop_sample_select_ETC = [COP_sample_ETC.Var1(ETC_select_ind) COP_sample_ETC.Var2(ETC_select_ind)];

Cop_sample_Combined = [cop_sample_select_TC; cop_sample_select_ETC];


%% Plotting using simulated KDE proportionate to the innitial threshold exceedences
figure


for i=1:length(rp)

    prd = ksdensity(Cop_sample_Combined,[cord_con_comb(i).X',cord_con_comb(i).Y']); % Calculating the prbability density
    
    pr=prd/max(prd); % normalize the relative probability
    fig4=scatter(cord_con_comb(i).X',cord_con_comb(i).Y',30,pr,'filled');hold on;

    % Plotting most lkely event
    [~,ind_MLI] = max(pr);
    MLE = scatter(cord_con_comb(i).X(ind_MLI)',cord_con_comb(i).Y(ind_MLI)','k','filled','^');hold on;

    str = ['RP = ',RP(i)];
    str=strcat(str(1,1),str(1,2));
    m=max(cord_con_comb(i).Y);
    text(0,m+5,str,"FontSize",F_Size-2,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
end

xlabel('NTR (m)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
ylabel('18 Hour Rainfall (mm)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
c = colorbar;c.Label.String = "Relative Probability";c.Label.Rotation = 270;c.Label.VerticalAlignment = "bottom"; c.Label.FontSize = F_Size;
colormap(flipud(c_map))
ax = gca;
ax.XTickLabelRotation= 0;
title('')
hold on
ylim([0 175]);xlim([0 3]);


%Plotting Extreme event sample: TC and ETC
TCC=scatter(Non_TC_ext_CON_ntr_yrs(:,2),Non_TC_ext_CON_ntr_yrs(:,4),15,'red','filled');hold on;
scatter(Non_TC_ext_CON_RF_yrs(:,2),Non_TC_ext_CON_RF_yrs(:,4),15,'red','filled');
ETCC=scatter(TC_ext_CON_ntr_yrs(:,2),TC_ext_CON_ntr_yrs(:,4),15,'blue','filled');hold on;
scatter(TC_ext_CON_RF_yrs(:,2),TC_ext_CON_RF_yrs(:,4),15,'blue','filled');


hh.Position=[0 0 width height];


legend([MLE,TCC,ETCC],'Most likely event','Sample Non-TC','Sample TC');
set(gcf,'units','inches','position',[x0,y0,width,height])
saveas(fig4,'Combined_ppl_with_simulated_kde')




%% Ploting with relative contribution

c_map = [ linspace(1,0,512)', linspace(0,1,512)', linspace(0,0,512)'];


figure
for i=1:length(rp)
    Zm_TC = griddata(TC.Var1, TC.Var2, An_ex_comb_TC, cord_con_comb(i).X,cord_con_comb(i).Y);
    Zm_ETC = griddata(ETC.Var1, ETC.Var2, An_ex_comb_ETC, cord_con_comb(i).X,cord_con_comb(i).Y);
    
    %Contribution to the isoline probabilities
    pr=Zm_TC./(Zm_ETC+Zm_TC)*100;

%     length(pr)

    fig10=scatter(cord_con_comb(i).X',cord_con_comb(i).Y',[],pr,'filled');hold on;
    str = ['RP = ',RP(i)];
    str=strcat(str(1,1),str(1,2));
    m=max(cord_con_comb(i).Y);
    text(-1,m+5,str,"FontSize",F_Size,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
    
    Isoline_data(i).cord_con_comb_X=cord_con_comb(i).X; % Creating a data file for saving the coridante and prob. data
    Isoline_data(i).cord_con_comb_Y=cord_con_comb(i).Y;
    Isoline_data(i).prob_TC=Zm_TC;
    Isoline_data(i).prob_ETC = Zm_ETC;

end
save('Isoline_data','Isoline_data') % Saving the isoline data file

xlabel('NTR (m)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
ylabel('18 Hour Rainfall (mm)');set(gca,'FontSize',F_Size,'FontName','Times New Roman');grid on
c = colorbar;c.Label.String = "Relative contribution for AEP from TC (%)";c.Label.Rotation = 270;c.Label.VerticalAlignment = "bottom"; c.Label.FontSize = F_Size;
%c.Ticks = linspace(0, 100, 5) ; %Create 8 ticks from zero to 1
%c.TickLabels = {'100','75','50','75','100'} ;

colormap(flipud(c_map))
clim([0 100])
ax = gca;
ax.XTickLabelRotation= 0;
title('')
hold on
ylim([0 175]);xlim([0 3]);



%Plotting Extreme event sample: TC and ETC
scatter(Non_TC_ext_CON_ntr_yrs(:,2),Non_TC_ext_CON_ntr_yrs(:,4),15,'red','filled');hold on;
scatter(Non_TC_ext_CON_RF_yrs(:,2),Non_TC_ext_CON_RF_yrs(:,4),15,'red','filled');
scatter(TC_ext_CON_ntr_yrs(:,2),TC_ext_CON_ntr_yrs(:,4),15,'blue','filled');hold on;
scatter(TC_ext_CON_RF_yrs(:,2),TC_ext_CON_RF_yrs(:,4),15,'blue','filled');

hh.Position=[0 0 width height];



set(gcf,'units','inches','position',[x0,y0,width,height])
legend('','','','','','','Sample Non-TC','','Sample TC');
saveas(fig10,'Relative_contribution_from_TC_for_combined_ppl')

