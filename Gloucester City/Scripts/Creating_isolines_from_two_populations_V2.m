clear 
clc
RP ={'5','10','20','50','100'};

TC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\RP_TC.csv");

% Here is for the last 30 years of data
% TC_ext_CON_ntr_yrs = load("C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_events_con_NTR_1991.txt");
% TC_ext_CON_RF_yrs = load("C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_events_con_RF_1991.txt");
% Non_TC_ext_CON_ntr_yrs = load("C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ETC_events_con_NTR_1991.txt");
% Non_TC_ext_CON_RF_yrs = load("C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ETC_events_con_RF_1991.txt");
% 
% Data_1901 = load("C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\NTR_and_18hr_RF_acc_data_30.mat");

% Here is for last 120 years of data
TC_ext_CON_ntr_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_events_conditioning_POT_NTR_18hr_Acc.txt");
TC_ext_CON_RF_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_events_conditioning_POT_RF_for_18hr_RF_acc.txt");
Non_TC_ext_CON_ntr_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ET_events_conditioning_POT_NTR_18hr_Acc.txt");
Non_TC_ext_CON_RF_yrs = load("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ET_events_conditioning_POT_RF_for_18hr_RF_acc.txt");

Data_1901 = load("C:\Users\pr109704\OneDrive - Knights - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\NTR_and_18hr_RF_acc_data.mat");
Data_1901.D=Data_1901.C;




COP_sample_TC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\TC_Cop_Sample.csv");
COP_sample_ETC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\ETC_Cop_Sample.csv");


event = {'Hazel','Unnamed','Sandy','Connie','Isabel','','Irene','','David','Floyd','Ida','','Ernesto','','Agnes','Able','','Groundhog',     'Danielle','Ivan','Bertha','Isaias','Jeanne'};
date ={'10/15/1954','08/23/1933','10/30/2012','08/12/1955','09/19/2003','', '08/29/2011', '','09/06/1979','09/17/1999','09/02/2021','','09/03/2006','','06/23/1972','09/01/1952','','02/02/1952',    '09/27/1992','09/21/2004','07/14/1996','08/15/2020','09/30/2004'};


c_map = [ linspace(0.8,1,512)', linspace(0,0.9,512)', linspace(0,0.2,512)'];
F_Size = 15;
Con_NTR_LW = 1;
Con_NTR_M_Size = 6;
Con_RF_LW = 1;
Con_RF_M_Size = 4;
width = 20;
height = 20;
resolution = 1000;
x0=0.2;
y0=0.2;

x_U_lim = 8;
x_L_lim = -1.2;
Y_U_ilm = 300;
Y_L_lim = 0;


%%
rp=[5 10 20 50 100];

xv = linspace(x_L_lim, x_U_lim, resolution);
yv = linspace(Y_L_lim, Y_U_ilm, resolution);
[Xm,Ym] = ndgrid(xv, yv);

Zm1 = griddata(TC.Var1, TC.Var2, TC.V3, Xm, Ym);

for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm1,[rp(i) rp(i)]); hold on;
    cord_con1(i).X = M(1,2:end);
    cord_con1(i).Y = M(2,2:end);
end


%% for conditioning parameter two

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

% Ploting
cord_con_comb_TC=[];
figure
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm3,[rp(i) rp(i)]); hold on;
    cord_con_comb_TC(i).X = M(1,2:end);
    cord_con_comb_TC(i).Y = M(2,2:end);
end
%%

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

% % Ploting the hourly events
% NTR = Data_1901.D(:,1);
% RF = Data_1901.D(:,2);
% fig2=plot(NTR,RF,Color=[0.7 0.7 0.7],Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');


%Plotting Extreme event sample: TC
ntr=plot(TC_ext_CON_ntr_yrs(:,2),TC_ext_CON_ntr_yrs(:,4),Color='r',Marker='o',MarkerSize=Con_NTR_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_NTR_LW);hold on;
rf=plot(TC_ext_CON_RF_yrs(:,2),TC_ext_CON_RF_yrs(:,4),Color='b',Marker='^',MarkerSize=Con_RF_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_RF_LW);
hh.Position=[0 0 width height];

%for TC events
loc1 = sortrows(TC_ext_CON_ntr_yrs,2,"descend");
% for m=1:11
%     str = [event(m),date(m)];
%     if m==1 
%         arrow([loc1(m,2)-0.3,loc1(m,4)+12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.6,loc1(m,4)+11,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
% %     elseif m==8
% %         arrow([loc1(m,2)-0.5,loc1(m,4)-20],[loc1(m,2),loc1(m,4)],6);
% %         text(loc1(m,2)-0.7,loc1(m,4)-19,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
% %     elseif  m==3
% %         arrow([loc1(m,2)-0.2,loc1(m,4)+6],[loc1(m,2),loc1(m,4)],6);
% %         text(loc1(m,2)-0.55,loc1(m,4)+13,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
%     else
%         arrow([loc1(m,2)+0.3,loc1(m,4)+12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)+0.3,loc1(m,4)+16,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
%     end
% end



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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading for ETC

ETC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\RP_ETC.csv");

%%

Zm1 = griddata(ETC.Var1, ETC.Var2, ETC.V3, Xm, Ym);
figure
for i=1:length(rp)
    [M,~] =contour(Xm, Ym, Zm1,[rp(i) rp(i)]); hold on;
    cord_con1(i).X = M(1,2:end);
    cord_con1(i).Y = M(2,2:end);
end


%%

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

% 
% % Ploting the hourly events
% NTR = Data_1901.D(:,1);
% RF = Data_1901.D(:,2);
% plot(NTR,RF,Color=[0.7 0.7 0.7],Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');
% 

%Plotting Extreme event sample: TC
ntr=plot(Non_TC_ext_CON_ntr_yrs(:,2),Non_TC_ext_CON_ntr_yrs(:,4),Color='r',Marker='o',MarkerSize=Con_NTR_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_NTR_LW);hold on;
rf=plot(Non_TC_ext_CON_RF_yrs(:,2),Non_TC_ext_CON_RF_yrs(:,4),Color='b',Marker='^',MarkerSize=Con_RF_M_Size,LineStyle='none',MarkerFaceColor='auto',LineWidth=Con_RF_LW);


hh.Position=[0 0 width height];
legend([fig3,MLE,ntr,rf],'JP Isoline','Most likely event','Sample Con NTR','Sample Con RF');
set(gcf,'units','inches','position',[x0,y0,width,height])
saveas(fig3,'Combined_Non_TC_with_Simulated_kde')

%% Combining two populations

% ETC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\RP_ETC.csv");
% TC = readtable("C:\Users\pr109704\OneDrive - University of Central Florida\Analysis_DW_5d_TC_350km_Philli_airport_AORC_Chelsa_V1\Copulas\RP_TC.csv");

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

% An_ex_comb_ppl = An_ex_comb_TC+An_ex_comb_ETC; This was the previous
% approach

% The new approach
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
%%

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


% % Ploting the hourly events
% NTR = Data_1901.D(:,1);
% RF = Data_1901.D(:,2);
% plot(NTR,RF,Color=[0.7 0.7 0.7],Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');


%Plotting Extreme event sample: TC and ETC
TCC=scatter(Non_TC_ext_CON_ntr_yrs(:,2),Non_TC_ext_CON_ntr_yrs(:,4),15,'red','filled');hold on;
scatter(Non_TC_ext_CON_RF_yrs(:,2),Non_TC_ext_CON_RF_yrs(:,4),15,'red','filled');
ETCC=scatter(TC_ext_CON_ntr_yrs(:,2),TC_ext_CON_ntr_yrs(:,4),15,'blue','filled');hold on;
scatter(TC_ext_CON_RF_yrs(:,2),TC_ext_CON_RF_yrs(:,4),15,'blue','filled');


hh.Position=[0 0 width height];


% %for TC events
% loc1 = sortrows(TC_ext_CON_ntr_yrs,2,"descend");
% for m=1:10
%     str = [event(m),date(m)];
%     if m==6 || m==7 
%         arrow([loc1(m,2)-0.3,loc1(m,4)-12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.6,loc1(m,4)-11,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
%     elseif m==8
%         arrow([loc1(m,2)-0.5,loc1(m,4)-20],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.7,loc1(m,4)-19,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
%     elseif m==2 || m==3
%         arrow([loc1(m,2)-0.2,loc1(m,4)+6],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.55,loc1(m,4)+13,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
%     else
%         arrow([loc1(m,2)+0.3,loc1(m,4)+12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)+0.3,loc1(m,4)+16,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
%     end
% end



legend([MLE,TCC,ETCC],'Most likely event','Sample Non-TC','Sample TC');
set(gcf,'units','inches','position',[x0,y0,width,height])
saveas(fig4,'Combined_ppl_with_simulated_kde')


%% plotting with Simulations and MDE
figure


for i=1:length(rp)

    prd = ksdensity(Cop_sample_Combined,[cord_con_comb(i).X',cord_con_comb(i).Y']); % Calculating the prbability density
    
    pr=prd/max(prd); % normalize the relative probability
    fig4=scatter(cord_con_comb(i).X',cord_con_comb(i).Y',[],pr,'filled');hold on;

    % Plotting most lkely event
    [~,ind_MLI] = max(pr);
    MLE = scatter(cord_con_comb(i).X(ind_MLI)',cord_con_comb(i).Y(ind_MLI)','k','filled','^');hold on;

    str = ['RP = ',RP(i)];
    str=strcat(str(1,1),str(1,2));
    m=max(cord_con_comb(i).Y);
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
ylim([0 250]);xlim([-1.2 7]);




% fig7=plot(COP_sample_TC.Var1,COP_sample_TC.Var2,Color='b',Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');
% 
% fig8=plot(COP_sample_ETC.Var1,COP_sample_ETC.Var2,Color='r',Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');

fig9=plot(Cop_sample_Combined(:,1),Cop_sample_Combined(:,2),Color=[0.7 0.7 0.7],Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');



%Plotting Extreme event sample: TC and ETC
% scatter(ETC_ext_CON_ntr_30yrs(:,2),ETC_ext_CON_ntr_30yrs(:,4),15,'red','filled');hold on;
% scatter(ETC_ext_CON_RF_30yrs(:,2),ETC_ext_CON_RF_30yrs(:,4),15,'red','filled');
% scatter(TC_ext_CON_ntr_30yrs(:,2),TC_ext_CON_ntr_30yrs(:,4),15,'blue','filled');hold on;
% scatter(TC_ext_CON_RF_30yrs(:,2),TC_ext_CON_RF_30yrs(:,4),15,'blue','filled');
% 

hh.Position=[0 0 width height];


%for TC events
% loc1 = sortrows(TC_ext_CON_ntr_30yrs,2,"descend");
% for m=1:10
%     str = [event(m),date(m)];
%     if m==6 || m==7 
%         arrow([loc1(m,2)-0.3,loc1(m,4)-12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.6,loc1(m,4)-11,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
%     elseif m==8
%         arrow([loc1(m,2)-0.5,loc1(m,4)-20],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.7,loc1(m,4)-19,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
%     elseif m==2 || m==3
%         arrow([loc1(m,2)-0.2,loc1(m,4)+6],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.55,loc1(m,4)+13,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
%     else
%         arrow([loc1(m,2)+0.3,loc1(m,4)+12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)+0.3,loc1(m,4)+16,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
%     end
% end

ind=MDA(Cop_sample_Combined,100);
fig3=plot(Cop_sample_Combined(ind,1),Cop_sample_Combined(ind,2),'k',Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');
set(gcf,'units','inches','position',[x0,y0,width,height])


legend([fig4,fig9,fig3],'Isolines','Simulated sample','MDA sample');
saveas(fig9,'Combined_ppl_with_simulated_kde_MDA_sample')


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


% % Ploting the hourly events
% NTR = Data_30.D(:,1);
% RF = Data_30.D(:,2);
% fig=plot(NTR,RF,Color=[0.7 0.7 0.7],Marker='.',MarkerSize=12,LineStyle='none',MarkerFaceColor='auto');


%Plotting Extreme event sample: TC and ETC
scatter(Non_TC_ext_CON_ntr_yrs(:,2),Non_TC_ext_CON_ntr_yrs(:,4),15,'red','filled');hold on;
scatter(Non_TC_ext_CON_RF_yrs(:,2),Non_TC_ext_CON_RF_yrs(:,4),15,'red','filled');
scatter(TC_ext_CON_ntr_yrs(:,2),TC_ext_CON_ntr_yrs(:,4),15,'blue','filled');hold on;
scatter(TC_ext_CON_RF_yrs(:,2),TC_ext_CON_RF_yrs(:,4),15,'blue','filled');

hh.Position=[0 0 width height];


% %for TC events
% loc1 = sortrows(TC_ext_CON_ntr_30yrs,2,"descend");
% for m=1:10
%     str = [event(m),date(m)];
%     if m==6 || m==7 
%         arrow([loc1(m,2)-0.3,loc1(m,4)-12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.6,loc1(m,4)-11,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
%     elseif m==8
%         arrow([loc1(m,2)-0.5,loc1(m,4)-20],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.7,loc1(m,4)-19,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
%     elseif m==2 || m==3
%         arrow([loc1(m,2)-0.2,loc1(m,4)+6],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)-0.55,loc1(m,4)+13,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');    
%     else
%         arrow([loc1(m,2)+0.3,loc1(m,4)+12],[loc1(m,2),loc1(m,4)],6);
%         text(loc1(m,2)+0.3,loc1(m,4)+16,str,"FontSize",14,"FontName",'Times New Roman', 'VerticalAlignment','top', 'HorizontalAlignment','left');
%     end
% end


set(gcf,'units','inches','position',[x0,y0,width,height])
legend('','','','','','','Sample Non-TC','','Sample TC');
saveas(fig10,'Relative_contribution_from_TC_for_combined_ppl')

%%
save("Combined_Isolines","cord_con_comb")