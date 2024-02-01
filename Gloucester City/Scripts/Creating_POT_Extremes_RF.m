

clear; close all; clc
load("Hourly_accumulation_Bias_Corrected_RF_data_Philli_airport.mat") % Please get the Bias Corrected Data set
RF_Data = Data(1:end,:); % Taking the Rainfall Data only Belongs to the Tide guage data period


%%  Events over threshold

ts = RF_Data;
events= struct;


th = ones(1,48)*11; % initial Threshold values
for ii = 1:48
    ii
    POT=nan(10000,2); %This is just to recognize for the loop
    while length(POT) > 580
        % events over threshold events
        EOT= ts(:,ii+1);
        EOT(EOT< th(ii))= NaN;
        EOT= [ts(:,1) EOT];
        EOT(isnan(EOT(:,2)),:)= [];
        
        % Plot the RF values over Threshold
        
        %hh= figure;
        %set(hh,'units','centimeters','Position',[0.2 1.5 28 12],...
        %    'color','w');
        %hold all
        %h1= plot(ts(:,1),ts(:,2),'.-','LineStyle', 'none');
        %ylabel('Rainfall (mm)');
        %set(gca,'FontSize',16,'FontName','Bell MT')
        %grid minor
        %ax = gca;
        %ax.XTickLabelRotation= 20;
        %title('Philadelphia')
        %axs= axis;
        %hold all; 
        %plot(EOT(:,1),EOT(:,2),'.-k','LineStyle', 'none' )
        
        %%%%%%%%%%%% Declustering %%%%%%%%%
        % decluster time
        dec_tim= 2.5; % in days
        xnan= sum(isnan(EOT(:,2)));
        POT= nan(length(EOT),2);
        
        while xnan < size(EOT,1)
            
            for i= 1: size(POT,2)
                
                [~,fmax]= max(EOT(:,2));
                
                POT(fmax,:)= EOT(fmax,:);
                %hold all; plot(EOT(fmax,1),EOT(fmax,2),'or','LineWidth',2);
                
                dec_wind= find(EOT(:,1)>= EOT(fmax,1)-dec_tim & EOT(:,1)<= EOT(fmax,1)+dec_tim);
                
                %hold all; hi= plot(EOT(dec_wind,1),EOT(dec_wind,2),...
                %    '.-','Color','r');
                
                EOT(dec_wind,2)= nan;
                
                xnan= sum(isnan(EOT(:,2)));
            end
        
        end
        
        POT(isnan(POT(:,2)),:)= [];
        
        
        
        % Plot
        
        %hh= figure;
        %set(hh,'units','centimeters','Position',[0.2 1.5 28 12],...
        %    'color','w');
        %hold all
        %plot(POT(:,1),POT(:,2),'.-')
        %ylabel('Rainfall (mm)');
        %set(gca,'FontSize',16,'FontName','Bell MT')
        %grid minor
        %ax = gca;
        %ax.XTickLabelRotation= 20;
        %title('Peaks over threshold (POT) Gloucester City RF')

        th(1,ii)=th(1,ii)+0.005;
        
    end
    events(ii).POT=POT;
    events(ii).Threshold=th(1,ii)-0.005;

    th(ii+1) = th(1,ii);
    th(1,ii)

end
%%

save("POT_Events_for_each_RF_Accumulation_time","events");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Peaks over the threshold (POT)

[x,is]= sort(POT(:,2));
xtime = POT(is,1);

x= x -th;
x(x<=0) = 0.001;

hold all; plot(xtime,x,'.');

dates= datevec(ts(:,1));
years= unique(dates(:,1)); 
n= length(years);
        
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
set(gca,'FontSize',16,'FontName','Bell MT')
grid minor

%% Return level

% return period
% RP= [2 5 10 25 50 100 250 1000];
RP= (2:1:1000);

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

ylabel('Return level (mm)'); xlabel('Return period (years)')

set(gca,'FontSize',16,'FontName','Bell MT')
grid minor

legend([h1,h2,h3],'Fitted Pareto','Weibull','Gringorten','location','Best');
set(gca, 'XScale', 'log','XTick',[2 5 10 25 50 100 250 1000])
title(['Pareto q= ' num2str(q)])

axs= axis;
xlim([1 axs(2)])


%% Save figures

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%     FigHandle = FigList(iFig);
%     FigName   = ['Fig ' num2str(get(FigHandle, 'Number')) '.png'];
%     saveas(FigHandle, fullfile(FSave, FigName));
% end

