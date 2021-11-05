% Comparing MC iterations: deterministic controller, ms controller with fewer scenarios and with multiple scenarios

clear
close all
clc

% loading data
noCont = load('HAC_no_controller'); 
det = load('HAC_deterministic_nt_MC_30_10'); 
ms = load('HAC_multiscenario_nt_MC_30_10_fewer'); 

%% Preparing data
nMC = 100;
par.T = 60; % simulation sampling time[s]
nFinal = 35*par.T; %[sampling time] - arbitrarily chosen
par.x_threshold = 0.316; % probe diameter [cm]

% minimum equivalent diameter
par.dMin = 0.3;
par.dMax = 0.3181;

markers = {'o','x','>'};
cc = {'b','k','g'};
leg = {'w_1','w_2','w_3'};

for dd = 1:3
    execTime{dd} = [];
    breakTime{dd} = [];
    profit{dd} = [];
    contFail{dd} = [];

end


for ii = 1:nMC
   % computing execution time
   execTime{1} = [execTime{1}, noCont.controlTimeArray{ii}];
   execTime{2} = [execTime{2}, det.controlTimeArray{ii}];
   execTime{3} = [execTime{3}, ms.controlTimeArray{ii}];
     
   % computing breakdown time
   breakTime{1} = [breakTime{1}, noCont.tgridArray{ii}(end)];
   breakTime{2} = [breakTime{2}, det.tgridArray{ii}(end)];
   breakTime{3} = [breakTime{3}, ms.tgridArray{ii}(end)];
   
   % computing total profit
   profit{1} = [profit{1}, sum(noCont.ofPlantArray{ii})];
   profit{2} = [profit{2}, sum(det.ofPlantArray{ii})];
   profit{3} = [profit{3}, sum(ms.ofPlantArray{ii})];
   
   contFail{1} = [contFail{1}, sum(noCont.flagArray{ii}(3,:))];
   contFail{2} = [contFail{2}, sum(det.flagArray{ii}(3,:))];
   contFail{3} = [contFail{3}, sum(ms.flagArray{ii}(3,:))];
end

%% plotting time distribution
figure(1)
for dd = 1:3    
    histogram(execTime{dd},'BinWidth',0.15)
    hold on
end

xlim([0, 10])
xlabel('time [s]')
ylabel('count')
legend({'det.','ms_{fewer}','ms'})
title('Execution time')

%% break down time distribution
figure(2)
boxplot([breakTime{1}',breakTime{2}',breakTime{3}'],'Labels',{'det.','ms_{fewer}','ms'})
grid on

%xlim([0, 10])
ylabel('time [min]')
title('Breakdown time')

%% failure percentage
figure(3)

temp1 = [];
temp2 = [];


for ii = 1:nMC
   temp1 = [temp1, [breakTime{1}(ii);breakTime{2}(ii);breakTime{3}(ii)]];
   temp2 = [temp2, [contFail{1}(ii)/(breakTime{1}(ii) - 1);contFail{2}(ii)/(breakTime{2}(ii) - 1);contFail{3}(ii)/(breakTime{3}(ii) - 1)]];
 
end

subplot(1,2,1)
    %plot(1:100, temp1,'r+')
    %plot(1:100, temp2(1,:),'bo');
    plot(1:100, 100*temp2(2,:),'b+');
    hold on
    plot(1:100, 100*temp2(3,:),'kd');

    grid on
    %yline(1,'r');
    
    ylim([80,100])
    
    xlabel('Monte Carlo Iteration [-]')
    ylabel('Percentage of failures [%]')
    
    title('Control Failure')

subplot(1,2,2)
    boxplot([100*temp2(2,:)',100*temp2(3,:)'],'Labels',{'det','ms'})
    grid on
    
    ylim([80 100])
    title('Control Failure')

%% profit distribution
figure(4)
boxplot([profit{1}',profit{2}',profit{3}'],'Labels',{'det.','ms_{fewer}','ms'})
grid on

%xlim([0, 10])
ylabel('time [min]')
title('Total profit')

%% Profit 
figure(5);

subplot(3,1,1)
    hold on
    for mm = 1:nMC
        plot(noCont.tgridArray{mm}, noCont.ofPlantArray{mm},'Color',[0,0,0,0.05],'Linewidth',1)
    end
    %ylim([0.05 0.95])
    
    xticks(0:(10*par.T/60):(1/60)*(nFinal))
    xlim([0 (1/60)*(nFinal)])
    
    xlabel('time [min]','FontSize',10)
    ylabel('J [$]','FontSize',10)
    
    title('Deterministic','FontSize',10)

subplot(3,1,2)
    hold on
    for mm = 1:nMC
        plot(det.tgridArray{mm}, det.ofPlantArray{mm},'Color',[0,0,0,0.05],'Linewidth',1)
    end
    %ylim([0.05 0.95])
    
    xticks(0:(10*par.T/60):(1/60)*(nFinal))
    xlim([0 (1/60)*(nFinal)])
    
    xlabel('time [min]','FontSize',10)
    ylabel('J [$]','FontSize',10)
    
    title('Multiscenario Fewer','FontSize',10)

subplot(3,1,3)
    hold on
    for mm = 1:nMC
        plot(ms.tgridArray{mm}, ms.ofPlantArray{mm},'Color',[0,0,0,0.05],'Linewidth',1)
    end
    %ylim([0.05 0.95])
    
    xticks(0:(10*par.T/60):(1/60)*(nFinal))
    xlim([0 (1/60)*(nFinal)])
    
    xlabel('time [min]','FontSize',10)
    ylabel('J [$]','FontSize',10)
    
    title('Multiscenario','FontSize',10)
    
%    

%% Probe diameter
figure(6);
for well = 1:3
    subplot(3,3,well)
        hold on
        for mm = 1:nMC
            plot(noCont.tgridArray{mm}, noCont.probeOrificeArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
        end
        yline(par.x_threshold,'r:','Linewidth',1);

        ylim([par.dMin 0.32])
        yticks(par.dMin:0.0036:par.dMax)
        lab = ['Well ',num2str(well),': d [cm]'];
        ylabel(lab,'FontSize',10)

        xticks(0:(10*par.T/60):(1/60)*(nFinal))
        xlim([0 (1/60)*(nFinal)])
        xlabel('time [min]','FontSize',10)
        if well == 2
            title('Deterministic','FontSize',10)
        end
        
        
    subplot(3,3,well + 3)
        hold on
        for mm = 1:nMC
            plot(det.tgridArray{mm}, det.probeOrificeArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
        end
        yline(par.x_threshold,'r:','Linewidth',1);

        ylim([par.dMin 0.32])
        yticks(par.dMin:0.0036:par.dMax)
        lab = ['Well ',num2str(well),': d [cm]'];
        ylabel(lab,'FontSize',10)

        xticks(0:(10*par.T/60):(1/60)*(nFinal))
        xlim([0 (1/60)*(nFinal)])
        xlabel('time [min]','FontSize',10)
        if well == 2
            title('Multiscenario Fewer','FontSize',10)
        end        
        
     subplot(3,3,well + 6)
        hold on
        for mm = 1:nMC
            plot(ms.tgridArray{mm}, ms.probeOrificeArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
        end
        yline(par.x_threshold,'r:','Linewidth',1);

        ylim([par.dMin 0.32])
        yticks(par.dMin:0.0036:par.dMax)
        lab = ['Well ',num2str(well),': d [cm]'];
        ylabel(lab,'FontSize',10)

        xticks(0:(10*par.T/60):(1/60)*(nFinal))
        xlim([0 (1/60)*(nFinal)])
        xlabel('time [min]','FontSize',10)
        if well == 2
            title('Multiscenario','FontSize',10)
        end        
        
end

%% Inputs
figure(7);
for well = 1:3
    subplot(3,3,well)
        hold on
        for mm = 1:nMC
            plot(noCont.tgridArray{mm}(2:end), noCont.uImpArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
        end
        ylim([0.05 0.95])
        lab = ['Well ',num2str(well),': v [-]'];
        ylabel(lab,'FontSize',10)

        xticks(0:(10*par.T/60):(1/60)*(nFinal))
        xlim([0 (1/60)*(nFinal)])
        xlabel('time [min]','FontSize',10)
        if well == 2
            title('Deterministic','FontSize',10)
        end
        
        
    subplot(3,3,well + 3)
        hold on
        for mm = 1:nMC
            plot(det.tgridArray{mm}(2:end), det.uImpArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
        end
        ylim([0.05 0.95])
        lab = ['Well ',num2str(well),': v [-]'];
        ylabel(lab,'FontSize',10)

        xticks(0:(10*par.T/60):(1/60)*(nFinal))
        xlim([0 (1/60)*(nFinal)])
        xlabel('time [min]','FontSize',10)
        if well == 2
            title('Multiscenario Fewer','FontSize',10)
        end        
        
     subplot(3,3,well + 6)
        hold on
        for mm = 1:nMC
            plot(ms.tgridArray{mm}(2:end), ms.uImpArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
        end
        ylim([0.05 0.95])
        lab = ['Well ',num2str(well),': v [-]'];
        ylabel(lab,'FontSize',10)

        xticks(0:(10*par.T/60):(1/60)*(nFinal))
        xlim([0 (1/60)*(nFinal)])
        xlabel('time [min]','FontSize',10)
        if well == 2
            title('Multiscenario','FontSize',10)
        end        
        
end

