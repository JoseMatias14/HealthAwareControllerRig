% This code computes the change in the probe area along the experiment and compares it with the pressure drop. 


% Author: Jose Otavio Matias
% Work address
% email: jose.o.a.matias@ntnu.no
% Novemebr 2020; Last revision: 
clear 
close all
clc
% 

%% Parameters tuning 
%number of wells being used in the experiment
par.nw = 3;

% Sampling times 
par.dtMeasu = 1; %[s]

%coefficient of discharge
par.Cd_hat = 0.00963481031516700;

% minimum equivalent diameter
par.dMin = 0.3*0.5;
par.dMax = 0.3181*1.5;

%% Data 
% exp 1:9
% cc (number of wells)
% data{exp,cc}.flowrate [L/min]
% data{exp,cc}.deltaP = [mbar]
% data{exp,cc}.temperature = [oC};
% data{exp,cc}.ptop = [mbar g];
% data{exp,cc}.ppump = [bar g];
% data{exp,cc}.pumpRotation = [0 - 1]
% data{exp,cc}.valveOpen = [0 - 1]
% data{exp,cc}.equivDiame = [cm]

par.currentdirectory = pwd;
par.nfolder = dir([par.currentdirectory,'\*.txt']);%pictures taken with png format
par.ni = size(par.nfolder,1); %number of experiments

%load('Experiments')
load('Experiments_dp_no_scale')

sampTime = 30; % [s]

% experimental table
%       | well 1 | well 2 | well 3 |
% exp 1 |   13   |   05   |   09   |
% exp 2 |   09   |   09   |   13   |
% exp 3 |   05   |   13   |   05   |
% exp 4 |   09   |   09   |   09   |
% exp 7 |   XX   |   11   |   XX   |
% exp 8 |   XX   |   07   |   XX   |
% exp 9 |   XX   |   XX   |   11   |

        % experiment, well, flowrate
expTable = [1,1,13;
            2,3,13;
            3,2,13;
            1,3,9;
            2,1,9;
            2,2,9;
            1,2,5;
            3,1,5;
            3,3,5];
           
%% Diameter
% new plot

figure(1)
for exp = 1:size(expTable,1)

    expTime = size(data{expTable(exp,1),1}.flowrate,1);
   
    tEval = 1:sampTime:expTime;
        
    if rem(tEval(end),sampTime) ~= 0
        tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
    else
        tgrid = 0:(sampTime/60):(1/60)*expTime;
    end
    
    filtTemp = medfilt1(data{expTable(exp,1),expTable(exp,2)}.equivDiame,sampTime);

        
    subplot(3,3,exp)
        plot(tgrid,data{expTable(exp,1),expTable(exp,2)}.equivDiame(tEval),'ko','Linewidth',1.5)
        hold on
        plot(tgrid,filtTemp(tEval),'b-','Linewidth',1.5)

        yline(0.3,'k:');
        yline(0.3181,'k:');
        
        xticks(0:5:tgrid(end))
        xlim([0, tgrid(end)])

        %             yticks(0:0.25:1)
        ylim([0.26, 0.33])

        if exp < 4
            lab = ['Q = 13; Well = ',num2str(expTable(exp,2))];
        elseif exp < 7
            lab = ['Q = 9; Well = ',num2str(expTable(exp,2))];
        else
            lab = ['Q = 5; Well = ',num2str(expTable(exp,2))];
        end
        
        title(lab)
        ylabel('d [cm]','FontSize',10)
        xlabel('time [min]','FontSize',10);

end

% %% Flowrate
% % new plot
% 
% figure(2)
% for exp = 1:size(expTable,1)
% 
%     expTime = size(data{expTable(exp,1),1}.flowrate,1);
%        
%     tEval = 1:sampTime:expTime;
%         
%     if rem(tEval(end),sampTime) ~= 0
%         tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
%     else
%         tgrid = 0:(sampTime/60):(1/60)*expTime;
%     end
%     
%     subplot(3,3,exp)
%         plot(tgrid,data{expTable(exp,1),expTable(exp,2)}.flowrate(tEval),'ko-','Linewidth',1.5)
% 
%         
%         xticks(0:5:tgrid(end))
%         xlim([0, tgrid(end)])
%         
%         if exp < 4
%             lab = ['Q = 13; Well = ',num2str(expTable(exp,2))];
%         elseif exp < 7
%             lab = ['Q = 9; Well = ',num2str(expTable(exp,2))];
%         else
%             lab = ['Q = 5; Well = ',num2str(expTable(exp,2))];
%         end
%                 
%         title(lab)
%         ylabel('Q [L/min]','FontSize',10)
%         xlabel('time [min]','FontSize',10);
% 
% end

%% Delta Pressure
figure(3)
for exp = 1:size(expTable,1)

    expTime = size(data{expTable(exp,1),1}.flowrate,1);
   
    tEval = 1:sampTime:expTime;
        
    if rem(tEval(end),sampTime) ~= 0
        tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
    else
        tgrid = 0:(sampTime/60):(1/60)*expTime;
    end
    
    Qtest = expTable(exp,3);
    
    if expTable(exp,2) == 1
        DPmax = 0.3081*Qtest^2 + 0.185*Qtest - 0.7528;
        DPmin = 0.1828*Qtest^2 + 1.0093*Qtest - 4.6835;
    elseif expTable(exp,2) == 2
        DPmax = 0.2934*Qtest^2 - 0.0592*Qtest - 0.7532;
        DPmin = 0.1896*Qtest^2 + 0.8415*Qtest - 4.9335;    
    else
        DPmax = 0.2902*Qtest^2 + 0.9971*Qtest - 1.5097;
        DPmin =  0.19*Qtest^2 + 1.1067*Qtest - 3.685;    
    end
    
    filtTemp = medfilt1(data{expTable(exp,1),expTable(exp,2)}.deltaP,sampTime);

    
    subplot(3,3,exp)
        plot(tgrid,data{expTable(exp,1),expTable(exp,2)}.deltaP(tEval),'ko','Linewidth',1.5)
        hold on
        plot(tgrid,filtTemp(tEval),'b-','Linewidth',1.5)
 
        yline(DPmin,'k:');
        yline(DPmax,'k:');
        
        xticks(0:5:tgrid(end))
        xlim([0, tgrid(end)])

        %ylim([0.26, 0.33])

        if exp < 4
            lab = ['Q = 13; Well = ',num2str(expTable(exp,2))];
        elseif exp < 7
            lab = ['Q = 9; Well = ',num2str(expTable(exp,2))];
        else
            lab = ['Q = 5; Well = ',num2str(expTable(exp,2))];
        end
        
        title(lab)
        ylabel('DP [mbar]','FontSize',10)
        xlabel('time [min]','FontSize',10);

end

%% Relative Delta Pressure
figure(4)
for exp = 1:size(expTable,1)

    expTime = size(data{expTable(exp,1),1}.flowrate,1);
       
    tEval = 1:sampTime:expTime;
        
    if rem(tEval(end),sampTime) ~= 0
        tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
    else
        tgrid = 0:(sampTime/60):(1/60)*expTime;
    end
    
    Qtest = expTable(exp,3);
    
    if expTable(exp,2) == 1
        DPmax = 0.3081*Qtest^2 + 0.185*Qtest - 0.7528;
        DPmin = 0.1828*Qtest^2 + 1.0093*Qtest - 4.6835;
    elseif expTable(exp,2) == 2
        DPmax = 0.2934*Qtest^2 - 0.0592*Qtest - 0.7532;
        DPmin = 0.1896*Qtest^2 + 0.8415*Qtest - 4.9335;    
    else
        DPmax = 0.2902*Qtest^2 + 0.9971*Qtest - 1.5097;
        DPmin =  0.19*Qtest^2 + 1.1067*Qtest - 3.685;    
    end
    
    DPNorm = (data{expTable(exp,1),expTable(exp,2)}.deltaP - DPmin)./(DPmax - DPmin);
    filtTemp = medfilt1(DPNorm,sampTime);
    
    subplot(3,3,exp)
        plot(tgrid,DPNorm(tEval),'ko','Linewidth',1.5)
        hold on
        plot(tgrid,filtTemp(tEval),'b-','Linewidth',1.5)
        
        hold on
 
        yline(1,'k:');
        yline(0,'k:');
        
        xticks(0:5:tgrid(end))
        xlim([0, tgrid(end)])

        %ylim([0.26, 0.33])

        if exp < 4
            lab = ['Q = 13; Well = ',num2str(expTable(exp,2))];
        elseif exp < 7
            lab = ['Q = 9; Well = ',num2str(expTable(exp,2))];
        else
            lab = ['Q = 5; Well = ',num2str(expTable(exp,2))];
        end
        
        title(lab)
        ylabel('DP [mbar]','FontSize',10)
        xlabel('time [min]','FontSize',10);

end