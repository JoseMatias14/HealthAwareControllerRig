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

% analysis samp. time
sampTime = 60; % [s]

% count for plotting
count = 1;

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
load('ExperimentsNN')


%% Delta Pressure (prediction vs. measured)
% for expe = 1:par.ni 
% 
%     % evaluating experimental time
%     expTime = size(data{expe,1}.flowrate,1);
%    
%     tEval = 1:sampTime:expTime;
%         
%     if rem(tEval(end),sampTime) ~= 0
%         tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
%     else
%         tgrid = 0:(sampTime/60):(1/60)*expTime;
%     end
%     
%     for well = 1:3
%         % start analysis at 5 minutes
%         DP_hat_kArray = [];
%        
%         for kk = 5*60:sampTime:expTime
%                        
%             RegMatrix = [normalize(data{expe,well}.relTime(kk - 150:kk)),normalize(data{expe,well}.flowrate(kk - 150:kk)),...
%                         normalize(data{expe,well}.ptop(kk - 150:kk)),normalize(data{expe,well}.temperature(kk - 150:kk)),...
%                         normalize(data{expe,well}.valvCurrent(kk - 150:kk))];
% 
%             Reg_k = RegMatrix(end,:)';
%             %using NN to predict the current normalized delta pressure
%             DP_hat_norm_k = NLIO_fun(Reg_k);
%             
%             % de-normalizing
%             mu_DP = mean(data{expe,well}.deltaP(kk - 150:kk));
%             std_DP = std(data{expe,well}.deltaP(kk - 150:kk));
%             DP_hat_k  = DP_hat_norm_k*std_DP + mu_DP; 
%             DP_hat_kArray = [DP_hat_kArray, DP_hat_k];
%             
%         end
%        
%         % filtering the data
%         filtDP = movmean(data{expe,well}.deltaP,150);
%         
%         figure(count)
%         count = count + 1;
%         subplot(3,1,well)
%             %plot(tgrid(tEval),filtDP(tEval),'ro','Linewidth',1.5)
%             plot(tgrid,data{expe,well}.deltaP(tEval),'ko','Linewidth',1.5)
%             hold on 
%             plot(5:(sampTime/60):(1/60)*expTime,DP_hat_kArray,'bx','Linewidth',1.5)
%             xticks(0:5:tgrid(end))
%             xlim([0, tgrid(end)])
%         
%         %ylim([0.26, 0.33])
%         
%         lab = ['Experiment = ',num2str(expe)];
%         title(lab)
%         ylabel('DP [mbar]','FontSize',10)
%         xlabel('time [min]','FontSize',10);
%     end
% 
% end

%% Relative Delta Pressure
figure(4)
for expe = 1:size(expTable,1)

    expTime = size(data{expTable(expe,1),1}.flowrate,1);
       
    tEval = 1:sampTime:expTime;
        
    if rem(tEval(end),sampTime) ~= 0
        tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
    else
        tgrid = 0:(sampTime/60):(1/60)*expTime;
    end
    
    Qtest = expTable(expe,3);
    
    if expTable(expe,2) == 1
        DPmax = 0.3081*Qtest^2 + 0.185*Qtest - 0.7528;
        DPmin = 0.1828*Qtest^2 + 1.0093*Qtest - 4.6835;
    elseif expTable(expe,2) == 2
        DPmax = 0.2934*Qtest^2 - 0.0592*Qtest - 0.7532;
        DPmin = 0.1896*Qtest^2 + 0.8415*Qtest - 4.9335;    
    else
        DPmax = 0.2902*Qtest^2 + 0.9971*Qtest - 1.5097;
        DPmin =  0.19*Qtest^2 + 1.1067*Qtest - 3.685;    
    end
    
    DPNorm =  1 - (data{expTable(expe,1),expTable(expe,2)}.deltaP - DPmin)./(DPmax - DPmin);
    filtTemp = medfilt1(DPNorm,sampTime);
    
    subplot(3,3,expe)
        plot(tgrid,DPNorm(tEval),'ko','Linewidth',1.5)
        hold on
        plot(tgrid,filtTemp(tEval),'b-','Linewidth',1.5)
        
        hold on
 
        yline(1,'k:');
        yline(0,'k:');
        
        xticks(0:5:tgrid(end))
        xlim([0, tgrid(end)])

        %ylim([0.26, 0.33])

        if expe < 4
            lab = ['Q = 13; Well = ',num2str(expTable(expe,2))];
        elseif expe < 7
            lab = ['Q = 9; Well = ',num2str(expTable(expe,2))];
        else
            lab = ['Q = 5; Well = ',num2str(expTable(expe,2))];
        end
        
        title(lab)
        ylabel('erosion indicator \epsilon','FontSize',10)
        xlabel('time [min]','FontSize',10);

end

%% Delta Pressure
count = 10;
for expe = 1:par.ni 

    % evaluating experimental time
    expTime = size(data{expe,1}.flowrate,1);
   
    tEval = 1:sampTime:expTime;
        
    if rem(tEval(end),sampTime) ~= 0
        tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
    else
        tgrid = 0:(sampTime/60):(1/60)*expTime;
    end
    
    % delta pressure
    for well = 1:3
        
        %calculating max and min delta P
        DPminArray = [];
        DPmaxArray = [];
        epsArray = [];
        
        for kk = 1:expTime
            Qtemp = data{expe,well}.flowrate(kk);
            if well == 1
                DPmax = 0.3081*Qtemp^2 + 0.185*Qtemp - 0.7528;
                DPmin = 0.1828*Qtemp^2 + 1.0093*Qtemp - 4.6835;
            elseif well == 2
                DPmax = 0.2934*Qtemp^2 - 0.0592*Qtemp - 0.7532;
                DPmin = 0.1896*Qtemp^2 + 0.8415*Qtemp - 4.9335;
            else % well == 3
                DPmax = 0.2902*Qtemp^2 + 0.9971*Qtemp - 1.5097;
                DPmin =  0.19*Qtemp^2 + 1.1067*Qtemp - 3.685;
            end
            DPminArray = [DPminArray, DPmin];
            DPmaxArray = [DPmaxArray, DPmax];
            
            epsNorm =  1 - (data{expe,well}.deltaP(kk) - DPmin)./(DPmax - DPmin);
            epsArray = [epsArray, epsNorm];
        end
        
        % filtering the data
        filtEps = movmean(epsArray,150);
        
        figure(count)
            count = count + 1;
            subplot(3,1,well)
                plot(tgrid,epsArray(tEval),'ko','Linewidth',1.5)
                hold on
                plot(tgrid,filtEps(tEval),'b-','Linewidth',1.5)
                plot(tgrid,DPminArray(tEval),'k:')
                plot(tgrid,DPmaxArray(tEval),'k:')

                xticks(0:5:tgrid(end))
                xlim([0, tgrid(end)])

        %ylim([0.26, 0.33])
        
        lab = ['Experiment = ',num2str(expe)];
        title(lab)
        ylabel('DP [mbar]','FontSize',10)
        xlabel('time [min]','FontSize',10);
    end

end

