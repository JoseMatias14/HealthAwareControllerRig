% This code computes the change in the probe area along the experiment and compares it with the pressure drop. 


% Author: Jose Otavio Matias
% Work address
% email: jose.o.a.matias@ntnu.no
% Novemebr 2020; Last revision: 
clear 
close all
clc
% 

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

load('Experiments_dp_no_scale')

% controller execution
sampTime = 60; % [s]

% Buffer length
BufferLength = 60; 

%% Run configuration file
InitializationLabViewMain %here we use the same syntax as in the rig
           
%% Run mockup loop
expTime = size(data{1,1}.flowrate,1);
tEval = 1:sampTime:expTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Production Optimization Methods Data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we save the same data that is saved in the Rig's labview code
flagArray = []; % flag SSD, Model Adaptation and Economic Optimization [flag == 0 (failed), == 1 (success)]
ofArray = [];   % OF value computed bu the Economic Optimization
thetaHatArray = []; % estimated parameters
yEstArray = []; % model prediction with estimated states
yOptArray = []; % model prediction @ new optimum
uOptArray = []; % computed inputs (u_k^\star)
uImpArray = []; % filtered inputs to be implemented (u_{k+1})


%%%% 
for well = 1:3
    diameterTemp{well} = [];
end

for kk = 1:expTime
    
    if kk > BufferLength && rem(kk,sampTime) == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rearranging the data vectors and units here, so they can the be exactly %
        % the same as in the actual rig. The goal is that LabViewRTO.m can be     %
        % directly plug in the Labview interface and it will work                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measurement buffer (dim = 19 X BufferLength)
        I_vector = [data{1,1}.flowrate((kk - BufferLength + 1):kk)'; % FI-101 [L/min]
                    data{1,2}.flowrate((kk - BufferLength + 1):kk)'; % FI-102 [L/min]
                    data{1,3}.flowrate((kk - BufferLength + 1):kk)'; % FI-103 [L/min]
                    ones(3,BufferLength);                           % FI-104 [sL/min] & FI-105 [sL/min] & FI-106 [sL/min]. Not used here
                    ones(3,BufferLength); %dummy values --> in the actual rig, they will the pressure at injection point (PI105, PI106, PI107). Not used here
                    data{1,1}.temperature((kk - BufferLength + 1):kk)'; %d(TI101, TI102, TI103). Not used here
                    data{1,2}.temperature((kk - BufferLength + 1):kk)';
                    data{1,3}.temperature((kk - BufferLength + 1):kk)';
                    data{1,1}.deltaP((kk - BufferLength + 1):kk)'; % dPI-101 [mbar]
                    data{1,2}.deltaP((kk - BufferLength + 1):kk)'; % dPI-102 [mbar] 
                    data{1,3}.deltaP((kk - BufferLength + 1):kk)'; % dPI-103 [mbar]
                    data{1,1}.ptop((kk - BufferLength + 1):kk)'; % PI-101 [mbar g]
                    data{1,2}.ptop((kk - BufferLength + 1):kk)'; % PI-102 [mbar g] 
                    data{1,3}.ptop((kk - BufferLength + 1):kk)'; % PI-103 [mbar g]
                    data{1,1}.ppump((kk - BufferLength + 1):kk)'];     % PI-104 [bar g]
                
       
        % values of the input variables at the previus rig sampling time (dim = nu[7] X 1)
        P_vector = [data{1,1}.valveOpen(kk);  % CV101 opening [-]
                    data{1,2}.valveOpen(kk);  % CV102 opening [-]
                    data{1,3}.valveOpen(kk);  % CV103 opening [-]
                    data{1,1}.ppump(kk);      % dummy values: in the actual rig, they will the pump rotation. Not used here
                    1;      % dummy values: FI-104 [sL/min]
                    1;      % dummy values: FI-105 [sL/min]
                    1];     % dummy values: FI-106 [sL/min]

        % values of the inputs (valve opening) of the last optimization run (dim = nQg[3] X 1)
        %O_vector = uPlantArray(1:3,kk - nExec);
        
        tic 
        
        % Run Labview/Matlab interface file
        LabViewMain
        
        % computing execution time
        controlTime = toc;

        flagArray = [flagArray, [SS;Estimation;Optimization]]; 
        ofArray = [ofArray, Result]; 
        thetaHatArray = [thetaHatArray, Parameter_Estimation'];
        yEstArray = [yEstArray, State_Variables_Estimation'];
        yOptArray = [yOptArray, State_Variables_Optimization'];
        uOptArray = [uOptArray, Optimized_Air_Injection'];
        uImpArray = [uImpArray, O_vector]; 
        
    else
        % update with dummy values
        flagArray = [flagArray, [0;0;0]];
        ofArray = [ofArray, 0];
        thetaHatArray = [thetaHatArray, zeros(1,6)'];
        yEstArray = [yEstArray, zeros(1,6)'];
        yOptArray = [yOptArray, zeros(1,6)'];
        uOptArray = [uOptArray, zeros(1,3)'];
        uImpArray = [uImpArray, zeros(1,3)']; 
    end
    
end

%%
temp = [];
kFilt = 0.5;

for jj = 1:58
    if jj == 1
        t1 = diameterTemp{1, 1}(5,1);
        t2 = diameterTemp{1, 2}(5,1);
        t3 = diameterTemp{1, 3}(5,1);
    else
        t1 = (1 - kFilt)*diameterTemp{1, 1}(5,jj) + kFilt*t1; 
        t2 = (1 - kFilt)*diameterTemp{1, 2}(5,jj) + kFilt*t2; 
        t3 = (1 - kFilt)*diameterTemp{1, 3}(5,jj) + kFilt*t3; 
    end
    temp = [temp, [t1;t2;t3]];
    
end

for well = 1:3
   
    figure(5 + well) 
        %plot(2:60,diameterTemp{1, well}(1,:))
        hold on
        %plot(2:60,diameterTemp{1, well}(2,:))
        %plot(2:60,diameterTemp{1, well}(3,:))
        %plot(2:60,diameterTemp{1, well}(4,:))
        plot(2:60,diameterTemp{1, well}(5,:),'ko')
        plot(3:60,temp(well,:),'b')
     
end

% 
%         
% 
% 
% % %% Flowrate
% % % new plot
% % 
% % figure(2)
% % for exp = 1:size(expTable,1)
% % 
% %     expTime = size(data{expTable(exp,1),1}.flowrate,1);
% %        
% %     tEval = 1:sampTime:expTime;
% %         
% %     if rem(tEval(end),sampTime) ~= 0
% %         tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
% %     else
% %         tgrid = 0:(sampTime/60):(1/60)*expTime;
% %     end
% %     
% %     subplot(3,3,exp)
% %         plot(tgrid,data{expTable(exp,1),expTable(exp,2)}.flowrate(tEval),'ko-','Linewidth',1.5)
% % 
% %         
% %         xticks(0:5:tgrid(end))
% %         xlim([0, tgrid(end)])
% %         
% %         if exp < 4
% %             lab = ['Q = 13; Well = ',num2str(expTable(exp,2))];
% %         elseif exp < 7
% %             lab = ['Q = 9; Well = ',num2str(expTable(exp,2))];
% %         else
% %             lab = ['Q = 5; Well = ',num2str(expTable(exp,2))];
% %         end
% %                 
% %         title(lab)
% %         ylabel('Q [L/min]','FontSize',10)
% %         xlabel('time [min]','FontSize',10);
% % 
% % end
% 
% %% Delta Pressure
% figure(3)
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
%     Qtest = expTable(exp,3);
%     
%     if expTable(exp,2) == 1
%         DPmax = 0.3081*Qtest^2 + 0.185*Qtest - 0.7528;
%         DPmin = 0.1828*Qtest^2 + 1.0093*Qtest - 4.6835;
%     elseif expTable(exp,2) == 2
%         DPmax = 0.2934*Qtest^2 - 0.0592*Qtest - 0.7532;
%         DPmin = 0.1896*Qtest^2 + 0.8415*Qtest - 4.9335;    
%     else
%         DPmax = 0.2902*Qtest^2 + 0.9971*Qtest - 1.5097;
%         DPmin =  0.19*Qtest^2 + 1.1067*Qtest - 3.685;    
%     end
%     
%     filtTemp = medfilt1(data{expTable(exp,1),expTable(exp,2)}.deltaP,sampTime);
% 
%     
%     subplot(3,3,exp)
%         plot(tgrid,data{expTable(exp,1),expTable(exp,2)}.deltaP(tEval),'ko','Linewidth',1.5)
%         hold on
%         plot(tgrid,filtTemp(tEval),'b-','Linewidth',1.5)
%  
%         yline(DPmin,'k:');
%         yline(DPmax,'k:');
%         
%         xticks(0:5:tgrid(end))
%         xlim([0, tgrid(end)])
% 
%         %ylim([0.26, 0.33])
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
%         ylabel('DP [mbar]','FontSize',10)
%         xlabel('time [min]','FontSize',10);
% 
% end
% 
% %% Relative Delta Pressure
% figure(4)
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
%     Qtest = expTable(exp,3);
%     
%     if expTable(exp,2) == 1
%         DPmax = 0.3081*Qtest^2 + 0.185*Qtest - 0.7528;
%         DPmin = 0.1828*Qtest^2 + 1.0093*Qtest - 4.6835;
%     elseif expTable(exp,2) == 2
%         DPmax = 0.2934*Qtest^2 - 0.0592*Qtest - 0.7532;
%         DPmin = 0.1896*Qtest^2 + 0.8415*Qtest - 4.9335;    
%     else
%         DPmax = 0.2902*Qtest^2 + 0.9971*Qtest - 1.5097;
%         DPmin =  0.19*Qtest^2 + 1.1067*Qtest - 3.685;    
%     end
%     
%     DPNorm = (data{expTable(exp,1),expTable(exp,2)}.deltaP - DPmin)./(DPmax - DPmin);
%     filtTemp = medfilt1(DPNorm,sampTime);
%     
%     subplot(3,3,exp)
%         plot(tgrid,DPNorm(tEval),'ko','Linewidth',1.5)
%         hold on
%         plot(tgrid,filtTemp(tEval),'b-','Linewidth',1.5)
%         
%         hold on
%  
%         yline(1,'k:');
%         yline(0,'k:');
%         
%         xticks(0:5:tgrid(end))
%         xlim([0, tgrid(end)])
% 
%         %ylim([0.26, 0.33])
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
%         ylabel('DP [mbar]','FontSize',10)
%         xlabel('time [min]','FontSize',10);
% 
% end