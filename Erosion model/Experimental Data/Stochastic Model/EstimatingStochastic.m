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

% minimum equivalent diameter
par.dMin = 0.3;
par.dMax = 0.312;

par.dLim = (0.312 - 0.3)*0.2 + 0.3;

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
load('Experiments_newCV')

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
figure(1)

% finding the breakdown time
tBreak_0 = [];
tBreak_end = [];

countBreak = 0;

cdArray = [];

for expe = 1:size(expTable,1)

    expTime = size(data{expTable(expe,1),1}.flowrate,1);
    cdArray = [cdArray, data{expTable(expe,1),expTable(expe,2)}.Cd_temp];
   
    tEval = 1:sampTime:expTime;
        
    if rem(tEval(end),sampTime) ~= 0
        tgrid = 0:(sampTime/60):(1/60)*(expTime - 1);
    else
        tgrid = 0:(sampTime/60):(1/60)*expTime;
    end
    
    filtTemp2 = medfilt1(data{expTable(expe,1),expTable(expe,2)}.equivDiameCvNew,sampTime);
 
    for kk = 1:length(filtTemp2)
        if filtTemp2(kk) > par.dLim && countBreak == 0
            tBreak_0 = [tBreak_0, kk];
            countBreak = 1;
        end
        if filtTemp2(kk) > par.dMax && countBreak == 1
            tBreak_end = [tBreak_end, kk];
            countBreak = 2;
        end
    end
    countBreak = 0; 
    
    subplot(3,3,expe)
        plot(tgrid,data{expTable(expe,1),expTable(expe,2)}.equivDiameCvNew(tEval),'ko','Linewidth',1.5)
        hold on
        plot(tgrid,filtTemp2(tEval),'b-','Linewidth',1.5)

        yline(par.dMin,'k:');
        yline(par.dMax,'k:');
        
        xline(tBreak_end(end)/60,'k:');
        xline(tBreak_0(end)/60,'k:');
        
        xticks(0:5:60)
        xlim([0, 60])

        %             yticks(0:0.25:1)
        ylim([0.28, 0.33])

        if expe < 4
            lab = ['Q = 13; Well = ',num2str(expTable(expe,2))];
        elseif expe < 7
            lab = ['Q = 9; Well = ',num2str(expTable(expe,2))];
        else
            lab = ['Q = 5; Well = ',num2str(expTable(expe,2))];
        end
        
        title(lab)
        ylabel('d [cm]','FontSize',10)
        xlabel('time [min]','FontSize',10);

end

cdArray
