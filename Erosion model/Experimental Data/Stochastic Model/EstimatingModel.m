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
           

%% Diameter - estimating CV
% initial time for the "CV estimation" visually determined
timeTable = [5, 10]*60;

% Modeling new probe (DPmax)
%%%%%%%%%%%%%%%%%%%%%
% Model with Casadi %
%%%%%%%%%%%%%%%%%%%%%
import casadi.*

% Declare variables to the model (parameters and controls)
Q_var = MX.sym('Q_var');
Cd = MX.sym('Cd');
d = MX.sym('d');
rho = 1000;
D = 2; %[cm]
d0 = 0.3; %[cm] 

% model
DP_var = Cd*((Q_var/(60*1000))*rho)^2*(1 - (d/D)^4)/((pi()/4)^2*2*rho*(0.01*d)^4)/100;

% Create a function that computes the delta pressure
dP_fun = Function('dP_fun',{Q_var, Cd, d},{DP_var});

%%%%%%%
% Map %
%%%%%%%
%simulating the system for all the inputs at the same time
N = size(length(timeTable(1):timeTable(2)),1);%number of points used in the estimation procedure - based on timeTable
all_samples = dP_fun.map(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving the estimated CV %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
CdhatArray = [];

for exp = 1:size(expTable,1)
        % extracting the data
        if expTable(exp,2) == 1
            dP_kk = data{expTable(exp,1),expTable(exp,2)}.deltaP(timeTable) - 1.5 + 0.7;
        elseif expTable(exp,2) == 2
            dP_kk = data{expTable(exp,1),expTable(exp,2)}.deltaP(timeTable) + 1.7 + 0.24;
        else %expTable(exp,3) == 3
            dP_kk = data{expTable(exp,1),expTable(exp,2)}.deltaP(timeTable) + 0 + 0.20;
        end
        
        Q_kk = data{expTable(exp,1),expTable(exp,2)}.flowrate(timeTable);
        dP_symbolic = all_samples(Q_kk, repmat(Cd,1,N), repmat(d0,1,N));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimating model parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        e = dP_kk - dP_symbolic;
        
        %formalizing nlp structure and nlp problem
        % Create an NLP solver
        opts2 = struct;
        %opts2.ipopt.max_iter = nmpcPar.maxiter;
        opts2.ipopt.print_level = 0;
        %opts2.ipopt.tol = nmpcPar.tol;
        %opts2.ipopt.acceptable_tol = 100*nmpcPar.tol; % optimality convergence tolerance
        %opts2.ipopt.linear_solver = 'mumps';
        nlp = struct('x', Cd, 'f', 0.5*dot(e,e));
        solver = nlpsol('solver','ipopt', nlp,opts2);
                
        sol = solver('x0',0.7,'lbx',0,'ubx',10);
        Cd_temp = full(sol.x);
        CdhatArray = [CdhatArray; Cd_temp];
        data{expTable(exp,1),expTable(exp,2)}.Cd_temp = Cd_temp;
        
        
        expTime = length(data{expTable(exp,1),expTable(exp,2)}.deltaP);
        for kk = 1:expTime
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimating model parameters %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dP_sym = dP_fun(data{expTable(exp,1),expTable(exp,2)}.flowrate(kk),Cd_temp,d);
            
            if expTable(exp,2) == 1
                temp2 = data{expTable(exp,1),expTable(exp,2)}.deltaP(kk) - 1.5 + 0.7;
            elseif expTable(exp,2) == 2
                temp2 = data{expTable(exp,1),expTable(exp,2)}.deltaP(kk) + 1.7 + 0.24;
            else %expTable(exp,2) == 3
                temp2 = data{expTable(exp,1),expTable(exp,2)}.deltaP(kk) + 0 + 0.20;
            end
            e = temp2 - dP_sym;
            
            nlp2 = struct('x', d, 'f', 0.5*dot(e,e));
            solver = nlpsol('solver','ipopt', nlp2,opts2);
            
            % Solve
            if kk == 1
                xGuess = 0.3;
            else
                xGuess = temp;
            end
            
            sol = solver('x0',xGuess,'lbx',par.dMin,'ubx',par.dMax);
            temp = full(sol.x);
            data{expTable(exp,1),expTable(exp,2)}.equivDiameCvNew(kk) = temp;
            
        end
        
        
end

%% Diameter
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
    filtTemp2 = medfilt1(data{expTable(exp,1),expTable(exp,2)}.equivDiameCvNew,sampTime);

    subplot(3,3,exp)
        plot(tgrid,data{expTable(exp,1),expTable(exp,2)}.equivDiameCvNew(tEval),'ko','Linewidth',1.5)
        hold on
        plot(tgrid,filtTemp(tEval),'r-','Linewidth',1.5)
        plot(tgrid,filtTemp2(tEval),'b-','Linewidth',1.5)

        yline(0.3,'k:');
        yline(0.3181,'k:');
        
        xticks(0:5:60)
        xlim([0, 60])

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

%% Diameter - estimating the noise
% initial time for the "noise estimation" visually determined
timeTable2 = [35, 55;
            25, 40;
            35, 55;
            35, 55;
            25, 40;
            25, 40;
            35, 55;
            35, 55;
            35, 55];
        
% noise standard deviation
meanNoise = [];
stdNoise = [];

figure(3)
for exp = 1:size(expTable,1)
   
    tEval = timeTable2(exp,1)*60:sampTime:timeTable2(exp,2)*60;
    tgrid = timeTable2(exp,1):(sampTime/60):timeTable2(exp,2);
    
    filtTemp = medfilt1(data{expTable(exp,1),expTable(exp,2)}.equivDiameCvNew,sampTime);
    temp = data{expTable(exp,1),expTable(exp,2)}.equivDiameCvNew(tEval) - filtTemp(tEval);
    stdNoise = [stdNoise, std(temp)];
    meanNoise = [meanNoise, mean(temp)];
        
    data{expTable(exp,1),expTable(exp,2)}.stdNoise = stdNoise(end);
    data{expTable(exp,1),expTable(exp,2)}.meanNoise = meanNoise(end);
    
    subplot(3,3,exp)
        plot(tgrid,temp,'ko','Linewidth',1)
        hold on
        yline(meanNoise(end),'b:','Linewidth',1.5);
        yline(meanNoise(end) - stdNoise(end),'b-.','Linewidth',1.5);
        yline(meanNoise(end) + stdNoise(end),'b-.','Linewidth',1.5);
        %yline(0,'k:');
        
        xticks(0:5:60)
        xlim([0, 60])
        ylim([-0.03, 0.03])

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

save('Experiments_newCV','data')


