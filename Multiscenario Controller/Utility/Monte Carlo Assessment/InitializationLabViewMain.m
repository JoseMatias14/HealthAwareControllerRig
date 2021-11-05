%clear
%clc

%% Standard Health-aware Controller
%initial condition
[dx0,z0,u0,theta0] = InitialConditionGasLift_NMPC(parPlant);

% Building Dynamic Model
[~,diff,alg,x_var,z_var,p_var] = BuildingDynModel_NMPC_prob(parPlant); % perfect model

% control tuning
nmpcConfig.umax = 0.95; %valve opening [-]
nmpcConfig.umin = 0.05;
nmpcConfig.dumax = 0.05; % 0.01 | 0.1
nmpcConfig.x_threshold = 0.316; % probe diameter [cm]
nmpcConfig.x_healthy = parPlant.dMin; % probe diameter [cm]
nmpcConfig.nx = size(dx0,1);
nmpcConfig.nz = size(z0,1);
nmpcConfig.nu = size(u0,1);
nmpcConfig.ntheta = size(theta0,1);

nmpcConfig.nm = 30; % 10 | 60 | 120 | 240 x sampling time = 60[s]
nmpcConfig.np = 10;

nmpcConfig.tol = 1e-4;% 1e-08
nmpcConfig.maxiter = 2000;
% % regularization
% nmpcConfig.R = 1*eye(nmpcConfig.nu); % 0 | 0.01 | 1  

% levels of the cfd used
probH = 0.9; % 0.75 | 0.95
probL = 0.50;

% testing the probabilities of different scenarios
%            w1 |  w2  | w3 
nmpcConfig.scendist =  [probH, probH, probH; %S1:   HHH
                        probH, probH, probL; %S2:   HHL
                        probH, probL, probH; %S3:   HLH
                        probL, probH, probH; %S4:   LHH
                        probH, probL, probL; %S5:   HLL
                        probL, probH, probL; %S6:   LHL
                        probL, probL, probH; %S7:   LLH
                        probL, probL, probL];%S8:   LLL

nmpcConfig.Nr = 2; 
nmpcConfig.Nlevels = size(nmpcConfig.scendist,1); 
nmpcConfig.Nscenarios = nmpcConfig.Nlevels^nmpcConfig.Nr;
nmpcConfig.ConsideredScen = [40,48,56,61,62,63,64];
%nmpcConfig.ConsideredScen = [16,24,32,37,38,39,40,45,46,47,48,53,54,55,56,58,59,60,61,62,63,64];

load('scenario_pruning_cumprob')
temp = scen_avg_prob(nmpcConfig.ConsideredScen);
temp2 = sum(temp);
nmpcConfig.scenProb = temp/temp2;


% Building NMPC
solverNMPC_MS = BuildingNMPC_MS_prob(diff,alg,x_var,z_var,p_var,parPlant,nmpcConfig);

