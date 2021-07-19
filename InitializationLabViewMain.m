%clear
%clc

%% Standard Health-aware Controller
% Building Dynamic Model
[diff,alg,x_var,z_var,p_var] = BuildingDynModel(parPlant); % perfect model

% control tuning
nmpcConfig.umax = 0.6; %valve opening [-]
nmpcConfig.umin = 0.25;
nmpcConfig.dumax = 0.05; % 0.01 | 0.1
nmpcConfig.x_threshold = 0.316; % probe diameter [cm]
nmpcConfig.nx = size(dk,1);
nmpcConfig.nz = size([xPlant0;zPlant0],1);
nmpcConfig.nu = size(uPlant0(1:3),1);

nmpcConfig.nm = 120; % x sampling time = 60[s]
nmpcConfig.np = 120;
% regularization
nmpcConfig.R = 1*eye(nmpcConfig.nu); % 0 | 0.01 | 1  

% slack
nmpcConfig.rho = 1e5*eye(nmpcConfig.nx); % 1e3 | 999999

% Building NMPC
solverNMPC = BuildingNMPC(diff,alg,x_var,z_var,p_var,parPlant,nmpcConfig);
