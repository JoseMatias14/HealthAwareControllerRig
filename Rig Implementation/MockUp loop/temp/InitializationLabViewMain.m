%clear
%clc

%% Standard Health-aware Controller
%plant parameters
par = ParametersGasLiftModel;
par.T = 60; % simulation sampling time[s]

%states to measurement mapping function
par.nMeas = 9;
par.H = zeros(par.nMeas,18);
par.H(1,1) = 1e-2*60*1e3/par.rho_o(1); %wro-oil rate from reservoir, well 1 [1e-2 kg/s] --> [L/min]
par.H(2,2) = 1e-2*60*1e3/par.rho_o(2); %wro-oil rate from reservoir, well 2
par.H(3,3) = 1e-2*60*1e3/par.rho_o(3); %wro-oil rate from reservoir, well 3
par.H(4,7) = 1; %prh - riser head pressure well 1
par.H(5,8) = 1; %prh - riser head pressure well 2
par.H(6,9) = 1; %prh - riser head pressure well 3
par.H(7,13) = 1; %dp - well 1
par.H(8,14) = 1; %dp - well 2
par.H(9,15) = 1; %dp - well 3

par.m0 = [5.05;5.60;5.55]; % [g]

% control tuning
nmpcConfig.umax = 0.95; %valve opening [-]
nmpcConfig.umin = 0.05;
nmpcConfig.dumax = 0.05; % 0.01 | 0.1
nmpcConfig.x_threshold = 0.316; % probe diameter [cm]
nmpcConfig.x_healthy = par.dMin; % probe diameter [cm]
nmpcConfig.nx = 3;
nmpcConfig.nz = 18;
nmpcConfig.nu = 3;

nmpcConfig.nm = 30; % 60 | 120 | 240 x sampling time = 60[s]
nmpcConfig.np = 10; % 20

nmpcConfig.tol = 1e-4;% 1e-08
nmpcConfig.maxiter = 2000;
% % regularization
% nmpcConfig.R = 1*eye(nmpcConfig.nu); % 0 | 0.01 | 1  


