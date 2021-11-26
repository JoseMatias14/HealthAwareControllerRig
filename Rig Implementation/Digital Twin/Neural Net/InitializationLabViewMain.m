%clear
%clc

%% Model tuning
%paramters
par = ParametersGasLiftModel;
[x0,u0,theta0] = InitialCondition_HAC(par);

%states to measurement mapping function
par.nMeas = 6;
H = zeros(par.nMeas,length(x0));
H(1,4) = 1e-2*60*1e3/par.rho_o(1); %wro-oil rate from reservoir, well 1 [1e2 kg/s] --> [L/min]
H(2,5) = 1e-2*60*1e3/par.rho_o(2); %wro-oil rate from reservoir, well 2
H(3,6) = 1e-2*60*1e3/par.rho_o(3); %wro-oil rate from reservoir, well 3
H(4,10) = 1; %prh - riser head pressure well 1 [bar]
H(5,11) = 1; %prh - riser head pressure well 2
H(6,12) = 1; %prh - riser head pressure well 2
par.H = H;

%% Standard Health-aware Controller
% control tuning
nmpcConfig.umax = 0.75; % 0.95 valve opening [-]
nmpcConfig.umin = 0.25; % 0.05
nmpcConfig.dumax = 0.01; % 0.01 | 0.1

% degradation indicator [normalized delta pressure]
nmpcConfig.x_threshold = 0.8; % degraded fo4r the controller (upper bound)
nmpcConfig.x_healthy = 0; % healthy (lower bound)
nmpcConfig.x_broken = 1; % completely broken (unknown top the controller)

nmpcConfig.nx = size(x0,1);
nmpcConfig.nu = size(uPlant0(1:3),1);

%sampling time = 60[s]
nmpcConfig.st = 60;

nmpcConfig.nm = 30; % 60 | 120 | 240 x sampling time = 60[s]
nmpcConfig.np = 30; % 10 | nmpcConfig.nm
% % regularization
% nmpcConfig.R = 1*eye(nmpcConfig.nu); % 0 | 0.01 | 1  




