function  par = ParametersGasLiftModel_model

%number of wells
par.n_w = 3;
%gas constant
par.R = 8.314; %[m3 Pa/(K mol)]
%molecular weigth
par.Mw = 0.029; %[kg/mol] -- Attention: this unit is not usual

%% Properties
%density of oil - dim:  nwells x 1
par.rho_o = 996.57*ones(par.n_w,1); %[kg/m3] - water
%1cP oil viscosity
par.mu_oil = 1*0.000853*ones(par.n_w,1); %[Pa s or kg/(m s)]  - water
%Density Gas
par.rho_g = 1.2041*ones(par.n_w,1); %[kg/m3]
%riser temperature
par.T_r = 23+273; %[K]
%separator pressure
par.p_s = 1.01; %[bar]

%% Project
%well parameters - dim:  nwells x 1
%length
par.L_w = 1.8*ones(par.n_w,1); %[m]
%height
par.H_w  = 0*ones(par.n_w,1); %[m]
%diameter
par.D_w = 0.02*ones(par.n_w,1); %[m]
%well transversal area
par.A_w = pi.*(par.D_w/2).^2;%[m2]

%well below injection - [m]
par.L_bh = 0.4*ones(par.n_w,1);
par.H_bh = 0*ones(par.n_w,1);
par.D_bh = 0.02*ones(par.n_w,1);
par.A_bh = pi.*(par.D_bh/2).^2;%[m2]

%riser - [m]
par.L_r = 2.2*ones(par.n_w,1);
par.H_r = 2.2*ones(par.n_w,1);
par.D_r = 0.02*ones(par.n_w,1);
%riser areas
par.A_r = pi.*(par.D_r/2).^2;%[m2]

%% Upper limits 
%Max gas lift flow 
% par.qGLMax = 4; % [L/min] 
%Max wellhead gas production rate 
par.QgMax = 7.5; % [L/min] 

%% Dynamic model
par.T = 10; %[s]


