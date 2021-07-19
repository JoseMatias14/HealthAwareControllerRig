function [diff,alg,x_var,z_var,p_var] = BuildingDynModel(par)
%    Build model for the MPC - deterministic degradation using average of
%    the estimated remaining useful life)

% Inputs:
%    par = system parameters
%
% Outputs:
%   diff,alg = model differential and algebraic equations 
%   x_var,z_var,p_var = model states and parameters

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Matias
% email: jose.o.a.matias@ntnu.no 
% July 2021; Last revision: 
        
    import casadi.*

%% Parameters
%number of wells
n_w = par.n_w; %[]
%gas constant
R = par.R; %[m3 Pa K^?1 mol^?1]
%molecular weigth
Mw = par.Mw; %[kg/mol?]

%properties
%density of oil - dim:  nwells x 1
rho_o = par.rho_o; %[kg/m3]
%1cP oil viscosity
mu_oil = par.mu_oil;% [Pa s] 
%density of gas
rho_g = par.rho_g; %[kg/m3]

%project
%well parameters - dim:  nwells x 1
L_w = par.L_w; %[m]
H_w = par.H_w; %[m]
D_w = par.D_w; %[m]
A_w = par.A_w;%[m2]

%well below injection - [m]
L_bh = par.L_bh;
H_bh = par.H_bh;
D_bh = par.D_bh;
A_bh = par.A_bh;%[m2]

%riser - [m]
L_r = par.L_r;
H_r = par.H_r;
D_r = par.D_r;
A_r = par.A_r;%[m2]

Cd = par.Cd_hat;

%% System states

%orifice diameter
d = MX.sym('d',n_w); %[cm]      % 1:3 [cm]

%oil holdup 
m_o = MX.sym('m_or',n_w);       % 1:3 [kg]
%oil rate from reservoir
w_ro = MX.sym('w_ro',n_w);      % 4:6 [1e2 kg/s]
%riser head total production rate
w_pr = MX.sym('w_pr',n_w);      % 7:9 [1e2 kg/s]

%riser head pressure
p_rh = MX.sym('p_rh',n_w);      % 10:12 [bar]
%pressure - below injection point (bottom hole)
p_bh = MX.sym('p_bh',n_w);      % 13:15 [bar]
%pressure - below injection point (bottom hole)
dP = MX.sym('dP',n_w);      % 16:18 [mbar]

%% System input
%valve oppening
vo=MX.sym('vo',n_w);            % 1:3 [0-1]
%Pump outlet pressure
Ppump=MX.sym('Ppump',1);      % 4 [bar g]

%% parameters
% fixed
%riser temperature
T_r = MX.sym('T_r',1); %[oC]
%separator pressure
p_s = MX.sym('p_s',1); %[bar]


% estimable
%reservoir data-driven model
res_theta = MX.sym('res_theta',n_w);%[m2]
%riser valve characteristics
val_theta = MX.sym('val_theta',n_w);%[m2]

%% Modeling
%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 

%oil from reservoir flowrate [kg/s] - same pump outlet pressure
% original equation
%f1 = -Ppump*ones(n_w,1)*1e5 + (w_ro.*1e-2).^2./((1e-5*vo.*res_theta).^2.*(rho_o)) + p_bh.*1e5 ; 
% with re-arranged parameters
f1 = -Ppump*ones(n_w,1)*1e5 + (w_ro.*1e-2).^2.*(res_theta.*1e9)./(vo.^2.*rho_o) + p_bh.*1e5 ; 
%riser head pressure [Pa]
%original equation
%f2 = -p_rh.*1e5 + (w_pr.*1e-2).^2./((1e-4*val_theta).^2.*(rho_r.*1e2)) + p_s.*1e5 ;
f2 = -p_rh.*1e5 + (w_pr.*1e-2).^2.*(val_theta.*1e8)./rho_o + p_s.*1e5 ;
%bottom hole pressure [Pa]
f3 = -p_bh.*1e5 + (dP*1e2 + p_rh.*1e5 + rho_o.*9.81.*H_r + 128.*mu_oil.*(L_w+L_r).*(w_ro.*1e-2)./(3.14.*D_bh.^4.*rho_o));
%Oil Mass [kg]
f4 = -m_o + (rho_o).*((A_w.*L_w + A_r.*L_r));
% orifice pressure drop 
f5 = dP*1e2 - Cd*(w_ro.*1e-2).^2.*(1 - (d/2).^4)./((pi()/4)^2*2*rho_o.*(0.01*d).^4);
% mass conservation
f6 = w_ro - w_pr;

% differential equations - (all symbolic) - [cm/min]
% from ...\Erosion model\gammaIncrementEstimation.m
QQ = CR*(w_ro.*1e-2)./rho_o; % kg/s --> L/min
dd_dt = (0.3181 - 0.3)./(-3*QQ + 51.333); %[cm/min]
df1 = dd_dt/60; %[cm/min] --> [cm/s]

% Form the DAE system
diff = vertcat(df1);
alg = vertcat(f1,f2,f3,f4,f5,f6);

% give parameter values
alg = substitute(alg,p_s,par.p_s);
alg = substitute(alg,T_r,par.T_r);

% concatenate the differential and algebraic states
x_var = vertcat(d);
z_var = vertcat(m_o,w_ro,w_pr,p_rh,p_bh,dP);
p_var = vertcat(vo,Ppump,res_theta,val_theta);
    
end

