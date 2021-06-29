function [F,S_xx,S_zz,S_xz,S_xp,S_zp,x_var,z_var,u_var,p_var,diff,alg,L] = ErosionRigDynModel_model(par)
%    Simulates SS model of the rig with the data-driven reservoir
%    model 

% Inputs:
%    par = system parameters
%
% Outputs:
%   F: system integrator
%   S's: system sensitivities
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Julio Paez
% email: julio.paez.oliveira@usp.br 
% June 2020; Last revision: 
     
%addpath('<yourpath>/casadi-matlabR2014a-v3.5.1')
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


%% System states
%riser head total production rate
w_pr = MX.sym('w_pr',n_w);      % 1:3 [kg/s]

%oil rate from reservoir
w_ro = MX.sym('w_ro',n_w);      % 1:3 [kg/s]
%oil holdup 
m_o = MX.sym('m_or',n_w);       % 4:6 [kg]
%riser head pressure
p_rh = MX.sym('p_rh',n_w);      % 7:9 [bar]
%pressure - below injection point (bottom hole)
p_bh = MX.sym('p_bh',n_w);      % 10:12 [bar]

%% System input
%valve oppening
vo = MX.sym('vo',n_w);            % 1:3 [0-1]
%velocity pump
Ppump = MX.sym('Ppump',1);      % 4 [bar]

%% parameters
% fixed
%riser temperature
T_r = MX.sym('T_r',1); %[oC]
%separator pressure
p_s = MX.sym('p_s',1); %[bar]
%separator pressure
w_pr_ref = MX.sym('w_pr_ref',n_w); %[kg/s]
%separator pressure
w_ro_ref = MX.sym('w_ro_ref',n_w); %[kg/s]

%time transformation: CASADI integrates always from 0 to 1 and the USER does the time
%scaling with T.
T = MX.sym('T',1); %[s]

% estimable
%reservoir data-driven model
res_theta = MX.sym('res_theta',n_w);
%riser valve characteristics
val_theta = MX.sym('val_theta',n_w);%[m2]

%% Modeling
%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 

%oil from reservoir flowrate [kg/s]
f1 = -Ppump*ones(n_w,1)*1e5 + (w_ro.*1e-2).^2.*(res_theta.*1e9)./(vo.^2.*rho_o) + p_bh.*1e5 ; 
%riser head pressure [Pa]
f2 = -p_rh.*1e5 + (w_pr.*1e-2).^2.*(val_theta.*1e8)./rho_o + p_s.*1e5 ;
%bottom hole pressure [Pa]
f3 = -p_bh.*1e5 + (p_rh.*1e5 + rho_o.*9.81.*H_r + 128.*mu_oil.*(L_w+L_r).*(w_ro.*1e-2)./(3.14.*D_bh.^4.*rho_o));
%Oil Mass [kg]
f4 = -m_o + (rho_o).*((A_w.*L_w + A_r.*L_r));

%dynamic: change these to differential 
%liquid [kg] - modeled with a first order delay
df1 = 1e2*(-(w_pr - [11.4411829228142;9.71651938137975;10.8136777449560]).*1e-2 + (w_ro - [11.4320492854616;9.70749851497498;10.8045832186035]).*1e-2)./1; 
df1 = 1e2*(-(w_ro).*1e-2)./100; 

%df1 = 1e2*1e2*(-(w_pr.*1e-2) + (w_ro.*1e-2)); 
 

% Form the DAE system
diff = vertcat(df1);
alg = vertcat(f1,f2,f3,f4);

% give parameter values
alg = substitute(alg,p_s,par.p_s);
alg = substitute(alg,T_r,par.T_r);

% concatenate the differential and algebraic states
x_var = vertcat(w_pr);
z_var = vertcat(w_ro,m_o,p_rh,p_bh);
u_var = vertcat(vo,Ppump);
p_var = vertcat(res_theta,val_theta,T);

%objective function 
L = 10*((w_ro(1)*1e-2)*CR/rho_o(1)) + 20*((w_ro(2)*1e-2)*CR/rho_o(2)) + 30*((w_ro(3)*1e-2)*CR/rho_o(3));
%end modeling

%% Casadi commands
%declaring function in standard DAE form (scaled time)
dae = struct('x',x_var,'z',z_var,'p',vertcat(u_var,p_var),'ode',T*diff,'alg',alg);

%calling the integrator, the necessary inputs are: label; integrator; function with IO scheme of a DAE (formalized); struct (options)
F = integrator('F','idas',dae);

    % ================================================
    %       Calculating sensitivity matrix (theta)
    % ================================================

S_xx = F.factory('sensStaStates',{'x0','z0','p'},{'jac:xf:x0'});
S_zz = F.factory('sensStaStates',{'x0','z0','p'},{'jac:zf:z0'});
S_xz = F.factory('sensStaStates',{'x0','z0','p'},{'jac:xf:z0'});

S_xp = F.factory('sensParStates',{'x0','z0','p'},{'jac:xf:p'});
S_zp = F.factory('sensParStates',{'x0','z0','p'},{'jac:zf:p'});

end
