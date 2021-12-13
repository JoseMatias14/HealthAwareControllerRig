function [thetaHat,zHat,Psi,flagEst] = HAC_SSEstimation(zGuess,thetaGuess,uk,yk,par)

%    Estimates Parameters of the steady-state rig model 

% Inputs:
%    xGuess,zGuess = guess for the states 
%    thetaGuess = parameters guess
%    uk = current inputs (measured)
%    yk = current measurements
%    diff,alg,x_var,z_var,p_var = casadi model
%    par = system parameters
%
% Outputs:
%   thetaHat = estimated parameters
%   xHat, yHat = estimated states and model prediction
%   dpHat = delta pressure estimate
%   Psi = estimation OF value
%   flagEst  = estimation didn't converge == 0 | success == 1

% Other m-files required: OptimizationBoundsHAC
% MAT-files required: none

%addpath ('C:\Users\lab\Documents\casadi-windows-matlabR2016a-v3.4.5')                
import casadi.*

%% Model
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
f3 = -p_bh.*1e5 + (yk(7:9).*1e2 + p_rh.*1e5 + rho_o.*9.81.*H_r + 128.*mu_oil.*(L_w+L_r).*(w_ro.*1e-2)./(3.14.*D_bh.^4.*rho_o));
%Oil Mass [kg]
f4 = -m_o + (rho_o).*((A_w.*L_w + A_r.*L_r));
% mass conservation
f5 = w_ro - w_pr;

alg = vertcat(f1,f2,f3,f4,f5);

% give parameter values
alg = substitute(alg,p_s,par.p_s);
alg = substitute(alg,T_r,par.T_r);

% concatenate the differential and algebraic states
x_var = vertcat(m_o,w_ro,w_pr,p_rh,p_bh);
p_var = vertcat(vo,Ppump,res_theta,val_theta);

%% Estimation problem (part 1)
% ===================================
%     Declaring 
% ===================================
% decision variables
w = {x_var,p_var};
% constraints (model) -- differential equations dont enter here! Integrator
g = {alg};

%objective function
% computing residual
yModel = [CR*(w_ro{1}*1e-2)./rho_o(1),CR*(w_ro{2}*1e-2)./rho_o(2),CR*(w_ro{3}*1e-2)./rho_o(3),p_rh{1},p_rh{2},p_rh{3}];
eps = (yModel' - yk(1:6));

J = eps'*par.Sigma(1:6,1:6)*eps;
%J = eps'*eps;

% formalize it into an NLP problem
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

% Assign solver
options = struct;
options.ipopt.print_level = 5;
% options.ipopt.acceptable_constr_viol_tol = 0.01;
Est = nlpsol('solver','ipopt',nlp,options);

% ===================================
%     Solving
% ===================================
%Variable bounds
[lbx,lbz,~,lbtheta,ubx,ubz,~,ubtheta] = OptimizationBoundsHAC(par);

% x_var,z_var,p_var -- > (vo,Ppump,res_theta,val_theta);
wk = [zGuess(1:15);uk(1:3);uk(4);thetaGuess];
lbw =   [lbz(1:15);uk(1:3);uk(4);lbtheta];
ubw =   [ubz(1:15);uk(1:3);uk(4);ubtheta];

lbg = zeros(15,1); %SS: diff and alg == 0
ubg = zeros(15,1);

% Solve
sol = Est('x0',wk,'lbg',lbg,'ubg',ubg,'lbx',lbw,'ubx',ubw);

flagEst = Est.stats.success;
% if Est.stats.success ~=1
%     msg = 'Error estimating parameters, model';
%     error(msg);
% end

%clc %clean ipopt output

% Extract Solution
zHat = full(sol.x(1:15));
%yHat = par.H*zHat;%conversion [kg/s] --> [L/min]
thetaHat = full(sol.x(15 + length(uk) + 1:end));
Psi = full(sol.f);

% %% diameter estimation
% % orifice pressure drop 
% f5 = dP*1e2 - Cd.*(w_ro.*1e-2).^2.*(1 - (d/2).^4)./((pi()/4)^2*2*rho_o.*(0.01*d).^4);


end
