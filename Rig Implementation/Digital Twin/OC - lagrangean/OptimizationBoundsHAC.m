function [lbx,lbz,lbu,lbtheta,ubx,ubz,ubu,ubtheta] = OptimizationBoundsHAC(par)
%   Defines the lower and upper bounds for the model variables

% Inputs:
%    par = system parameters
%
% Outputs:
%   bounds on differential states (x)
%             algebraic states (z)
%             inputs (u)
%             parameters (p)

% Other m-files required: none
% MAT-files required: none
% concatenate the differential and algebraic states
    
%number of wells
n_w = par.n_w;

%% Differential States (x)
%diameter [cm]
d_lb = par.dMin*ones(n_w,1);
d_ub = par.dMax*ones(n_w,1);

%% Algebraic States (z)
%liquid holdup [kg]
m_o_lb = 0.*ones(n_w,1);
m_o_ub = 1e3.*ones(n_w,1);

%water rate from reservoir [1e-2 kg/s]
w_ro_lb = 0.*ones(n_w,1);
w_ro_ub = 20.*ones(n_w,1);

% total production rate [1e-2 kg/s]
w_pr_lb = 0.*ones(n_w,1);
w_pr_ub = 60.*ones(n_w,1); 

%riser head pressure [bar]
p_rh_lb = 0.*ones(n_w,1);
p_rh_ub = 1.500.*ones(n_w,1); %49

%pressure - before injection point (bottom hole) [bar]
p_bh_lb = 0.*ones(n_w,1);
p_bh_ub = 2.*ones(n_w,1);

%pressure - before injection point (bottom hole) [bar]
dP_lb = ones(n_w,1);
dP_ub = 100*ones(n_w,1);

%% Inputs
% valve opening [0-1]
vo_lb = 0.1*ones(n_w,1);
vo_ub = 0.995*ones(n_w,1); 

%% Parameters
%Theta
res_theta_lb= 0.*(ones(n_w,1));
res_theta_ub= 10.*(ones(n_w,1));

%Riser valve caracteristics
top_theta_lb=0.*ones(n_w,1);
top_theta_ub=10.*ones(n_w,1);

%% 
lbx = vertcat(d_lb);
lbz = vertcat(m_o_lb,w_ro_lb,w_pr_lb,p_rh_lb,p_bh_lb,dP_lb);
lbu = vertcat(vo_lb);
lbtheta = vertcat(res_theta_lb,top_theta_lb);

ubx = vertcat(d_ub);
ubz = vertcat(m_o_ub,w_ro_ub,w_pr_ub,p_rh_ub,p_bh_ub,dP_ub);
ubu = vertcat(vo_ub);
ubtheta = vertcat(res_theta_ub,top_theta_ub);