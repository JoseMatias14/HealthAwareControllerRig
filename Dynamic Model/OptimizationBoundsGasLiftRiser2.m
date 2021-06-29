function [lbx,lbz,lbu,ubx,ubz,ubu,lbp,ubp] = OptimizationBoundsGasLiftRiser2(par)

%number of wells
n_w = par.n_w;

%% States
%oil holdup [kg]
m_o_lb = 0.*ones(n_w,1);
m_o_ub = 1e3.*ones(n_w,1);

%oil rate from reservoir [1e2 kg/s]
w_ro_lb = 0.*ones(n_w,1);
w_ro_ub = 20.*ones(n_w,1);

% total production rate [1e2 kg/s]
w_pr_lb = 0.*ones(n_w,1);
w_pr_ub = 60.*ones(n_w,1); 

%riser head pressure [bar]
p_rh_lb = 0.*ones(n_w,1);
p_rh_ub = 1.500.*ones(n_w,1); %49

%pressure - below injection point (bottom hole) [bar]
p_bh_lb = 0.*ones(n_w,1);
p_bh_ub = 2.*ones(n_w,1);

%% Inputs
% valve opening [0-1]
vo_lb = 0.1*ones(n_w,1);
vo_ub = 0.9*ones(n_w,1); 

% pump rotation [bar]
Ppump_lb = 1.01325*ones(1,1); %atmosferic pressure
Ppump_ub = 2*1.01325.*ones(1,1); 

%% Parameters
%Theta
res_theta_lb= 0.*(ones(n_w,1));
res_theta_ub= 10*(ones(n_w,1));

%Riser valve caracteristics
val_theta_lb=0.*ones(n_w,1);
val_theta_ub=10.*ones(n_w,1);

%% 
lbx = vertcat(m_o_lb);
lbz = vertcat(w_ro_lb,w_pr_lb,p_rh_lb,p_bh_lb);
lbu = vertcat(vo_lb,Ppump_lb);
lbp = vertcat(res_theta_lb,val_theta_lb);

ubx = vertcat(m_o_ub);
ubz = vertcat(w_ro_ub,w_pr_ub,p_rh_ub,p_bh_ub);
ubu = vertcat(vo_ub,Ppump_ub);
ubp = vertcat(res_theta_ub,val_theta_ub);