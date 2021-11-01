%Initial Condition Simplified
function [x0,z0,u0,theta0] = InitialConditionGasLift_NMPC(par)
%   Defines the initial state
%   pre-computed from experimental data using a SS estimation routine

% Inputs:
%    par = system parameters
%
% Outputs:
%   initial value of differential states (x)
%                    algebraic states (z)
%                    inputs (u)
%                    parameters (p)

% Other m-files required: none
% MAT-files required: InitialState_2021-03-03_104459_NoRTO_test2_1.mat

% previously computed
temp = load('InitialState_2021-03-03_104459_NoRTO_test2_1');

x0 = 0.3*ones(3,1);
Ql = temp.w_ro_0*1e-2*60*1e3./par.rho_o;
dP0 = 0.2788*Ql.^2 + 1.143*Ql - 2.3831; % dP model previosly calculated
dd0 = 0.1*ones(3,1); % initial increment
prob0 = 0.05*ones(3,1); % initial increment

z0 = vertcat(temp.m_o_0,temp.w_ro_0,temp.w_pr_0,temp.p_rh_0,temp.p_bh_0,dP0,dd0,prob0);
u0 = vertcat(temp.vo_0);

%temp2 = (temp.w_ro_0.*1e-2).^2.*(1 - (0.01*0.3./par.D_w).^4) - 2*(dP0*1e2).*par.Cd_hat^2*(pi()/4)^2.*(0.001*0.3).^4.*par.rho_o;
probCum0 = [0.75;0.5;0.25];
theta0  = vertcat(temp.vpump_0,temp.res_theta_0,temp.val_theta_0,probCum0);