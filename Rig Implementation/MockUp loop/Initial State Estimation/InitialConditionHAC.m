%Initial Condition Simplified
function [x0,z0,u0,theta0] = InitialConditionHAC(par)
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

%% Differential states (x)
% initial diameter
x0 = 0.3*ones(par.n_w,1);

%% Algebraic states (z)
Ql = temp.w_ro_0*1e-2*60*1e3./par.rho_o;
dP0 = 0.2788*Ql.^2 + 1.143*Ql - 2.3831; % dP model previosly calculated
z0 = vertcat(temp.m_o_0,temp.w_ro_0,temp.w_pr_0,temp.p_rh_0,temp.p_bh_0,dP0);

%% Inputs
u0 = vertcat(temp.vo_0,temp.vpump_0);

%% Parameters
theta0  = vertcat(temp.res_theta_0,temp.val_theta_0);