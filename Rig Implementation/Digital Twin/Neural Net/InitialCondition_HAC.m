%Initial Condition Simplified
function [x0,u0,theta0] = InitialCondition_HAC(par)
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

% concatenate the differential and algebraic states
% x_var = vertcat(m_o,w_ro,w_pr,p_rh,p_bh,dP);
% p_var = vertcat(vo,Ppump,T_s,tCurr,res_theta,val_theta,mu_rg,std_rg,mu_rp,std_rp);


%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 
Q_pr = CR*(temp.w_pr_0*1e-2)./par.rho_o;

DP_0 = [];
for well = 1:3
    if well == 1
        DPmin = 0.1828*Q_pr(well)^2 + 1.0093*Q_pr(well) - 4.6835;
    elseif well == 2
        DPmin = 0.1896*Q_pr(well)^2 + 0.8415*Q_pr(well) - 4.9335;    
    else % well == 3
        DPmin =  0.19*Q_pr(well)^2 + 1.1067*Q_pr(well)- 3.685;    
    end

    DP_0 = [DP_0; DPmin];
end



x0 = vertcat(temp.m_o_0,temp.w_ro_0,temp.w_pr_0,temp.p_rh_0,temp.p_bh_0,DP_0);
u0 = vertcat(temp.vo_0,temp.vpump_0,[25;25;25],0);
theta0 = vertcat(temp.res_theta_0,temp.val_theta_0);
