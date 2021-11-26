function [thetaHat,dpHat,xHat,yHat,Psi,flagEst] = ErosionRigSSEstimation(xGuess,thetaGuess,uk,yk,tk,Reg_mu,Reg_std,Resp_mu,Resp_std,eq,x_var,p_var,par)

%    Estimates Parameters of the steady-state rig model 

% Inputs:
%    xGuess = guess for the states 
%    thetaGuess = parameters guess
%    uk = current inputs (measured)
%    yk = current measurements
%    eq,x_var,p_var = casadi model
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

%% Estimation problem
% ===================================
%     Declaring 
% ===================================
% decision variables
w = {x_var,p_var};
% constraints (model)
g = eq;

%objective function
% computing residual
eps = (par.H*x_var - yk);

J = eps'*par.Sigma*eps;
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
[lbx,~,lbtheta,ubx,~,ubtheta] = OptimizationBoundsHAC(par);

% p_var = vertcat(vo,tCurr,Ppump,T_s,res_theta,val_theta,mu_rg,std_rg,mu_rp,std_rp);
wk = [xGuess;uk(1:3);tk;uk(4);uk(5:7);thetaGuess;Reg_mu;Reg_std;Resp_mu;Resp_std];
lbw =   [lbx;uk(1:3);tk;uk(4);uk(5:7);lbtheta;Reg_mu;Reg_std;Resp_mu;Resp_std];
ubw =   [ubx;uk(1:3);tk;uk(4);uk(5:7);ubtheta;Reg_mu;Reg_std;Resp_mu;Resp_std];

lbg = zeros(length(xGuess),1); %SS: diff and alg == 0
ubg = zeros(length(xGuess),1);

% Solve
sol = Est('x0',wk,'lbg',lbg,'ubg',ubg,'lbx',lbw,'ubx',ubw);

flagEst = Est.stats.success;
% if Est.stats.success ~=1
%     msg = 'Error estimating parameters, model';
%     error(msg);
% end

%clc %clean ipopt output

% Extract Solution
xHat = full(sol.x(1:length(xGuess)));
dpHat = xHat(16:18);
yHat = par.H*xHat;%conversion [kg/s] --> [L/min]
thetaHat = full(sol.x(length(xGuess) + length(uk) + 1 + 1:length(xGuess) + length(uk) + 1 + 6));
Psi = full(sol.f);

end
