function dHat = DiameterEstimation(dguess,Qplant,DPplant,par)

% Modeling new probe
%%%%%%%%%%%%%%%%%%%%%
% Model with Casadi %
%%%%%%%%%%%%%%%%%%%%%
import casadi.*

% Declare variables to the model (parameters and controls)
Q_var = MX.sym('Q_var');
Cd = MX.sym('Cd');
d = MX.sym('d');
rho = 1000;
D = 2; %[cm]

% model
DP_var = Cd*((Q_var/(60*1000))*rho)^2*(1 - (d/D)^4)/((pi()/4)^2*2*rho*(0.01*d)^4)/100;

% Create a function that computes the delta pressure
dP_fun = Function('dP_fun',{Q_var, Cd, d},{DP_var});

%%%%%%%
% Map %
%%%%%%%
%simulating the system for all the inputs at the same time
N = size(Qplant,2);%number of points used in the estimation procedure
all_samples = dP_fun.map(N);

dHat = [];
for well = 1:3
    % previously estimated cd (depends on the well)
    Cd_temp = par.Cd_hat(well);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimating model parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dP_sym = dP_fun(Qplant(well,:),Cd_temp,d);
    
    if well == 1
        temp2 = DPplant(well,:) - 1.5 + 0.7;
    elseif well == 2
        temp2 = DPplant(well,:) + 1.7 + 0.24;
    else %expTable(exp,2) == 3
        temp2 = DPplant(well,:) + 0 + 0.20;
    end
    e = temp2 - dP_sym;
    
    % Create an NLP solver
    opts2 = struct;
    %opts2.ipopt.max_iter = nmpcPar.maxiter;
    opts2.ipopt.print_level = 0;
    %opts2.ipopt.tol = nmpcPar.tol;
    %opts2.ipopt.acceptable_tol = 100*nmpcPar.tol; % optimality convergence tolerance
    %opts2.ipopt.linear_solver = 'mumps';
    nlp2 = struct('x', d, 'f', 0.5*dot(e,e));
    solver = nlpsol('solver','ipopt', nlp2,opts2);

    sol = solver('x0',dguess(well),'lbx',par.dMin,'ubx',par.dMax);
    dHat = [dHat; full(sol.x)];
    
end

end

