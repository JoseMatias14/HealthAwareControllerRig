function solver = BuildingHAC(eq,x_var,p_var,par,nmpcPar)

import casadi.*

%% Defining system OF

% liquid production
%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 
Q_pr = CR*(x_var(7:9)*1e-2)./par.rho_o;

deltaDeg = {};
for well = 1:3
    if well == 1
        DPmax = 0.3081*Q_pr{well}^2 + 0.185*Q_pr{well} - 0.7528;
        DPmin = 0.1828*Q_pr{well}^2 + 1.0093*Q_pr{well} - 4.6835;
    elseif well == 2
        DPmax = 0.2934*Q_pr{well}^2 - 0.0592*Q_pr{well} - 0.7532;
        DPmin = 0.1896*Q_pr{well}^2 + 0.8415*Q_pr{well} - 4.9335;    
    else % well == 3
        DPmax = 0.2902*Q_pr{well}^2 + 0.9971*Q_pr{well} - 1.5097;
        DPmin =  0.19*Q_pr{well}^2 + 1.1067*Q_pr{well} - 3.685;    
    end

    DPNorm =  1 - (x_var(15 + well) - DPmin)./(DPmax - DPmin);
    
    % constraint violation
    deltaDeg = {deltaDeg{:}, (DPNorm - nmpcPar.x_healthy)./(nmpcPar.x_threshold - nmpcPar.x_healthy)};
end

% computing using utility function
L =  - ((1 - exp(-deltaDeg{1}*Q_pr{1}))/(deltaDeg{1}+0.001) + (1 - exp(-deltaDeg{2}*Q_pr{2}))/(deltaDeg{2}+0.001) + (1 - exp(-deltaDeg{3}*Q_pr{3}))/(deltaDeg{3}+0.001));
%L =  - sum(Q_pr);
%L =  0;

% creating system function (LHS of the dynamic equations)
f = Function('f',{x_var,p_var},{eq,L});

%% Defining empty nlp-problem
% objective function
J = 0;

% declare variables (bounds and initial guess)
w = {};
% w0 = [];
% lbw =[];
% ubw = [];

% declare constraints and its bounds
g = {};
% lbg = [];
% ubg = [];

%% declaring parameters
% Fixed during the simulation
% Ppump (1),T_s (3),res_theta (3),val_theta (3),mu_rg (15),std_rg (15),mu_rp (3),std_rp (3) = 46;
fixedP = MX.sym('p',46);

%% Lifting initial conditions

% initial time
tk = MX.sym('t0',1);
w = {w{:}, tk}; % 1
% w0 = [w0;uk_meas];
% lbw = [lbw;uk_meas];
% ubw = [ubw;uk_meas];

% initial input
Uk = MX.sym('U0',nmpcPar.nu);
w = {w{:}, Uk}; % 2-4
% w0 = [w0;uk_meas];
% lbw = [lbw;uk_meas];
% ubw = [ubw;uk_meas];

%% Looping through until time end
for k = 1:nmpcPar.np
    
    % storing the previous input
    uprev = Uk; 
    
    % creating current input
    Uk = MX.sym(['U_' num2str(k)],nmpcPar.nu);  
    w = {w{:}, Uk}; % 5-7
%     w0 = [w0;  uk_meas];
%     lbw = [lbw;nmpcPar.umin*ones(nmpcPar.nu,1)];
%     ubw = [ubw;nmpcPar.umax*ones(nmpcPar.nu,1)];
        
    % Adding constraint for delta_u
    duk = Uk - uprev;
    g = {g{:},duk};

%     if k > nmpcPar.nm
%         lbg = [lbg;zeros(nmpcPar.nu,1)];
%         ubg = [ubg; zeros(nmpcPar.nu,1)];
%     else
%         lbg = [lbg;-nmpcPar.dumax*ones(nmpcPar.nu,1)];
%         ubg = [ubg;nmpcPar.dumax*ones(nmpcPar.nu,1)];
%     end
    
    % creating states at the evaluation points
    Xk = MX.sym(['Xk_' num2str(k)],nmpcPar.nx); 
    w = {w{:}, Xk}; % 8-26
    
    % Calculating the residuals
    [fk,qk] = f(Xk,vertcat(Uk,tk,fixedP));% 1.7159
    g = {g{:},fk};   % residuals must be equal to zero 
%    lbg = [lbg;zeros(nmpcPar.nx,1)];
%    ubg = [ubg;zeros(nmpcPar.nx,1)];         
        
    J = J + qk;
end

% formalizing problem
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',fixedP);

% Create an NLP solver
opts = struct;
%opts.ipopt.max_iter = nmpcPar.maxiter;
opts.ipopt.print_level = 5;
%opts.print_time = 0;
%opts.ipopt.tol = nmpcPar.tol;
%opts.ipopt.acceptable_tol = 100*nmpcPar.tol; % optimality convergence tolerance
%opts.ipopt.linear_solver = 'mumps';

% Assigning solver (IPOPT)
solver = nlpsol('solver','ipopt',nlp,opts);


end

