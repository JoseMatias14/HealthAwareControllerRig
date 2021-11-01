function solver = BuildingNMPC_MS_prob(diff,alg,x_var,z_var,p_var,par,nmpcPar)

import casadi.*

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0 collocation_points(d, 'legendre')];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);
end

%% Control discretization
h = par.T;

%% Defining system OF
% for computing Du
U_1 = MX.sym('U_1',nmpcPar.nu);
U1 = MX.sym('U1',nmpcPar.nu);

% liquid production
%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 
Q_pr = CR*(z_var(7:9)*1e-2)./par.rho_o;

% constraint violation
constV = (x_var - nmpcPar.x_healthy)./(nmpcPar.x_threshold - nmpcPar.x_healthy);

% computing using utility function
L =  - ((1 - exp(-constV{1}*Q_pr{1}))/(constV{1}+0.001) + (1 - exp(-constV{2}*Q_pr{2}))/(constV{2}+0.001) + (1 - exp(-constV{3}*Q_pr{3}))/(constV{3}+0.001));% + 1/2 * ((U1 - U_1)'*nmpcPar.R*(U1 - U_1));

% creating system function (LHS of the dynamic equations)
f = Function('f',{x_var,z_var,p_var,U1,U_1},{diff,alg,L});

%% Defining empty nlp-problem

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
xk_meas = MX.sym('xk_meas',nmpcPar.nx);
zk_meas = MX.sym('zk_meas',nmpcPar.nz);
uk_meas = MX.sym('uk_meas',nmpcPar.nu);

% original p: pPump, res_theta_{1,2,3}, val_theta_{1,2,3}, probCum_{1,2,3}
% but probCum will change depending on the scenario
p = MX.sym('p',nmpcPar.ntheta - 3);

% total number of scenarios
S = nmpcPar.Nlevels^nmpcPar.Nr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating scenario combinations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transforming the scenarios (decimal) into base Number of levels
temp = dec2base(0:S - 1,nmpcPar.Nlevels);
%  generate all tuples
%  temp is 'char'
% we need to substitute the zeros (Matlab has 1-based index) by 1 
% and by adding the double 1, we obtain the indexes that we want
combs = temp-'0'+1;  %// 

% non-antecipative constraints
U_ant = {};

% probability of each scenario
prob_scen = {};

% objective function
% calculated individually for each scenario
J_scen = {};

for ll = nmpcPar.ConsideredScen
    %% Lifting initial conditions

    J_temp = 0;
%     prob_temp = 1;
    
    % lifiting initial state
    Xk = MX.sym(['X0_',num2str(ll)],nmpcPar.nx);
    w = {w{:}, Xk}; 

    % initial input
    Uk = MX.sym(['U0_',num2str(ll)],nmpcPar.nu);
    w = {w{:}, Uk}; 

    %% Looping through until timeend
    for k = 1:nmpcPar.np

        % checking which of the uncertainty realizations (pre-defined) is going to used in the
        % current scenario
        if k < nmpcPar.Nr
            idx = combs(ll,k);
        else
            idx = combs(ll,nmpcPar.Nr);
        end
        
        % storing the previous input
        uprev = Uk; 

        % creating current input
        Uk = MX.sym(['U_',num2str(k),'_' num2str(ll)],nmpcPar.nu);  
        w = {w{:}, Uk};
        
        % saving the antecipative constraints
        U_ant = {U_ant{:}, Uk};

        % Adding constraint for delta_u
        duk = Uk - uprev;
        g = {g{:},duk};

        % State at collocation points
        Xkj = {};
        Zkj = {};
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j),'_' num2str(ll)],nmpcPar.nx);
            Zkj{j} = MX.sym(['Z_' num2str(k),'_',num2str(d),'_' num2str(ll)],nmpcPar.nz);
            w = {w{:}, Xkj{j},Zkj{j}};
        end
        
        % Loop over collocation points
        Xk_end = D(1)*Xk;
        for j=1:d
            % Expression for the state derivative at the collocation point
            xp = C(1,j+1)*Xk;
            for r=1:d
                xp = xp + C(r+1,j+1)*Xkj{r};
            end
            
            % Append collocation equations
            [fj,gj,qj] = f(Xkj{j},Zkj{j},vertcat(Uk,p,nmpcPar.scendist(idx,:)'),Uk,uprev);
            g = {g{:}, h*fj - xp};
            g = {g{:}, gj};
            
            % Add contribution to the end state
            Xk_end = Xk_end + D(j+1)*Xkj{j};
            
            % updating objective function
            % Add contribution to quadrature function
            J_temp = J_temp + B(j+1)*qj*h;
            % N.B. we apply the weights to compute the expectation
            % afterwards
        end       
                  
        % New NLP variable for state at end
        Xk  = MX.sym(['X_',num2str(k+1),'_' num2str(ll)],nmpcPar.nx); 
        w = {w{:}, Xk}; 

        % Gap
        g = {g{:},Xk_end-Xk};

    end
    
%     prob_scen = {prob_scen{:}, prob_temp};
    J_scen = {J_scen{:}, J_temp};
    
end

a = 1;
ac = 1;

for i = 1:nmpcPar.np
    for j = 1:S
        nonant(i,j) = a;
        if i == 1
            a = 1;
        else if i <= nmpcPar.Nr
                ac = ac+1;
                if ac > nmpcPar.Nlevels
                    a = a+1;
                    ac = 1;
                end
            else
                a = NaN;
            end
        end
    end
end
nonant = nonant';

nonant = nonant(nmpcPar.ConsideredScen,:);

% Computing the number of considered scenarios
Scons = length(nmpcPar.ConsideredScen);

% Construct input matrix and E matrix
U_ant = reshape(U_ant, nmpcPar.np, Scons)';

% Add Non-anticipativity constraints
for k = 1:nmpcPar.Nr
    for js = 1:Scons-1
        if ~isnan(nonant(js,k))
            if nonant(js,k) == nonant(js+1,k)
                g = {g{:},(U_ant{js,k} - U_ant{js+1,k})};
%                 lbg = [lbg; zeros(par.nu,1)];
%                 ubg = [ubg; zeros(par.nu,1)];
            end
        end
    end
end

% calculating the objective function (weighted by the probability of each
% scenario)
J = vertcat(J_scen{:})'*nmpcPar.scenProb;%/(sum(vertcat(prob_scen{:})));


% Formalizing problem 
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',vertcat(xk_meas,zk_meas,uk_meas,p));

% Create an NLP solver
opts = struct;
opts.ipopt.max_iter = nmpcPar.maxiter;
%opts.ipopt.print_level = 0;
%opts.print_time = 0;
opts.ipopt.tol = nmpcPar.tol;
opts.ipopt.acceptable_tol = 100*nmpcPar.tol; % optimality convergence tolerance
%opts.ipopt.linear_solver = 'mumps';

% Assigning solver (IPOPT)
solver = nlpsol('solver','ipopt',nlp,opts);


end

