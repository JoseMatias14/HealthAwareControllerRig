function solver = BuildingNMPC_MS_prob(diff,alg,x_var,z_var,p_var,par,nmpcPar)

import casadi.*

%% Using 3 collocation points:
h = par.T;
% Radau
t = [collocation_points(3, 'radau')];
% Finding M 
M = [t',0.5*t'.^2,1/3*t'.^3]*inv([[1;1;1],t',t'.^2]);

%% Defining system OF
% for computing Du
U_1 = MX.sym('U_1',nmpcPar.nu);
U1 = MX.sym('U1',nmpcPar.nu);

% slack variable
S1 = MX.sym('s',nmpcPar.nx);

% objective function
w_po = z_var(4:6);
L =  -sum(w_po)+ S1'*nmpcPar.rho*S1  + 1/2 * ((U1 - U_1)'*nmpcPar.R*(U1 - U_1));

% creating system function (LHS of the dynamic equations)
f = Function('f',{x_var,z_var,p_var,U1,U_1,S1},{diff,alg,L});

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

for ll = 1:S
    %% Lifting initial conditions

    J_temp = 0;
    prob_temp = 1;
    
    % initial state
    x_prev = MX.sym(['X0_',num2str(ll)],nmpcPar.nx);
    w = {w{:},x_prev}; 

    % initial input
    uk = MX.sym(['uk_init_',num2str(ll)],nmpcPar.nu);
    w = {w{:}, uk}; 

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
        uprev = uk; 

        % creating current input
        uk = MX.sym(['uk_',num2str(k),'_' num2str(ll)],nmpcPar.nu);  
        w = {w{:}, uk};
        
        % saving the antecipative constraints
        U_ant = {U_ant{:}, uk};

        % creating current slack variables
        s = MX.sym(['s_',num2str(k),'_' num2str(ll)],nmpcPar.nx);
        w = {w{:},s}; 

        % Adding constraint for delta_u
        duk = uk - uprev;
        g = {g{:},duk};

        % Collocation points
        fk = [];
        Xk1 = [];
        gk = [];
        quad = [];

        for d = 1:3
            % creating states at collocation points
            Xk = MX.sym(['Xk_' num2str(k),'_',num2str(d),'_' num2str(ll)],nmpcPar.nx);
            Zk = MX.sym(['Zk_' num2str(k),'_',num2str(d),'_' num2str(ll)],nmpcPar.nz);
            w = {w{:}, Xk, Zk}; 

            % for continuinity
            Xk1 = [Xk1,Xk];

            % Calculating xdot and objective function
            [fk1,gk1,qj] = f(Xk,Zk,vertcat(uk,p,nmpcPar.scendist(idx,:)'),uk,uprev,s);

            fk = [fk, fk1];
            gk = [gk, gk1];
            quad = [quad;qj];

        end

        % integrating the system
        x_next1 = [];
        for d = 1:3 
            % Calculating M*xdot for each collocation point
            Mfk = M(d,1)*fk(:,1) + M(d,2)*fk(:,2) + M(d,3)*fk(:,3);

            % Calculating x
            x_next = x_prev+h*Mfk;
            x_next1 = [x_next1,x_next];

            % Adding xk and Xk1 as constrains as they must be equal - in
            % collocation intervals
            % algebraic constraints are set to zero in the collocation point
            g = {g{:},x_next-Xk1(:,d),gk(:,d)};
        end

        if k <= nmpcPar.Nr
            % calculating the probability of a given scenario
            prob_temp = prob_temp*Zk{22}*Zk{23}*Zk{24};
        end

        % updating objective function
        %ML = M*L1;
        J_temp = J_temp + h*M(end,:)*quad; 
        % N.B. we are using same weight for the scenarios
                  
        % New NLP variable for state at end
        x_prev = MX.sym(['x_init_',num2str(k),'_' num2str(ll)],nmpcPar.nx); 
        w = {w{:}, x_prev}; 

        % Gap
        g = {g{:},x_next-x_prev};

        % Constraint on erosion
        g = {g{:},x_prev-s};

    end
    
    prob_scen = {prob_scen{:}, prob_temp};
    J_scen = {J_scen{:}, J_temp};
    
end

% Construct input matrix and E matrix
U_ant = reshape(U_ant, nmpcPar.np, S)';

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

% Add Non-anticipativity constraints
for k = 1:nmpcPar.Nr
    for js = 1:S-1
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
J = vertcat(J_scen{:})'*vertcat(prob_scen{:});%/(sum(vertcat(prob_scen{:})));


% Formalizing problem 
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',vertcat(xk_meas,zk_meas,uk_meas,p));

% Create an NLP solver
opts = struct;
%opts.ipopt.max_iter = par.maxiter;
%opts.ipopt.print_level = 0;
%opts.print_time = 0;
%opts.ipopt.tol = par.tol;
%opts.ipopt.acceptable_tol = 100*par.tol; % optimality convergence tolerance
opts.ipopt.linear_solver = 'mumps';

% Assigning solver (IPOPT)
solver = nlpsol('solver','ipopt',nlp,opts);


end

