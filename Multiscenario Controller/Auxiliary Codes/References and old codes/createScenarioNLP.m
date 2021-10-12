function [solver, par] = createScenarioNLP(f, par, Nlevels, Nr, options)

import casadi.*

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0 collocation_points(d, 'radau')]; %can be 'legendre'

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
    % Construct Lagrange polynomials to get the polynomial basis
    % at the collocation point
    coeff = 1;
    for r=1:d+1
        if r ~= j
            coeff = conv(coeff, [1, -tau_root(r)]);
            coeff = coeff / (tau_root(j)-tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the
    % coefficients of the continuity equation
    D(j) = polyval(coeff, 1.0);
    
    % Evaluate the time derivative of the polynomial at all collocation
    % points to get the coefficients of the continuity equation
    pder = polyder(coeff);
    for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
    end
    
    % Evaluate the integral of the polynomial to get the coefficients
    % of the quadrature function
    pint = polyint(coeff);
    B(j) = polyval(pint, 1.0);
end

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

S = Nlevels^Nr;
U_ant = {};

% "Lift" initial conditions by setting initial value as a parameter p
p = MX.sym('p', par.nx);

% Formulate the scenario NLP
for l = 1:S
    
    % States (at stage k and scenario l)
    Xkl = MX.sym('X0', par.nx);
    w = {w{:}, Xkl};
    lbw = [lbw; par.xlb];
    ubw = [ubw; par.xub];
    w0 = [w0; par.x_init];
    
    g = {g{:}, Xkl - p};
    lbg = [lbg; zeros(par.nx,1)];
    ubg = [ubg; zeros(par.nx,1)];
    
    for k = 0:par.Nt-1
        
        if mod(l,Nlevels) ~= 0
            idx = mod(l,Nlevels);
        else
            idx = Nlevels;
        end
        
        % New NLP variable for the control
        Ukl = MX.sym(['U_' num2str(k) '_' num2str(l)], par.nu);
        w = {w{:}, Ukl};
        U_ant = {U_ant{:}, Ukl};
        lbw = [lbw; par.ulb];
        ubw = [ubw; par.uub];
        w0 = [w0; par.u_init];
        
        % State at collocation points
        Xklj = {};
        
        for j = 1:d
            % Differential states
            Xklj{j} = MX.sym(['X_' num2str(k) '_' num2str(l) '_' num2str(j)], par.nx);
            w = {w{:}, Xklj{j}};
            lbw = [lbw; par.xlb];
            ubw = [ubw; par.xub];
            w0 = [w0; par.x_init];
        end
        
        % Loop over collocation points
        Xkl_end = D(1)*Xkl;
        
        for j = 1:d
            % Expression for the state derivative at the collocation point
            xp = C(1,j+1)*Xkl;
            for r=1:d
                xp = xp + C(r+1,j+1)*Xklj{r};
            end
            
            % Append collocation equations
            [fj, qj] = f(Xklj{j}, Ukl, par.scendist(:,idx));
            g = {g{:}, par.h*fj - xp};
            lbg = [lbg; zeros(par.nx,1)];
            ubg = [ubg; zeros(par.nx,1)];
            
            % Add contribution to the end state
            Xkl_end = Xkl_end + D(j+1)*Xklj{j};
            
            % Add contribution to quadrature function
            J = J + B(j+1)*qj*par.h;
        end
        
        % Add regularization terms to economic cost objective
        switch options.reg_terms
            case 'on'
                Phi_reg = par.qA/2*((Xkl(1)-par.xreg(1,idx))^2) ...
                    + par.qTk/2*((Xkl(4)-par.xreg(4,idx))^2) + par.qTr/2*((Xkl(3)-par.xreg(3,idx))^2);  
            case 'off'
                Phi_reg = 0;
        end
        
        if k>0
            switch options.control_penalty
                case 'on'                    
                    Phi_control = par.r1*(Ukl(1) - Ukl_prev(1))^2 + par.r2*(Ukl(2) - Ukl_prev(2))^2;                    
                case 'off'                    
                    Phi_control = 0;
            end
        else
            Phi_control = 0;
        end
        
        J = J + Phi_reg + Phi_control;
        Ukl_prev = Ukl;
        
        % New NLP variable for state at end of interval
        Xkl = MX.sym(['X_' num2str(k+1) '_' num2str(l)], par.nx);
        w = {w{:}, Xkl};
        lbw = [lbw; par.xlb];
        ubw = [ubw; par.xub];
        w0 = [w0; par.x_init];
        
        % Shooting gap constraint (continuity for differential states)
        g = {g{:}, Xkl_end-Xkl};
        lbg = [lbg; zeros(par.nx,1)];
        ubg = [ubg; zeros(par.nx,1)];
        
        % Add contribution of terminal cost and constraints for each scenario
        if k == par.Nt-1
            switch options.termin_con
                case 'on'
                    Vn = (Xkl - par.x_final(:,l))'*par.P(:,:,l)*(Xkl - par.x_final(:,l));                    
                    % Add terminal constraints
                    g = {g{:}, (Xkl - par.x_final(:,l))'*(Xkl - par.x_final(:,l))};
                    lbg = [lbg; 0];
                    ubg = [ubg; par.Cf(idx)];                    
                case 'off'
                    Vn = 0;
            end
            J = (J + Vn)/S;
        end
    end
end


% Construct input matrix and E matrix
U_ant = reshape(U_ant, par.Nt, S)';

a = 1;
ac = 1;

for i = 1:par.Nt
    for j = 1:S
        nonant(i,j) = a;
        if i == 1
            a = 1;
        else if i <= Nr
                ac = ac+1;
                if ac > Nlevels
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
for k = 1:Nr
    for js = 1:S-1
        if ~isnan(nonant(js,k))
            if nonant(js,k) == nonant(js+1,k)
                g = {g{:},(U_ant{js,k} - U_ant{js+1,k})};
                lbg = [lbg; zeros(par.nu,1)];
                ubg = [ubg; zeros(par.nu,1)];
            end
        end
    end
end

par.w0 = w0;
par.lbw = lbw;
par.ubw = ubw;
par.lbg = lbg;
par.ubg = ubg;

% Create an NLP solver
opts = struct;
opts.ipopt.max_iter = par.maxiter;
opts.ipopt.print_level = 0;
opts.print_time = 0;
%opts.ipopt.tol = par.tol;
%opts.ipopt.acceptable_tol = 100*par.tol; % optimality convergence tolerance
opts.ipopt.linear_solver = 'mumps';

prob = struct('f', J, 'x', vertcat(w{:}), 'p', p, 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob, opts);

par.prob = prob;
par.d = d;
par.S = S;
par.Nlevels = Nlevels;
par.Nr = Nr;
end
