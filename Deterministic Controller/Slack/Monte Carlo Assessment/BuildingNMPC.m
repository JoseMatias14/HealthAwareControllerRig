function solver = BuildingNMPC(diff,alg,x_var,z_var,p_var,par,nmpcPar)

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
w_pr = z_var(7:9);
%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 

L =  -sum(CR*(w_pr*1e-2)./par.rho_o)+ S1'*nmpcPar.rho*S1  + 1/2 * ((U1 - U_1)'*nmpcPar.R*(U1 - U_1));

% creating system function (LHS of the dynamic equations)
f = Function('f',{x_var,z_var,p_var,U1,U_1,S1},{diff,alg,L});

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
xk_meas = MX.sym('xk_meas',nmpcPar.nx);
zk_meas = MX.sym('zk_meas',nmpcPar.nz);
uk_meas = MX.sym('uk_meas',nmpcPar.nu);

% Ppump,res_theta,val_theta - which are fixed
p = MX.sym('p',7);

%% Lifting initial conditions

% initial state
x_prev = MX.sym('X0',nmpcPar.nx);
w = {w{:},x_prev}; % 1-3
% lbw = [lbw,xk_meas]; 
% ubw = [ubw,xk_meas];
% w0 = [w0;xk_meas];

% initial input
uk = MX.sym('uk_init',nmpcPar.nu);
w = {w{:}, uk}; % 4-6
% w0 = [w0;uk_meas];
% lbw = [lbw;uk_meas];
% ubw = [ubw;uk_meas];



%% Looping through until timeend
for k = 1:nmpcPar.np
    
    % storing the previous input
    uprev = uk; 
    
    % creating current input
    uk = MX.sym(['uk_' num2str(k)],nmpcPar.nu);  
    w = {w{:}, uk}; % 7-9
%     w0 = [w0;  uk_meas];
%     lbw = [lbw;nmpcPar.umin*ones(nmpcPar.nu,1)];
%     ubw = [ubw;nmpcPar.umax*ones(nmpcPar.nu,1)];
    
    % creating current slack variables
    s = MX.sym(['s_',num2str(k)],nmpcPar.nx);
    w = {w{:},s}; % 10-12
%     w0 = [w0;0*ones(nmpcPar.nu,1)];
%     lbw = [lbw;0*ones(nmpcPar.nu,1)];
%     ubw = [ubw;1*ones(nmpcPar.nu,1)];
    
    % Adding constraint for delta_u
    duk = uk - uprev;
    g = {g{:},duk};

%     if k > nmpcPar.nm
%         lbg = [lbg;zeros(nmpcPar.nu,1)];
%         ubg = [ubg; zeros(nmpcPar.nu,1)];
%     else
%         lbg = [lbg;-nmpcPar.dumax*ones(nmpcPar.nu,1)];
%         ubg = [ubg;nmpcPar.dumax*ones(nmpcPar.nu,1)];
%     end
    
    % Collocation points
    fk = [];
    Xk1 = [];
    gk = [];
    quad = [];
    
    for d = 1:3
        % creating states at collocation points
        Xk = MX.sym(['Xk_' num2str(k),'_',num2str(d)],nmpcPar.nx);
        Zk = MX.sym(['Zk_' num2str(k),'_',num2str(d)],nmpcPar.nz);
        w = {w{:}, Xk, Zk}; % 13-15 | 16 - 77
%         w0 = [w0;xk_meas;zk_meas];
%         lbw = [lbw;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];
%         ubw = [ubw;inf*ones(nmpcPar.nx,1);inf*ones(nmpcPar.nz,1)];
     
        % for continuinity
        Xk1 = [Xk1,Xk];
        
        % Calculating xdot and objective function
        [fk1,gk1,qj] = f(Xk,Zk,vertcat(uk,p),uk,uprev,s);
        
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
%         lbg = [lbg;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];
%         ubg = [ubg;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];  
    end
    
    
    % updating objective function
    %ML = M*L1;
    J = J + h*M(end,:)*quad;
    
    % New NLP variable for state at end
    x_prev = MX.sym(['x_init_' num2str(k)],nmpcPar.nx); 
    w = {w{:}, x_prev}; % 78 - 80
%     w0 = [w0;xk_meas];
%     lbw = [lbw;zeros(nmpcPar.nx,1)];
%     ubw = [ubw;inf*ones(nmpcPar.nx,1)];
%     
    % Gap
    g = {g{:},x_next-x_prev};
%     lbg = [lbg;zeros(nmpcPar.nx,1)];
%     ubg = [ubg;zeros(nmpcPar.nx,1)];
    
    % Constraint on erosion
    g = {g{:},x_prev-s};
%     lbg = [lbg;zeros(nmpcPar.nx,1)];
%     ubg = [ubg;nmpcPar.x_threshold*ones(nmpcPar.nx,1)];
    
end

% Formalizing problem 
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',vertcat(xk_meas,zk_meas,uk_meas,p));

% Assigning solver (IPOPT)
solver = nlpsol('solver','ipopt',nlp);


end

