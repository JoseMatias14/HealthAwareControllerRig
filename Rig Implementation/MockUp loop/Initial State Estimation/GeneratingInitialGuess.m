function flag = GeneratingInitialGuess(thetaHat,xEstHat,zEstHat,uPlant,par,name)
% Function generates initial guess for the models (dynamic and
% steady-state)

%% Inputs
%valve oppening
vo_0 = 0.05*ones(3,1); % uPlant(1:3); %[0-1]
%velocity pump
vpump_0 = uPlant(4);% [bar]

%% Degradation
d_0 = xEstHat;

%% States
%oil holdup 
m_o_0 = zEstHat(1:3);%[kg] 
%oil rate from reservoir %%%%%%%%%%%%%%%%%%%%
w_ro_0 = zEstHat(4:6); %[1e2 kg/s]
%riser head total production rate %%%%%%%%%%%%%%%%%%%%%%
w_pr_0 = zEstHat(7:9);%[1e2 kg/s] 
%riser head pressure
p_rh_0 = zEstHat(10:12);%[bar]
%pressure - below injection point (bottom hole)
p_bh_0 = zEstHat(13:15);%[bar]
% Delta Pressure
Ql = w_ro_0*1e-2*60*1e3./par.rho_o;
dP_0 = 0.2788*Ql.^2 + 1.143*Ql - 2.3831; % dP model previosly calculated

%% Parameters
res_theta_0 = thetaHat(1:3); % rearranged 1e9./(1e-5*theta0(1:3)).^2 [1e-4 m2]
val_theta_0 = thetaHat(4:6); % rearranged 1e8./(1e-4*theta0(4:6)).^2 [1e-5 m2]

%%
results = ['InitialState_',name]; 
save(results,'d_0','m_o_0','w_ro_0','w_pr_0','p_rh_0','p_bh_0','dP_0','vo_0','vpump_0','res_theta_0','val_theta_0');

%%
flag = 1; %finished

end

