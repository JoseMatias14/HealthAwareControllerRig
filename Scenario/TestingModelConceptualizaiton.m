clear
close all
clc

import casadi.*
%% Building the model
%plant parameters
par = ParametersGasLiftModel;
par.T = 60; % simulation sampling time[s]

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
%number of wells
n_w = par.n_w; %[]
%gas constant
R = par.R; %[m3 Pa K^?1 mol^?1]
%molecular weigth
Mw = par.Mw; %[kg/mol?]

%properties
%density of oil - dim:  nwells x 1
rho_o = par.rho_o; %[kg/m3]
%1cP oil viscosity
mu_oil = par.mu_oil;% [Pa s] 
%density of gas
rho_g = par.rho_g; %[kg/m3]

%project
%well parameters - dim:  nwells x 1
L_w = par.L_w; %[m]
H_w = par.H_w; %[m]
D_w = par.D_w; %[m]
A_w = par.A_w;%[m2]

%well below injection - [m]
L_bh = par.L_bh;
H_bh = par.H_bh;
D_bh = par.D_bh;
A_bh = par.A_bh;%[m2]

%riser - [m]
L_r = par.L_r;
H_r = par.H_r;
D_r = par.D_r;
A_r = par.A_r;%[m2]

Cd = par.Cd_hat;

%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
% System states %
%%%%%%%%%%%%%%%%%
% --> orifice diameter 
d = MX.sym('d',n_w); %[cm]      % 1:3 [cm]

%oil holdup 
m_o = MX.sym('m_or',n_w);       % 1:3 [kg]
%oil rate from reservoir
w_ro = MX.sym('w_ro',n_w);      % 4:6 [1e2 kg/s]
%riser head total production rate
w_pr = MX.sym('w_pr',n_w);      % 7:9 [1e2 kg/s]
%riser head pressure
p_rh = MX.sym('p_rh',n_w);      % 10:12 [bar]
%pressure - below injection point (bottom hole)
p_bh = MX.sym('p_bh',n_w);      % 13:15 [bar]
%pressure - below injection point (bottom hole)
dP = MX.sym('dP',n_w);          % 16:18 [mbar]
%cummulative probability of a given increment
probCum = MX.sym('probCum',n_w);      % 19:21 [-]
%probability of a given increment
prob = MX.sym('prob',n_w);      % 19:21 [-]
%increment in the diameter due to erosion
dD = MX.sym('dD',n_w);          % 22:24 [cm]

%%%%%%%%%%%%%%%%
% System input %
%%%%%%%%%%%%%%%%
%valve oppening
vo=MX.sym('vo',n_w);            % 1:3 [0-1]
%Pump outlet pressure
Ppump=MX.sym('Ppump',1);      % 4 [bar g]

%%%%%%%%%%%%%%
% parameters %
%%%%%%%%%%%%%%
% fixed
%riser temperature
T_r = MX.sym('T_r',1); %[oC]
%separator pressure
p_s = MX.sym('p_s',1); %[bar]

%time transformation: CASADI integrates always from 0 to 1 and the USER does the time
%scaling with T --> sampling time
t_samp = MX.sym('t_samp',1); %[s]

% estimable
%reservoir data-driven model
res_theta = MX.sym('res_theta',n_w);%[m2]
%riser valve characteristics
val_theta = MX.sym('val_theta',n_w);%[m2]


%%%%%%%%%%%%
% Modeling %
%%%%%%%%%%%%
%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 

%oil from reservoir flowrate [kg/s] - same pump outlet pressure --> with re-arranged parameters
f1 = -Ppump*ones(n_w,1)*1e5 + (w_ro.*1e-2).^2.*(res_theta.*1e9)./(vo.^2.*rho_o) + p_bh.*1e5 ; 
%riser head pressure [Pa] --> with re-arranged parameters
f2 = -p_rh.*1e5 + (w_pr.*1e-2).^2.*(val_theta.*1e8)./rho_o + p_s.*1e5 ;
%bottom hole pressure [Pa]
f3 = -p_bh.*1e5 + (dP*1e2 + p_rh.*1e5 + rho_o.*9.81.*H_r + 128.*mu_oil.*(L_w+L_r).*(w_ro.*1e-2)./(3.14.*D_bh.^4.*rho_o));
%Oil Mass [kg]
f4 = -m_o + (rho_o).*((A_w.*L_w + A_r.*L_r));
% orifice pressure drop 
f5 = dP*1e2 - Cd*(w_ro.*1e-2).^2.*(1 - (d/2).^4)./((pi()/4)^2*2*rho_o.*(0.01*d).^4);
% mass conservation
f6 = w_ro - w_pr;
% gamma distribution cfd (of a given increment)
QQ = CR*(w_ro.*1e-2)./rho_o; % kg/s --> L/min
theta = (0.0043*QQ.^3 - 0.0949*QQ.^2 + 0.7305.*QQ - 1.32).^-1; %empirical
alpha = 2;
f7 = - probCum + (1 - exp(-theta.*dD) - (theta.*dD).*exp(-theta.*dD));
f8 = - prob + (theta.^alpha).*(dD.^(alpha - 1)).*exp(-theta.*dD)/1;


% d_{k+1} = d_{k} + dD
df1 = 0.0005*dD/60; %% [min] to [s]

% Form the DAE system
diff = vertcat(df1);
alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8);

% give parameter values
alg = substitute(alg,p_s,par.p_s);
alg = substitute(alg,T_r,par.T_r);
alg = substitute(alg,probCum,[0.75;0.5;0.25]);

% concatenate the differential and algebraic states
x_var = vertcat(d);
z_var = vertcat(m_o,w_ro,w_pr,p_rh,p_bh,dP,dD,prob);
u_var = vertcat(vo,Ppump);
p_var = vertcat(res_theta,val_theta);

% end modeling

% Casadi commands
%declaring function in standard DAE form (scaled time)
dae = struct('x',x_var,'z',z_var,'p',vertcat(u_var,p_var,t_samp),'ode',t_samp*diff,'alg',alg);

%calling the integrator, the necessary inputs are: label; integrator; function with IO scheme of a DAE (formalized); struct (options)
F = integrator('F','idas',dae);

%% Simulation tuning
%initial condition
[dxk,zk,uk,thetak] = InitialConditionGasLift(par);

xArray = [];
zArray = [];

for kk = 1:40

    % integrating the "plant" (the model sampling time is defined by parPlant.T)
    Fend = F('x0',dxk,'z0',zk,'p',[uk;thetak;par.T]);

        %extracting the results (from Casadi symbolic to numerical)
        dxk = full(Fend.xf);
        zk = full(Fend.zf);

        xArray = [xArray, dxk];
        zArray = [zArray, zk];

end

figure
subplot(2,1,1)
    plot(1:40, xArray(1,:),'ro')
    hold on 
    plot(1:40, xArray(2,:),'kx')
    plot(1:40, xArray(3,:),'bd')

    ylabel('d [cm]')
    xlabel('t [min]')
    legend('1','2','3')
    
subplot(2,1,2)
    plot(1:40, zArray(22,:),'ro')
    hold on 
    plot(1:40, zArray(23,:),'kx')
    plot(1:40, zArray(24,:),'bd')

    ylabel('prob [-]')
    xlabel('t [min]')
    legend('1','2','3')

% %%%%%%%%%%%%
% % Plotting %
% %%%%%%%%%%%%
% %% plotting the data 
% % checking sampling rate
% markers = {'o','x','>'};
% cc = {'b','k','g'};
% leg = {'w_1','w_2','w_3'};
%  
% %% Inputs (valve opening)
% f = figure(1);
% for well = 1:3
%     subplot(3,1,well)
%         plot(tgrid, uPlantArray(well,:),'bo','Linewidth',1)
%         
%         ylim([0.05 0.95])
%         
%         xticks(0:(10*parPlant.T/60):(1/60)*(nFinal))
%         xlim([0 (1/60)*(nFinal)])
% 
%         xlabel('time [min]','FontSize',10)
%         ylabel('v_o [-]','FontSize',10)
%         
%         name = ['Valve opening - Well ',num2str(well)];
%         title(name,'FontSize',10)   
% end
% 
% %% Outputs (valve opening)
% f = figure(2);
% for well = 1:3
%      % Reservoir Valve Parameters
%     subplot(3,1,well)
%     
%         yyaxis left
%         plot(tgrid, measPlantArray(well,:),'bd-','Linewidth',1)
%         ylim([2 11])
%         ylabel('Q_l [L/min]','FontSize',10)
%         
%         yyaxis right
%         plot(tgrid, measPlantArray(3 + well,:),'rx-','Linewidth',1)
%         ylim([0.9 1.1])
%         ylabel('P_{top} [mbar G]','FontSize',10)
% 
%         xticks(0:(10*parPlant.T/60):(1/60)*(nFinal))
%         xlim([0 (1/60)*(nFinal)])
%         xlabel('time [min]','FontSize',10)
% 
%         name = ['Measurements - Well ',num2str(well)];
%         title(name,'FontSize',10)  
%         
% end
% 
% %% Outputs (dP)
% f = figure(3);
% for well = 1:3
%      % Reservoir Valve Parameters
%     subplot(3,1,well)
%     
%         yyaxis left
%         plot(tgrid, measPlantArray(6 + well,:),'bx-','Linewidth',1.5)
%         %ylim([2 11])
%         ylabel('dP [mbar]','FontSize',10)
%         
%         yyaxis right
%         stairs(tgrid, probeOrificeArray(well,:),'k-','Linewidth',1.5)
%         hold on 
%         yline(nmpcConfig.x_threshold,'r:','Linewidth',1);
%         
%         ylim([parPlant.dMin 0.32])
%         yticks(parPlant.dMin:0.0036:parPlant.dMax)
%         ylabel('Orifice diameter','FontSize',10)
% 
%         xticks(0:(10*parPlant.T/60):(1/60)*(nFinal))
%         xlim([0 (1/60)*(nFinal)])
%         xlabel('time [min]','FontSize',10)
% 
%         name = ['Measurements - Well ',num2str(well)];
%         title(name,'FontSize',10)  
%         
% end
 