clear
close all
clc

import casadi.*
%% Building the model
%plant parameters
par = ParametersGasLiftModel;
par.T = 60; % simulation sampling time[s]

%initial condition
[dx0,z0,u0,theta0] = InitialConditionGasLift(par);

% control tuning
nmpcConfig.umax = 0.9;
nmpcConfig.umin = 0.1;
nmpcConfig.dumax = 0.01; 
nmpcConfig.x_threshold = 0.305;  %0.305 | 0.3181
nmpcConfig.nx = size(dx0,1);
nmpcConfig.nz = size(z0,1);
nmpcConfig.nu = size(u0,1);
nmpcConfig.ntheta = size(theta0,1);

nmpcConfig.nm = 10;
nmpcConfig.np = 10;
nmpcConfig.rho = 1e9*eye(nmpcConfig.nx); % 1e3 | 999999
nmpcConfig.R = 1*eye(nmpcConfig.nu); % 0 | 0.01 | 1 

% Building Dynamic Model
[F,diff,alg,x_var,z_var,p_var] = BuildingDynModel(par);

% Building NMPC
solverNMPC = BuildingNMPC(diff,alg,x_var,z_var,p_var,par,nmpcConfig);

% Calculating input with NMPC
[uk,s,erosionHat,inputSeq,OF,solFlag] = SolvingNMPC(solverNMPC,dx0,z0,u0,theta0,nmpcConfig);

% figure
subplot(3,1,1)
    plot(0:nmpcConfig.np, inputSeq(1,:),'ro')
    hold on 
    plot(0:nmpcConfig.np, inputSeq(2,:),'kx')
    plot(0:nmpcConfig.np, inputSeq(3,:),'bd')

    xlim([0, nmpcConfig.np])
    
    ylabel('vo [-]')
    xlabel('k [-]')
    legend('1','2','3')
    
subplot(3,1,2)
    plot(0:nmpcConfig.np, erosionHat(1,:),'ro')
    hold on 
    plot(0:nmpcConfig.np, erosionHat(2,:),'kx')
    plot(0:nmpcConfig.np, erosionHat(3,:),'bd')
    yline(nmpcConfig.x_threshold,'k:');
    
    xlim([0, nmpcConfig.np])
    
    ylabel('d_{hat} [cm]')
    xlabel('k [-]')
    legend('1','2','3')
    
subplot(3,1,3)
    plot(0:nmpcConfig.np, s(1,:),'ro')
    hold on 
    plot(0:nmpcConfig.np, s(2,:),'kx')
    plot(0:nmpcConfig.np, s(3,:),'bd')
    
    xlim([0, nmpcConfig.np])

    ylabel('slack [cm]')
    xlabel('k [-]')
    legend('1','2','3')

