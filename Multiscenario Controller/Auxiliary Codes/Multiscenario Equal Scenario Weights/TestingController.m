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
nmpcConfig.x_threshold = 0.305;  % %0.3181
nmpcConfig.nx = size(dx0,1);
nmpcConfig.nz = size(z0,1);
nmpcConfig.nu = size(u0,1);
nmpcConfig.ntheta = size(theta0,1);

nmpcConfig.nm = 10;
nmpcConfig.np = 10;
nmpcConfig.rho = 1e9*eye(nmpcConfig.nx); % 1e3 | 999999
nmpcConfig.R = 1*eye(nmpcConfig.nu); % 0 | 0.01 | 1 

% number of ramifications
nmpcConfig.Nr = 2; 
nmpcConfig.Nlevels = 4; 
% levels of the cfd used
%                       w1 |  w2  | w3 
nmpcConfig.scendist = [0.95, 0.95, 0.95; %S1
                       0.50, 0.95, 0.95; %S2
                       0.50, 0.50, 0.95; %S3
                       0.50, 0.50, 0.50];%S4

% How the scenarios are built
% with 16 scenarios
%S1: 1     1
%S2: 1     2
%S3: 1     3
%S4: 1     4
%S5: 2     1
%S6: 2     2
%S7: 2     3
%S8: 2     4
%S9: 3     1
%S10:3     2
%S11:3     3
%S12:3     4
%S13:4     1
%S14:4     2
%S15:4     3
%S16:4     4

% Building Dynamic Model
[F,diff,alg,x_var,z_var,p_var] = BuildingDynModel(par);

% Building NMPC
solverNMPC_MS = BuildingNMPC_MS(diff,alg,x_var,z_var,p_var,par,nmpcConfig);

% Calculating input with NMPC
[uk,s,erosionHat,inputSeq,OF,solFlag] = SolvingNMPC_MS(solverNMPC_MS,dx0,z0,u0,theta0,nmpcConfig);

%% For plotting
S = nmpcConfig.Nlevels^nmpcConfig.Nr;
palette = jet(S);

% lab = {'vo [-]','d_{hat} [cm]','slack [cm]'};
% 
% leg = {};
% for l = 1:S
%     leg = {leg{:},['S_',num2str(l)]};
% 
% end    
    
for l = 1:S % for scenarios
    figure(l)
    
   
    for well = 1:3
        subplot(3,1,well)
        
        % INPUTS
        yyaxis left
        plot(0:nmpcConfig.np, inputSeq{l}(well,:),'bd-','linewidth',1.5)
        hold on
        ylabel('vo [-]','FontSize',8)
        
        % EROSION
        yyaxis right
        plot(0:nmpcConfig.np, erosionHat{l}(well,:),'rx-','linewidth',1.5)
        yline(nmpcConfig.x_threshold,'r:','linewidth',1);
        
        xlim([0, nmpcConfig.np])
        ylabel('d_{hat} [cm]','FontSize',8)
        
        title(['Well: ',num2str(well),' - Scenario: ',num2str(l)])        
        xlabel('k [-]')
        
    end
end

%% SLACK
figure(S + 1)
for l = 1:S % for scenarios
    
    for well = 1:3
    subplot(S,3,well + (l - 1))    
        plot(0:nmpcConfig.np, s{l}(well,:),'Color',palette(l,:),'linewidth',1.5)
        hold on
        
        xlim([0, nmpcConfig.np])

        if l == 1
            title(['Well: ',num2str(well)])
        else
            if l == S
                xlabel('k [-]')
            end
        end
        ylabel('slack [cm]')
        
    end
end

