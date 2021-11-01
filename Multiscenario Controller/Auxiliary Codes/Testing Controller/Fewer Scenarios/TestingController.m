clear
close all
clc

import casadi.*
%% Building the model
%plant parameters
par = ParametersGasLiftModel;
par.T = 60; % simulation sampling time[s]

%initial condition
[dx0,z0,u0,theta0] = InitialConditionGasLift_model(par);

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

% levels of the cfd used
probH = 0.9; % 0.75 | 0.95
probL = 0.50;

% testing the probabilities of different scenarios
%            w1 |  w2  | w3 
nmpcConfig.scendist =  [probH, probH, probH; %S1:   HHH
                        probH, probH, probL; %S2:   HHL
                        probH, probL, probH; %S3:   HLH
                        probL, probH, probH; %S4:   LHH
                        probH, probL, probL; %S5:   HLL
                        probL, probH, probL; %S6:   LHL
                        probL, probL, probH; %S7:   LLH
                        probL, probL, probL];%S8:   LLL

nmpcConfig.Nr = 2; 
nmpcConfig.Nlevels = size(nmpcConfig.scendist,1); 
nmpcConfig.Nscenarios = nmpcConfig.Nlevels^nmpcConfig.Nr;
nmpcConfig.ConsideredScen = [40,48,56,61,62,63,64];

load('scenario_pruning_cumprob')
temp = scen_avg_prob(nmpcConfig.ConsideredScen);
temp2 = sum(temp);
nmpcConfig.scenProb = temp/temp2;

% Building Dynamic Model
[F,diff,alg,x_var,z_var,p_var] = BuildingDynModel(par);

% Building NMPC
solverNMPC_MS = BuildingNMPC_MS_prob(diff,alg,x_var,z_var,p_var,par,nmpcConfig);

% Calculating input with NMPC
[uk,s,erosionHat,inputSeq,OF,solFlag] = SolvingNMPC_MS(solverNMPC_MS,dx0,z0,u0,theta0,nmpcConfig);

%% For plotting
Scons = length(nmpcConfig.ConsideredScen);
tit = {};
for ll = 1:Scons
    tit = {tit{:},['Scen.: ',num2str(nmpcConfig.ConsideredScen(ll))]};
end

palette = jet(Scons);

for ll = 1:Scons % for scenarios
    figure(ll)
    
    for well = 1:3
        subplot(3,1,well)
        
        % INPUTS
        yyaxis left
        plot(0:nmpcConfig.np, inputSeq{ll}(well,:),'bd-','linewidth',1.5)
        hold on
        ylabel('vo [-]','FontSize',8)
        
        % EROSION
        yyaxis right
        plot(0:nmpcConfig.np, erosionHat{ll}(well,:),'rx-','linewidth',1.5)
        yline(nmpcConfig.x_threshold,'r:','linewidth',1);
        
        xlim([0, nmpcConfig.np])
        ylabel('d_{hat} [cm]','FontSize',8)
        
        title(['Well: ',num2str(well),' - ',tit{ll}])        
        xlabel('k [-]')
        
    end
end

%% SLACK
figure(Scons + 1)
for ll = 1:Scons % for scenarios
    
    for well = 1:3
    subplot(Scons,3,well + 3*(ll - 1))    
        plot(0:nmpcConfig.np, s{ll}(well,:),'Color',palette(ll,:),'linewidth',1.5)
        hold on
        
        xlim([0, nmpcConfig.np])

        if ll == 1
            title(['Well: ',num2str(well)])
        else
            if ll == Scons
                xlabel('k [-]')
            end
        end
        
        if well == 1
            ylabel(tit{ll})
        else
            if well == 2
                ylabel('slack [cm]')
            end
        end
    end
end

