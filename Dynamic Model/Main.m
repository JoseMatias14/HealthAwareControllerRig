clear
close all
clc

%% Simulation settings
% 
%paramters
par = ParametersGasLiftModel_model;
%initial condition
[x0,z0,u0,theta0] = InitialConditionGasLift_model_SS(par);
%states to measurement mapping function
par.nMeas = 6;
H = zeros(par.nMeas,length(z0));
H(1,1) = 1e-2*60*1e3/par.rho_o(1); %wro-oil rate from reservoir, well 1 [1e2 kg/s] --> [L/min]
H(2,2) = 1e-2*60*1e3/par.rho_o(2); %wro-oil rate from reservoir, well 2
H(3,3) = 1e-2*60*1e3/par.rho_o(3); %wro-oil rate from reservoir, well 3
H(4,7) = 1; %prh - riser head pressure well 1
H(5,8) = 1; %prh - riser head pressure well 2
H(6,9) = 1; %prh - riser head pressure well 2
par.H = H;

uArray = [0.50, 0.50, 0.50;
          0.50, 0.25, 0.50;
          0.50, 0.25, 0.75;
          0.75, 0.50, 0.75]';

%% Steady state model 
nSS = size(uArray,2);

ySSArray = [];

for ii = 1:nSS
    
    % choose inputs
    uii = [uArray(:,ii);u0(4:end)];
    
    % simulate the SSS model
    [xSS,zSS,ySS,phiSS,flag] = ErosionRigSSSimulation(x0,z0,theta0,uii,par);
    
    if ii == 1
        %intialization
        xkk = xSS;
        zkk = zSS;
    end
    
    %save the data
    ySSArray = [ySSArray, ySS];

end

clc
%% Dynamics Model 
dt = 10;
tem = 60/dt; 

nDyn = nSS*tem; %1 minute for each SS

[F_model,S_xx,S_zz,S_xz,S_xp,S_zp,x_var,z_var,u_var,p_var,difEq,algEq,L] = ErosionRigDynModel_model(par);

%for plotting
yPlot{1} = []; % SS measurements
yPlot{2} = []; % Dyn measurements

for kk = 1:nDyn

    ssCurrent = ceil(kk/tem);
    % choose inputs
    ukk = [uArray(:,ssCurrent);u0(4:end)];
    
    Fend = F_model('x0',xkk,'z0',zkk,'p',[ukk;theta0;dt]);
    
    % %extracting the results (from symbolic to numerical)
    xkk = full(Fend.xf);
    zkk = full(Fend.zf);
    
    yPlot{1} = [yPlot{1}, ySSArray(:,ssCurrent)]; % SS measurements
    yPlot{2} = [yPlot{2}, par.H*zkk]; % Dyn measurements
    
end


%% Plotting
figure(1)

% time array
time = dt:dt:dt*nDyn; 

% top pressure
tag = {'Q_{l,1}','Q_{l,2}','Q_{l,3}','P_{top,1}','P_{top,2}','P_{top,3}'};

for ii = 1:6 %ny
    subplot(2,3,ii)
    hold on
        % Steady-state value
        stairs(time, yPlot{1}(ii,:),'k:','Linewidth',1.5)

        % Dynamic value
        plot(time, yPlot{2}(ii,:),'k','Linewidth',1.5)
    hold off
    if ii == 1
        legend({'SS','Dyn'},'Location','best')
    end
    
    xlim([0, par.T*nDyn])
    xlabel('time [s]')
    ylabel(tag{ii})
    
end

