% Simulation emulating the degradation of a PVA probe in the experimental
% rig

% Other m-files required: SmoothPlantParameters2.mat; Disturbances2.mat; RigNoise.mat; InitialState_2021-03-03_104459_NoRTO_test2_1.mat
% MAT-files required: 
%       For Mockup Labview Interface:
%           1. InitializationLabViewMain.m
%           2. LabViewMain.m

%       For "Plant":
%           1. ErosionRigDynModel.m

%       Bounds, Parameters and Initial conditon:
%           1. InitialConditionGasLift.m
%           2. OptimizationBoundsGasLiftRiser2.m
%           3. ParametersGasLiftModel.m
clear
close all
clc

%saving data
name = 'HAC_multiscenario_MC_30_10';
%noise seed
rng('default')

%% Loading .mat files
% previously defined disturbance profiles
% disturbances = load('Disturbances3'); 

% previously computed system parameter profiles - used as "plant" true
% parameters
parProfile = load('SmoothPlantParameters3'); 

% rig noise characteristics - compute previously from actual rig data
noise = load('RigNoise'); 

%% Simulation tuning
%plant parameters
parPlant = ParametersGasLiftModel;
parPlant.T = 60; % simulation sampling time[s]

% Buffer length
BufferLength = 60/parPlant.T; % measurement buffer contains 60 measurements (one per second)

% Controller sampling time
nExec = 60/parPlant.T; %--> executes every 60s 

%simulation parameters
nInit = 0; %[s]
nFinal = 50*parPlant.T; %[sampling time] - arbitrarily chosen
tgrid = (nInit:parPlant.T:nFinal)/parPlant.T; %[min] one measurements per second

%initial condition
[xPlant0,zPlant0,uPlant0,thetaPlant0] = InitialConditionGasLift(parPlant);

%states to measurement mapping function
parPlant.nMeas = 6;
parPlant.H = zeros(6,length(zPlant0));
parPlant.H(1,1) = 1e-2*60*1e3/parPlant.rho_o(1); %wro-oil rate from reservoir, well 1 [1e-2 kg/s] --> [L/min]
parPlant.H(2,2) = 1e-2*60*1e3/parPlant.rho_o(2); %wro-oil rate from reservoir, well 2
parPlant.H(3,3) = 1e-2*60*1e3/parPlant.rho_o(3); %wro-oil rate from reservoir, well 3
parPlant.H(4,7) = 1; %prh - riser head pressure well 1
parPlant.H(5,8) = 1; %prh - riser head pressure well 2
parPlant.H(6,9) = 1; %prh - riser head pressure well 3
parPlant.H(7,13) = 1; %dp - well 1
parPlant.H(8,14) = 1; %dp - well 2
parPlant.H(9,15) = 1; %dp - well 3

%% Run configuration file
InitializationLabViewMain %here we use the same syntax as in the rig

%% 
% number of Monte Carlo iterations
nMC = 100;

% Monte Carlo assessment
for mm = 1:nMC
    
    %% Initializing simulation
    % Plant states
    xk = xPlant0;
    zk = zPlant0;
    probeStatusk = [0;0;0]; % flag --> 0 = healthy | 1 = degraded
    dk = parPlant.dMin*ones(3,1); %initial probe diameter
    
    % Inputs
    uk = [0.5;     %CV-101 opening [-]
        0.5;     %CV-102 opening [-]
        0.5;     %CV-103 opening [-]
        1.29833006756757];    %PI-104 [bar]
    
    % "Setpoints" for the gas lift controllers
    O_vector = [uk(1); uk(2); uk(3)];
    
    % "Plant" parameters
    thetak = thetaPlant0;
    
    %% Run mock-up loop
    % arrays for plotting
    %%%%%%%%%%%%%%%%
    % "Plant" Data %
    %%%%%%%%%%%%%%%%
    xPlantArray{mm} = xk;
    zPlantArray{mm} = zk;
    uPlantArray{mm} = uk;
    measPlantArray{mm} = [parPlant.H*zk];
    thetaPlantArray{mm} = thetak;
    ofPlantArray{mm} = 20*(measPlantArray{mm}(1)) + 10*(measPlantArray{mm}(2)) + 30*(measPlantArray{mm}(3));
    probeStatusArray{mm} = probeStatusk;
    probeOrificeArray{mm} = dk;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Production Optimization Methods Data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % we save the same data that is saved in the Rig's labview code
    flagArray{mm} = []; % flag SSD, Model Adaptation and Economic Optimization [flag == 0 (failed), == 1 (success)]
    ofArray{mm} = [];   % OF value computed bu the Economic Optimization
    thetaHatArray{mm} = []; % estimated parameters
    yEstArray{mm} = []; % model prediction with estimated states
    yOptArray{mm} = []; % model prediction @ new optimum
    uOptArray{mm} = []; % computed inputs (u_k^\star)
    uImpArray{mm} = []; % filtered inputs to be implemented (u_{k+1})
    controlTimeArray{mm} = []; % execution time
    
    for kk = 1:nFinal/parPlant.T
        
        % printing the loop evolution in minutes
        fprintf('     kk >>> %6.4f [min]\n',tgrid(kk + 1))
        
        % simulate the SS model
        [xk,zk,yk,~,~] = ErosionRigSSPlant(xk,zk,thetak,uk,dk,parPlant);
        
        % evolving probe degradation
        [dk,probeStatusk] = ProbeErosionModel(dk,yk(1:3),probeStatusk,parPlant);
        
        % saving the results
        xPlantArray{mm} = [xPlantArray{mm}, xk];
        zPlantArray{mm} = [zPlantArray{mm}, zk];
        measPlantArray{mm} = [measPlantArray{mm}, yk]; %adding artificial noise to the measurements  + noise.output*randn(6,1)
        ofPlantArray{mm} = [ofPlantArray{mm}, 20*(measPlantArray{mm}(1,end)) + 10*(measPlantArray{mm}(2,end)) + 30*(measPlantArray{mm}(3,end));];
        probeStatusArray{mm} = [probeStatusArray{mm}, probeStatusk];
        probeOrificeArray{mm} = [probeOrificeArray{mm}, dk];
        
        % we execute the production optimization:
        % a. after the initial [BufferLength]-second buffer
        % b. every [nExec] seconds
        if kk > BufferLength && rem(kk,nExec) == 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Rearranging the data vectors and units here, so they can the be exactly %
            % the same as in the actual rig. The goal is that LabViewRTO.m can be     %
            % directly plug in the Labview interface and it will work                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % measurement buffer (dim = 19 X BufferLength)
            I_vector = [measPlantArray{mm}(1,(kk - BufferLength + 2):kk + 1); % FI-101 [L/min]
                measPlantArray{mm}(2,(kk - BufferLength + 2):kk + 1); % FI-102 [L/min]
                measPlantArray{mm}(3,(kk - BufferLength + 2):kk + 1); % FI-103 [L/min]
                ones(3,BufferLength);                             % FI-104 [sL/min] & FI-105 [sL/min] & FI-106 [sL/min]. Not used here
                ones(3,BufferLength); %dummy values --> in the actual rig, they will the pressure at injection point (PI105, PI106, PI107). Not used here
                ones(3,BufferLength); %dummy values --> in the actual rig, they will the well temperature (TI101, TI102, TI103). Not used here
                (measPlantArray{mm}(7,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % dPI-101 [mbar]
                (measPlantArray{mm}(8,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % dPI-102 [mbar]
                (measPlantArray{mm}(9,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % dPI-103 [mbar]
                (measPlantArray{mm}(4,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % PI-101 [mbar g]
                (measPlantArray{mm}(5,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % PI-102 [mbar g]
                (measPlantArray{mm}(6,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % PI-103 [mbar g]
                uPlantArray{mm}(4,(kk - BufferLength + 1):kk) - 1.01325];     % PI-104 [bar g]
            
            
            % values of the input variables at the previus rig sampling time (dim = nu[7] X 1)
            P_vector = [uk(1);  % CV101 opening [-]
                uk(2);  % CV102 opening [-]
                uk(3);  % CV103 opening [-]
                1;      % dummy values: in the actual rig, they will the pump rotation. Not used here
                1;      % dummy values: FI-104 [sL/min]
                1;      % dummy values: FI-105 [sL/min]
                1];     % dummy values: FI-106 [sL/min]
            
            % values of the inputs (valve opening) of the last optimization run (dim = nQg[3] X 1)
            %O_vector = uPlantArray(1:3,kk - nExec);
            
            tic
            
            % Run Labview/Matlab interface file
            LabViewMain
            
            % computing execution time
            controlTime = toc;
            
            flagArray{mm} = [flagArray{mm}, [SS;Estimation;Optimization]];
            ofArray{mm} = [ofArray{mm}, Result];
            thetaHatArray{mm} = [thetaHatArray{mm}, Parameter_Estimation'];
            yEstArray{mm} = [yEstArray{mm}, State_Variables_Estimation'];
            yOptArray{mm} = [yOptArray{mm}, State_Variables_Optimization'];
            uOptArray{mm} = [uOptArray{mm}, Optimized_Air_Injection'];
            
            uImpArray{mm} = [uImpArray{mm}, O_vector];
            controlTimeArray{mm} = [controlTimeArray{mm}, controlTime];
            
        else
            % update with dummy values
            flagArray{mm} = [flagArray{mm}, [0;0;0]];
            ofArray{mm} = [ofArray{mm}, 0];
            thetaHatArray{mm} = [thetaHatArray{mm}, zeros(1,6)'];
            yEstArray{mm} = [yEstArray{mm}, zeros(1,6)'];
            yOptArray{mm} = [yOptArray{mm}, zeros(1,6)'];
            uOptArray{mm} = [uOptArray{mm}, zeros(1,3)'];
            
            uImpArray{mm} = [uImpArray{mm}, zeros(1,3)'];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % updating input and parameter vectors %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Saving setpoints
        uk = [O_vector(1);  %CV-101 opening [-]
            O_vector(2);  %CV-102 opening [-]
            O_vector(3);  %CV-103 opening [-]
            uPlantArray{mm}(4,end)]; %PI-104 [bar]
        
        uPlantArray{mm} = [uPlantArray{mm}, uk];
        
        % Updating plant parameters according to pre-computed array
        %   the values are updated every 10s
        thetak = parProfile.thetaPlant(:,kk*6);
        thetaPlantArray{mm} = [thetaPlantArray{mm}, thetak];
        
        
        % emulating system breaking if threshold is achieved
        if ~all(dk < nmpcConfig.x_threshold)
            
            % saving breakdown time
            tBreak = kk;
            % updating tgrid
            tgridArray{mm} = (nInit:parPlant.T:(tBreak*parPlant.T))/parPlant.T; %[min] one measurements per second
            
            %breaking loop
            break
        end
        
    end
end

%%
save(name,'flagArray','ofArray','thetaHatArray','yEstArray','yOptArray','uOptArray','uImpArray','xPlantArray','zPlantArray','measPlantArray','ofPlantArray','probeStatusArray','probeOrificeArray','tgridArray','controlTimeArray'); 

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% plotting the data 
% checking sampling rate
markers = {'o','x','>'};
cc = {'b','k','g'};
leg = {'w_1','w_2','w_3'};
 
%% Inputs (valve opening)
f = figure(1);

for well = 1:3
    subplot(3,1,well)
    hold on
    for mm = 1:nMC
        plot(tgridArray{mm}, uPlantArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
    end
    ylim([0.05 0.95])
    
    xticks(0:(10*parPlant.T/60):(1/60)*(nFinal))
    xlim([0 (1/60)*(nFinal)])
    
    xlabel('time [min]','FontSize',10)
    ylabel('v_o [-]','FontSize',10)
    
    name = ['Valve opening - Well ',num2str(well)];
    title(name,'FontSize',10)
end


%% Outputs (Liquid Production)
f = figure(2);
for well = 1:3
     % Reservoir Valve Parameters
    subplot(3,1,well)
    hold on
    for mm = 1:nMC
        plot(tgridArray{mm}, measPlantArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
    end
        ylim([2 11])
        ylabel('Q_l [L/min]','FontSize',10)

        xticks(0:(10*parPlant.T/60):(1/60)*(nFinal))
        xlim([0 (1/60)*(nFinal)])
        xlabel('time [min]','FontSize',10)

        name = ['Measurements - Well ',num2str(well)];
        title(name,'FontSize',10)  
        
end


%% Outputs (dP)
f = figure(3);
for well = 1:3
    % Reservoir Valve Parameters
    subplot(3,1,well)
    hold on
    for mm = 1:nMC
        plot(tgridArray{mm}, measPlantArray{mm}(6 + well,:),'Color',[0,0,0,0.05],'Linewidth',1)
    end
    %ylim([2 11])
    ylabel('dP [mbar]','FontSize',10)
    
    xticks(0:(10*parPlant.T/60):(1/60)*(nFinal))
    xlim([0 (1/60)*(nFinal)])
    xlabel('time [min]','FontSize',10)
    
    name = ['Measurements - Well ',num2str(well)];
    title(name,'FontSize',10)
        
end

%% Outputs (dP)
f = figure(4);
for well = 1:3
    % Reservoir Valve Parameters
    subplot(3,1,well)
    hold on
    for mm = 1:nMC
        plot(tgridArray{mm}, probeOrificeArray{mm}(well,:),'Color',[0,0,0,0.05],'Linewidth',1)
    end
    yline(nmpcConfig.x_threshold,'r:','Linewidth',1);
    
    ylim([parPlant.dMin 0.32])
    yticks(parPlant.dMin:0.0036:parPlant.dMax)
    ylabel('Orifice diameter','FontSize',10)
    
    xticks(0:(10*parPlant.T/60):(1/60)*(nFinal))
    xlim([0 (1/60)*(nFinal)])
    xlabel('time [min]','FontSize',10)
    
    name = ['Measurements - Well ',num2str(well)];
    title(name,'FontSize',10)
        
end
 