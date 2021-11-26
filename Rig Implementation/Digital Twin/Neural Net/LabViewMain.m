% Main program
% Run Initialization file first

%%%%%%%%%%%%%%%%
% Get Variables
%%%%%%%%%%%%%%%%
% disturbances 
%valve opening [-] 
cv101 = P_vector(1);
cv102 = P_vector(2);
cv103 = P_vector(3);
% if value is [A] from 0.004 to 0.020
% if you want to convert to 0 (fully closed) to 1 (fully open)
% vo_n = (vo - 0.004)./(0.02 - 0.004);

%pump rotation [%]
pRate = P_vector(4);
% if value is [A] from 0.004 to 0.020
% if you want to convert to (min speed - max speed)
% goes from 12% of the max speed to 92% of the max speed
% pRate = 12 + (92 - 12)*(P_vector(4) - 0.004)./(0.02 - 0.004);

% always maintain the inputs greater than 0.5
% inputs computed in the previous MPC iteration
% Note that the inputs are the setpoints to the gas flowrate PID's
fic104sp = P_vector(5);
fic105sp = P_vector(6); 
fic106sp = P_vector(7);
%current inputs of the plant
u0old=[P_vector(5);P_vector(6);P_vector(7)];

% cropping the data vector
nd = size(I_vector,2);
dataCrop = (nd - BufferLength + 1):nd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOST RECENT VALUE IS THE LAST ONE! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% liquid flowrates [L/min]
fi101 = I_vector(1,dataCrop);
fi102 = I_vector(2,dataCrop);
fi103 = I_vector(3,dataCrop);

% actual gas flowrates [sL/min]
fic104 = I_vector(4,dataCrop);
fic105 = I_vector(5,dataCrop);
fic106 = I_vector(6,dataCrop);

% pressure @ injection point [mbar g]
pi105 = I_vector(7,dataCrop);
pi106 = I_vector(8,dataCrop);
pi107 = I_vector(9,dataCrop);

% reservoir outlet temperature [oC]
ti101 = I_vector(10,dataCrop);
ti102 = I_vector(11,dataCrop);
ti103 = I_vector(12,dataCrop);

% DP @ erosion boxes [mbar]
dp101 = I_vector(13,dataCrop);
dp102 = I_vector(14,dataCrop);
dp103 = I_vector(15,dataCrop);

% top pressure [mbar g]
% for conversion [bar a]-->[mbar g]
% ptop_n = ptop*10^-3 + 1.01325;
pi101 = I_vector(16,dataCrop);
pi102 = I_vector(17,dataCrop);
pi103 = I_vector(18,dataCrop);

% reservoir pressure [bar g]
% for conversion [bar g]-->[bar a]
% ptop_n = ptop + 1.01325;
pi104 = I_vector(19,dataCrop);

% number of measurements in the data window
dss = size(pi104,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE GOES HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. prepare the initial run, in which the controller is built

% check for first iteration
if ~exist('kTime','var') 
    %initial condition
    [x0,u0,theta0] = InitialCondition_HAC(par);

    xEstHat = x0;
    thetaHat = theta0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regressors and Normalization %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Building Dynamic Model
    [eq,x_var,p_var] = BuildingHACModel(par); % perfect model
        
    % Building NMPC
    solverHAC = BuildingHAC(eq,x_var,p_var,par,nmpcConfig);

    % we start at 180s to make the code cleaner --> regressor of the NN that estimates the current degradation level    
    kTime = 3; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regressors and Normalization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WELL 1:
RegMatrix{1} = [(kTime*60 - BufferLength + 1):(kTime*60); % time in seconds (kTime*60 - BufferLength + 1):(kTime*60)
             fi101(end - BufferLength + 1:end);
             pi101(end - BufferLength + 1:end);
             ti101(end - BufferLength + 1:end);
             cv101*ones(1,BufferLength).*(0.02 - 0.004) + 0.004]; % in [A] --> cv doesnt change
         
RespMatrix{1} = dp101(end - BufferLength + 1:end); % in [A]         
         
% WELL 2:
RegMatrix{2} = [(kTime*60 - BufferLength + 1):(kTime*60); % time in seconds
             fi102(end - BufferLength + 1:end);
             pi102(end - BufferLength + 1:end);
             ti102(end - BufferLength + 1:end);
             cv102*ones(1,BufferLength).*(0.02 - 0.004) + 0.004]; % in [A]
         
RespMatrix{2} = dp102(end - BufferLength + 1:end); % in [A]   
         
% WELL 3:
RegMatrix{3} = [(kTime*60 - BufferLength + 1):(kTime*60); % time in seconds
             fi103(end - BufferLength + 1:end);
             pi103(end - BufferLength + 1:end);
             ti103(end - BufferLength + 1:end);
             cv103*ones(1,BufferLength).*(0.02 - 0.004) + 0.004]; % in [A]
         
RespMatrix{3} = dp103(end - BufferLength + 1:end); % in [A]   
     
% finding the normalization constants to center the data to have mean 0 and scales it to have standard 
%    deviation 1

%regressor normalization 
mu_Reg  = [mean(RegMatrix{1},2); mean(RegMatrix{2},2); mean(RegMatrix{3},2)]; 
std_Reg = [std(RegMatrix{1},[],2); std(RegMatrix{2},[],2); std(RegMatrix{3},[],2)]; 

%response normalization 
mu_Resp  = [mean(RespMatrix{1},2); mean(RespMatrix{2},2); mean(RespMatrix{3},2)]; 
std_Resp = [std(RespMatrix{1},[],2); std(RespMatrix{2},[],2); std(RespMatrix{3},[],2)]; 

% estimating the model parameters and states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating SS model parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yPlant = [fi101;
    fi102;
    fi103;
    1.01325 + 10^-3*pi101; %[mbarg]-->[bar a];
    1.01325 + 10^-3*pi102;
    1.01325 + 10^-3*pi103];

uPlant = [cv101*ones(1,dss); %workaround - i just have the last measurement here. Since it is the disturbance, it doesn't really matter;
    cv102*ones(1,dss);
    cv103*ones(1,dss);
    pi104 + 1.01325;
    ti101;
    ti102;
    ti103];

% since we run the controller with a low frequency, we assume that the last
% 10s are always at steady-state
uEst = mean(uPlant(:,end - 10:end),2);
yEst = mean(yPlant(:,end - 10:end),2);

% running casadi ErosionRigSSEstimation(xGuess,thetaGuess,uk,yk,tk,RegN,RespN,eq,x_var,p_var,par)
[thetaHat,dpHat,~,yEstHat,~,flagEst] = ErosionRigSSEstimation(xEstHat,thetaHat,uEst,yEst,1.7159,mu_Reg,std_Reg,mu_Resp,std_Resp,eq,x_var,p_var,par);

% estimation flag --> indicating if it has failed or not
if flagEst
    Estimation = 1;
else
    Estimation = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the controller %
%%%%%%%%%%%%%%%%%%%%%%%%%%
[uk,erosionSeq,inputSeq,OF,solFlag] = SolvingHAC(solverHAC,xEstHat,uEst(1:3),1.7159,uEst(4:7),thetak(1:3),thetak(4:6),mu_Reg,std_Reg,mu_Resp,std_Resp,par,nmpcConfig);

% compute new values for the valve opening setpoints
% using random values for testing
%o_new = 0.1*ones(3,1) + (0.9*ones(3,1) - 0.1*ones(3,1)).*rand(3,1);

%%%% updating time
kTime = kTime + 1;

if solFlag
    O_vector = uk';
    Optimization = 1;
else % controller failed
    O_vector = O_vector'; %dummy - doesn't do anything
    Optimization = 0;
end

SS = kTime;

Result = 0;
Parameter_Estimation = thetaHat';
State_Variables_Estimation = (par.H*xEstHat)';
State_Variables_Optimization = [0,0,0,0,0,0];
Optimized_Air_Injection = [0,0,0];

