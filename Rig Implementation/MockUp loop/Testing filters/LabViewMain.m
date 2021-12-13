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

%%%%%%%%%%%%%%%%%%
% CODE GOES HERE %
%%%%%%%%%%%%%%%%%%
% check for first iteration
if ~exist('kCheck','var') 
    %initial condition
    [xHat,zHat,u0,thetaHat] = InitialConditionHAC(par);
    O_vector = u0(1:3);

    % Building Dynamic Model
    [diff,alg,x_var,z_var,p_var] = BuildingDynModel(par); 

    % Building NMPC
    solverHAC = BuildingHAC(diff,alg,x_var,z_var,p_var,par,nmpcConfig);

    % first iteration indication    
    kCheck = 1; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diagnostics - estimating the states and the current degradation level %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yPlant = [fi101;
    fi102;
    fi103;
    1.01325 + 10^-3*pi101; %[mbarg]-->[bar a];
    1.01325 + 10^-3*pi102;
    1.01325 + 10^-3*pi103;
    dp101;
    dp102;
    dp103];

uPlant = [cv101*ones(1,dss); %workaround - i just have the last measurement here. Since it is the disturbance, it doesn't really matter;
    cv102*ones(1,dss);
    cv103*ones(1,dss);
    pi104 + 1.01325];

% since we run the controller with a low frequency, we assume that the last
% 10s are always at steady-state
uEst = mean(uPlant,2);
yEst = mean(yPlant,2);
% uEst = uPlant;
% yEst = yPlant;

[thetaHat,zHat,Psi,flagEst] = HAC_SSEstimation(zHat,thetaHat,uEst,yEst,par);

xHatArray = [];
for jj = 1:dss
    xHat = DiameterEstimation(xHat,yPlant(1:3,jj),yPlant(7:9,jj),par);
    xHatArray = [xHatArray, xHat];
end
xHatFiltered = mean(xHatArray,2);

xHat_MI = [];
for well = 1:3
    % fitting monotonic increasing line
    n = length(1:dss);
    C = eye(n);
    D = xHatArray(well,:);
    A = diag(ones(n,1),0) - diag(ones(n-1,1),1);
    A(end,:) = [];
    b = zeros(n-1,1);
    
    opts = optimset('lsqlin');
    opts.LargeScale = 'off';
    opts.Display = 'none';
    xhat_temp = lsqlin(C,D,A,b,[],[],[],[],[],opts);
    
    xHat_MI = [xHat_MI; xhat_temp'];

end

xHat_mean = DiameterEstimation(xHat,yEst(1:3),yEst(7:9),par);

%figure(1)
%clf;
for well = 1:3
    %subplot(3,1,well)
        % raw data
        %plot(1:dss,xHatArray(well,:),'ok');
        hold on
        
        % median filter
        temp = medfilt1(xHatArray(well,:),20);
        %plot(20:dss,temp(20:end),'b');
        
        % monotonic increasing function
        %plot(1:dss,xHat_MI(well,:),'k');
        
        % mean of the estimates
        temp2 = mean(xHatArray(well),2);
        %yline(temp2,'r');
        
        % estimates using mean data
        %yline(xHat_mean(well),'g');
        
        % saving data
        diameterTemp{well} = [diameterTemp{well}, [xHatArray(well,end);temp(end);xHat_MI(well,end);temp2;xHat_mean(well)]];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving Health-aware controller %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[uk,solFlag,uHacSeq,dHacSeq,~,~] = SolvingHAC(solverHAC,xHatFiltered,[zHat;yEst(7:9)],uEst(1:3),uEst(4),thetaHat(1:3),thetaHat(4:6),nmpcConfig);

% compute new values for the valve opening setpoints
% using random values for testing
%o_new = 0.1*ones(3,1) + (0.9*ones(3,1) - 0.1*ones(3,1)).*rand(3,1);
if solFlag
    O_vector = uk;
    Optimization = 1;
else % controller failed
    %O_vector = O_vector'; %dummy - doesn't do anything
    Optimization = 0;
end

SS = 0;
Estimation = 0;
Result = 0;
Parameter_Estimation = [thetaHat(1),thetaHat(2),thetaHat(3),thetaHat(4),thetaHat(5),thetaHat(6)];
State_Variables_Estimation = [xHat(1),xHat(2),xHat(3),0,0,0];
State_Variables_Optimization = [0,0,0,0,0,0];
Optimized_Air_Injection = [0,0,0];

