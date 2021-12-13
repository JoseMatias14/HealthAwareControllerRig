clear 
clc
close all

%% Import the data
name = [pwd,'/2021-07-09_101557_Exp5.txt'];

%% IMPORT LABVIEW DATA
opts = delimitedTextImportOptions("NumVariables", 36);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

opts.VariableNames = ["Date", "Time", "Relativetimes", "FIC104setptslmin", "FIC104slmin", "FIC105setptslmin", "FIC105slmin", "FIC106setptslmin", "FIC106slmin", "CV101setptlmin", "FI101lmin", "CV102setptlmin", "FI102lmin", "CV103setptlmin", "FI103lmin", "dP101mbarD", "dP102mbarD", "dP103mbarD", "CV107setptmbarG", "PI101mbarG", "CV108setptmbarG", "PI102mbarG", "CV109setptmbarG", "PI103mbarG", "PumpoutputpressuresetptbarG", "PI104barG", "TI101C", "TI102C", "TI103C", "CV101currentA", "CV102currentA", "CV103currentA", "CV107currentA", "CV108currentA", "CV109currentA", "PumpctrlcurrentA"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Date", "Time"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["Date", "Time", "Relativetimes", "FIC104setptslmin", "FIC104slmin", "FIC105setptslmin", "FIC105slmin", "FIC106setptslmin", "FIC106slmin", "CV101setptlmin", "FI101lmin", "CV102setptlmin", "FI102lmin", "CV103setptlmin", "FI103lmin", "dP101mbarD", "dP102mbarD", "dP103mbarD", "CV107setptmbarG", "PI101mbarG", "CV108setptmbarG", "PI102mbarG", "CV109setptmbarG", "PI103mbarG", "PumpoutputpressuresetptbarG", "PI104barG", "TI101C", "TI102C", "TI103C", "CV101currentA", "CV102currentA", "CV103currentA", "CV107currentA", "CV108currentA", "CV109currentA", "PumpctrlcurrentA"], "DecimalSeparator", ",");
opts = setvaropts(opts, ["Date", "Time"], "ThousandsSeparator", ".");

% Import the data
expData = readtable(name, opts);

% Convert to output type
labviewData = table2array(expData);

% Clear temporary variables
clear opts

%% preparing data

% number of data points used in the estimation
nData = 60;

% liquid flowrates [L/min]
fi101 = labviewData(1:nData,11)';
fi102 = labviewData(1:nData,13)';
fi103 = labviewData(1:nData,15)';

% pressure @ injection point [mbar g]
pi101 = labviewData(1:nData,20)';
pi102 = labviewData(1:nData,22)';
pi103 = labviewData(1:nData,24)';

% DP @ erosion boxes [mbar]
dp101 = labviewData(1:nData,16)';
dp102 = labviewData(1:nData,17)';
dp103 = labviewData(1:nData,18)';

% reservoir pressure [bar g]
% for conversion [bar g]-->[bar a]
% ptop_n = ptop + 1.01325;
pi104 = labviewData(1:nData,26)';

% disturbances 
%valve opening [-] 
cv101 = (labviewData(1:nData,30)' - 0.004)./(0.02 - 0.004);
cv102 = (labviewData(1:nData,31)' - 0.004)./(0.02 - 0.004);
cv103 = (labviewData(1:nData,32)' - 0.004)./(0.02 - 0.004);

%% calculating the initial condition
%plant parameters
par = ParametersGasLiftModel;
par.T = 60; % simulation sampling time[s]

%initial condition
[xHat,zHat,u0,thetaHat] = InitialConditionHAC(par);
O_vector = u0(1:3);

%% estimating intial values
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

uPlant = [cv101; %workaround - i just have the last measurement here. Since it is the disturbance, it doesn't really matter;
    cv102;
    cv103;
    pi104 + 1.01325];

% since we run the controller with a low frequency, we assume that the last
% 10s are always at steady-state
uEst = mean(uPlant,2);
yEst = mean(yPlant,2);

[thetaHat,zHat,Psi,flagEst] = HAC_SSEstimation(zHat,thetaHat,uEst,yEst,par);

% initial diameter
xEstHat = 0.3*ones(par.n_w,1);
fileName = '2021-07-09_101557_Exp5';

GeneratingInitialGuess(thetaHat,xEstHat,zHat,uEst,par,fileName);
