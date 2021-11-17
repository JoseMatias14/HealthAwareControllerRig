% This code computes the change in the probe area along the experiment and compares it with the pressure drop. 


% Author: Jose Otavio Matias
% Work address
% email: jose.o.a.matias@ntnu.no
% Novemebr 2020; Last revision: 
clear 
close all
clc
% 

%% Parameters tuning 
%number of wells being used in the experiment
par.nw = 3;

% Sampling times 
par.dtMeasu = 1; %[s]

%%  getting the files
par.currentdirectory = pwd;
par.nfolder = dir([par.currentdirectory,'\*.txt']);%pictures taken with png format
par.ni = size(par.nfolder,1); %number of experiments

% reading data
opts = delimitedTextImportOptions("NumVariables", 36);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Date", "Time", "Relativetimes", "FIC104setptslmin", "FIC104slmin", "FIC105setptslmin", "FIC105slmin", "FIC106setptslmin", "FIC106slmin", "CV101setptlmin", "FI101lmin", "CV102setptlmin", "FI102lmin", "CV103setptlmin", "FI103lmin", "dP101mbarD", "dP102mbarD", "dP103mbarD", "CV107setptmbarG", "PI101mbarG", "CV108setptmbarG", "PI102mbarG", "CV109setptmbarG", "PI103mbarG", "PumpoutputpressuresetptbarG", "PI104barG", "TI101C", "TI102C", "TI103C", "CV101currentA", "CV102currentA", "CV103currentA", "CV107currentA", "CV108currentA", "CV109currentA", "PumpctrlcurrentA"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Date", "Time"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["Date", "Time", "Relativetimes", "FIC104setptslmin", "FIC104slmin", "FIC105setptslmin", "FIC105slmin", "FIC106setptslmin", "FIC106slmin", "CV101setptlmin", "FI101lmin", "CV102setptlmin", "FI102lmin", "CV103setptlmin", "FI103lmin", "dP101mbarD", "dP102mbarD", "dP103mbarD", "CV107setptmbarG", "PI101mbarG", "CV108setptmbarG", "PI102mbarG", "CV109setptmbarG", "PI103mbarG", "PumpoutputpressuresetptbarG", "PI104barG", "TI101C", "TI102C", "TI103C", "CV101currentA", "CV102currentA", "CV103currentA", "CV107currentA", "CV108currentA", "CV109currentA", "PumpctrlcurrentA"], "DecimalSeparator", ",");
opts = setvaropts(opts, ["Date", "Time"], "ThousandsSeparator", ".");

for expe = 1:par.ni 

    expDataTemp = table2array(readtable(par.nfolder(expe).name, opts));
    expTime = size(expDataTemp,1);

    for cc = 1:par.nw
        data{expe,cc}.relTime = expDataTemp(:, 3);
        data{expe,cc}.flowrate = expDataTemp(:, 11 + 2*(cc - 1));
        data{expe,cc}.deltaP = expDataTemp(:, 16 + (cc - 1));
        data{expe,cc}.ptop = expDataTemp(:, 20 + 2*(cc - 1));
        data{expe,cc}.temperature = expDataTemp(:, 27 + (cc - 1));
        data{expe,cc}.ppump = expDataTemp(:, 26);
        data{expe,cc}.valvCurrent = expDataTemp(:, 30 + (cc - 1));
        data{expe,cc}.pumpCurrent = expDataTemp(:, 36);

        data{cc}.units = {'t [s]','Q [L/min] ','dP [mbar]','T [oC]','P_{top} [mbar g]','P_{pump} [bar g]','v_{open} [mA]','P_{rot} [mA]'}; 
    end
end

save('ExperimentsNN','data')





