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

%coefficient of discharge
par.Cd_hat = 0.00963481031516700;

% minimum equivalent diameter
par.dMin = 0.3*0.5;
par.dMax = 0.3181*1.5;

%% Modeling new probe (DPmax)
%%%%%%%%%%%%%%%%%%%%%
% Model with Casadi %
%%%%%%%%%%%%%%%%%%%%%
import casadi.*

% Declare variables to the model (parameters and controls)
Q_var = MX.sym('Q_var');
d = MX.sym('d');
rho = 1000;
D = 2; %[cm]

DP_var = par.Cd_hat*((Q_var/(60*1000))*rho)^2*(1 - (d/D)^4)/((pi()/4)^2*2*rho*(0.01*d)^4)/100;

% Create a function that computes the delta pressure
dP_fun = Function('dP_fun',{Q_var, d},{DP_var});

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

for exp = 1:par.ni 

    expDataTemp = table2array(readtable(par.nfolder(exp).name, opts));
    expTime = size(expDataTemp,1);

    for cc = 1:par.nw
        data{exp,cc}.flowrate = expDataTemp(:, 11 + 2*(cc - 1));
        data{exp,cc}.deltaP = expDataTemp(:, 16 + (cc - 1));
        data{exp,cc}.temperature = expDataTemp(:, 27 + (cc - 1));
        data{exp,cc}.ptop = expDataTemp(:, 20 + 2*(cc - 1));
        data{exp,cc}.ppump = expDataTemp(:, 26);
        data{exp,cc}.pumpRotation = (12 + (92 - 12)*(expDataTemp(:, 36) - 0.004)./(0.02 - 0.004))./100; % converting from mA to [0 - 1]
        data{exp,cc}.valveOpen = (expDataTemp(:, 30 + (cc - 1)) - 0.004)./(0.02 - 0.004); % converting from mA to [0 - 1]

        data{exp,cc}.equivDiame = zeros(expTime,1);
        
        for kk = 1:expTime
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimating model parameters % 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dP_sym = dP_fun(data{exp,cc}.flowrate(kk),d);
            
            if cc == 1
                temp2 = data{exp,cc}.deltaP(kk) - 1.5 + 0.7;
            elseif cc == 2
                temp2 = data{exp,cc}.deltaP(kk) + 1.7 + 0.24;
            else
                temp2 = data{exp,cc}.deltaP(kk) + 0 + 0.20;
            end
            e = temp2 - dP_sym;
            
            %formalizing nlp structure and nlp problem
            % Create an NLP solver
            opts2 = struct;
            %opts2.ipopt.max_iter = nmpcPar.maxiter;
            opts2.ipopt.print_level = 0;
            %opts2.ipopt.tol = nmpcPar.tol;
            %opts2.ipopt.acceptable_tol = 100*nmpcPar.tol; % optimality convergence tolerance
            %opts2.ipopt.linear_solver = 'mumps';
            nlp = struct('x', d, 'f', 0.5*dot(e,e));
            solver = nlpsol('solver','ipopt', nlp,opts2);
            
            % Solve
            if kk == 1
                xGuess = 0.3;
            else
                xGuess = temp; 
            end
            
            sol = solver('x0',xGuess,'lbx',par.dMin,'ubx',par.dMax);
            temp = full(sol.x);
            data{exp,cc}.equivDiame(kk) = temp;
                     
        end
        
        data{cc}.units = {'Q [L/min] ','dP [mbar]','T [oC]','P_{top} [mbar g]','P_{pump} [bar g]','v_{open} [0-1]','P_{rot} [0-1]','d_{eq} [cm]'}; 
       
        % new plot
        figure(exp)

        tgrid = 0:(1/60):(1/60)*(expTime - 1);
        subplot(3,1,cc)
            plot(tgrid,data{exp,cc}.equivDiame,'k','Linewidth',1.5)
            grid on

            xticks(0:5:tgrid(end))
            xlim([0, tgrid(end)])
            
%             yticks(0:0.25:1)
             ylim([par.dMin/0.5, par.dMax/1.5])
            
            title(['Well ',num2str(cc)])
            ylabel('d [cm]','FontSize',10)
            xlabel('time [min]','FontSize',10);
        
    end
   
end

save('Experiments','data')
