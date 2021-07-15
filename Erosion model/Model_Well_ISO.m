close all
clear all
clc

%liquid flowrate - pilling up the three wells
Ql = [2:10;2:10;2:10]; % [L/min]

% DP [mbar]
DPmin = [1.36, 3.20, 4.65, 6.76, 9.77, 13.13, 16.91, 21.56, 26.86;
         0.82, 1.89, 3.53, 5.80, 9.23, 12.59, 16.42, 20.64, 26.40;
         1.17, 2.99, 4.54, 7.42, 9.44, 13.57, 18.26, 22.53, 27.34];

DPmax = [2.33, 4.00, 6.07, 9.13, 13.22, 16.52, 21.66, 26.01, 32.31;
         0.61, 2.79, 4.51, 7.76, 11.81, 14.48, 19.67, 25.17, 30.17;
         1.40, 3.52, 6.22, 9.67, 14.42, 20.31, 24.76, 30.33, 36.69];

  
Ql_data = reshape(Ql,1,27);
DPmax_data = reshape(DPmax,1,27);
DPmin_data = reshape(DPmin,1,27);

%% Modeling new probe (DPmax)
%%%%%%%%%%%%%%%%%%%%%
% Model with Casadi %
%%%%%%%%%%%%%%%%%%%%%
import casadi.*

% Declare variables to the model (parameters and controls)
Q_var = MX.sym('Q_var');
Cd = MX.sym('Cd');
d = MX.sym('d');

rho = 1000;
D = 2; %[cm]
DP_var = Cd*((Q_var/(60*1000))*rho)^2*(1 - (d/D)^4)/((pi()/4)^2*2*rho*(0.01*d)^4)/100;

% Create a function that computes the delta pressure
quadFun = Function('quad',{Q_var, Cd, d},{DP_var});

%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%% %
%simulating the system for all the inputs at the same time
N = size(Ql_data,1);%number of points used in the estimation procedure - changes
all_samples = quadFun.map(N);

Y_symbolic = all_samples(Ql_data, repmat(Cd,1,N), repmat(0.3,1,N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating model parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = DPmax_data - Y_symbolic;

%formalizing nlp structure and nlp problem
nlp = struct('x', Cd, 'f', 0.5*dot(e,e));
solver = nlpsol('solver','ipopt', nlp);

% Solve
sol = solver('x0',0.7,'lbx',0,'ubx',10);
Cd_hat = full(sol.x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting to check the functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
for well = 1:3
    subplot(3,1,well)
        X = Ql(well,:);
        Ydata = DPmax(well,:);
        Ymodel = full(quadFun(X,Cd_hat,0.3));
        
        plot(X,Ydata,'ro')
        hold on 
        plot(X,Ymodel,'kx')
        legend({'Data','Model'})
        xlabel('Ql [L/min]') ; ylabel('dP [mbar]'); title('Estimating Cd')
end


%% Finding the orifice size relative to dPmin
%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%% %
Yd_symbolic = all_samples(Ql_data, repmat(Cd_hat,1,N), repmat(d,1,N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating model parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = DPmin_data - Yd_symbolic;

%formalizing nlp structure and nlp problem
nlp = struct('x', d, 'f', 0.5*dot(e,e));
solver = nlpsol('solver','ipopt', nlp);

% Solve
sol = solver('x0',30,'lbx',0,'ubx',30);
d_max_hat = full(sol.x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting to check the functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
for well = 1:3
    subplot(3,1,well)
        X = Ql(well,:);
        Ydata = DPmax(well,:);
        Ymodel = full(quadFun(X,Cd_hat,d_max_hat));
        
        plot(X,Ydata,'ro')
        hold on 
        plot(X,Ymodel,'kx')
        legend({'Data','Model'})
        xlabel('Ql [L/min]') ; ylabel('dP [mbar]'); title('Estimating d_max')
end

