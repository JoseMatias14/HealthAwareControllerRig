function [dk_1,probeStatusk] = ProbeErosionModel(dk,Qk,probeStatusk,par)
%    Simulates the degradation of the PVA probe in the lab rig - the models
%    were obtained with rig data. See Erosion Model folder

% Inputs:
%    dk = current orifice diameter [cm]
%    Qk = current system flowrate [L/min]
%    probeStatusk_1 = current probe status
%    par = system parameters

% Outputs:
%   dk_1 = future orifice diameter [cm]
%   probeStatusk_1 = probe healht [ 0 = healthy | 1 = degraded ]

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Matias
% email: jose.o.a.matias@ntnu.no
% June 2021; Last revision: 
     
import casadi.*

%% Models



for well = 1:3
    % calculating random increment
    %state dependent time evolution
    if Qk(well) < 5
        alpha = 0.5;
    else
        alpha = 0.0043*Qk(well)^3 - 0.0949*Qk(well)^2 + 0.7305*Qk(well) - 1.32;
    end
    
    dk_1(well,1) = dk(well) - 0.0005*gamrnd(alpha,2);

    % cannot be more degraded than 100%
    if dk_1(well,1) < par.dMin % all need to be lower than 1
        dk_1(well,1) = par.dMin;
        probeStatusk(well) = 1;
    end
end

end
