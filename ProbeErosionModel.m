function [dk_1,probeStatusk] = ProbeErosionModel(dk,probeStatusk,par)
%    Simulates the degradation of the PVA probe in the lab rig - the models
%    were obtained with rig data. See Erosion Model folder

% Inputs:
%    dk = current orifice diameter [cm]
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
% calculating random increment
dk_1 = dk - 0.0005*gamrnd(0.8,2,[3 1]);

% cannot be more degraded than 100%
for well = 1:3
    if dk_1(well) < par.dMin % all need to be lower than 1
        dk_1(well) = par.dMin;
        probeStatusk(well) = 1;
    end
end

end
