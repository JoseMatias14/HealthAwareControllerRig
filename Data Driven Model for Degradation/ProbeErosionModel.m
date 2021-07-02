function [dP_k_1,eProbe_k_1,status] = ProbeErosionModel(Ql_k,eProbe_k,DeltaT,erosionEvolution)
%    Simulates the degradation of the PVA probe in the lab rig - the models
%    were obtained with rig data. See Erosion Model folder

% Inputs:
%    Ql_k = current well liquid flowrate [L/min]
%    t_k = current experiment time [min]
%    DeltaT = system sampling time [min] 
%    dP_k = past delta pressure along the erosion chamber [mbar]
%    eProbe_k = past probe degradation stage [0-1]
%    erosionEvolution = choosing type of model used in the erosion

% Outputs:
%   dP_k_1 = future delta pressure along the erosion chamber [mbar]
%   eProbe_k_1 = future probe degradation stage [0-1]
%   status = probe healht [ 0 = healthy | 1 = degraded ]

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Matias
% email: jose.o.a.matias@ntnu.no
% June 2021; Last revision: 
     
import casadi.*


% dPModel2 = -4.45057218353704 + 0.150265436491970*tData + 1.88054326803212*QlData -0.0503621574689086*tData.*QlData + 0.213170980385915*QlData.^2;

%% Models
switch erosionEvolution
    case 'deterministic'
        
        % calculating min and max delta pressure value
        dPmin = 0.2788*Ql_k.^2 + 1.143*Ql_k - 2.3831;
        dPmax = 0.2693*Ql_k.^2 + 0.0504*Ql_k + 0.1256;
      
        % calculating DP change in the current point
        % dPModel = a + b*tData + c*QlData + d*tData.*QlData + e*QlData.^2;
        % ddPModel_dt = b + d*QlData;
        ddPModel_dt = 0.150265436491970 - 0.0503621574689086*Ql_k;
        
        % calculating the current dPvalue based on the probe degradation
        dP_true = ddPModel_dt*DeltaT + dPmin + eProbe_k.*(dPmax - dPmin);
       
        % calculating the future dPvalue based on the probe degradation
        eProbe_k_1 = (dP_true - dPmin)./(dPmax - dPmin);
        
        % cannot be more degraded than 100%
        for well = 1:3
            if eProbe_k_1(well) > 1 % all need to be lower than 1
                eProbe_k_1(well) = 1;
                dP_k_1(well) = dPmax(well);
            end
        end
        
        % adding noise to the dP measurement
        dP_k_1 = dP_true + 0.0001*randn(3,1);
       
        status = [];
        for well = 1:3
            if eProbe_k_1(well) < 1 % all need to be lower than 1
                status = [status; 0];
            else
                status = [status; 1];
            end
        end
        
    case 'deterministicWithBreak'
        disp('zero')
        
    case 'randomIncrements'
        disp('positive one')
end


end
