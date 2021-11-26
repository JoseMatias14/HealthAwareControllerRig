function [u_,ehatArray,inputArray,ofValue,solFlag] = SolvingHAC(solver,x_k,u_k,t_k,u_k_fixed,res_theta,val_theta,Reg_mu,Reg_std,Resp_mu,Resp_std,par,nmpcPar)

% declare variables (bounds and initial guess)
w0 = [];
lbw =[];
ubw = [];

% declare constraints and its bounds
lbg = [];
ubg = [];

%Variable bounds
[lbx,~,~,ubx,~,~] = OptimizationBoundsHAC(par);

%% Lifting initial conditions

% initial time
w0 = [w0;t_k];
lbw = [lbw;t_k];
ubw = [ubw;t_k];

% initial input
w0 = [w0;u_k];
lbw = [lbw;u_k];
ubw = [ubw;u_k];

%% Looping through until time end
for k = 1:nmpcPar.np
    w0 = [w0; u_k];
    lbw = [lbw;nmpcPar.umin*ones(nmpcPar.nu,1)];
    ubw = [ubw;nmpcPar.umax*ones(nmpcPar.nu,1)];
%     lbw = [lbw;u_k];
%     ubw = [ubw;u_k];

        
    if k > nmpcPar.nm
        lbg = [lbg;zeros(nmpcPar.nu,1)];
        ubg = [ubg; zeros(nmpcPar.nu,1)];
    else
        lbg = [lbg;-nmpcPar.dumax*ones(nmpcPar.nu,1)];
        ubg = [ubg;nmpcPar.dumax*ones(nmpcPar.nu,1)];
    end
    
    % creating states at the evaluation points
    w0 = [w0;x_k];
    lbw = [lbw;lbx];
    ubw = [ubw;ubx];%inf
%     lbw = [lbw;-inf*ones(nmpcPar.nx,1)];
%     ubw = [ubw;inf*ones(nmpcPar.nx,1)];%inf
     
    % Calculating the residuals
    lbg = [lbg;zeros(nmpcPar.nx,1)];
    ubg = [ubg;zeros(nmpcPar.nx,1)];
         
end

% Solving the problem
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,'p',[u_k_fixed(1);u_k_fixed(2:4);res_theta;val_theta;Reg_mu;Reg_std;Resp_mu;Resp_std]);

%% Extracting solution
w_opt = full(sol.x);

% catch error
if solver.stats.success ~=1
    % solution failed 
    solFlag = 0;
    
    u_ = u_k; % fix inputs
    ofValue = -1;

    %dummy
    ehatArray = zeros(3,nmpcPar.np);
    inputArray = zeros(nmpcPar.nu,nmpcPar.np);

    
else
    % solution succeeded 
    solFlag = 1;
    
    % variable order  
    % 1: t0
    % 2-4: u0
    
    % 5-7: u1
    % 8-25: x1 --> 11:13 = w_ro_1, 23:25 dp_1
    
    % 26-28: u2
    % 29-46: x2 --> 32:34 = w_ro_2, 44:46 dp_1
    
    % 47-49: u3
    % 50-67: x4
    
    % total variables in one iteration = 3(u) + 18(xk) = 21
    
    u_ = w_opt(5:7);
    ofValue = full(sol.f);
   
    ehatArray = [];
    inputArray = [];
    for ii = 1:nmpcPar.np
        temp1 = 5 + (ii - 1)*21;
        inputArray = [inputArray, w_opt(temp1:temp1 + 2)];
                
        temp2 = 23 + (ii - 1)*21;
        dP_temp = w_opt(temp2:temp2 + 2);
        
        temp3 = 11 + (ii - 1)*21;
        w_ro_temp = w_opt(temp3:temp3 + 2);
        
        % liquid production
        %conversion
        CR = 60*10^3; % [L/min] -> [m3/s] 
        Q_pr = CR*(w_ro_temp*1e-2)./par.rho_o;

        DPNorm = [];
        for well = 1:3
            if well == 1
                DPmax = 0.3081*Q_pr(well)^2 + 0.185*Q_pr(well) - 0.7528;
                DPmin = 0.1828*Q_pr(well)^2 + 1.0093*Q_pr(well) - 4.6835;
            elseif well == 2
                DPmax = 0.2934*Q_pr(well)^2 - 0.0592*Q_pr(well) - 0.7532;
                DPmin = 0.1896*Q_pr(well)^2 + 0.8415*Q_pr(well) - 4.9335;    
            else % well == 3
                DPmax = 0.2902*Q_pr(well)^2 + 0.9971*Q_pr(well) - 1.5097;
                DPmin =  0.19*Q_pr(well)^2 + 1.1067*Q_pr(well) - 3.685;    
            end

            DPNorm =  [DPNorm; 1 - (dP_temp(well) - DPmin)./(DPmax - DPmin)];
        end
        
        ehatArray = [ehatArray, DPNorm];
               
    end

end

end

