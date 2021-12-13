function [u_,solFlag,inputArray,ehatArray,xColArray,zColArray] = SolvingHAC(solver,x_next,z_next,u_k,Ppump,res_theta,val_theta,nmpcPar)

% declare variables (bounds and initial guess)
w0 = [];
lbw =[];
ubw = [];

% declare constraints and its bounds
lbg = [];
ubg = [];

% initial state
lbw = [lbw,x_next]; 
ubw = [ubw,x_next];
w0 = [w0;x_next];

% initial input
w0 = [w0;u_k];
lbw = [lbw;u_k];
ubw = [ubw;u_k];

%% Looping through until timeend
for k = 1:nmpcPar.np
    w0 = [w0; u_k];
    lbw = [lbw;nmpcPar.umin*ones(nmpcPar.nu,1)];
    ubw = [ubw;nmpcPar.umax*ones(nmpcPar.nu,1)];
        
    if k > nmpcPar.nm
        lbg = [lbg;zeros(nmpcPar.nu,1)];
        ubg = [ubg; zeros(nmpcPar.nu,1)];
    else
        lbg = [lbg;-nmpcPar.dumax*ones(nmpcPar.nu,1)];
        ubg = [ubg;nmpcPar.dumax*ones(nmpcPar.nu,1)];
    end
    
    for d = 1:3
        % creating states at collocation points
        w0 = [w0;x_next;z_next];
        lbw = [lbw;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];
        ubw = [ubw;10*ones(nmpcPar.nx,1);inf*ones(nmpcPar.nz,1)];%inf
     
    end
    
    % integrating the system
    for d = 1:3 
        % Adding xk and Xk1 as constrains as they must be equal - in
        % collocation intervals
        % algebraic constraints are set to zero in the collocation point
        lbg = [lbg;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];
        ubg = [ubg;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];  
    end
    
       
    % New NLP variable for state at end
    w0 = [w0;x_next];
    lbw = [lbw;zeros(nmpcPar.nx,1)];
    ubw = [ubw;10*ones(nmpcPar.nx,1)]; %inf
    
    % Gap
    lbg = [lbg;zeros(nmpcPar.nx,1)];
    ubg = [ubg;zeros(nmpcPar.nx,1)];
        
end


% Solving the problem
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,'p',[x_next;z_next;u_k;Ppump;res_theta;val_theta]);

%% Extracting solution
w_opt = full(sol.x);

% catch error
if solver.stats.success ~=1
    % solution failed
    solFlag = 0;
    
    % new inputs
    u_ = u_k; % fix inputs
    
    %dummy values
    ofValue = 0;
    
     inputArray = [];
     xColArray  = [];
     zColArray  = [];
     ehatArray  = [];
            
else
    % solution succeeded
    solFlag = 1;
    
    % new inputs
    u_ = w_opt(7:9);
    
    % computed values
    ofValue = full(sol.f);
    
    % variable order
    % 1-3: x0
    % 4-6: u0
    
    % total variables in one iteration = 3(u) + 9(xkd) +
    % 54(zkd) + 3(xprev) = 69
    
    % first
    % 7-9: u1
    % 10-12 | 31-33 | 52 -54: x11, x12, x13
    % 13-30 | 34-51 | 55 -72: z11, z12, z13
    % 73-75: xprev1
    
    % 76-78: u2
    % 79-81   | 100-102 | 121 - 123: x21, x22, x23
    % 82-99   | 103-120 | 124 - 141: z21, z22, z23
    % 142-144: xprev2
    
    % 145-147: u3
    % 148-150 | 169-171 | 190 - 192: x31, x32, x33
    % 151-168 | 172-189 | 193 - 210: z31, z32, z33
    % 211-213: xprev3
    
    % for k = 1:nmpcPar.np - 1
    % temp = (k - 1)*69
    
    % (temp + 7):(temp + 9): uk
    % (temp + 10):(temp + 12): xk_1
    % (temp + 31):(temp + 33): xk_2
    % (temp + 52):(temp + 54): xk_3
    % (temp + 73):(temp + 75): xk_next
    
    % (temp + 13):(temp + 30): zk_1
    % (temp + 34):(temp + 51): zk_2
    % (temp + 55):(temp + 72): zk_3
    
    % variables p/ iteration
    % = 69*nmpcPar.np + 3(u0) + 3(x0)
    
     inputArray = [];
     xColArray  = [];
     zColArray  = [];
     ehatArray  = [];
   
    for ii = 1:nmpcPar.np - 1
        temp = (ii - 1)*69;
        inputArray = [inputArray, w_opt((temp + 7):(temp + 9))];
        xColArray     = [xColArray, ...
            w_opt((temp + 10):(temp + 12)),...
            w_opt((temp + 31):(temp + 33)),...
            w_opt((temp + 52):(temp + 54))];
        zColArray    = [zColArray, ...
            w_opt((temp + 13):(temp + 30)),...
            w_opt((temp + 34):(temp + 51)),...
            w_opt((temp + 55):(temp + 72))];
        ehatArray  = [ehatArray, w_opt((temp + 73):(temp + 75))];
    end
end


end

