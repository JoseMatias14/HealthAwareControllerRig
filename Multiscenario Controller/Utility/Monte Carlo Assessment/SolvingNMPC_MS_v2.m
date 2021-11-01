function [u_,ehatArray,inputArray,xColArray,zColArray,ofValue,solFlag] = SolvingNMPC_MS_v2(solver,x_next,z_next,u_k,theta_k,u_k_1Array,x_k_1Array,xCol_k_1Array,zCol_k_1Array,nmpcPar,par)

% declare variables (bounds and initial guess)
w0 = [];
lbw =[];
ubw = [];

% declare constraints and its bounds
lbg = [];
ubg = [];

% total number of scenarios
S = nmpcPar.Nlevels^nmpcPar.Nr;

% Computing the number of considered scenarios
Scons = length(nmpcPar.ConsideredScen);

%Variable bounds
[lbx,lbz,~,ubx,ubz,~,~,~] = OptimizationBounds_MPC(par);

for l = 1:Scons
    % initial state
    lbw = [lbw;x_next]; 
    ubw = [ubw;x_next];
    w0 = [w0;x_next];

    % initial input
    w0 = [w0;u_k];
    lbw = [lbw;u_k];
    ubw = [ubw;u_k];

    %% Looping through until timeend
    for k = 1:nmpcPar.np
        w0 = [w0; u_k_1Array{l}(:,k)];
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
            w0 = [w0;xCol_k_1Array{l}(:,(k - 1)*3 + d);zCol_k_1Array{l}(:,(k - 1)*3 + d)];
            lbw = [lbw;lbx;lbz];
            ubw = [ubw;ubx;ubz];

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
        w0 = [w0;x_k_1Array{l}(:,k + 1)];
        lbw = [lbw;lbx];
        ubw = [ubw;ubx]; %inf

        % Gap
        lbg = [lbg;zeros(nmpcPar.nx,1)];
        ubg = [ubg;zeros(nmpcPar.nx,1)];

    end
     
end

% building no-nantecipative matrices
a = 1;
ac = 1;
for i = 1:nmpcPar.np
    for j = 1:S
        nonant(i,j) = a;
        if i == 1
            a = 1;
        else if i <= nmpcPar.Nr
                ac = ac+1;
                if ac > nmpcPar.Nlevels
                    a = a+1;
                    ac = 1;
                end
            else
                a = NaN;
            end
        end
    end
end
nonant = nonant';
nonant = nonant(nmpcPar.ConsideredScen,:);

% Add Non-anticipativity constraints
for k = 1:nmpcPar.Nr
    for js = 1:Scons-1
        if ~isnan(nonant(js,k))
            if nonant(js,k) == nonant(js+1,k)
                 lbg = [lbg; zeros(nmpcPar.nu,1)];
                 ubg = [ubg; zeros(nmpcPar.nu,1)];
            end
        end
    end
end

    % Solving the problem
    sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,'p',[x_next;z_next;u_k;theta_k(1:7)]);

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
        for l = 1:Scons
            ehatArray{l} = NaN;
            inputArray{l} = NaN;
            xColArray{l} = NaN;
            zColArray{l} = NaN;
            
        end

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
        % 72(zkd) + 3(xprev) = 87
        
        % first 
        % 7-9: u1
        % 10-12 | 37-39 | 64 - 66: x11, x12, x13
        % 13-36 | 40-63 | 67 - 90: z11, z12, z13
        % 91-93: xprev1

        % 94-96: u2
        % 97-99   | 124-126 | 151 - 153: x21, x22, x23
        % 100-123 | 127-150 | 154 - 177: z21, z22, z23
        % 178-180: xprev2
       
        % 181-183: u3
        % 184-186 | 211-213 | 238 - 240: x31, x32, x33
        % 187-210 | 214-237 | 241 - 264: z31, z32, z33
        % 265-267: xprev3
        
        % for k = 1:nmpcPar.np - 1
        % temp = (k - 1)*87 
        
        % (temp + 7):(temp + 9): uk
        % (temp + 10):(temp + 12): xk_1
        % (temp + 37):(temp + 39): xk_2
        % (temp + 64):(temp + 66): xk_3
        % (temp + 91):(temp + 93): xk_next
        
        % (temp + 13):(temp + 26): zk_1
        % (temp + 40):(temp + 63): zk_2
        % (temp + 67):(temp + 90): zk_3
        
        % variables p/ iteration 
        % = 90*nmpcPar.np + 3(u0) + 3(x0)
        
        for l = 1:Scons
            
            % initial values of the scenario
            s_idx = (87*nmpcPar.np + 3 + 3)*(l - 1);
            
            ehatArray{l} = w_opt(s_idx + 1:s_idx + 3); % !x_prev
            inputArray{l} = w_opt(s_idx + 4:s_idx + 6);
            xColArray{l} = [];
            zColArray{l} = [];
            
            for ii = 1:nmpcPar.np
                temp = (ii - 1)*87;
                inputArray{l} = [inputArray{l}, w_opt((s_idx + temp + 7):(s_idx + temp + 9))];
                xColArray{l}     = [xColArray{l}, ...
                                    w_opt((s_idx + temp + 10):(s_idx + temp + 12)),...
                                    w_opt((s_idx + temp + 37):(s_idx + temp + 39)),...
                                    w_opt((s_idx + temp + 64):(s_idx + temp + 66))];
                 zColArray{l}     = [zColArray{l}, ...
                                    w_opt((s_idx + temp + 13):(s_idx + temp + 36)),...
                                    w_opt((s_idx + temp + 40):(s_idx + temp + 63)),...
                                    w_opt((s_idx + temp + 67):(s_idx + temp + 90))];   
                ehatArray{l}  = [ehatArray{l}, w_opt((s_idx + temp + 91):(s_idx + temp + 93))];   
            end
        end
    end

end

