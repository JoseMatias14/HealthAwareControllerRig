function [u_,sArray,ehatArray,inputArray,ofValue,solFlag] = SolvingNMPC_MS(solver,x_next,z_next,u_k,theta_k,nmpcPar)

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
        w0 = [w0; u_k];
        lbw = [lbw;nmpcPar.umin*ones(nmpcPar.nu,1)];
        ubw = [ubw;nmpcPar.umax*ones(nmpcPar.nu,1)];

        % creating current slack variables
        w0 = [w0;0*ones(nmpcPar.nx,1)];
        lbw = [lbw;0*ones(nmpcPar.nx,1)];
        ubw = [ubw;10000*ones(nmpcPar.nx,1)]; 

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

        % Constraint on erosion
        lbg = [lbg;-inf*ones(nmpcPar.nx,1)]; %zeros
        ubg = [ubg;nmpcPar.x_threshold*ones(nmpcPar.nx,1)];

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
            sArray{l} = zeros(nmpcPar.nx,nmpcPar.np);
            ehatArray{l} = zeros(nmpcPar.nx,nmpcPar.np + 1);
            inputArray{l} = zeros(nmpcPar.nx,nmpcPar.np);
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
        
        % total variables in one iteration = 3(u) + 3(s) + 9(xkd) +
        % 72(zkd) + 3(xprev) = 90
        
        % first 
        % 7-9: u1
        % 10-12: s1
        % 13-15 | 40-42 | 67 - 69: x11, x12, x13
        % 16-39 | 43-66 | 70 - 93: z11, z12, z13
        % 94-96: xprev1

        % 97-99: u2
        % 100-102: s2
        % 103-105 | 130-132 | 157 - 159: x21, x22, x23
        % 106-129 | 133-156 | 160 - 183: z21, z22, z23
        % 184-186: xprev2
       
        % 187-189: u3
        % 190-192: s3
        % 193-195 | 220-222 | 247 - 249: x31, x32, x33
        % 196-219 | 223-246 | 250 - 273: z31, z32, z33
        % 274-276: xprev3
        
        % for k = 1:nmpcPar.np - 1
        % temp = (k - 1)*90 
        % (temp + 7):(temp + 9): uk
        % (temp + 10):(temp + 12): sk
        % (temp + 94):(temp + 96): xk
        
        % variables p/ iteration 
        % = 90*nmpcPar.np + 3(u0) + 3(x0)
        
        for l = 1:Scons
            
            % initial values of the scenario
            s_idx = (90*nmpcPar.np + 3 + 3)*(l - 1);
            
            ehatArray{l} = w_opt(s_idx + 1:s_idx + 3); % !x_prev
            inputArray{l} = w_opt(s_idx + 4:s_idx + 6);
            sArray{l} = zeros(3,1); %s1

            for ii = 1:nmpcPar.np
                temp = (ii - 1)*90;
                inputArray{l} = [inputArray{l}, w_opt((s_idx + temp + 7):(s_idx + temp + 9))];
                sArray{l}     = [sArray{l}, w_opt((s_idx + temp + 10):(s_idx + temp + 12))];
                ehatArray{l}  = [ehatArray{l}, w_opt((s_idx + temp + 94):(s_idx + temp + 96))];   
            end
        end
    end

end

