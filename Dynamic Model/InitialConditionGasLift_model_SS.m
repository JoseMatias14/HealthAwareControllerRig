%Initial Condition Simplified
function [dx0,z0,u0,theta0] = InitialConditionGasLift_model_SS(par)


%% Inputs
%valve oppening
vo_0 = [0.799999999999999;0.400000000000000;1]; %[0-1]
%velocity pump
vpump_0 = [1.27717459016393];% [bar]

%% States
%riser head total production rate %%%%%%%%%%%%%%%%%%%%%%
w_pr_0 = [11.4411829228142;9.71651938137975;10.8136777449560];%[1e2 kg/s] 

%oil rate from reservoir %%%%%%%%%%%%%%%%%%%%
w_ro_0 = [11.4320492854616;9.70749851497498;10.8045832186035]; %[1e2 kg/s]
%oil holdup 
m_o_0 = [0.791996795274916;0.743820368926861;0.775391489111477];%[kg] 
%riser head pressure
p_rh_0 = [1.02477465749425;1.02025016819727;1.02249599492403];%[bar]
%pressure - below injection point (bottom hole)
p_bh_0 = [1.16247966051223;1.14953998299781;1.15729790026840];%[bar]

%% Parameters
res_theta_0 = [0.559736756017292;0.215964104321753;1.02335546722807]; % rearranged 1e9./(1e-5*theta0(1:3)).^2 
val_theta_0 = [0.711927821010912;0.643237182252521;0.659934091089500]; % rearranged 1e8./(1e-4*theta0(4:6)).^2 

dx0 = vertcat(w_pr_0);
z0 = vertcat(w_ro_0,m_o_0,p_rh_0,p_bh_0);
u0 = vertcat(vo_0,vpump_0);
theta0  = vertcat(res_theta_0,val_theta_0);