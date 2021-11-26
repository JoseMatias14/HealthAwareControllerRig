function [eq,x_var,p_var] = BuildingHACModel(par)
%    Build model for the MPC - deterministic degradation using a NN to
%    compute the static value of DP

% Inputs:
%    par = system parameters
%
% Outputs:
%   eq = model equations 
%   x_var,p_var = model states and parameters/inputs

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Matias
% email: jose.o.a.matias@ntnu.no 
% November 2021; Last revision: 
        
    import casadi.*

%% Parameters
%number of wells
n_w = par.n_w; %[]

%properties
%density of oil - dim:  nwells x 1
rho_o = par.rho_o; %[kg/m3]
%1cP oil viscosity
mu_oil = par.mu_oil;% [Pa s] 

%project
%well parameters - dim:  nwells x 1
L_w = par.L_w; %[m]
A_w = par.A_w;%[m2]

%well below injection - [m]
D_bh = par.D_bh;

%riser - [m]
L_r = par.L_r;
H_r = par.H_r;
A_r = par.A_r;%[m2]

%% System states

%liquid holdup 
m_o = MX.sym('m_or',n_w);       % 1:3 [kg]
%liquid rate from reservoir
w_ro = MX.sym('w_ro',n_w);      % 4:6 [1e2 kg/s]
%riser head total production rate
w_pr = MX.sym('w_pr',n_w);      % 7:9 [1e2 kg/s]

%riser head pressure
p_rh = MX.sym('p_rh',n_w);      % 10:12 [bar]
%pressure - below injection point (bottom hole)
p_bh = MX.sym('p_bh',n_w);      % 13:15 [bar]

%delta pressure (calculated using the NN)
dP = MX.sym('dP',n_w);          % 16:18 [mbar]

%% System input
%valve oppening
vo = MX.sym('vo',n_w);            % 1:3 [0-1]
%Pump outlet pressure
Ppump = MX.sym('Ppump',1);      % 4 [bar g]
%current time
tCurr = MX.sym('tCurr',1);      % 5 [s]

%% parameters
%%%%%%%%%%% fixed %%%%%%%%%%%%%%%
%separator pressure (open to atmosphere)
p_s = MX.sym('p_s',1); %[bar]

%%%%%%%%%%%%%%% variable %%%%%%%%%%%%%%%
%system temperature
T_s = MX.sym('T_r',n_w); %[oC]

%regressor normalization 
mu_rg = MX.sym('mu_rg',15); % NN has 15 regressors (5 for each well)
std_rg = MX.sym('std_rg',15); % NN has 15 regressors (5 for each well)

%response normalization 
mu_rp = MX.sym('mu_rp',3); % NN has 3 outputs (1 for each well)
std_rp = MX.sym('std_rp',3); % NN has 3 outputs (1 for each well)

%%%%%%%%%%%%%%% estimable %%%%%%%%%%%%%%%
%reservoir data-driven model
res_theta = MX.sym('res_theta',n_w);%[m2]
%riser valve characteristics
val_theta = MX.sym('val_theta',n_w);%[m2]

%% Modeling
%conversion
% CR = 60*10^3; % [L/min] -> [m3/s] 

%oil from reservoir flowrate [kg/s] - same pump outlet pressure
% original equation
%f1 = -Ppump*ones(n_w,1)*1e5 + (w_ro.*1e-2).^2./((1e-5*vo.*res_theta).^2.*(rho_o)) + p_bh.*1e5 ; 
% with re-arranged parameters
f1 = -Ppump*ones(n_w,1)*1e5 + (w_ro.*1e-2).^2.*(res_theta.*1e9)./(vo.^2.*rho_o) + p_bh.*1e5 ; 
%riser head pressure [Pa]
%original equation
%f2 = -p_rh.*1e5 + (w_pr.*1e-2).^2./((1e-4*val_theta).^2.*(rho_r.*1e2)) + p_s.*1e5 ;
f2 = -p_rh.*1e5 + (w_pr.*1e-2).^2.*(val_theta.*1e8)./rho_o + p_s.*1e5 ;
%bottom hole pressure [Pa]
f3 = -p_bh.*1e5 + (dP*1e2 + p_rh.*1e5 + rho_o.*9.81.*H_r + 128.*mu_oil.*(L_w+L_r).*(w_ro.*1e-2)./(3.14.*D_bh.^4.*rho_o));
%liquid mass [kg]
f4 = -m_o + (rho_o).*((A_w.*L_w + A_r.*L_r));
% mass conservation
f5 = w_ro - w_pr;

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Use neural net model %
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Input 1
    x1_step1.xoffset = [-1.73110675005468;-4.46948475875232;-4.45396066189363;-2.7215822742063;-3.6952483850769];
    x1_step1.gain = [0.577665179480862;0.208395588336176;0.234668224625366;0.402885386891299;0.266822575542922];
    x1_step1.ymin = -1;

    % Layer 1
    b1 = [-3.6579783652621227219;-10.953904718353749459;0.57531937824306511597;1.7626972192892471636;0.69903768126504994829;25.346048275816599471;-6.7231855818824470816;-7.291408508150973411;-7.8671658911697530669;0.40720020299180348378;0.22090074002400741926;-0.95235090921732035163;-0.55245597575828986336;9.7930861930456742925;-2.6789648091887778847;-0.23647299623944248448;-1.3513719480995785016;0.90772138693394754938;0.23607765548607728689;-3.0651687430551128877;-0.19952869571892181688;11.38761553520214953;-1.5743053052916036183;7.844286740368879407;-11.861108049576342793];
    IW1_1 = [-5.2386844165460786371 -0.048265052628785790056 -0.098286748332031048658 -0.44056296704068942383 -1.4272682239958736083;6.8826368299586597743 -0.087797603139359414737 0.22295934625200908297 2.09087299208620081 -8.583371597512776674;-1.4064125697962837069 -0.22868191103792939733 -0.23280978573811950127 -1.1563412041735039448 -8.0561647493957870836;-7.4795172575455826092 0.076098425296620364477 0.094367120544386645653 6.2261527259452336125 1.3657091037126252164;-2.6491960241421783628 -0.29976966190534171108 -0.15916879089077609621 -1.6380822990633150393 -10.729203395871408588;-27.401359632682598999 -0.090924971801460183207 0.26389415171592317533 2.7962262460129423403 -1.3155730505492340132;-10.121087035459824577 -0.17469957185008316847 -0.10710776050485093513 0.099792077668577228167 -2.9965120693333049395;2.6999918796083215433 -0.040202914689202458309 0.30299938449891178349 7.575448174770666121 -3.354394957447159431;2.9969039816325526715 -0.042809616660182042702 0.32865709167356405018 8.1046426821212556035 -3.827698881142777676;-2.3960837123894211942 0.038736635999167642885 -0.073089546582151804222 0.1420293460258429119 -0.27840178433655388135;1.8795441452028753826 0.15734795933744896712 0.18166478390823834199 -4.5623331356572647266 2.0462292101888261975;12.238054226817290271 -0.47586662480100361261 0.26853306150843825328 -4.4504677081375971781 -2.2693977736592438887;1.819498022618531996 0.25465694933241050935 0.18973759234667586515 1.2317237544606951527 8.4398447113596315461;11.005641950053147937 -0.043338580783989347212 0.013438780101932114655 -0.48675136999504442503 -0.019095025488966802657;-3.581515879686110182 -0.02482636754947721916 -0.10057910949607175299 -0.37951005170530777155 -0.96338065880999945723;4.525471298823998012 0.7482981607681268299 0.16226423829832464785 -5.4211636331715791215 3.5569425965106633569;4.5926482260770162824 -0.31939919604044192392 0.046429729624900530072 -6.4482658126901339912 -3.3959484902642680382;-3.4228302789565510089 -0.33523790620529397843 -0.14594016532857445601 -2.0894633943015339383 -13.393378587886912712;-4.8397717698951243648 -0.77182326358874542294 -0.20267658711909869451 5.8358217761573039084 -3.5333039292525332264;-3.0321895363589317363 0.20070804797772015537 0.013348553984252860619 8.2277316198581011975 0.41409048819568289312;-1.701128492942367787 -0.15387947270973209579 -0.16605013213713756826 4.0825369117096244054 -2.1339575265091390577;17.510216985296462866 0.42573641425593639065 0.24171746004520300311 -0.94251212439309484115 6.3173946532318119296;-7.0562914809814483164 -0.06075008370310496647 -0.26685942231627890475 3.5712353144062478627 -0.13046919136274601203;8.9043874723859293141 0.0089125458429335953359 0.031168723584626772949 -0.50127503872802692264 0.15074996122465558757;-11.559915647522753801 0.25192858012718055694 0.042395930804416073756 -1.7166098619250025337 0.071443551077339836897];

    % Layer 2
    b2 = 0.33963056254321927208;
    LW2_1 = [-5.3016346069287303422 0.31640685838966930987 2.3058539048244233172 -0.20855038986123569655 6.5632467016327717957 -0.42816363990229522329 2.3592487303063069959 -2.2993537450439629488 2.1046130686702184676 0.82507831277984322593 -1.625114281953945472 0.15179407965154054216 6.2063485262685480492 -4.4986565620625853512 5.0609015464736346601 -1.9356640491594723219 -0.1927748079886754673 -2.8382307653280438586 -1.8463510081818608199 0.15382343138343787525 -1.6603149098766967207 0.70341288494909004481 -0.29147998618512771518 4.9252751507043557311 -0.96015302825757076199];

    % Output 1
    y1_step1.ymin = -1;
    y1_step1.gain = 0.692218221781655;
    y1_step1.xoffset = -1.48750400518539;

    % ===== COMPUTING NN  REGRESSORS ========
    regr = [tCurr*ones(1,3);
            w_ro'.*1e-2*60*1e3/par.rho_o(1); % changing units from [1e2 kg/s] --> [L/min]
            (p_rh' - 1.01325).*10^3; %[bar a]-->[mbar g];
            T_s';
            vo'.*(0.02 - 0.004) + 0.004]; %[0-1]-->[A];

    % normalizing regressors with data from (ExpDataAnalysis.m)
    temp1 = [mu_rg{1:5},mu_rg{6:10},mu_rg{11:15}];
    temp2 = [std_rg{1:5},std_rg{6:10},std_rg{11:15}];
    regrN = (regr - temp1)./temp2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WORKAROUND: we always want to evaluate at t_curr        %
    % and we normalize t as (t_curr*60 - 150 + 1):(t_curr*60) %
    % hence, t_curr_normm is always equal to 1.7159           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    regrN(1,:) = tCurr*ones(1,3);

    % ===== COMPUTING NN OUTPUTS ========
    xp1 = (regrN - x1_step1.xoffset).*x1_step1.gain + x1_step1.ymin;

    % Layer 1
    temp3 = b1 + IW1_1*xp1;
    a1 = 2 ./ (1 + exp(-2*temp3)) - 1;

    % Layer 2
    a2 = b2 + LW2_1*a1;

    % Output 1
    dPN = (a2 - y1_step1.ymin)/y1_step1.gain + y1_step1.xoffset;

%rescaling the outputs: calculating current DP
f6 = dP - (dPN'.*std_rp + mu_rp);

%f6 = dP - ([1;1;1].*std_rp + mu_rp);

% Form the DAE system
eq = vertcat(f1,f2,f3,f4,f5,f6);

% give parameter values
eq = substitute(eq,p_s,par.p_s);

% concatenate the differential and algebraic states
x_var = vertcat(m_o,w_ro,w_pr,p_rh,p_bh,dP);
p_var = vertcat(vo,tCurr,Ppump,T_s,res_theta,val_theta,mu_rg,std_rg,mu_rp,std_rp);
    
end

