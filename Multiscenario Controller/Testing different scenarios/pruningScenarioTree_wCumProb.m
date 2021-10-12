clear all
close all
clc

probH = 0.9; % 0.75 | 0.95
probL = 0.50;

% testing the probabilities of different scenarios
%            w1 |  w2  | w3 
scendist = [probH, probH, probH; %S1:   HHH
            probH, probH, probL; %S2:   HHL
            probH, probL, probH; %S3:   HLH
            probL, probH, probH; %S4:   LHH
            probH, probL, probL; %S5:   HLL
            probL, probH, probL; %S6:   LHL
            probL, probL, probH; %S7:   LLH
            probL, probL, probL];%S8:   LLL

Nr = 2; 
Nlevels = size(scendist,1); 
Nscenarios = Nlevels^Nr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating scenario combinations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transforming the scenarios (decimal) into base Number of levels
temp = dec2base(0:Nscenarios - 1,Nlevels);
%  generate all tuples
%  temp is 'char'
% we need to substitute the zeros (Matlab has 1-based index) by 1
% and by adding the double 1, we obtain the indexes that we want
combs = temp-'0'+1;  %//


% number of MonteCarlo iterations
nMC = 100;

% creating arrays for saving data                   
Q_array = [];
prob_total_array = [];

for ss = 1:Nscenarios
    incr_array{ss} = [];
    prob_array{ss} = [];
end

%%
for ii = 1:nMC
    
    ii
    
    % sampling the input distribution (continuous uniform distribution)
    %   - pick one random flowrate value for each ramification and for each
    %   well
    QQ = unifrnd(5,15,[3,Nr]);
    Q_array = [Q_array, QQ];

    % calculating gamma distribution based on the sampled flowrate
    % parameter 1 (equation derived "empirically")
    theta = (0.0043*QQ.^3 - 0.0949*QQ.^2 + 0.7305.*QQ - 1.32).^-1; 
    % parameter 2 (fixed)
    alpha = 2;
    
    for ss = 1:Nscenarios % looping ramifications
        scen_prob = ones(3,1);  % for calculating probability of a given scenario
        d_inc = zeros(3,1);     % for calculating total increment of a given scenario
        
        for kk = 1:Nr % looping ramifications
            % chosen cummulative probability
            probCum = scendist(combs(ss,kk),:); % probL | probH
                
            for well = 1:3 % looping wells

                % creating cdf function (probCum is the cummulative
                % probability that we want. Function finds the value of the
                % equivalent diameter that corresponds to that probability
                fun1 = @(x) cdf_gamma(x,theta(well,kk),probCum(well));
                dD_0 = 0.9;

                % evaluating equivalent diameter values
                x_dD = fsolve(fun1,dD_0);
                d_inc(well) = d_inc(well) + x_dD;
                
                % creating pdf function
                fun2 = @(x) pdf_gamma(x,alpha,theta(well,kk),x_dD);
                p_0 = 0.99;
                
                %evaluating probability
                x_p = fsolve(fun2,p_0);
                scen_prob(well) = scen_prob(well)*x_p;
                clc
                
            end
        end
        
        % saving
        incr_array{ss} = [incr_array{ss}, d_inc];
        prob_array{ss} = [prob_array{ss}, scen_prob];
        
    end
    
    % normalizing probability
    probTemp = zeros(Nscenarios,1);
    for ss = 1:Nscenarios
        probTemp(ss) = prod(prob_array{ss}(:,end)); % normalizing for the current MC iteration
    end
    probTotal = sum(probTemp);
    
  
    prob_total_array = [prob_total_array, probTemp/probTotal];
end

%% preparing for plotting
scen_avg_prob = mean(prob_total_array,2);
scen_avg_std = std(prob_total_array,[],2);

% sorting increments by well
inc_avg = zeros(Nscenarios,3);
inc_std = zeros(Nscenarios,3);

for ss = 1:Nscenarios
    for well = 1:3
        inc_avg(ss,well) = mean(incr_array{ss}(well,:),2);
        inc_std(ss,well) = std(incr_array{ss}(well,:),[],2);
    end
end

save('scenario_pruning_cumprob','prob_total_array','scen_avg_prob','scen_avg_std','inc_avg','inc_std','incr_array','prob_array')

%%

%load('scenario_pruning_cumprob')
h = figure(1);
    plot(1:64,100*scen_avg_prob,'ko','markersize',10)
    hold on 
    plot(1:64,100*scen_avg_prob + 2*scen_avg_std,'rx')
    plot(1:64,100*scen_avg_prob - 2*scen_avg_std,'rx')

    grid on

    xlim([0.5, 64.5])
    xticks(1:64)

    xlabel('Scenarios')
    ylabel('Probabilities [%]')
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'prob_comparison','-dpdf','-r0')


for well = 1:3
   figure(1 + well)
        plot(1:64,inc_avg(:,well),'ko','markersize',10)
        hold on 
        plot(1:64,inc_avg(:,well) + 2*inc_std(:,well),'rx')
        plot(1:64,inc_avg(:,well) - 2*inc_std(:,well),'rx')
        
        %expected value of the incre
        plot(1:64,inc_avg(:,well).*scen_avg_prob,'bd')

        grid on

        xlim([0.5, 64.5])
        xticks(1:64)

        xlabel('Scenarios')
        ylabel('increment')

        title(['Well :',num2str(well)])
   
end
    
  
%%
function F = cdf_gamma(x,theta,probCum)
% alpha is hardcoded here
F = probCum - (1 - exp(-theta.*x) - (theta.*x).*exp(-theta.*x));
end

function f = pdf_gamma(x,alpha,theta,dD)

f = x - (theta.^alpha).*(dD.^(alpha - 1)).*exp(-theta.*dD)/1;
end