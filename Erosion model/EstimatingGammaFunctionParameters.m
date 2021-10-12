clear 
%close all
clc

% x = 0:0.1:50;
% y1 = gampdf(x,1,2);
% y2 = gampdf(x,3,2);
% y3 = gampdf(x,6,2);
% 
% %Plot the pdfs.
% figure(1);
% plot(x,y1)
% hold on
% plot(x,y2)
% plot(x,y3)
% hold off
% xlabel('Observation')
% ylabel('Probability Density')
% legend('a = 1, b = 10','a = 3, b = 5','a = 6, b = 4')

%noise seed
rng('default')    
  
trajectories = [];

figure(1)
count = 1;

for QQ = [5, 7.5, 10, 12.5, 15]

    theta = 0.0043*QQ^3 - 0.0949*QQ^2 + 0.7305*QQ - 1.32;

    tEnd = -3*QQ + 51.333;
    
    data = [0.3, 0;
        0.3181, ceil(tEnd)];
    
    tMax = ceil(tEnd);

    %%%%%%%%%%%%%%%
    % Monte Carlo %
    %%%%%%%%%%%%%%%
    trajectories = [];
    for mc = 1:500  % samples

        d0 = 0.3;
        traj = d0;
        dk = d0;
        for kk = 1:tMax %time
            
            dk = dk + 0.0005*gamrnd(2,theta);
            
            traj = [traj, dk];
        end
        trajectories = [trajectories; traj];

    end

    meanTrajectory = mean(trajectories,1);
    stdTrajectory = std(trajectories,[],1);

    subplot(3,2,count)
        hold on
        for ii = 1:500
            plot(0:tMax,trajectories(ii,:),'Color',[0,0,0,0.05])
        end

        plot(0:tMax,meanTrajectory,'k','Linewidth',1.5)
        plot(0:tMax,meanTrajectory + 2*stdTrajectory,'k:','Linewidth',1.5)
        plot(0:tMax,meanTrajectory - 2*stdTrajectory,'k:','Linewidth',1.5)

        plot(data(:,2),data(:,1),'r','Linewidth',1.5)
        xlim([0, 40])

        xlabel('time [min]')
        ylabel('d [cm]')

        title(['Q: ',num2str(QQ)])

        count = count + 1;
end

%% plotting gamma distributions
count = 2;

lab = {'low Q', 'nanan', 'lalala', 'medium Q','high Q'};

for QQ = [5, 7.5, 10, 12.5, 15]

    theta = 0.0043*QQ^3 - 0.0949*QQ^2 + 0.7305*QQ - 1.32;

    tEnd = -3*QQ + 51.333;
    data = [0.3, 0;
        0.3181, ceil(tEnd)];
    
    tMax = ceil(tEnd);
    
    alpha = 2;
    
    % closed form pdf
    X = 0:0.01:20;

    % shape parameter
    theta1 = 1/theta;
    f = (theta1.^alpha).*(X.^(alpha - 1)).*exp(-theta1.*X)/gamma(alpha);
   
    figure(count)
    plot(X,f,'k','LineWidth',1.5)
    grid on 
    
    %set(gca,'xtick',[])
    %set(gca,'xticklabel',[])
    %set(gca,'ytick',[])
    %set(gca,'yticklabel',[])
    
    %xlabel('\Delta X_{k + \tau |k}','Fontsize',16)
    xlabel('d_{inc}','Fontsize',16)
    ylabel('pdf','Fontsize',16)
    
    title(['Q: ',num2str(QQ)],'Fontsize',16)
    
%     title(lab{count - 1},'Fontsize',16)
     count = count + 1;

end

% found manually
% Q = 5, alpha = 0.5
% Q = 7.5, alpha = 0.6
% Q = 10, alpha = 0.8
% Q = 12.5, alpha = 1.3
% Q = 15, alpha = 2.7
% function = y = 0.0043x3 - 0.0949x2 + 0.7305x - 1.32 --> Excel


