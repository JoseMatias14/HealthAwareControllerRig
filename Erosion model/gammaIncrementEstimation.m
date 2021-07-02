clear 
close all
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

%plotting behavior
Q = 12.5;
tEnd = -3*Q + 51.333;
data = [0.3, 0;
        0.2829, tEnd];

tMax = ceil(tEnd);
    
%noise seed
rng('default')    
  
trajectories = [];
for mc = 1:500  % samples
    
    t0 = 0.3;
    traj = t0;
    for kk = 1:tMax %time
        t0 = t0 - 0.0005*gamrnd(1.3,2);
        traj = [traj, t0];
    end
    trajectories = [trajectories; traj];

end

meanTrajectory = mean(trajectories,1);
stdTrajectory = std(trajectories,[],1);

figure(2)
hold on
for ii = 1:500
    plot(0:tMax,trajectories(ii,:),'Color',[0,0,0,0.05])
end

plot(0:tMax,meanTrajectory,'k','Linewidth',1.5)
plot(0:tMax,meanTrajectory + 2*stdTrajectory,'k:','Linewidth',1.5)
plot(0:tMax,meanTrajectory - 2*stdTrajectory,'k:','Linewidth',1.5)

plot(data(:,2),data(:,1),'r','Linewidth',1.5)
xlabel('time [min]')
ylabel('d [cm]')

% found manually
% Q = 5, alpha = 0.5
% Q = 7.5, alpha = 0.6
% Q = 10, alpha = 0.8
% Q = 12.5, alpha = 1.3
% Q = 15, alpha = 2.7
% function = y = 0.0043x3 - 0.0949x2 + 0.7305x - 1.32 --> Excel


