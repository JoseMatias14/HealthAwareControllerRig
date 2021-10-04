clear
clc
close all

% testing gamma approximation 

alpha = 2;

% for plotting
count = 1; 

for QQ = [5, 7.5, 10, 12.5, 15]

    beta = 0.0043*QQ^3 - 0.0949*QQ^2 + 0.7305*QQ - 1.32;

    samples = [];

    for ii = 1:10000
        samples = [samples, gamrnd(alpha,beta)];
    end

    % closed form pdf
    X = 0:0.01:40;

    % shape parameter
    theta = 1/beta;
    f = (theta.^alpha).*(X.^(alpha - 1)).*exp(-theta.*X)/gamma(alpha);
    F = 1 - exp(-theta.*X) - (theta.*X).*exp(-theta.*X);

    subplot(2,6,count)
        histogram(samples,'Normalization','pdf')
        hold on
        plot(X,f,'r','LineWidth',1.5)
        
        
    subplot(2,6,6 + count)
        histogram(samples,'Normalization','cdf')
        hold on
        plot(X,F,'r','LineWidth',1.5)
        
   count = count + 1;
end