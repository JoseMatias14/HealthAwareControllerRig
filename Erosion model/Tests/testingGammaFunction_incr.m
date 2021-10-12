clear all
clc


incr_array = [];
prob_array = [];

for QQ = 5:15

    % parameter 1
    theta = (0.0043*QQ.^3 - 0.0949*QQ.^2 + 0.7305.*QQ - 1.32).^-1; %empirical

    % parameter 2
    alpha = 2;

    % chosen cummulative probability
    x_dD = 10; % 0.5 | 0.95

    x_p = (theta.^alpha).*(x_dD.^(alpha - 1)).*exp(-theta.*x_dD)/1;

    incr_array = [incr_array, x_dD];
    prob_array = [prob_array, x_p];
    
    clc

end

figure
subplot(2,1,1), plot(5:15,incr_array), title('Increment size')
subplot(2,1,2), plot(5:15,prob_array), title('Probability')
