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
    probCum = 0.95; % 0.5 | 0.95

    % creating cdf function
    fun1 = @(x) cdf_gamma(x,theta,probCum);
    dD_0 = 0.9;

    x_dD = fsolve(fun1,dD_0);

    % creating pdf function
    fun2 = @(x) pdf_gamma(x,alpha,theta,x_dD);
    p_0 = 0.99;

    x_p = fsolve(fun2,p_0);
    
    incr_array = [incr_array, x_dD];
    prob_array = [prob_array, x_p];
    
    clc

end

figure
subplot(2,1,1), plot(5:15,incr_array), title('Increment size')
subplot(2,1,2), plot(5:15,prob_array), title('Probability')


function F = cdf_gamma(x,theta,probCum)
% alpha is hardcoded here
F = probCum - (1 - exp(-theta.*x) - (theta.*x).*exp(-theta.*x));
end

function f = pdf_gamma(x,alpha,theta,dD)

f = x - (theta.^alpha).*(dD.^(alpha - 1)).*exp(-theta.*dD)/1;
end