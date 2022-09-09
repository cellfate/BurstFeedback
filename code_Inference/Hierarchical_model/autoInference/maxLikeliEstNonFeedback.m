function result = maxLikeliEstNonFeedback(data,lb,ub)
% -- Remark
% This function uses maximum likelihood estimation to estimate the parameters
% based on the data, and calculates the theoretical mean and noise determined
% by the estimated parameters.
%
% -- Input
% data: Gene expression data.
% lb: Lower bound of parameters.
% ub: Upper bound of parameters.
%
% -- Output
% result: [theta_est,x_mean,x_noise,fval].
% theta_est: [a_est,b_est,0,0].
% x_mean: theoretical mean.
% x_noise: theoretical noise.
% fval: The value of negative log likelihood.
% -------------------------------------------------------------------------

% Fmincon solve optimization problem
A = [];
b = [];
Aeq = [];
beq = [];
options = optimoptions('fmincon','Display','off','Algorithm','interior-point','HessianApproximation','lbfgs');

exitflag = 0;
while exitflag <= 0
    % Initial point
    theta0 = initialPoint(data,lb,ub,'non-feedback');
    [theta_est,fval,exitflag,~] = fmincon(@(theta) negLogLikeNonFeedback(theta,data),theta0,A,b,Aeq,beq,lb,ub,[],options);
    [x_mean,x_noise] = computeTheoNoiseNonFeedback(theta_est,data);
end

result = [theta_est,0,0,x_mean,x_noise,fval];
end