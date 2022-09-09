function p = computePoissonFeedbackProb(theta,m_max)
% -- Remark
% This function computes the probability value of the Poisson-Feedback
% distribution.
%
% -- Input 
% theta : [a,b,k,h], parameters of distribution with feedback.
% m_max : The range [0,m_max] of probability values.
%
% -- Output
% p: Unnormalized probability.
% -------------------------------------------------------------------------

% Parameter setting
a = theta(1);
b = theta(2);
k = theta(3);
h = theta(4);
e = 0.05;
r = 0.5;

% Compute normalized constant
x_range = linspace(eps,1e7,1e4);
max_1  = max(max((a*(1+e)-1)*log(x_range) - 1/b*x_range - a/h*log(1+(x_range/k).^h)),1);
A = integral(@(x) unnormFeedback(x,theta,max_1),0,Inf);

% Generalized Gauss¨CLaguerre quadrature
alpha = 0.9 * a * e;
beta = 1/b + r;
n_points = 200;
[x,w] = genLaguerreRule(n_points,alpha,0,beta,[]);

% Compute probability 
p_x = 0:m_max;
p = sum( exp(repmat((a*(1+e)+p_x-1-alpha),n_points,1) .* log(repmat(x,1,length(p_x)))+repmat(p_x*log(r)-gammaln(p_x+1),n_points,1)-...
    repmat(a/h*log(1+(x/k).^h),1,length(p_x))-max_1 + repmat(log(w),1,length(p_x))),1)/A;
p(1) = max(eps,1-sum(p(2:end)));
p = p./sum(p);
p(p<eps) = eps;

end