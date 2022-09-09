function LogLik = negLogLikeFeedback(theta,data)
% -- Remark
% This function generates the initial point of the optimization problem.
% a b: The initial point of a, b is estimated by the moment estimation of 
% Poisson-Gamma distribution.
% k: Log uniform distribution randomly generates the initial point of k.
% H: The integer [-5,5] is randomly generated as the initial point of H.
%
% -- Input 
% theta: [a,b,k,h], parameters of distribution with feedback.
% data: Gene expression data.
% 
% -- Output
% LogLik: The value of negative log likelihood.
% -------------------------------------------------------------------------

% Parameter setting
a = theta(1);
b = theta(2);
e = 0.05;
r = 0.5;

% Compute probability
m_max = max(round(r*a*b + 10*sqrt(r*a*b+r^2*a*b^2)),max(data)+100);
% m_max = max(data)+100;
p = computePoissonFeedbackProb(theta,m_max);

% Compute negative log likelihood  
p_vals = p(data+1);
LogLik = -sum(log(p_vals));

end