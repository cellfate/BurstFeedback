function [vals,p] = generateSample(theta,m_max,n_sample)
% -- Remark
% This code generates simulation data.
%
% -- Input
% theta : [a,b,k,h], parameters for generating data.
% m_max : The range [0,m_max] of probability values.
% n_sample : Sample size of simulated data.
%
% -- Output
% vals: simulation data.
% p:the probability of the true parameter.
% -------------------------------------------------------------------------

% Parameter setting
a = theta(1);
b = theta(2);
k = theta(3);
h = theta(4);
e = 0.05;
r = 0.5;
m = 0:m_max;

if h == 0
    %  Case 1: No feedback (poisson-gamma)
    X = gamrnd(a,b,1,n_sample);
    vals = poissrnd(r*X);
    p = 1./( (m+a) .* beta(m+1,a) .* (1/(b*r)+1).^(m+a) .* (b*r)^a );
    
else
    %  Case 2: feedback (poisson-Feedback)
    p = computePoissonFeedbackProb(theta,m_max);
    vals = discreteSample(p,n_sample)-1;
end

end


