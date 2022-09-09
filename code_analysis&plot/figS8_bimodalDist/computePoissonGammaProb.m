function p = computePoissonGammaProb(theta,m_max)
% -- Remark
% This function computes the probability value of the Poisson-Gamma
% distribution.
%
% -- Input 
% theta : [a,b], parameters of Gamma distribution.
% m_max : The range [0,m_max] of probability values.
%
% -- Output
% p: Poisson-Gamma probability.
% -------------------------------------------------------------------------

% Parameter setting
a = theta(1);
b = theta(2);
r = 0.5;

% Compute probability 
p_x = 0:m_max;
p = 1./( (p_x+a) .* beta(p_x+1,a) .* (1/(b*r)+1).^(p_x+a) .* (b*r)^a );
p(p<eps) = eps;

end