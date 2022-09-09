function [x_mean,x_noise] = computeTheoNoiseFeedback(theta,data)
% -- Remark
% This function computes the theoretical mean and noise for a given parameter.
%
% -- Input 
% theta: [a,b,k,h], parameters of distribution with feedback.
% data: Gene expression data.
% 
% -- Output
% x_mean: theoretical mean.
% x_noise: theoretical noise.
% -------------------------------------------------------------------------

% Parameter setting
a = theta(1);
b = theta(2);
k = theta(3);
h = theta(4);
e = 0.05;
r = 0.5;

% Compute probability
m_max = max(round(r*a*b + 10*sqrt(r*a*b+r^2*a*b^2)),max(data)+100);
p = computePoissonFeedbackProb(theta,m_max);

% Compute theoretical mean and noise
x_mean = sum(p.*(0:length(p)-1));
x_2m = sum(p.*(0:length(p)-1).^2);
x_noise = (x_2m - x_mean^2)/x_mean^2;

end

