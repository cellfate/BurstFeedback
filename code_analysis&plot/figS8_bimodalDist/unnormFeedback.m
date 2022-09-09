function p = unnormFeedback(x,theta,log_max)
% -- Remark
% This function computes the unnormalized probability value of the Feedback
% distribution (/exp^log_max for numerical overflow).
%
% -- Input 
% x: Protein concentration.
% theta : [a,b,k,h], parameters of Feedback distribution.
% log_max : The maximum value of the Log Feedback distribution.
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

% Compute probability
p = exp( (a*(1+e)-1)*log(x) - x./b - a/h*log(1+(x/k).^h) - log_max );

end