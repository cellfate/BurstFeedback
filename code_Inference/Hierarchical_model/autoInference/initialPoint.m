function theta0 = initialPoint(data,lb,ub,varargin)
% -- Remark
% This function generates the initial point of the optimization problem.
% a b: The initial point of a, b is estimated by the moment estimation of
% Poisson-Gamma distribution.
% k: Log uniform distribution randomly generates the initial point of k.
% H: The integer [-5,5] is randomly generated as the initial point of H.
%
% -- Input
% data: Gene expression data.
% lb: Lower bound of parameters.
% ub: Upper bound of parameters.
% varargin£º'non-feedback' corresponding to non-feedback initial point
%           'feedback' corresponding to feedback initial point.

% -- Output
% theta0: initial point.
% -------------------------------------------------------------------------

r = 0.5;

if ~isempty(varargin)
    switch varargin{1}
        % Case 1: No feedback
        case 'non-feedback'
            theta0 = [mean(data)^2/(var(data)-mean(data)),(var(data)/mean(data)-1)/r];
            if ~(sum(theta0 > lb) == length(lb) && sum(theta0 < ub) == length(ub))
                theta0 = [lb(1) + (ub(1)-lb(1))*rand, lb(2) + (ub(2)-lb(2))*rand];
            end
          
        %  Case 2: feedback    
        case 'feedback'
            theta0 = [mean(data)^2/(var(data)-mean(data)),(var(data)/mean(data)-1)/r,...
                logunif(0,2),sign(lb(4))*randsample([1,2,3,4,5],1)];
            if ~(sum(theta0 > lb) == length(lb) && sum(theta0 < ub) == length(ub))
                theta0 = [lb(1) + (ub(1)-lb(1))*rand, lb(2) + (ub(2)-lb(2))*rand,...
                    logunif(0,2),sign(lb(4))*randsample([1,2,3,4,5],1)];
            end
    end
end

end

function y = logunif(a,b)

x = unifrnd(a,b);
y = 10^x;

end