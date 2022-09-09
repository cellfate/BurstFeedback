function [theta_est,fval,exitflag] = burstInference(vals)
theta0 = momentInference(vals);
if isempty(theta0) || any(theta0(1:2)<=0) || any(theta0(1:2)>1e3) || theta0(3)<1 || theta0(3)>1e4
    theta0 = [10,10,20];
end
[theta_est,fval,exitflag] = MLInference(vals,theta0);
end
