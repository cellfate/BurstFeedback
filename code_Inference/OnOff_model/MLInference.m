function [theta_est,fval,exitflag] = MLInference(vals,theta0)
% Parameters range
lb = [1e-3;1e-3;1];
ub = [1e3;1e3;1e4];
A = [];
b = [];
Aeq = [];
beq = [];
options = optimoptions('fmincon','Display','off','Algorithm','interior-point','HessianApproximation','lbfgs');
[theta_est,fval,exitflag] = fmincon(@(theta) LogLikelihood(theta,vals),theta0,A,b,Aeq,beq,lb,ub,[],options);
end

function LogLik = LogLikelihood(theta,vals)
vals = vals(:);
kon = theta(1);
koff = theta(2);
ksyn = theta(3);
[bp,wf] = gaussJacob(50,kon-1,koff-1);
A = 1/beta(kon,koff) * 2^(1-kon-koff);
p = exp( repmat(vals,1,50).*log(repmat(ksyn*(1+bp')/2,length(vals),1)+eps)...
    - repmat(gammaln(vals+1),1,50)...
    - repmat(ksyn*(1+bp')/2,length(vals),1) );
LogLik = -sum(log(A*p*wf + eps));
end
