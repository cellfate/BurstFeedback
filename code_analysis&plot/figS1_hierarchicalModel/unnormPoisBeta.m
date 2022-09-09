function p = unnormPoisBeta(x,n,param)
kon = param.kon/param.delta;
koff = param.koff/param.delta;
ksyn = param.mu/param.delta;
k = param.K_A/ksyn;
h = param.h;
e = param.epslion;

p = exp( n*log(ksyn.*x) - gammaln(n+1) - ksyn*x +...
    (kon*(1+e)-1)*log(x) + (koff-1)*log(1-x) - kon/h*log(1+(x/k).^h));