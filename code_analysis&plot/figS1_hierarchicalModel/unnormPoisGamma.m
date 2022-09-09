function p = unnormPoisGamma(x,n,param)
a = param.kon/param.delta;
b = 1/(param.koff/param.delta);
ksyn = param.mu/param.delta;
k = param.K_A/ksyn;
h = param.h;
e = param.epslion;

p = exp( n*log(ksyn.*x) - gammaln(n+1) - ksyn*x +...
    (a*(1+e)-1)*log(x) - x./b - a/h*log(1+(x/k).^h));