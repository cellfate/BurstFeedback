function ic = informationCriterion(logL,numParam,numObs)

aic = -2.*logL + 2.*numParam;
bic = -2*logL + log(numObs).*numParam;
aicc = aic + (2.*numParam.*(numParam + 1))./(numObs - numParam - 1);
caic = -2.*logL + (log(numObs) + 1).*numParam;
hqc = -2.*logL + 2.*log(log(numObs)).*numParam;
ic = [aic,bic,aicc,caic,hqc];

end