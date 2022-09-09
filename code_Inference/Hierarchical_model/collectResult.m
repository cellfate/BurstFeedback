%% collect result
clear;clc;
genefiles = dir('results\results_MEF\*.mat');
results = zeros(length(genefiles),8);
for i = 1:length(genefiles)
    s = genefiles(i).name;
    s = str2num(s(isstrprop(s,'digit')));
    filename = sprintf('results\\results_MEF\\%s',genefiles(i).name);
    load(filename);
    results(s,:) = [theta_est,ifgood];
    disp(i)
end
csvwrite('results_MEF.csv',results);