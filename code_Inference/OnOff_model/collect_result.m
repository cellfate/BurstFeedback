clear;clc;
genefiles = dir('results\results_MEF\*.mat');
results = zeros(length(genefiles),8);
for i = 1:length(genefiles)
    s = genefiles(i).name;
    s = str2num(s(isstrprop(s,'digit')));
    filename = sprintf('results\\results_MEF\\%s',genefiles(i).name);
    load(filename);
    [~,min_index] = min(data(:,4));
    results(s,:) = data(min_index,:);
    disp(i)
end

csvwrite('results_onoff_MEF.csv',results)
