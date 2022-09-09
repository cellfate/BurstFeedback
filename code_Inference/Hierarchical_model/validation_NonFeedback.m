%% validation：non-feedback model
clear;clc;
addpath(genpath(pwd))

tic;
if isempty(gcp('nocreate'))
    parpool(4);
end
toc;

%% generate the simulation data
% Parameter setting
a_set = [0.2, 0.5, 1:2:29];
b_set = [2:2:18];
k_set = 0;
h_set = 0;
theta_true_grid = zeros(length(a_set)*length(b_set)*length(k_set),4);

iter_grid = 1;
for h = h_set
    for k = k_set
        for a = a_set
            for b = b_set
                theta_true_grid(iter_grid,:) = [a,b,k,h];
                iter_grid = iter_grid+1;
            end
        end
    end
end

save_folder  = fullfile(pwd,'results/validationNonFeedback1');

parfor iter = 1:size(theta_true_grid,1)
    e = 0.05;
    r = 0.5;
    theta_true = theta_true_grid(iter,:);
    
    a = theta_true(1);
    b = theta_true(2);
    k = theta_true(3);
    h = theta_true(4);
    
    % Generating simulated data
    m_max = round(r*a*b + 10*sqrt(r*a*b+r^2*a*b^2));
    [data,~] = generateSample(theta_true,m_max,50000);
    data_mean = mean(data);
    data_var = var(data);
    data_noise = data_var/data_mean^2;
    
    %% Inference
    % Inference for no feedback
    results_non_total= [];
    lb_non = [1e-1, 1];% Lower bound
    ub_non = [30 , 20 ];% Upper bound
    if data_mean < data_var
        results_non_total = inferenceKinetic(data,lb_non,ub_non,'non-feedback');
    end
    
    
    % Inference for positive feedback (H<0)
    lb_positive = [1e-1, 1,    1, -10];% Lower bound
    ub_positive = [30 , 20 ,  1e3,  -1];% Upper bound
    results_positive_total = inferenceKinetic(data,lb_positive,ub_positive,'feedback');
    
    
    % Inference for negative feedback (H>0)
    lb_negative = [1e-1, 1,    1,   1];% Lower bound
    ub_negative = [30 , 20 ,  1e3,  10];% Upper bound
    results_negative_total = inferenceKinetic(data,lb_negative,ub_negative,'feedback');
    
    % Model selection
    results_total = [results_non_total;results_positive_total;results_negative_total];
    theta_est = modelSelect(data,results_total,lb_non,ub_non);
    
    % error
    error1 = -1;
    error2 = -1;
    if sum(theta_est) ~= 0
        error1 = sum(abs(theta_est(1:2) - theta_true(1:2))./theta_true(1:2));
        error2 = sum(log(theta_est(1:2)./theta_true(1:2)).^2);
    end
    
    % collect
    result_validation = [theta_true,theta_est,error1,error2];
    
    filename = sprintf('//validation_%d',iter);
    parsaveValidation([save_folder,filename],result_validation);
    fprintf('validation of %d iters complete!\n',iter);
end

%% collect result
clear;clc;
genefiles = dir('results\validation_NonFeedback\*.mat');
results = zeros(length(genefiles),13);
for i = 1:length(genefiles)
    s = genefiles(i).name;
    s = str2num(s(isstrprop(s,'digit')));
    filename = sprintf('results\\validation_NonFeedback\\%s',genefiles(i).name);
    load(filename);
    results(s,:) = result_validation;
    disp(i)
end
csvwrite('results\validation_Nonfeedback_results.csv',results)

results = csvread('results\validation_Nonfeedback_results.csv');
f = figure;
validationPlot(results,[-15.5,-10]);
set(f,'position',[100 100 300 250]);
colorbar