clear;clc;
addpath(genpath(pwd))

tic;
if isempty(gcp('nocreate'))
    parpool(4);
end
toc;

%% generate the simulation data
% Parameter setting
e = 0.05;
r = 0.5;
% theta_true_set = [3 4 0 0;
%     5 10 0 0;
%     2 15 2 -3;
%     4 10 10 -3;
%     3 18 2 3;
%     5 10 10 3];
theta_true_set = [2 10 10 3];
sensitive_set = [0.1,0.3,0.5,0.7,0.9];

param_set = [];
for theta_index = 1:size(theta_true_set,1)
    for sensitive = sensitive_set
        mat_temp = [sensitive,theta_true_set(theta_index,:)];
        param_set = [param_set;[[1:50]',repmat(mat_temp,50,1)]];
    end
end

data_set = [];
for theta_index = 1:size(theta_true_set,1)
    theta_true = theta_true_set(theta_index,:);
    a = theta_true(1);
    b = theta_true(2);
    k = theta_true(3);
    h = theta_true(4);
    m_max = round(r*a*b + 10*sqrt(r*a*b+r^2*a*b^2));
    [data,~] = generateSample(theta_true,m_max,2000);
    data_set = [data_set;data];
end

save_folder  = fullfile(pwd,'results/validation_Sensitive_negative/');

parfor iter = 1:size(param_set,1)
    times = param_set(iter,1);
    theta_true = param_set(iter,3:6);
    a = theta_true(1);
    b = theta_true(2);
    k = theta_true(3);
    h = theta_true(4);
    sensitive = param_set(iter,2);
    data_iter = find(sum(theta_true_set == theta_true,2) == 4);
    data_all = data_set(data_iter,:);
    data = data_all(rand(1,2000) < sensitive);
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
    
    % collect
    result_validation = [times,sensitive,theta_true,theta_est];
    
    filename = sprintf('validation_sensitive=%.1f_a=%d_b=%d_k=%d_h=%d_%d.mat',sensitive,a,b,k,h,times);
    parsaveValidation([save_folder,filename],results_total,result_validation);
    fprintf('validation of %d iters complete!\n',iter);
end

clear;clc;
genefiles = dir('results\validation_Sensitive\*.mat');
results = [];
for i = 1:length(genefiles)
    s = genefiles(i).name;
    filename = sprintf('results\\validation_Sensitive\\%s',s);
    load(filename);
    results = [results;result_validation];
    disp(i)
end

csvwrite('results\validation_Sensitive_results.csv',results)