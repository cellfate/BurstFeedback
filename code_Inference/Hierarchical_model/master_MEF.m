%% master program for MEFs scRNA-seq data
clear;clc;
addpath(genpath(pwd))

tic;
if isempty(gcp('nocreate'))
    parpool(4);
end
toc;

%% load preprocessed data
data_all = csvread('data/MEF_QC_all.csv');
e = 0.05;
r = 0.5;
save_folder  = fullfile(pwd,'results/results_MEF');

% check the uninferred gene
gene_number_wait = [];
for gene_number = 1:size(data_all,1)
    filename = sprintf('//result_gene_%d.mat',gene_number);
    if exist([save_folder,filename])==0
        gene_number_wait = [gene_number_wait,gene_number];
    end
end
fprintf('Begining\n');

%% main program
parfor infer_index = 1:length(gene_number_wait)
    % Delete 5% of the tail data
    gene = gene_number_wait(infer_index);
    data = data_all(gene,:);
    data = data(data>=0);
    [~,inter] = mink(data,floor(0.95*length(data)));
    data = data(inter);
    data_mean = mean(data);
    data_var = var(data);
    data_noise = data_var/data_mean^2;
    
    % Inference for no feedback
    results_non_total= [];
    lb_non = [1e-1, 1 ];% Lower bound
    ub_non = [30 , 20 ];% Upper bound
    if data_mean < data_var
        results_non_total = inferenceKinetic(data,lb_non,ub_non,'non-feedback');
    end
    
    % Inference for positive feedback (H<0)
    lb_positive = [1e-1, 1   ,   1, -10];% Lower bound
    ub_positive = [30 , 20 ,  1e3,  -1];% Upper bound
    results_positive_total = inferenceKinetic(data,lb_positive,ub_positive,'feedback');
    
    % Inference for negative feedback (H>0)
    lb_negative = [1e-1, 1   ,   1,   1];% Lower bound
    ub_negative = [30 , 20 ,  1e3,  10];% Upper bound
    results_negative_total = inferenceKinetic(data,lb_negative,ub_negative,'feedback');
    
    % Model selection
    results_total = [results_non_total;results_positive_total;results_negative_total];
    theta_est = modelSelect(data,results_total,lb_non,ub_non);
    
    % goodness-Fit test for final result
    ifgood = goodnessFit(data,theta_est);
    
    % save result
    filename = sprintf('//result_gene_%d',gene);
    parsave([save_folder,filename],results_non_total,results_positive_total,results_negative_total,...
        theta_est,ifgood,data,gene);
    fprintf('The simulation accepted the results of gene %d and goodness-fit %d\n',gene,ifgood);
end