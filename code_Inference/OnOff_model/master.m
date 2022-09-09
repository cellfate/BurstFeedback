clear;clc;
tic;
if isempty(gcp('nocreate'))
    parpool(4);
end
toc;

data_all = csvread('data\MEF_QC_all.csv');
save_folder  = fullfile(pwd, 'results\results_MEF');

gene_number_wait = [];
for gene_number = 1:size(data_all,1)
    filename = sprintf('//result_gene_%d.mat', gene_number);
    if exist([save_folder,filename]) == 0
        gene_number_wait = [gene_number_wait, gene_number];
    end
end
fprintf('Begining\n');

for infer_index = 1:length(gene_number_wait)
    gene = gene_number_wait(infer_index);
    vals = data_all(gene,:);
    vals = vals(vals >= 0);
    [~,inter] = mink(vals,floor(0.95*length(vals)));
    vals = vals(inter);
    result_gene = [];
    while size(result_gene,1) < 10
        exitflag = 0;
        while exitflag ~= 2
            [theta_est,fval,exitflag] = burstInference(vals);
        end
        result_gene = [result_gene;[theta_est,fval]];
    end
    % compute statics
    kon_est = result_gene(:,1);
    koff_est = result_gene(:,2);
    ksyn_est = result_gene(:,3);
    bs_est = ksyn_est./koff_est;
    bf_est = kon_est;
    mean_est = ksyn_est.*kon_est./(kon_est+koff_est);
    cv2_est = (kon_est+koff_est)./(ksyn_est.*kon_est)+...
        koff_est./(kon_est.*(1+kon_est+koff_est));
    result_gene = [result_gene,bf_est,bs_est,mean_est,cv2_est];
    filename = sprintf('//result_gene_%d',gene);
    parsave([save_folder,filename],result_gene);
    fprintf('Simulation results of accepted gene %d\n',gene);
end