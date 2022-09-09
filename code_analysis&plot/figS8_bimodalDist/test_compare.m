clear;clc;
data_all = csvread('..\..\code_Inference\Hierarchical_model\data\MEF_QC_all.csv');
result_onoff = csvread('..\..\code_Inference\OnOff_model\results_onoff_MEF.csv');
result_Hierarchical = csvread('..\..\code_Inference\Hierarchical_model\results_MEF.csv');

index_posi = find(result_Hierarchical(:,4) < 0);
num_figure = ceil(length(index_posi)/20);
numParam_onoff = 3;
numParam_Hiera_Posi = 4;
numObs = zeros(size(data_all,1),1);
for iter = 1:size(data_all,1)
    index = iter;
    % Data
    data = data_all(index,:);
    data = data(data >= 0);
    [~,inter] = mink(data,floor(0.95*length(data)));
    data = data(inter);
    numObs(index) = length(data);
end

aicc_onoff = informationCriterion(-result_onoff(index_posi,4), numParam_onoff, numObs(index_posi));
aicc_Hiera = informationCriterion(-result_Hierarchical(index_posi,7), numParam_Hiera_Posi, numObs(index_posi));

genename = importdata('..\genename.txt');
example_gene = [37,73,113,152,162,174,211,212,241,247,249,261,306,339,514,561,595,690,714,748,960,986,1066,1104];
example_genename = genename(example_gene);

f = figure;
for ind = 1:24
    index = example_gene(ind);
    name = example_genename{ind,1};
    
    data = data_all(index,:);
    data = data(data >= 0);
    [~,inter] = mink(data,floor(0.95*length(data)));
    data = data(inter);
    
    vals = 0:100;
    vals = vals(:);
    
    % Distribution of ONOFF model
    kon = result_onoff(index,1);
    koff = result_onoff(index,2);
    ksyn = result_onoff(index,3);
    [bp,wf] = gaussJacob(50,kon-1,koff-1);
    p_onoff = exp( repmat(vals,1,50).*log(repmat(ksyn*(1+bp')/2,length(vals),1)+eps)...
        - repmat(gammaln(vals+1),1,50)...
        - repmat(ksyn*(1+bp')/2,length(vals),1));
    A = 1/beta(kon,koff) * 2^(1-kon-koff);
    p_onoff = A*p_onoff*wf + eps;
    
    % Distribution of Hierarchical model
    theta = result_Hierarchical(index,1:4);
    if theta(3) == 0
        p_Hierarchical = computePoissonGammaProb(theta,max(vals));
    else
        p_Hierarchical = computePoissonFeedbackProb(theta,max(vals));
    end
    
    subplot(6,4,ind)
    histogram(data,'normalization','pdf','BinEdges',0:max(data),'EdgeColor','none','FaceColor','#CDCFD3')
    hold on
    p1 = plot(vals+0.5,p_onoff,'Color','#38827E','LineWidth',1);
    p2 = plot(vals+0.5,p_Hierarchical,'Color','#DC3A2E','LineWidth',1);
    xlim([0,max(data)+5])
    title(name)
    legend([p1,p2],{sprintf('%.1f',aicc_onoff(index_posi == index,3)),sprintf('%.1f',aicc_Hiera(index_posi == index,3))},'FontSize',4)
end
set(gcf,'position',[200 200 800 1600])
