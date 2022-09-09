clear;clc;
addpath(genpath(pwd))
% -- Remark
% This code generates simulation data and tests the inferences.
% a : burst frequency
% b : burst size
% k : equilibrium constant
% H : feedback coefficient.
% -------------------------------------------------------------------------

%% generate the simulation data
% Parameter setting
a = 2;
b = 10;
k = 10;
h = 3;
e = 0.05;
r = 0.5;
theta_true = [a,b,k,h];

% Generating simulated data
m_max = round(r*a*b + 10*sqrt(r*a*b+r^2*a*b^2));
[data,p_true] = generateSample(theta_true,m_max,50000);
data_mean = mean(data);
data_var = var(data);
data_noise = data_var/data_mean^2;

% Figure the data and p_true
figure
histogram(data,max(data),'normalization','pdf')
hold on
plot(0:length(p_true)-1,p_true)
hold off

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

%% check distribution

mRNA = 0:m_max;
[~,p_true] = generateSample(theta_true,m_max,50000);
[~,p_est] = generateSample(theta_est,m_max,50000);
f = figure;
histogram(data,max(data),'normalization','pdf','EdgeColor','none','FaceColor',[220,220,221]/255)
hold on
plot(mRNA,p_true,'-','LineWidth',1,'color',[230,4,48]/255)
plot(mRNA,p_est,'o','MarkerEdgeColor',[0,114,119]/255,'MarkerFaceColor',[197,229,250]/255,'MarkerIndices',1:1:length(p_est))
hold off
ylim([0,1.2*max(p_true)])
set(f,'position',[100 100 130/0.618 130]);
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

%% figure model selection
ic = informationCriterion(-results_total(:,7),sum(abs(results_total(:,1:4))>0,2),length(data));
results_ic = [results_total,ic(:,3)];
min0 = min(results_ic(results_ic(:,4) == 0,8));
minpo = min(results_ic(results_ic(:,4) < 0,8));
minne = min(results_ic(results_ic(:,4) > 0,8));

X = categorical({'positive','non','negative'});
X = reordercats(X,{'positive','non','negative'});
Y = [minpo,min0,minne];
Y = Y - mean(Y);
f = figure;
barf = bar(X,Y,0.5);
xtips1 = barf(1).XEndPoints;
ytips1 = barf(1).YEndPoints;
labels1 = string(barf(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',6)
ylim([-1.4*max(abs(Y)),1.4*max(abs(Y))])
set(f,'position',[100 100 150 150]);
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

save(sprintf('results\\validation_example\\validation_example_a=%d_b=%d_k=%d_h=%d.mat',a,b,k,h))