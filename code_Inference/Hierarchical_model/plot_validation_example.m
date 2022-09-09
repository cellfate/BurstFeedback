clear;clc;
%% fitting figure
f1 = figure;
load('results\validation_example\validation_example_a=3_b=4_k=0_h=0.mat')
subplot(3,2,1)
mRNA = 0:m_max;
[~,p_true] = generateSample(theta_true,m_max,50000);
[~,p_est] = generateSample(theta_est,m_max,50000);
histogram(data,max(data),'normalization','pdf','EdgeColor','none','FaceColor',[220,220,221]/255)
hold on
plot(mRNA,p_true,'-','LineWidth',1,'color',[230,4,48]/255)
plot(mRNA,p_est,'o','MarkerEdgeColor',[0,114,119]/255,'MarkerFaceColor',[197,229,250]/255,'MarkerIndices',[1:2:length(p_est)])
xlim([-0.95,40+0.95])
ylim([-0.006,1.2*max(p_true)])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,2)
load('results\validation_example\validation_example_a=5_b=10_k=0_h=0.mat')
mRNA = 0:m_max;
[~,p_true] = generateSample(theta_true,m_max,50000);
[~,p_est] = generateSample(theta_est,m_max,50000);
histogram(data,max(data),'normalization','pdf','EdgeColor','none','FaceColor',[220,220,221]/255)
hold on
plot(mRNA,p_true,'-','LineWidth',1,'color',[230,4,48]/255)
plot(mRNA,p_est,'o','MarkerEdgeColor',[0,114,119]/255,'MarkerFaceColor',[197,229,250]/255,'MarkerIndices',[1:5:length(p_est)])
xlim([-2.5,100+4])
ylim([-0.0017,1.2*max(p_true)])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);
hold off

subplot(3,2,3)
load('results\validation_example\validation_example_a=2_b=15_k=2_h=-3.mat')
mRNA = 0:m_max;
[~,p_true] = generateSample(theta_true,m_max,50000);
[~,p_est] = generateSample(theta_est,m_max,50000);
histogram(data,max(data),'normalization','pdf','EdgeColor','none','FaceColor',[220,220,221]/255)
hold on
plot(mRNA,p_true,'-','LineWidth',1,'color',[230,4,48]/255)
plot(mRNA,p_est,'o','MarkerEdgeColor',[0,114,119]/255,'MarkerFaceColor',[197,229,250]/255,'MarkerIndices',[1,2:5:length(p_est)])
hold off
xlim([-2,50+0.8])
ylim([-0.005,1.2*max(p_true)])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,4)
load('results\validation_example\validation_example_a=4_b=10_k=10_h=-3.mat')
mRNA = 0:m_max;
[~,p_true] = generateSample(theta_true,m_max,50000);
[~,p_est] = generateSample(theta_est,m_max,50000);
histogram(data,max(data),'normalization','pdf','EdgeColor','none','FaceColor',[220,220,221]/255)
hold on
plot(mRNA,p_true,'-','LineWidth',1,'color',[230,4,48]/255)
plot(mRNA,p_est,'o','MarkerEdgeColor',[0,114,119]/255,'MarkerFaceColor',[197,229,250]/255,'MarkerIndices',[1,2:5:length(p_est)])
hold off
xlim([-0.95,40+0.8])
ylim([-0.005,1.1*max(p_true)])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,5)
load('results\validation_example\validation_example_a=3_b=18_k=2_h=3.mat')
mRNA = 0:m_max;
[~,p_true] = generateSample(theta_true,m_max,50000);
[~,p_est] = generateSample(theta_est,m_max,50000);
histogram(data,max(data),'normalization','pdf','EdgeColor','none','FaceColor',[220,220,221]/255)
hold on
plot(mRNA,p_true,'-','LineWidth',1,'color',[230,4,48]/255)
plot(mRNA,p_est,'o','MarkerEdgeColor',[0,114,119]/255,'MarkerFaceColor',[197,229,250]/255,'MarkerIndices',[1,2:3:length(p_est)])
hold off
xlim([-0.95,30+0.8])
ylim([-0.006,1.1*max(p_true)])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,6)
load('results\validation_example\validation_example_a=2_b=10_k=10_h=3.mat')
mRNA = 0:m_max;
[~,p_true] = generateSample(theta_true,m_max,50000);
[~,p_est] = generateSample(theta_est,m_max,50000);
histogram(data,max(data),'normalization','pdf','EdgeColor','none','FaceColor',[220,220,221]/255)
hold on
plot(mRNA,p_true,'-','LineWidth',1,'color',[230,4,48]/255)
plot(mRNA,p_est,'o','MarkerEdgeColor',[0,114,119]/255,'MarkerFaceColor',[197,229,250]/255,'MarkerIndices',[1:2:length(p_est)])
hold off
xlim([-0.95,40+1.1])
ylim([-0.005,1.1*max(p_true)])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

set(f1,'position',[100 100 320/0.618 480]);

%% model selection
f2 = figure;
load('results\validation_example\validation_example_a=3_b=4_k=0_h=0.mat')
subplot(3,2,1)
ic = informationCriterion(-results_total(:,7),sum(abs(results_total(:,1:4))>0,2),length(data));
results_ic = [results_total,ic(:,3)];
min0 = min(results_ic(results_ic(:,4) == 0,8));
minpo = min(results_ic(results_ic(:,4) < 0,8));
minne = min(results_ic(results_ic(:,4) > 0,8));

X = categorical({'positive','non','negative'});
X = reordercats(X,{'positive','non','negative'});
Y = [minpo,min0,minne];
Y = Y - mean(Y);
barf = bar(X,Y,0.5);
xtips1 = barf(1).XEndPoints;
ytips1 = barf(1).YEndPoints;
labels1 = string(barf(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',6)
ylim([-1.4*max(abs(Y)),1.4*max(abs(Y))])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,2)
load('results\validation_example\validation_example_a=5_b=10_k=0_h=0.mat')
ic = informationCriterion(-results_total(:,7),sum(abs(results_total(:,1:4))>0,2),length(data));
results_ic = [results_total,ic(:,3)];
min0 = min(results_ic(results_ic(:,4) == 0,8));
minpo = min(results_ic(results_ic(:,4) < 0,8));
minne = min(results_ic(results_ic(:,4) > 0,8));

X = categorical({'positive','non','negative'});
X = reordercats(X,{'positive','non','negative'});
Y = [minpo,min0,minne];
Y = Y - mean(Y);
barf = bar(X,Y,0.5);
xtips1 = barf(1).XEndPoints;
ytips1 = barf(1).YEndPoints;
labels1 = string(barf(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',6)
ylim([-1.4*max(abs(Y)),1.4*max(abs(Y))])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,3)
load('results\validation_example\validation_example_a=2_b=15_k=2_h=-3.mat')
ic = informationCriterion(-results_total(:,7),sum(abs(results_total(:,1:4))>0,2),length(data));
results_ic = [results_total,ic(:,3)];
min0 = min(results_ic(results_ic(:,4) == 0,8));
minpo = min(results_ic(results_ic(:,4) < 0,8));
minne = min(results_ic(results_ic(:,4) > 0,8));

X = categorical({'positive','non','negative'});
X = reordercats(X,{'positive','non','negative'});
Y = [minpo,min0,minne];
Y = Y - mean(Y);
barf = bar(X,Y,0.5);
xtips1 = barf(1).XEndPoints;
ytips1 = barf(1).YEndPoints;
labels1 = string(barf(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',6)
ylim([-1.4*max(abs(Y)),1.4*max(abs(Y))])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,4)
load('results\validation_example\validation_example_a=4_b=10_k=10_h=-3.mat')
ic = informationCriterion(-results_total(:,7),sum(abs(results_total(:,1:4))>0,2),length(data));
results_ic = [results_total,ic(:,3)];
min0 = min(results_ic(results_ic(:,4) == 0,8));
minpo = min(results_ic(results_ic(:,4) < 0,8));
minne = min(results_ic(results_ic(:,4) > 0,8));

X = categorical({'positive','non','negative'});
X = reordercats(X,{'positive','non','negative'});
Y = [minpo,min0,minne];
Y = Y - mean(Y);
barf = bar(X,Y,0.5);
xtips1 = barf(1).XEndPoints;
ytips1 = barf(1).YEndPoints;
labels1 = string(barf(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',6)
ylim([-1.4*max(abs(Y)),1.4*max(abs(Y))])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,5)
load('results\validation_example\validation_example_a=3_b=18_k=2_h=3.mat')
ic = informationCriterion(-results_total(:,7),sum(abs(results_total(:,1:4))>0,2),length(data));
results_ic = [results_total,ic(:,3)];
min0 = min(results_ic(results_ic(:,4) == 0,8));
minpo = min(results_ic(results_ic(:,4) < 0,8));
minne = min(results_ic(results_ic(:,4) > 0,8));

X = categorical({'positive','non','negative'});
X = reordercats(X,{'positive','non','negative'});
Y = [minpo,min0,minne];
Y = Y - mean(Y);
barf = bar(X,Y,0.5);
xtips1 = barf(1).XEndPoints;
ytips1 = barf(1).YEndPoints;
labels1 = string(barf(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',6)
ylim([-1.4*max(abs(Y)),1.4*max(abs(Y))])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

subplot(3,2,6)
load('results\validation_example\validation_example_a=2_b=10_k=10_h=3.mat')
ic = informationCriterion(-results_total(:,7),sum(abs(results_total(:,1:4))>0,2),length(data));
results_ic = [results_total,ic(:,3)];
min0 = min(results_ic(results_ic(:,4) == 0,8));
minpo = min(results_ic(results_ic(:,4) < 0,8));
minne = min(results_ic(results_ic(:,4) > 0,8));

X = categorical({'positive','non','negative'});
X = reordercats(X,{'positive','non','negative'});
Y = [minpo,min0,minne];
Y = Y - mean(Y);
barf = bar(X,Y,0.5);
xtips1 = barf(1).XEndPoints;
ytips1 = barf(1).YEndPoints;
labels1 = string(barf(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',6)
ylim([-1.4*max(abs(Y)),1.4*max(abs(Y))])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);

set(f2,'position',[100 100 200 400]);
