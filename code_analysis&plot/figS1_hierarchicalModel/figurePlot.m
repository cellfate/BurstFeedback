%% plot
clear;clc;

figure
load('example1.mat')
subplot(2,2,1)
h = histogram(vals,'BinEdges',0:max(vals),'normalization','probability','EdgeColor','none','FaceColor','#989898'); % 
hold on
plot((0:n_max)+0.5,p_beta/A_Beta,'Color','#2F62A0')
plot((0:n_max)+0.5,p_gamma/A_Gamma,'Color','#EC7232')
xlim([0,150])
ylim([0,max(h.Values)])
axis square

load('example2.mat')
subplot(2,2,2)
h = histogram(vals,'BinEdges',0:max(vals),'normalization','probability','EdgeColor','none','FaceColor','#989898'); % 
hold on
plot((0:n_max)+0.5,p_beta/A_Beta,'Color','#2F62A0')
plot((0:n_max)+0.5,p_gamma/A_Gamma,'Color','#EC7232')
xlim([0,150])
ylim([0,max(h.Values)])
axis square

load('example3.mat')
subplot(2,2,3)
h = histogram(vals,'BinEdges',0:max(vals),'normalization','probability','EdgeColor','none','FaceColor','#989898'); % 
hold on
plot((0:n_max)+0.5,p_beta/A_Beta,'Color','#2F62A0')
plot((0:n_max)+0.5,p_gamma/A_Gamma,'Color','#EC7232')
xlim([0,100])
ylim([0,max(h.Values)])
axis square

load('example4.mat')
subplot(2,2,4)
h = histogram(vals,'BinEdges',0:max(vals),'normalization','probability','EdgeColor','none','FaceColor','#989898'); % 
hold on
plot((0:n_max)+0.5,p_beta/A_Beta,'Color','#2F62A0')
plot((0:n_max)+0.5,p_gamma/A_Gamma,'Color','#EC7232')
xlim([0,200])
ylim([0,max(h.Values)])
axis square