clear;clc;
% parameters
param.h = 4;
param.kon = 10;
param.epslion = 0.05; 
param.mu = 1000;
param.koff = 50;
param.K_A = 50;
param.delta = 1;

n_max = 500;
p_gamma = zeros(1,n_max+1);
p_beta = zeros(1,n_max+1);
for n_x = 0:n_max
    p_gamma(n_x+1) = integral(@(x) unnormPoisGamma(x,n_x,param),0,1); 
end

for n_x = 0:n_max
    p_beta(n_x+1) = integral(@(x) unnormPoisBeta(x,n_x,param),0,1); 
end

A_Gamma = integral(@(x) unnormGamma(x,param),0,1); 
A_Beta = integral(@(x) unnormBeta(x,param),0,1); 

param.x0 = [1 0 0];
param.t_c = 500;
[t,x] = gillespieOnOffFeedback(param);
t_q = linspace(param.t_c/4,param.t_c,10000);
x_q = interp1(t,x(:,3),t_q,'previous');
vals = x_q;

figure
h = histogram(vals,'BinEdges',0:max(vals),'normalization','probability','EdgeColor','none','FaceColor','#989898'); % 
hold on
plot((0:n_max)+0.5,p_beta/A_Beta,'Color','#2F62A0')
plot((0:n_max)+0.5,p_gamma/A_Gamma,'Color','#EC7232')
xlim([0,200])

save('example4.mat')