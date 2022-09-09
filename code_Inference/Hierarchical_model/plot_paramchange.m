%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter sensitivity analysis of hierarchical model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
addpath(genpath(pwd))

%% a change: a = 0.1:0.1:1, a = 3:3:30, h = -2, h = 2, h = 0 (non feedback)
%  e = 0.05, r = 0.5, b = 10, k = 5.
figure
subplot(2,3,1)
result = [];
for a = 0.1:0.1:1
    h = -2;
    b = 10;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\a_1.csv',result)

subplot(2,3,2)
result = [];
for a = 0.1:0.1:1
    h = 2;
    b = 10;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:20;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\a_2.csv',result)

subplot(2,3,3)
result = [];
for a = 0.1:0.1:1
    h = 0;
    b = 10;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:20;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\a_3.csv',result)

subplot(2,3,4)
result = [];
for a = 3:3:30
    h = -2;
    b = 10;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:300;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\a_4.csv',result)

subplot(2,3,5)
result = [];
for a = 3:3:30
    h = 2;
    b = 10;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:50;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\a_5.csv',result)

subplot(2,3,6)
result = [];
for a = 3:3:30
    h = 0;
    b = 10;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:300;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\a_6.csv',result)

%% b change:  b = 2:2:20, h = -2, h = 2, h = 0 (non feedback)
%  fix a = 0.5 or a = 5, e = 0.05, r = 0.5, k = 5.
figure
subplot(2,3,1)
result = [];
for b = 2:2:20
    h = -2;
    a = 0.5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\b_1.csv',result)

subplot(2,3,2)
result = [];
for b = 2:2:20
    h = 2;
    a = 0.5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\b_2.csv',result)

subplot(2,3,3)
result = [];
for b = 2:2:20
    h = 0;
    a = 0.5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\b_3.csv',result)

subplot(2,3,4)
result = [];
for b = 2:2:20
    h = -2;
    a = 5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:100;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\b_4.csv',result)

subplot(2,3,5)
result = [];
for b = 2:2:20
    h = 2;
    a = 5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:50;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\b_5.csv',result)

subplot(2,3,6)
result = [];
for b = 2:2:20
    h = 0;
    a = 5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:150;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\b_6.csv',result)

%% H change:  a = 0.5 or a = 5, h = 1:10, h = -10:-1
figure
subplot(2,2,1)
result = [];
for h = 1:10
    b = 10;
    a = 0.5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\h_1.csv',result)

subplot(2,2,2)
result = [];
for h = -10:-1
    b = 10;
    a = 0.5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\h_2.csv',result)

subplot(2,2,3)
result = [];
for h = 1:10
    b = 10;
    a = 5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:30;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\h_3.csv',result)

subplot(2,2,4)
result = [];
for h = -10:-1
    b = 10;
    a = 5;
    k = 5;
    theta = [a,b,k,h];
    xvals = 0:100;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\h_4.csv',result)

%% k change:  a = 0.5 or a = 5, h = 2, h = -2, k = [1, 5, 10, 50, 100]
figure
subplot(2,2,1)
result = [];
for k = [1, 5, 10, 50, 100]
    b = 10;
    a = 0.5;
    h = -2;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\k_1.csv',result)

subplot(2,2,2)
result = [];
for k = [1, 5, 10, 50, 100]
    b = 10;
    a = 0.5;
    h = 2;
    theta = [a,b,k,h];
    xvals = 0:10;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\k_2.csv',result)

subplot(2,2,3)
result = [];
for k = [1, 5, 10, 50, 100]
    b = 10;
    a = 5;
    h = -2;
    theta = [a,b,k,h];
    xvals = 0:80;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
    ylim([0,0.25])
end
csvwrite('misc\data_paramchange\k_3.csv',result)

subplot(2,2,4)
result = [];
for k = [1, 5, 10, 50, 100]
    b = 10;
    a = 5;
    h = 2;
    theta = [a,b,k,h];
    xvals = 0:70;
    [~,p] = generateSample(theta,max(xvals),5000);
    result = [result;p];
    hold on
    plot(xvals,p)
    hold off
end
csvwrite('misc\data_paramchange\k_4.csv',result)