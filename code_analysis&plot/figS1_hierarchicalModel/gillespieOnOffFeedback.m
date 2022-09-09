function [t,x] = gillespieOnOffFeedback(param)
% defines the reactions and the stochastic rate constants
%         OFF ON mRNA
r_mu  = [ -1  1  0;
           1 -1  0;
           0  0  1;
           0  0 -1] ;
[m,~] = size(r_mu); %number of reactions
t_c = param.t_c;
rate = [param.kon param.koff param.mu param.delta];
x0 = param.x0;

t = 0;  %time
x = x0;  %initial number of reactant molecules

% Simulation the reaction
while (t(end) < t_c)
    % step 1: Calculate a_mu & a_0
    a_mu = computeh_mu(x(end,:),rate,param);
    a_0  = sum(a_mu);
    
    %Step 2: calculate tau and mu using random number generators
    r1  = rand;
    tau = (1/a_0)*log(1/r1);
    
    r2 = rand;
    for iter = 1:m
        if (sum(a_mu(1:iter)) >= r2*a_0)
            next_mu = iter;
            break;
        end
    end
    
    %Step 3: carry out the reaction mu_next
    t =[t; t(end) + tau];
    %carry out reaction next_mu
    gain = r_mu(next_mu,:); %extract substrate stoichiometry
%     if next_mu == 1
%         gain = 1;
%    end
    last = x(end,:) + gain;
    x = [x; last];
    disp(t(end))
end %end of while (t(end) < totalTime)

% function retunring vector of M h_mu which are functions of N X
function h_mu = computeh_mu(x,rate,param)
rate_temp = rate;
rate_temp(1) = rate(1) * (param.epslion + param.K_A^param.h/(param.K_A^param.h + x(end,3)^param.h));
if x(end,1) == 1
    h_mu = rate_temp.*[1 0 0 x(end,3)];
elseif x(end,2) == 1
    h_mu = rate_temp.*[0 1 1 x(end,3)];
end

