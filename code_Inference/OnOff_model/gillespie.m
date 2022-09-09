function [t,x] = gillespie(reactionmatrix,rate,tottime,x0)

[nrow,ncol] = size(reactionmatrix);
m = nrow;     % number of reactions  
n = fix(ncol/2);     % number of reactants;fix:prevent something is
                     %created out of thin air

rmu = reactionmatrix(:,n+1:end) - reactionmatrix(:,1:n);
p = cell(1,m);
for iter = 1:m
%order of the first coefficient that greater than zero in the left matrix
    p{iter} = find(reactionmatrix(iter,1:n)>0);                                              
end

x    = x0; % initial number of reactant molecules
cnt  = 0;  % reaction counter
t    = 0;  % time 

% Simulation the reaction
while (t(end)<tottime)
    
    % step 1: Calculate a_mu & a_0
    amu = hmu(x(end,:),m,rate,reactionmatrix,p);
    %amu = hmu(x,m,rate,reactionmatrix,p);
    a0  = sum(amu);
        
    %Step 2: calculate tau and mu using random number generators
    r1  = rand; 
    tau = (1/a0)*log(1/r1);
        
    r2 = rand;
    for iter = 1:m
        if (sum(amu(1:iter)) >= r2*a0)
            next_mu = iter;
            break;
        end
    end
        
    %Step 3: carry out the reaction mu_next
    t =[t; t(end) + tau];         %add a new time in the vector
    cnt = cnt + 1;
        
    %carry out reaction next_mu
    gain = rmu(next_mu,:);        %extract substrate stoichiometry
    last = x(end,:) + gain;       %add a new row
    x = [x; last];
    disp(t(end))
    
end %end of while (t(end) < totalTime)

% function retunring vector of M h_mu which are functions of N X
function hmu = hmu(x,m,rate,reactionmatrix,p)
hmu = zeros(1,m);
for iter = 1:m
    if x(p{iter})
        %construction of propensity function matrix
        hmu(iter) = rate(iter)*prod(matnchoosek(x(p{iter}),reactionmatrix(iter,p{iter})));  
    elseif isempty(p{iter})
        hmu(iter) = rate(iter);
    else
        hmu(iter) = 0;  
    end
end

function c = matnchoosek(n,k)

g = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);

c = floor(exp(g)+.5);%ceil

