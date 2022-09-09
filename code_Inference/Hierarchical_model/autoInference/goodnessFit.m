function  ifgood = goodnessFit(data,result)
% -- Remark
% This function checks that the inferred result fits the data distribution.
%
% -- Input
% data: Gene expression data.
% result : Inferred result.
%
% -- Output
% ifgood: 1 means good fit; 0 means bad fit.
% -------------------------------------------------------------------------

% Parameter setting
if sum(result) == 0
    ifgood = 0;
else
    theta_est = result(1:4);
    a = theta_est(1);
    b = theta_est(2);
    k = theta_est(3);
    h = theta_est(4);
    e = 0.05;
    r = 0.5;
    m_max = max(round(r*a*b + 10*sqrt(r*a*b+r^2*a*b^2)),max(data));
    
    tic
    [~,p_est] = generateSample(theta_est,m_max,length(data));
    escapetime = toc;
    fprintf('%f',escapetime)
    
    expect = p_est * length(data);
    actual = zeros(1,length(p_est));
    h = tabulate(data);
    actual(h(:,1)+1) = h(:,2);
    chi_score_val = sum((expect-actual).^2./expect);
    chi_score = [];
    
    % Simulation
    if escapetime < 30
        for i = 1:1000
            simuls = [];
            [simuls,~] = generateSample(theta_est,m_max,length(data));
            actual = zeros(1,length(p_est));
            h = tabulate(simuls);
            h = h(h(:,1)+1 < length(actual),:);
            actual(h(:,1)+1) = h(:,2);
            chi_score = [chi_score,sum((expect-actual).^2./expect)];
            warning off
        end
        ifgood = sum(chi_score_val<chi_score) > 50;
    else
        ifgood = 0;
    end
end

end

