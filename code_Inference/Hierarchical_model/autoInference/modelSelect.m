function result = modelSelect(data,results_total,lb_non,ub_non)
% -- Remark
% This function selects models for different models.
%
% -- Input
% data: Gene expression data.
% results_non_total: All results without feedback.
% results_positive_total:The result of all positive feedback.
% results_negative_total:The result of all negative feedback.
%
% -- Output
% result: Final inference result.
% -------------------------------------------------------------------------
index_reliable = (results_total(:,1) > lb_non(1)+1e-1 & results_total(:,1) < ub_non(1)-1e-1) &...
    (results_total(:,2) > lb_non(2)+1e-1 & results_total(:,2) < ub_non(2)-1e-1);

if sum(index_reliable) == 0
    result = zeros(1,7);
else
    results_reliable = results_total(index_reliable,:);
    ic = informationCriterion(-results_reliable(:,7),sum(abs(results_reliable(:,1:4))>0,2),length(data));
    [~,min_index] = min(ic(:,3));
    result = results_reliable(min_index,:);
end

end