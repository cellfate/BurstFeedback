function results = inferenceKinetic(data,lb,ub,varargin)
% -- Remark
% This function estimates kinetic parameter of scRNA-seq data.
%
% -- Input
% data: Gene expression data.
% lb: Lower bound of parameters.
% ub: Upper bound of parameters.
% varargin£º'non-feedback' corresponding to non-feedback inference
%           'feedback' corresponding to feedback inference.
%
% -- Output
% results: Estimated parameter.
% -------------------------------------------------------------------------

results = [];

if ~isempty(varargin)
    switch varargin{1}
        case 'non-feedback'
            % Case 1: No feedback
            while size(results,1) < 30 % accept 10 results
                result_1 = maxLikeliEstNonFeedback(data,lb,ub);
                results = [results;result_1];
            end
            
            
        case 'feedback'
            %  Case 2: feedback
            while size(results,1) < 30 % accept 10 results
                result_1 = maxLikeliEstFeedback(data,lb,ub);
                results = [results;result_1];
            end
            
    end
end
    
end
