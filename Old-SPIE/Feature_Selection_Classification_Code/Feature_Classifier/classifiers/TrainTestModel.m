function [stats, methodstring, Mdl] = TrainTestModel(method, training_set, testing_set, training_labels, testing_labels, threshold, options)
% Wrapper function to call any one of CCIPD's implemented classifier
% methods in a training and validation fashion

% THRESHOLD: ideally passed in, otherwise will compute from training dataset

if ~ismember(method,{'LDA','QDA','SVM','RANDOMFOREST'})
        error('[TrainTestModel:] Unsupported classification method');
end

if nargin < 7
    options = {};
end

if length(unique(training_labels)) > 2 || length(unique(testing_labels)) > 2
    error('[TrainTestModel:] only meant for binary classification.');
end

[nan_rows,~] = ind2sub(size(training_set),find(isnan(training_set)));
if ~isempty(nan_rows)
    disp('     Removing nans from training set');
    training_set(unique(nan_rows),:) = [];
    training_labels(unique(nan_rows)) = [];
end

%1 - identify the optimal threshold based on training dataset
if nargout == 3 %desire to return Mdl
    if strcmp(method,'RANDOMFOREST')
        [~,stats, Mdl]=feval(method, training_set, training_set, training_labels, training_labels, options);
        threshold = stats.optThresh;
    else
        error('Output argument Mdl only supported for RANDOMFOREST currently.');
    end
end
        
if nargin < 6 || (exist('threshold','var') && isempty(threshold))
    [~,stats]=feval(method, training_set, training_set, training_labels, training_labels, options);
    threshold = stats.optThresh;
end
 
%2 - get predictions on testing dataset
[methodstring,stats]=feval(method, training_set, testing_set, training_labels, testing_labels, options);

%3 - now apply training threshold
targetclass_name = max(testing_labels);
stats.optThresh = threshold;
stats.decisions = stats.predictions >= stats.optThresh;
stats.tp = length(find(stats.decisions(testing_labels==targetclass_name)==targetclass_name)); 
stats.tn = length(find(stats.decisions(testing_labels~=targetclass_name)~=targetclass_name)); 
stats.fp = length(find(stats.decisions(testing_labels~=targetclass_name)==targetclass_name)); 
stats.fn = length(find(stats.decisions(testing_labels==targetclass_name)~=targetclass_name)); 
stats.acc = (stats.tp+stats.tn)/(stats.tp+stats.tn+stats.fp+stats.fn);
stats.ppv = stats.tp/(stats.tp+stats.fp);
stats.sens = stats.tp/(stats.tp+stats.fn);
stats.spec = stats.tn/(stats.fp+stats.tn);
Pre = ((stats.tp+stats.fp)*(stats.tp+stats.fn) + (stats.tn+stats.fn)*(stats.tn+stats.fp)) / (stats.tp+stats.tn+stats.fp+stats.fn)^2;
stats.kappa = (stats.acc - Pre) / (1 - Pre);        
stats.Fscore = 2*stats.tp/(2*stats.tp+stats.fp+stats.fn);
stats.MCC = (stats.tp*stats.tn - stats.fp*stats.fn) / (((stats.tp+stats.fp)*(stats.fn+stats.tn)*(stats.fp+stats.tn)*(stats.tp+stats.fn))^0.5); %matthew correlation coefficient: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
[FPR,TPR,T,~] = perfcurve(testing_labels,stats.predictions,1);
stats.FPR = FPR;
stats.TPR = TPR;
stats.T = T;
% stats.OPTROCPT = OPTROCPT;

end
