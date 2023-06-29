function [stats, methodstring,Mdl] = TrainTestModel(method, training_set, testing_set, training_labels, testing_labels, threshold, options)
% Wrapper function to call any one of CCIPD's implemented classifier
% methods in a training and validation fashion

% THRESHOLD: ideally passed in, otherwise will compute from training dataset

if ~ismember(method,{'LDA','QDA','SVM','RANDOMFOREST'})
        error('Unsupported classification method');
end

if nargin < 7
    options = {};
end

if nargin < 6
    [methodstring,stats]=feval(method, training_set, training_set, training_labels, training_labels, options);
    threshold = stats.optThresh;
end

%get predictions on testing dataset
[methodstring,stats]=feval(method, training_set, testing_set, training_labels, testing_labels, options);

%now apply training threshold
targetclass_name = max(testing_labels);
stats.optThresh = threshold;
stats.decision = stats.prediction >= stats.optThresh;
stats.tp = length(find(stats.decision(testing_labels==targetclass_name)==targetclass_name)); 
stats.tn = length(find(stats.decision(testing_labels~=targetclass_name)~=targetclass_name)); 
stats.fp = length(find(stats.decision(testing_labels~=targetclass_name)==targetclass_name)); 
stats.fn = length(find(stats.decision(testing_labels==targetclass_name)~=targetclass_name)); 
stats.acc = (stats.tp+stats.tn)/(stats.tp+stats.tn+stats.fp+stats.fn);
stats.ppv = stats.tp/(stats.tp+stats.fp);
stats.sens = stats.tp/(stats.tp+stats.fn);
stats.spec = stats.tn/(stats.fp+stats.tn);
Pre = ((stats.tp+stats.fp)*(stats.tp+stats.fn) + (stats.tn+stats.fn)*(stats.tn+stats.fp)) / (stats.tp+stats.tn+stats.fp+stats.fn)^2;
stats.kappa = (stats.acc - Pre) / (1 - Pre);        
stats.Fscore = 2*stats.tp/(2*stats.tp+stats.fp+stats.fn);
stats.OPTROCPT = stats.OPTROCPT;

end
