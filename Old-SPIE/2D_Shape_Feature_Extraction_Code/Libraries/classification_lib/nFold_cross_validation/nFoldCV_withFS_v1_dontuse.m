function [stats]= nFoldCV_withFS_v1(data_set,data_labels,params)
error('DO NOT USE. NEEDS TO BE UPDATED');

%***THE ORIGINAL CODE***%

% Using k runs of n-fold cross validation embedded with a fresh feature selection each fold of all k runs
%
% INPUTS:
%      data_set: m x n matrix of m observations and n variables (e.g. observations = patients, variables = features)
%   data_labels: m x 1 column of class labels corresponding to each observation (1 = positive class, -1 = negative class)
%        params: (optional) struct containing the fields:
%                   params.classifier -- options: 'LDA','QDA','SVM','RANDOMFOREST' (Default: 'LDA');
%                   params.svm_kernel (needed for 'SVM' only) --  options: 'linear','quadratic','polynomial','rbf','mlp')
%                   params.fsname -- options: 'mrmr','ttest','wilcoxon','entropy','bhattacharyya','roc') (Default: 'wilcoxon')
%                   params.num_top_feats -- the number of top features to select each iteration
%                   params.feature_idxs: vector of pre-selected set of feature indices (i.e. from pruning or anti-correlation) (Default: include all variables in data_set)
%                   params.shuffle: 1 for random, 0 for non-random partition (Default: 1)
%                   params.n: Number of folds to your cross-validation (Default: 3)
%                   params.nIter: Number of cross-validation iterations (Default: 25)
%                   params.subsets: training and testing subsets for each iteration (Default: computer will generate using 'nFold')
%
% OUTPUTS:
%         stats: struct containing TP, FP, TN, FN, etc. & topfeatinds
%
% The template for this function was originally written by George Lee @2009
% Function has been last modified by Jacob Antunes @2018
%
% ============ example call:====================
% load fisheriris
% 
% data_set = meas;
% data_labels = species;
% params.classifier='LDA';
% params.fsname='wilcoxon';
% params.shuffle = 1;
% params.n = 3;
% params.nIter = 25;
% params.num_top_feats = 1;
% params.subsets = {};
% stats = nFold_AnyClassifier_withFeatureselection(data_set,data_labels,params);
%===============================================

%% Check for valid inputs
if isa(data_labels,'cell'), error('DATA_LABELS must be a non-cell array with labels 1 and -1.'); end
if any(~xor(data_labels == 1, data_labels == -1)), error('DATA_LABELS must be 1 and -1'); end
if size(data_set,1)~=length(data_labels), error('DATA_LABELS must contain same number of labels as observations in DATA_SET'); end

if nargin==3 && ~isstruct(params), error('PARAMS must be a struct'); end

%%  Set default inputs

% if you didn't define any parameters
if nargin < 3
    params.classifier = 'LDA';
    params.fsname = 'wilcoxon';
%     params.feature_list = [];
    params.shuffle = 1;
    params.n = 3;
    params.nIter = 25;
    params.num_top_feats = 1;
    params.feature_idxs = 1:size(data_set,2);
    params.subsets = {};
end

% if you forgot certain parameters. i'll add them in for you
if nargin == 3
    if ~isfield(params,'classifier'),      params.classifier = 'LDA';                      end
    if ~isfield(params,'fsname'),          params.fsname = 'wilcoxon';                     end
    if ~isfield(params,'feature_idxs'),    params.feature_idxs = 1:size(data_set,2);       end
    if ~isfield(params,'shuffle'),         params.shuffle = 1;                             end
    if ~isfield(params,'n'),               params.n = 3;                                   end
    if ~isfield(params,'nIter'),           params.nIter = 25;                              end
    if ~isfield(params,'num_top_feats'),   params.num_top_feats = 1;                       end
    if ~isfield(params,'subsets'),         params.subsets = {};                            end
end
        
%% update dataset based on feature_idxs
data_set = data_set(:,params.feature_idxs);

% % make into 3D struct
% for i = 1:1:length(pts_used)
%     j = 1;
%     while j <= 3
%         new_data_set=(i,data_set(j,:),j);
%         j = j + 1;
%     end
% end

%% nFold CV with FS
stats = struct; %in the end: contains (params.nIter)-rows of information ;
fprintf('\nExecuting %i runs of %i-fold cross-validation using %s with updated feature selection on %i observations\n',params.nIter,params.n,params.classifier,length(data_labels))
pause(0.5);
for j=1:params.nIter
    fprintf('Iteration: %i\n',j);
    
    % reset total statistics
    Ttp = 0; Ttn = 0; Tfp = 0; Tfn = 0; decision=zeros(size(data_labels)); prediction=zeros(size(data_labels));
    
    if isempty(params.subsets)
        [tra, tes]=GenerateSubsets('nFold',data_set,data_labels,params.shuffle,params.n);
    else
        tra = params.subsets(j).training;
        tes = params.subsets(j).testing;
    end
    
    for i=1:params.n
        fprintf(['Fold # ' num2str(i) '\n']);
        
        training_set = data_set(tra{i},:);
        testing_set = data_set(tes{i},:);
        training_labels = data_labels(tra{i});
        testing_labels = data_labels(tes{i});
        
        %adding all three slices to one training set
%         training_set = data_set(

        % *********STEP 1: FEATURE SELECTION*************%
        if strcmp(params.fsname, 'mrmr')
            fea = mrmr_mid_d(training_set,training_labels,params.num_top_feats); % MRMR
            fea = fea'; 
        else
            [IDX, ~] = rankfeatures(training_set',training_labels,'criterion',params.fsname,'CCWeighting', 1); %%CCWeighting is the parameter that weights by correlation w/ previously selected features
            fea(i,:) = IDX(1:params.num_top_feats);
        end    
        
        % *********STEP 2: CLASSIFICATION****************%
        [temp_stats,~] = Classify(params.classifier, training_set(:,fea) , testing_set(:,fea), training_labels(:), testing_labels(:));
       
        Ttp = Ttp + temp_stats.tp;
        Ttn = Ttn + temp_stats.tn;
        Tfp = Tfp + temp_stats.fp;
        Tfn = Tfn + temp_stats.fn;
        decision(tes{i}) = temp_stats.decision;
        prediction(tes{i}) = temp_stats.prediction;
    end
    decision(decision==0) = -1;
    
    % output statistics
    if numel(unique(data_labels))>1 %numel(unique(testing_labels))>1
        if params.n == 1
            [FPR,TPR,~,AUC,~,~,~] = perfcurve(data_labels(tes{i}),prediction(tes{i}),1);
        else
            [FPR,TPR,~,AUC,~,~,~] = perfcurve(data_labels,prediction,1);
        end
        stats(j).AUC = AUC;
        stats(j).TPR = TPR;
        stats(j).FPR = FPR;
    else
        stats(j).AUC = [];
        stats(j).TPR = [];
        stats(j).FPR = [];
    end
    
    stats(j).tp = Ttp;
    stats(j).tn = Ttn;
    stats(j).fp = Tfp;
    stats(j).fn = Tfn;
    stats(j).acc = (Ttp+Ttn)/(Ttp+Ttn+Tfp+Tfn);
    stats(j).ppv = Ttp/(Ttp+Tfp);
    stats(j).sens = Ttp/(Ttp+Tfn);
    stats(j).spec = Ttn/(Tfp+Ttn);
    stats(j).subsets.training = tra;
    stats(j).subsets.testing = tes;
    stats(j).decision = decision;
    stats(j).prediction = prediction;
    Pre = ((Ttp+Tfp)*(Ttp+Tfn) + (Ttn+Tfn)*(Ttn+Tfp)) / (Ttp+Ttn+Tfp+Tfn)^2;
    stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);
    stats(j).topfeatinds = params.feature_idxs(fea);
end
            
fprintf('\nDONE.\n');