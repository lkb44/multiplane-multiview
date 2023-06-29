function [stats]= nFoldCV_withFS_v4(data_set,data_labels,params)

error('DO NOT USE. NEEDS TO BE UPDATED');

%***PERFORMS PERFCURVE ACROSS ALL FOLDS SEPARATELY. VOTES ACROSS SLICES PER PATIENT***%

% Using k runs of n-fold cross validation embedded with a fresh feature selection each fold of all k runs
%
% INPUTS:
%      data_set: M x N matrix of m observations and n variables (e.g. observations = patients, variables = features)
%   data_labels: M x 1 column of class labels corresponding to each observation (e.g. 1 = positive class, -1 = negative class)[any non-positive value is acceptable]
%        params: (optional) struct containing the fields:
%                   params.classifier -- options: 'LDA','QDA','SVM','RANDOMFOREST' (Default: 'LDA');
%                   params.svm_kernel (needed for 'SVM' only) --  options: 'linear','quadratic','polynomial','rbf','mlp')
%                   params.fsname -- options: 'mrmr','ttest','wilcoxon','entropy','bhattacharyya','roc') (Default: 'wilcoxon')
%                   params.num_top_feats -- the number of top features to select each iteration
%                   params.feature_idxs: vector of pre-selected set of feature indices (i.e. from pruning or anti-correlation) (Default: include all features in data_set)
%                   params.shuffle: 1 for random, 0 for non-random partition (Default: 1)
%                   params.n: Number of folds to your cross-validation (Default: 3)
%                   params.classbal: (optional) for class balancing on training set. 0=leave default ratios,1=undersample larger class,2=oversample smaller class. (Default: 0)
%                   params.nIter: Number of cross-validation iterations (Default: 25)
%                   params.subsets: training and testing subsets for each iteration (Default: computer will generate using 'nFold')
%                   params.patient_ids: (needed for augmented data, i.e. multiple slices per patient) -- M x 1 vector of unique patient identifier (an integer) from which each data_set and data_label came from 
%                   
% OUTPUTS:
%         stats: struct containing TP, FP, TN, FN, etc. & topfeatinds
%
% The template for this function was originally written by George Lee @2009
% Function has been modified by Jacob Antunes @2018
% Updates provided by Prathyush Chirra, Nitya Talasila @2018
%
% UPDATES:
% 01-01-2009 - Initial Development
% 03-01-2018 - Modified to include updated feature selection each iteration
% 04-25-2018 - Modified to include option for handling multiple slices per patient ("augmented" data)
% 05-01-2018 - Modified to include option for class balancing (see param.classbal field and nFold.m)

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
if isa(data_labels,'cell'), error('DATA_LABELS must be a non-cell array of positive and non-positive labels.'); end
if any(~xor(data_labels > 0, data_labels <= 0)), error('DATA_LABELS must contain least one positive and non-positive label'); end
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
    params.classbal = 0;
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
    if ~isfield(params,'classbal'),        params.classbal = 0;                            end
end
        

%% nFold CV with FS
stats = struct; %in the end: contains (params.nIter)-rows of information ;
fprintf('\nExecuting %i runs of %i-fold cross-validation using %s with updated feature selection on %i observations\n',params.nIter,params.n,params.classifier,length(data_labels))
pause(0.5);
if isfield(params,'patient_ids')
    [pts,idxs] = unique(params.patient_ids);
end
 
for j=1:params.nIter
    fprintf('Iteration: %i\n',j);

    if isempty(params.subsets)
        if isfield(params,'patient_ids')
            [tra, tes]=GenerateSubsets('nFold',[],data_labels(idxs),params.shuffle,params.n);
        else
            [tra, tes]=GenerateSubsets('nFold',[],data_labels,params.shuffle,params.n);
        end
    else
        tra = params.subsets(j).training;
        tes = params.subsets(j).testing;
    end
    
    for i=1:params.n
        fprintf(['Fold # ' num2str(i) '\n']);

        if isfield(params,'patient_ids')
            training_set = data_set(ismember(params.patient_ids,pts(tra{i})),:);
            testing_set = data_set(ismember(params.patient_ids,pts(tes{i})),:);
            training_labels = data_labels(ismember(params.patient_ids,pts(tra{i})));
            testing_labels = data_labels(ismember(params.patient_ids,pts(tes{i})));
            testing_patient_ids = params.patient_ids(ismember(params.patient_ids,pts(tes{i})));
        else
            training_set = data_set(tra{i},:);
            testing_set = data_set(tes{i},:);
            training_labels = data_labels(tra{i});
            testing_labels = data_labels(tes{i});  
        end     

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
       
        p = []; d = []; l = [];
        if isfield(params,'patient_ids')
            for k = pts(tes{i})'
                p = [p; double(median(temp_stats.prediction(ismember(testing_patient_ids,k))))];
                if p(end)>=0.5,
                    d = [d; max(unique(data_labels))];
                else
                    d = [d; min(unique(data_labels))];
                end
                l = [l; data_labels(find(ismember(params.patient_ids,k),1))];
            end
            decisions{i} = d;
            predictions{i} = p;
            labels{i} = l;
        else
            decisions(tes{i}) = temp_stats.decision;
            predictions{i} = temp_stats.prediction;
            labels{i} = testing_labels;
        end
    end
    [~,~,~,AUCs] = perfcurve(labels, predictions, 1,'XVals', [0:0.02:1]);
    
    % output statistics    
    stats(j).AUCs = AUCs;
    stats(j).subsets.training = tra;
    stats(j).subsets.testing = tes;
    stats(j).predictions = predictions;
    stats(j).decisions = decisions;
    stats(j).labels = labels;
    stats(j).topfeatinds = params.feature_idxs(fea);
    
end
            
fprintf('\nDONE.\n');