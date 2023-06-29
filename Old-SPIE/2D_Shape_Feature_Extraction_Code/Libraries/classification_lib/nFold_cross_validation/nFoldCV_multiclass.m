function [stats]= nFoldCV_multiclass(data_set,data_labels,params)

% Using k runs of n-fold cross validation embedded with a fresh feature selection each fold of all k runs
%***PERFORMS PERFCURVE ACROSS ALL FOLDS SEPARATELY***%
%
% INPUTS:
%      data_set: M x N matrix of m observations and n variables (e.g. observations = patients, variables = features)
%   data_labels: M x 1 column of class labels corresponding to each observation (1 = positive class, -1 = negative class)
%        params: (optional) struct containing the fields:
%                   params.classifier: -- 'RANDOMFOREST' (default),'SVM' are only implementations so far
%                   params.classifieroptions: (optional) A dditional {...,'Name','Value'...'} parameters for your classifier. Pass in as a cell array of strings
%                   params.num_top_feats: -- N/A for now
%                   params.feature_idxs: (optional) vector of pre-selected set of feature indices (i.e. from pruning or anti-correlation) (Default: include all variables in data_set)
%                   params.shuffle: 1 for random, 0 for non-random partition (Default: 1)
%                   params.n: Number of folds to your cross-validation (Default: 3)
%                   params.nIter: Number of cross-validation iterations (Default: 25)
%                   params.subsets: (optional) training and testing subsets for each iteration (Default: computer will generate using 'nFold')
%                   params.classbalancetestfold: (option) how to class balance the testing folds -- options: 'none' (default),'equal'
%                   params.patient_ids: (needed for augmented data, i.e. multiple slices per patient) -- M x 1 vector of unique patient identifier (an integer) from which each data_set and data_label came from 
%                   params.osname: (optional) string for oversampling method) -- options: 'SMOTE','ADASYN','UNDERSAMPLE'
%                   params.remouts: (optional) string for removing outlier method -- options: 'median','mean','quartile'
%                   params.decsmeth: (optional) string for how to assign estimated class label on testing observation based on scores -- options: 'max' (Default), 'matlab'
%                   xxx params.threshmeth: (optional) string of optimal ROC threshold method -- options: 'euclidean', 'matlab' (Default: 'matlab')
%
% OUTPUTS:
%         stats: struct containing TP, FP, TN, FN, etc. & topfeatinds
%
% The template for this function was originally written by George Lee @2009
% Function has been modified by Jacob Antunes @2018
% Updates provided by Prathyush Chirra, Nitya Talasila, Neil Sehgal @2018
%
% UPDATES:
% 04-021-2019 - Initial Development

%% Check for valid inputs
% if isa(data_labels,'cell'), error('DATA_LABELS must be a non-cell array with labels 1 and -1.'); end
% if all([any(~xor(data_labels == 1, data_labels == -1)), any(~xor(data_labels == 1, data_labels == 0))]), error('DATA_LABELS must be 1 and -1 (or 1 and 0)'); end
if size(data_set,1)~=length(data_labels), error('DATA_LABELS must contain same number of labels as observations in DATA_SET (equal # of rows)'); end
if length(unique(data_labels))~=3, error('%s.m requires exactly 3 distinct classes',mfilename); end
if nargin==3 && ~isstruct(params), error('PARAMS must be a struct'); end    

%%  Set default inputs

%add in all necessary subfolders
funcpath = which(mfilename);
funcdir = funcpath(1:end-length(mfilename));
idx1 = strfind(funcdir,'\'); idx2 = strfind(funcdir,'/');
if isempty(idx1)
    codepath = funcdir(1:idx2(end-1));
else
    codepath = funcdir(1:idx1(end-1));
end
addpath(genpath(codepath));


% if you didn't define any parameters
if nargin < 3
    params.classifier = 'RANDOMFOREST';
    params.classifieroptions = {};
%     params.feature_list = [];
    params.shuffle = 1;
    params.n = 3;
    params.nIter = 25;
    params.num_top_feats = 1;
    params.feature_idxs = 1:size(data_set,2);
%     params.threshmeth = 'euclidean';
%     params.classbalancetestfold = 'none';
    params.decsmeth = 'max';
end

% if you forgot certain parameters. i'll add them in for you
if nargin == 3
    if ~isfield(params,'classifier')           || isempty(params.classifier),      params.classifier = 'RANDOMFOREST';                      end
    if ~isfield(params,'classifieroptions'),                                       params.classifieroptions = {};                  end
    if ~isfield(params,'shuffle')              || isempty(params.shuffle),         params.shuffle = 1;                             end
    if ~isfield(params,'n')                    || isempty(params.n),               params.n = 3;                                   end
    if ~isfield(params,'nIter')                || isempty(params.nIter),           params.nIter = 25;                              end
    if ~isfield(params,'num_top_feats')        || isempty(params.num_top_feats),   params.num_top_feats = 1;                       end
    if ~isfield(params,'feature_idxs')         || isempty(params.feature_idxs),    params.feature_idxs = 1:size(data_set,2);       end   
   % if ~isfield(params,'threshmeth')           || isempty(params.threshmeth),      params.threshmeth = 'euclidean';                end
%     if ~isfield(params,'classbalancetestfold') || isempty(params.classbalancetestfold), params.classbalancetestfold ='none';       end
    if ~isfield(params,'decsmethds')           || isempty(params.decsmeth),        params.decsmeth = 'max';                        end
end  

%% Final checks
if params.n == size(data_set,1), error('LOOCV not implemented for %s',mfilename); end

%% nFold CV with FS
stats = struct; %in the end: contains (params.nIter)-rows of information ;

data_set = data_set(:,params.feature_idxs);
% data_labels = double(data_labels);

fprintf('\nExecuting %i run(s) of %i-fold cross-validation for %s to compute Gini Indices',params.nIter,params.n,params.classifier)
if isfield(params,'osname') && ~isempty(params.osname)
    fprintf('\t Oversampling of training data to be performed using %s method\n',params.osname);
end
if isfield(params,'remouts') && ~isempty(params.remouts)
    fprintf('\t Training data will have outliers removed using %s method\n',params.remouts);
%     warning('Not recommended. Classifier will probably fail due to NaNs in data');
end

pause(0.5);
if isfield(params,'patient_ids') && ~isempty(params.patient_ids)
    [pts,idxs] = unique(params.patient_ids,'stable');
end

for j=1:params.nIter
    fprintf('Iteration: %i\n',j);
    
    % reset total statistics
    Ttp = 0; Ttn = 0; Tfp = 0; Tfn = 0; feats = [];
    all_est_labels  = cell(length(data_labels),1);%-inf*ones(length(data_labels),1);
    all_est_scores = -inf*ones(length(data_labels),3);
    all_labels = cell(length(data_labels),1);;%-inf*ones(length(data_labels),1);
    all_cmat = zeros(length(unique(data_labels)));
    foldpredictions = cell(1,params.n);
    foldlabels = cell(1,params.n);
    
    if ~isfield(params,'subsets') || isempty(params.subsets)
        if isfield(params,'patient_ids') && ~isempty(params.patient_ids)
            [tra, tes]=GenerateSubsets('nFold_multiclass',[],data_labels(idxs),params.shuffle,params.n); %for augmented data
        else
            [tra, tes] = GenerateSubsets('nFold_multiclass',[],data_labels,params.shuffle,params.n);
        end
    else
        tra = params.subsets(j).training;
        tes = params.subsets(j).testing;
    end
    
    for i=1:params.n
        fprintf(['Fold # ' num2str(i) '\n']);

        if isfield(params,'patient_ids') && ~isempty(params.patient_ids)
            training_set = data_set(ismember(params.patient_ids,pts(tra{i})),:);
            testing_set = data_set(ismember(params.patient_ids,pts(tes{i})),:);
            training_labels = data_labels(ismember(params.patient_ids,pts(tra{i})));
            testing_labels = data_labels(ismember(params.patient_ids,pts(tes{i})));
        else
            training_set = data_set(tra{i},:);
            testing_set = data_set(tes{i},:);
            training_labels = data_labels(tra{i});
            testing_labels = data_labels(tes{i}); 
        end 
        
        %oversampling?
        if isfield(params,'osmethod') && ~isempty(params.osmethod)
            [training_set, training_labels] = Oversample(params.osname,training_set,training_labels);
        end
        
        
%         %class balance testing fold?
%         if isfield(params,'classbalancetestfold') && ~isempty(params.classbalancetestfold) && strcmp(params.classbalancetestfold,'equal')
%             max_label = max(testing_labels);
%             min_label = min(testing_labels);
%             [~, max_idx] = max([nnz(testing_labels==max_label),nnz(testing_labels==min_label)]);
%             if max_idx==1
%                 major_label = max_label;
%                 minor_label = min_label;
%             elseif max_idx==2
%                 major_label = min_label;
%                 minor_label = max_label;
%             end
%             max_idx = find(testing_labels==major_label);
%             min_idx = find(testing_labels==minor_label);
%             max_idx2keep = datasample(max_idx,length(min_idx),'Replace',false);
%             testing_set = testing_set(sort([min_idx; max_idx2keep]),:);
%             testing_labels = testing_labels(sort([min_idx; max_idx2keep]));
%             tes{i} = tes{i}(sort([min_idx; max_idx2keep]));
%         end
        
        %remove outliers?
        training_set_copy = training_set;
        if isfield(params,'remouts') && ~isempty(params.remouts)
            TF = isoutlier(training_set,params.remouts);
            training_set(TF) = nan;
        end
   
        % *********STEP: CLASSIFICATION****************%
%         options = {'threshmeth',params.threshmeth};
%         options = cat(2,options,params.classifieroptions);

        training_set = training_set_copy;

        if strcmp(params.classifier,'RANDOMFOREST')
            Mdl = TreeBagger(50,training_set,training_labels,'OOBPrediction','on','OOBPredictorImportance','on','SplitCriterion','gdi','PredictorSelection','allsplits');
            imp = Mdl.OOBPermutedPredictorDeltaMeanMargin;
            [est_labels,est_scores,~] = predict(Mdl,testing_data);
            
            foldimp{i} = imp;
          
       
        elseif strcmp(params.classifier,'SVM')
            Mdl = fitcecoc(training_set,training_labels,'Learners',templateSVM('KernelFunction','polynomial'));
            [est_labels,est_scores,~] = predict(Mdl,testing_set);
            if strcmp(params.decsmeth,'max')
                u = unique(data_labels);
                est_labels = u(maxthresholding(est_scores));
            end %otherwise just use matlab's predict decisions
        else
           error('invalid classifier defined. Only RANDOMFOREST and SVM supported for multiclass nfold cv'); 
        end
        
        if isfield(params,'patient_ids')
            all_est_scores(ismember(params.patient_ids,pts(tes{i})),:) = est_scores;
            all_est_labels(ismember(params.patient_ids,pts(tes{i})),:) = est_labels;
            all_labels(ismember(params.patient_ids,pts(tes{i}))) = testing_labels;
        else
            all_scores(tes{i},:) = est_scores;
            all_est_labels(tes{i}) = est_labels;
            all_labels(tes{i}) = testing_labels;
        end
        cmat = confusionmat(testing_labels,est_labels);
        all_cmat = all_cmat + cmat;
        
        %for generating upper and lower AUC bounds
        foldcmat{i} = confusionmat(testing_labels,est_labels);
        foldpreds{i} = est_scores;
        foldlabels{i} = testing_labels;
        folddecisions{i} = est_labels;
    end
    
    % output statistics          
        %compute perfcurve across each testing fold per iteration
% %         [FPR,TPR,T,AUC,OPTROCPT] = perfcurve(foldlabels, foldpredictions, 1,'XVals', [0:0.02:1]);
%         
%         if strcmp(params.threshmeth,'euclidean')
%             [~, optim_idx] = determine_threshold_euclidean_distance(FPR,TPR,T);
%         elseif strcmp(params.threshmeth,'matlab')
%             optim_idx = find(FPR == OPTROCPT(1) & TPR == OPTROCPT(2));
%         else
%             error('Invalid "threshmeth" field option in params struct.');
%         end
    
%         stats(j).AUCs = AUC;
%         stats(j).optThresh = T(optim_idx(1));
%         stats(j).thresholdmethod = params.threshmeth;
        stats(j).predictions = all_est_scores;
        stats(j).decisions = all_est_labels;
        stats(j).confusionmat = all_cmat;
        stats(j).acc = sum(diag(all_cmat))/sum(all_cmat(:));
       
        stats(j).subsets.training = tra;
        stats(j).subsets.testing = tes;
    
end
            
fprintf('\nDONE.\n');