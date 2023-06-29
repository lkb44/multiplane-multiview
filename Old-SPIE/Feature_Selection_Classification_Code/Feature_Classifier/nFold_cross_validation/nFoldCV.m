function [stats]= nFoldCV(data_set,data_labels,params)

% Using k runs of n-fold cross validation
%***PERFORMS PERFCURVE ACROSS ALL FOLDS SEPARATELY***%
%
% INPUTS:
%      data_set: M x N matrix of m observations and n variables (e.g. observations = patients, variables = features)
%   data_labels: M x 1 column of class labels corresponding to each observation (1 = positive class, -1 = negative class)
%        params: (optional) struct containing the fields:
%                   params.classifier: -- options: 'LDA','QDA','SVM','RANDOMFOREST' (Default: 'LDA');
%                   params.classifieroptions: (optional) A dditional {...,'Name','Value'...'} parameters for your classifier. Pass in as a cell array of strings
%                   params.feature_idxs: (optional) vector of pre-selected set of feature indices (i.e. from pruning or anti-correlation) (Default: include all variables in data_set)
%                   params.shuffle: 1 for random, 0 for non-random partition (Default: 1)
%                   params.n: Number of folds to your cross-validation (Default: 3)
%                   params.nIter: Number of cross-validation iterations (Default: 25)
%                   params.subsets: (optional) training and testing subsets for each iteration (Default: computer will generate using 'nFold')
%                   params.classbalancetestfold: (option) how to class balance the testing folds -- options: 'none' (default),'equal'
%                   params.patient_ids: (needed for augmented data, i.e. multiple slices per patient) -- M x 1 vector of unique patient identifier (an integer) from which each data_set and data_label came from 
%                   params.osname: (optional) string for oversampling method) -- options: 'SMOTE','ADASYN','UNDERSAMPLE'
%                   params.remouts: (optional) string for removing outlier method -- options: 'median','mean','quartile'
%                   params.threshmeth: (optional) string of optimal ROC threshold method -- options: 'euclidean', 'matlab' (Default: 'euclidean')
%
% OUTPUTS:
%         stats: struct containing TP, FP, TN, FN, etc.
%
% The template for this function was originally written by George Lee @2009
% Function has been modified by Jacob Antunes @2018
% Updates provided by Prathyush Chirra, Nitya Talasila, Neil Sehgal @2018
%
% UPDATES:
% 01-01-2009 - Initial Development
% 03-01-2018 - Modified to include updated feature selection each iteration
% 04-25-2018 - Modified to include option for handling multiple slices per patient ("augmented" data)
% 05-02-2018 - Modified to include option for oversampling data
% 05-22-2018 - Modified to include option for removing outliers from training data
%            - Can handle 0s and 1s as well as -1s and 1s class labels
% 07-03-2018 - Small fixes to handle empty fields in params
% 07-18-2018 - Added in stats for LOOCV, can now pass in classifier options, automatically adds all necessary subfolders into path.
% 10-08-2018 - Added params field "threshmeth." Reworked stats to compute once thru for all testing patients across all folds (verses summing them up).
% 02-06-2019 - Added params field "classbalancetestfold."
% 11-13-2019 - Removed feature selection (as this is nFoldCV without FS)
% 05-20-2020 - Added in AUPRC calculations

% ============ example call:====================
% load fisheriris
% 
% addpath('.')
% 
% data_set = meas;
% data_labels = (contains(species,'virginica'));
% params.classifier='QDA';
% params.shuffle = 1;
% params.n = 3;
% params.nIter = 25;
% params.threshmeth = 'matlab';
% stats = nFoldCV(data_set,data_labels,params);
%===============================================

%% Check for valid inputs
if isa(data_labels,'cell'), error('[nFoldCV:] DATA_LABELS must be a non-cell array with labels 1 and -1.'); end
if all([any(~xor(data_labels == 1, data_labels == -1)), any(~xor(data_labels == 1, data_labels == 0))]), error('[nFoldCV:] DATA_LABELS must be 1 and -1 (or 1 and 0)'); end
if size(data_set,1)~=length(data_labels), error('[nFoldCV:] DATA_LABELS must contain same number of labels as observations in DATA_SET (equal # of rows)'); end
if nargin==3 && ~isstruct(params), error('[nFoldCV:] PARAMS must be a struct'); end    

%%  Set default inputs

%add in all necessary subfolders
funcname = 'nFoldCV.m';
funcpath = which(funcname);
funcdir = funcpath(1:end-length(funcname));
idx1 = strfind(funcdir,'\'); idx2 = strfind(funcdir,'/');
if isempty(idx1)
    codepath = funcdir(1:idx2(end-1));
else
    codepath = funcdir(1:idx1(end-1));
end
addpath(genpath(codepath));


% if you didn't define any parameters
if nargin < 3
    params.classifier = 'LDA';
    params.classifieroptions = {};
%     params.feature_list = [];
    params.shuffle = 1;
    params.n = 3;
    params.nIter = 25;
    params.feature_idxs = 1:size(data_set,2);
    params.threshmeth = 'euclidean';
end

% if you forgot certain parameters. i'll add them in for you
if nargin == 3
    if ~isfield(params,'classifier')           || isempty(params.classifier),      params.classifier = 'LDA';                      end
    if ~isfield(params,'classifieroptions'),                                       params.classifieroptions = {};                  end
    if ~isfield(params,'shuffle')              || isempty(params.shuffle),         params.shuffle = 1;                             end
    if ~isfield(params,'n')                    || isempty(params.n),               params.n = 3;                                   end
    if ~isfield(params,'nIter')                || isempty(params.nIter),           params.nIter = 25;                              end
    if ~isfield(params,'feature_idxs')         || isempty(params.feature_idxs),    params.feature_idxs = 1:size(data_set,2);       end   
    if ~isfield(params,'threshmeth')           || isempty(params.threshmeth),      params.threshmeth = 'euclidean';                end
    if ~isfield(params,'classbalancetestfold') || isempty(params.classbalancetestfold), params.classbalancetestfold ='none';       end
end  

%% nFold CV 
stats = struct; %in the end: contains (params.nIter)-rows of information ;

data_set = data_set(:,params.feature_idxs);
data_labels = double(data_labels);

if params.n == size(data_set,1)
    fprintf('\nExecuting leave-one-out cross-validation using %s on %i observations\n',params.classifier,length(data_labels))
else
    fprintf('\nExecuting %i run(s) of %i-fold cross-validation using %s on %i observations\n',params.nIter,params.n,params.classifier,length(data_labels))
end
if isfield(params,'osname') && ~isempty(params.osname) && ~strcmp(params.osname,'none')
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
    decisions  = -inf*ones(length(data_labels),1);
    predictions = -inf*ones(length(data_labels),1);
    labels = -inf*ones(length(data_labels),1);
    foldpredictions = cell(1,params.n);
    foldlabels = cell(1,params.n);
    
    if ~isfield(params,'subsets') || isempty(params.subsets)
        if isfield(params,'patient_ids') && ~isempty(params.patient_ids)
            [tra, tes]=GenerateSubsets('nFold',[],data_labels(idxs),params.shuffle,params.n); %for augmented data
        else
            [tra, tes] = GenerateSubsets('nFold',[],data_labels,params.shuffle,params.n);
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
        if isfield(params,'osname') && ~isempty(params.osname) && ~strcmp(params.osname,'none')
            [training_set, training_labels] = Oversample(params.osname,training_set,training_labels);
        end
        
        
        %class balance testing fold?
        if isfield(params,'classbalancetestfold') && ~isempty(params.classbalancetestfold) && strcmp(params.classbalancetestfold,'equal')
            max_label = max(testing_labels);
            min_label = min(testing_labels);
            [~, max_idx] = max([nnz(testing_labels==max_label),nnz(testing_labels==min_label)]);
            if max_idx==1
                major_label = max_label;
                minor_label = min_label;
            elseif max_idx==2
                major_label = min_label;
                minor_label = max_label;
            end
            max_idx = find(testing_labels==major_label);
            min_idx = find(testing_labels==minor_label);
            max_idx2keep = datasample(max_idx,length(min_idx),'Replace',false);
            testing_set = testing_set(sort([min_idx; max_idx2keep]),:);
            testing_labels = testing_labels(sort([min_idx; max_idx2keep]));
            tes{i} = tes{i}(sort([min_idx; max_idx2keep]));
        end
        
        %remove outliers?
        training_set_copy = training_set;
        if isfield(params,'remouts') && ~isempty(params.remouts)
            TF = isoutlier(training_set,params.remouts);
            training_set(TF) = nan;
        end
        
        % *********STEP 2: CLASSIFICATION****************%
        options = {'threshmeth',params.threshmeth};
        options = cat(2,options,params.classifieroptions);
        % Go back to using the original set of observations for building model
        training_set = training_set_copy;
        [temp_stats,~] = Classify_wrapper(params.classifier, training_set , testing_set, training_labels(:), testing_labels(:), options);

        if numel(unique(testing_labels(:))) > 1 %won't work for Leave One Out CV
            
            if isfield(params,'patient_ids')
                decisions(ismember(params.patient_ids,pts(tes{i}))) = temp_stats.decisions;
            else
                decisions(tes{i}) = temp_stats.decisions;
            end
            
            if any(ismember(data_labels,-1))
                decisions(decisions==0) = -1;
            end
        end
        
        if isfield(params,'patient_ids')
            predictions(ismember(params.patient_ids,pts(tes{i}))) = temp_stats.predictions;
            labels(ismember(params.patient_ids,pts(tes{i}))) = testing_labels;
        else
            predictions(tes{i}) = temp_stats.predictions;
            labels(tes{i}) = testing_labels;
        end

        %for generating upper and lower AUC bounds
        foldpredictions{i} =  temp_stats.predictions;
        foldlabels{i} = testing_labels;
        
        %***for testing: call this script at this point. Otherwise, leave commented***%
%         temp_analysis_nFoldCV_withFS;

    end
    
    % output statistics  
    if numel(unique(testing_labels)) > 1 %won't work for Leave One Out CV
        
        %compute perfcurve across each testing fold per iteration
        [FPR,TPR,T,AUC,OPTROCPT] = perfcurve(foldlabels, foldpredictions, 1,'XVals', [0:0.02:1]);
        [RECALL,PRECISION,~,AUPRC,~] = perfcurve(foldlabels, foldpredictions, 1,'XVals', [0:0.02:1],'XCrit','reca','YCrit','prec');

        if strcmp(params.threshmeth,'euclidean')
            [~, optim_idx] = determine_threshold_euclidean_distance(FPR,TPR,T);
        elseif strcmp(params.threshmeth,'matlab')
            optim_idx = find(FPR == OPTROCPT(1) & TPR == OPTROCPT(2));
        else
            error('[nFoldCV_withFS:] Invalid "threshmeth" field option in params struct.');
        end
    
        stats(j).AUCs = AUC;
        stats(j).AUPRCs = AUPRC;
        stats(j).PRECISION = PRECISION;
        stats(j).RECALL = RECALL;
        stats(j).optThresh = T(optim_idx(1));
        stats(j).thresholdmethod = params.threshmeth;
        stats(j).predictions = predictions;
        stats(j).decisions = stats(j).predictions >= stats(j).optThresh;
        stats(j).tp = length(find(stats(j).decisions(labels==1)==1)); 
        stats(j).tn = length(find(stats(j).decisions(labels~=1)~=1)); 
        stats(j).fp = length(find(stats(j).decisions(labels~=1)==1)); 
        stats(j).fn = length(find(stats(j).decisions(labels==1)~=1)); 
        stats(j).acc = (stats(j).tp+stats(j).tn)/(stats(j).tp+stats(j).tn+stats(j).fp+stats(j).fn);
        stats(j).ppv = stats(j).tp/(stats(j).tp+stats(j).fp);
        stats(j).sens = stats(j).tp/(stats(j).tp+stats(j).fn);
        stats(j).spec = stats(j).tn/(stats(j).fp+stats(j).tn);
        Pre = ((stats(j).tp+stats(j).fp)*(stats(j).tp+stats(j).fn) + (stats(j).tn+stats(j).fn)*(stats(j).tn+stats(j).fp)) / (stats(j).tp+stats(j).tn+stats(j).fp+stats(j).fn)^2;
        stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);        
        stats(j).Fscore = 2*stats(j).tp/(2*stats(j).tp+stats(j).fp+stats(j).fn);
        stats(j).MCC = (stats(j).tp*stats(j).tn - stats(j).fp*stats(j).fn) / (((stats(j).tp+stats(j).fp)*(stats(j).fn+stats(j).tn)*(stats(j).fp+stats(j).tn)*(stats(j).tp+stats(j).fn))^0.5); %matthew correlation coefficient: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
        stats(j).OPTROCPT = OPTROCPT;
        stats(j).FPR = FPR;
        stats(j).TPR = TPR;
        stats(j).T = T;
        Pre = ((stats(j).tp+stats(j).fp)*(stats(j).tp+stats(j).fn) + (stats(j).tn+stats(j).fn)*(stats(j).tn+stats(j).fp)) ...
            / (stats(j).tp+stats(j).tn+stats(j).fp+stats(j).fn)^2;
        stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);
        stats(j).Fscore = 2*stats(j).tp/(2*stats(j).tp+stats(j).fp+stats(j).fn);
        stats(j).subsets.training = tra;
        stats(j).subsets.testing = tes;
    else % for leave one out
        try
            [FPR,TPR,T,AUC,OPTROCPT] = perfcurve(foldlabels, foldpredictions, 1,'XVals', [0:0.02:1]);
            [RECALL,PRECISION,~,AUPRC,~] = perfcurve(foldlabels, foldpredictions, 1,'XVals', [0:0.02:1],'XCrit','reca','YCrit','prec');
        catch %unresolved error
            if iscolumn(foldlabels{1})
                l = cell2mat(foldlabels');
                p = cell2mat(foldpredictions');
            else
                l = cell2mat(foldlabels);
                p = cell2mat(foldpredictions);
            end
            [FPR,TPR,T,AUC,OPTROCPT] = perfcurve(l(:), p(:), 1,'XVals', [0:0.02:1]);
            [RECALL,PRECISION,~,AUPRC,~] = perfcurve(l(:), p(:), 1,'XVals', [0:0.02:1],'XCrit','reca','YCrit','prec');
        end
        
        if strcmp(params.threshmeth,'euclidean')
            [~, optim_idx] = determine_threshold_euclidean_distance(FPR,TPR,T);
        elseif strcmp(params.threshmeth,'matlab')
            optim_idx = find(FPR == OPTROCPT(1) & TPR == OPTROCPT(2));
        else
            error('[nFoldCV_withFS:] Invalid "threshmeth" field option in params struct.');
        end
        

        stats(j).AUCs = AUC;
        stats(j).AUPRCs = AUPRC;
        stats(j).PRECISION = PRECISION;
        stats(j).RECALL = RECALL;
        stats(j).optThresh = T(optim_idx);
        stats(j).thresholdmethod = params.threshmeth;
        stats(j).predictions = predictions;
        stats(j).decisions = stats(j).predictions >= stats(j).optThresh(1);
        stats(j).tp = length(find(stats(j).decisions(labels==1)==1)); 
        stats(j).tn = length(find(stats(j).decisions(labels~=1)~=1)); 
        stats(j).fp = length(find(stats(j).decisions(labels~=1)==1)); 
        stats(j).fn = length(find(stats(j).decisions(labels==1)~=1)); 
        stats(j).acc = (stats(j).tp+stats(j).tn)/(stats(j).tp+stats(j).tn+stats(j).fp+stats(j).fn);
        stats(j).ppv = stats(j).tp/(stats(j).tp+stats(j).fp);
        stats(j).sens = stats(j).tp/(stats(j).tp+stats(j).fn);
        stats(j).spec = stats(j).tn/(stats(j).fp+stats(j).tn);
        Pre = ((stats(j).tp+stats(j).fp)*(stats(j).tp+stats(j).fn) + (stats(j).tn+stats(j).fn)*(stats(j).tn+stats(j).fp)) / (stats(j).tp+stats(j).tn+stats(j).fp+stats(j).fn)^2;
        stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);        
        stats(j).Fscore = 2*stats(j).tp/(2*stats(j).tp+stats(j).fp+stats(j).fn);
        stats(j).MCC = (stats(j).tp*stats(j).tn - stats(j).fp*stats(j).fn) / (((stats(j).tp+stats(j).fp)*(stats(j).fn+stats(j).tn)*(stats(j).fp+stats(j).tn)*(stats(j).tp+stats(j).fn))^0.5); %matthew correlation coefficient: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
        stats(j).OPTROCPT = OPTROCPT;
        stats(j).FPR = FPR;
        stats(j).TPR = TPR;
        stats(j).T = T;
        Pre = ((stats(j).tp+stats(j).fp)*(stats(j).tp+stats(j).fn) + (stats(j).tn+stats(j).fn)*(stats(j).tn+stats(j).fp)) ...
            / (stats(j).tp+stats(j).tn+stats(j).fp+stats(j).fn)^2;
        stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);
        stats(j).Fscore = 2*stats(j).tp/(2*stats(j).tp+stats(j).fp+stats(j).fn);
        stats(j).subsets.training = tra;
        stats(j).subsets.testing = tes;
    end
end
            
fprintf('\nDONE.\n');
