function [stats]= nFoldCV_withFS_v3(data_set,data_labels,params)

% Using k runs of n-fold cross validation embedded with a fresh feature selection each fold of all k runs
%***PERFORMS PERFCURVE ACROSS ALL FOLDS SEPARATELY***%
%
% INPUTS:
%      data_set: M x N matrix of m observations and n variables (e.g. observations = patients, variables = features)
%   data_labels: M x 1 column of class labels corresponding to each observation (1 = positive class, -1 = negative class)
%        params: (optional) struct containing the fields:
%                   params.classifier: -- options: 'LDA','QDA','SVM','RANDOMFOREST' (Default: 'LDA');
%                   params.classifieroptions: (optional) Additional {...,'Name','Value'...'} parameters for your classifier. Pass in as a cell array of strings
%                   params.fsname -- options: 'mrmr','ttest','wilcoxon','entropy','bhattacharyya','roc') (Default: 'wilcoxon')
%                   params.num_top_feats: -- the number of top features to select each iteration
%                   params.feature_idxs: (optional) vector of pre-selected set of feature indices (i.e. from pruning or anti-correlation) (Default: include all variables in data_set)
%                   params.shuffle: 1 for random, 0 for non-random partition (Default: 1)
%                   params.n: Number of folds to your cross-validation (Default: 3)
%                   params.nIter: Number of cross-validation iterations (Default: 25)
%                   params.subsets: (optional) training and testing subsets for each iteration (Default: computer will generate using 'nFold')
%                   params.patient_ids: (needed for augmented data, i.e. multiple slices per patient) -- M x 1 vector of unique patient identifier (an integer) from which each data_set and data_label came from 
%                   params.osname: (optional) string for oversampling method) -- options: 'SMOTE','ADASYN','UNDERSAMPLE'
%                   params.remouts: (optional) string for removing outlier method -- options: 'median','mean','quartile'
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
% 05-02-2018 - Modified to include option for oversampling data
% 05-22-2018 - Modified to include option for removing outliers from training data
%            - Can handle 0s and 1s as well as -1s and 1s class labels
% 07-03-2018 - Small fixes to handle empty fields in params
% 07-18-2018 - Added in stats for LOOCV, can now pass in classifier options, automatically adds all necessary subfolders into path.

% ============ example call:====================
% load fisheriris
% 
% addpath('G:\Team Drives\INVent\jacob\Code\classification\nFold_cross_validation')
%
% data_set = meas;
% data_labels = (contains(species,'virginica'));
% params.classifier='LDA';
% params.fsname='wilcoxon';
% params.shuffle = 1;
% params.n = 3;
% params.nIter = 25;
% params.num_top_feats = 1;
% stats = nFoldCV_withFS_v3(data_set,data_labels,params);
%===============================================

%% Check for valid inputs
if isa(data_labels,'cell'), error('DATA_LABELS must be a non-cell array with labels 1 and -1.'); end
if all([any(~xor(data_labels == 1, data_labels == -1)), any(~xor(data_labels == 1, data_labels == 0))]), error('DATA_LABELS must be 1 and -1 (or 1 and 0)'); end
if size(data_set,1)~=length(data_labels), error('DATA_LABELS must contain same number of labels as observations in DATA_SET (equal # of rows)'); end
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
    params.classifier = 'LDA';
    params.classifieroptions = {};
    params.fsname = 'wilcoxon';
%     params.feature_list = [];
    params.shuffle = 1;
    params.n = 3;
    params.nIter = 25;
    params.num_top_feats = 1;
    params.feature_idxs = 1:size(data_set,2);
end

% if you forgot certain parameters. i'll add them in for you
if nargin == 3
    if ~isfield(params,'classifier')    || isempty(params.classifier),      params.classifier = 'LDA';                      end
    if ~isfield(params,'classifieroptions'),                                params.classifieroptions = {};                  end
    if ~isfield(params,'fsname')        || isempty(params.fsname),          params.fsname = 'wilcoxon';                     end
    if ~isfield(params,'shuffle')       || isempty(params.shuffle),         params.shuffle = 1;                             end
    if ~isfield(params,'n')             || isempty(params.n),               params.n = 3;                                   end
    if ~isfield(params,'nIter')         || isempty(params.nIter),           params.nIter = 25;                              end
    if ~isfield(params,'num_top_feats') || isempty(params.num_top_feats),   params.num_top_feats = 1;                       end
    if ~isfield(params,'feature_idxs')  || isempty(params.feature_idxs),    params.feature_idxs = 1:size(data_set,2);       end       
end


%% nFold CV with FS
stats = struct; %in the end: contains (params.nIter)-rows of information ;

data_set = data_set(:,params.feature_idxs);
data_labels = double(data_labels);

fprintf('\nExecuting %i run(s) of %i-fold cross-validation using %s with updated feature selection on %i observations\n',params.nIter,params.n,params.classifier,length(data_labels))
if isfield(params,'osname') && ~isempty(params.osname)
    fprintf('\t Oversampling of training data to be performed using %s method\n',params.osname);
end
if isfield(params,'remouts') && ~isempty(params.remouts)
    fprintf('\t Training data will have outliers removed using %s method\n',params.osname);
    warning('Not recommended. Classifier will probably fail due to NaNs in data');
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
        if isfield(params,'osname') && ~isempty(params.osname)
            [training_set, training_labels] = Oversample(params.osname,training_set,training_labels);
        end
        
        %remove outliers?
        if isfield(params,'remouts') && ~isempty(params.remouts)
            TF = isoutlier(training_set,params.remouts);
            training_set(TF) = nan;
        end
        
  % *********STEP 1: FEATURE SELECTION*************%
        if strcmp(params.fsname, 'mrmr')
            % descritizing to help mrmr -- will fail if data hasn't been normalized!!
            if ~any(abs(mean(training_set))<.2) || ~any(abs(std(training_set))<10)
                warning('Data not normalized to mean=0, std=1. Discretiztion for mrmr will fail.');
            end
            training_set_d = discretize(training_set, [-inf -1 1 inf]);
            fea = mrmr_mid_d(training_set_d,training_labels,params.num_top_feats); % MRMR
            fea = fea'; 
        else
            [IDX, ~] = rankfeatures(training_set',training_labels,'criterion',params.fsname,'CCWeighting', 1); %%CCWeighting is the parameter that weights by correlation w/ previously selected features
%                 [IDX,~] = myrankfeatures(training_set,training_labels,params.num_top_feats);
            
            if length(find(IDX==1))>1 && size(data_set,2)>20 %occurs if CCWeighting is too high and not many features to choose from or if too few (<20) features in data_set
                warning('Too many correlated features found. Only %i features chosen. Reduce CCWeighting in call to rankfeatures() to correct.',length(IDX)-length(find(IDX==1))); 
                fea = IDX(1:length(IDX)-length(find(IDX==1)));
            else
                fea = IDX(1:params.num_top_feats);
            end    
        end
        
        % *********STEP 2: CLASSIFICATION****************%
        options = params.classifieroptions;
        [temp_stats,~] = Classify_wrapper(params.classifier, training_set(:,fea) , testing_set(:,fea), training_labels(:), testing_labels(:), options);

        feats = [feats; fea']; %may throw error later b/c the issue with CCWeighting and not full set of features being selected... will handle later 

        if numel(unique(testing_labels(:))) > 1 %won't work for Leave One Out CV

            Ttp = Ttp + temp_stats.tp;
            Ttn = Ttn + temp_stats.tn;
            Tfp = Tfp + temp_stats.fp;
            Tfn = Tfn + temp_stats.fn;
            
            if isfield(params,'patient_ids')
                decisions(ismember(params.patient_ids,pts(tes{i}))) = temp_stats.decision;
            else
                decisions(tes{i}) = temp_stats.decision;
            end
            
            if any(ismember(data_labels,-1))
                decisions(decisions==0) = -1;
            end
        end
        
        if isfield(params,'patient_ids')
            predictions(ismember(params.patient_ids,pts(tes{i}))) = temp_stats.prediction;
            labels(ismember(params.patient_ids,pts(tes{i}))) = testing_labels;
        else
            predictions(tes{i}) = temp_stats.prediction;
            labels(tes{i}) = testing_labels;
        end

        %for generating upper and lower AUC bounds
        foldpredictions{i} =  temp_stats.prediction;
        foldlabels{i} = testing_labels;
        
        %***for testing: call this script at this point. Otherwise, leave commented***%
%         temp_analysis_nFoldCV_withFS;

    end
    
    % output statistics  
    if numel(unique(testing_labels)) > 1 %won't work for Leave One Out CV
        
        %compute perfcurve across each testing fold per iteration
        [X,Y,T,AUC,OPTROCPT] = perfcurve(foldlabels, foldpredictions, 1,'XVals', [0:0.02:1]);

        stats(j).AUCs = AUC;
        stats(j).optThresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2))); %MATLAB'S OLD THRESHOLD COMPUTATION
        stats(j).OPTROCPT = OPTROCPT;
        stats(j).predictions = predictions;
        stats(j).decisions = decisions;
        stats(j).labels = labels;
        stats(j).tp = Ttp;
        stats(j).tn = Ttn;
        stats(j).fp  = Tfp;
        stats(j).fn = Tfn;
        stats(j).acc = (Ttp+Ttn)/(Ttp+Ttn+Tfp+Tfn);
        stats(j).ppv = Ttp/(Ttp+Tfp);
        stats(j).sens = Ttp/(Ttp+Tfn);
        stats(j).spec = Ttn/(Tfp+Ttn);
        stats(j).X = X;
        stats(j).Y = Y;
        stats(j).T = T;
        Pre = ((Ttp+Tfp)*(Ttp+Tfn) + (Ttn+Tfn)*(Ttn+Tfp)) / (Ttp+Ttn+Tfp+Tfn)^2;
        stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);
        stats(j).Fscore = 2*Ttp/(2*Ttp+Tfp+Tfn);
        stats(j).subsets.training = tra;
        stats(j).subsets.testing = tes;
        stats(j).topfeatinds = params.feature_idxs(feats);
    else
         %compute perfcurve across each testing fold per iteration        
        try
            [X,Y,T,AUC,OPTROCPT] = perfcurve(foldlabels, foldpredictions, 1, 'XVals', [0:0.02:1]);
        catch %unresolved error
            if iscolumn(foldlabels{1})
                l = cell2mat(foldlabels');
                p = cell2mat(foldpredictions');
            else
                l = cell2mat(foldlabels);
                p = cell2mat(foldpredictions);
            end
            [X,Y,T,AUC,OPTROCPT] = perfcurve(l(:), p(:), 1,'XVals', [0:0.02:1]);
        end
%     [TP FN] = perfcurve(labels,predictions,1,'xCrit','TP','yCrit','FN');
%     [FP TN] = perfcurve(labels,predictions,1,'xCrit','FP','yCrit','TN');
%     [~,ACC] = perfcurve(labels,predictions,1,'xCrit','TP','yCrit','accu');
%     [~,PPV] = perfcurve(labels,predictions,1,'xCrit','TP','yCrit','PPV');
%     
%     optim_idx = find(FPR == OPTROCPT(1) & TPR == OPTROCPT(2));
%     stats.tp = TP(optim_idx);
%     stats.fn = FN(optim_idx);
%     stats.fp = FP(optim_idx);
%     stats.tn = TN(optim_idx);
%     stats.auc = AUC;
%     stats.spec = 1-FPR(optim_idx);
%     stats.sens = TPR(optim_idx);
%     stats.acc = ACC(optim_idx);
%     stats.ppv = PPV(optim_idx);
%     stats.threshold = T(optim_idx);
%     stats.decision = stats.prediction >= stats.threshold;
        stats(j).AUCs = AUC;
        stats(j).X = X;
        stats(j).Y = Y;
        stats(j).T = T;
        stats(j).OPTROCPT = OPTROCPT;    
        stats(j).predictions = predictions;
        stats(j).labels = labels;
        stats(j).subsets.training = tra;
        stats(j).subsets.testing = tes;
        stats(j).topfeatinds = params.feature_idxs(feats);
        try
            stats(j).optThresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));
            
            guessLabels = (predictions>=stats(j).optThresh);
            trueLabels = labels;
                        
            truePositiveVector = guessLabels==1 & trueLabels==1;
            Ttp = (sum(truePositiveVector(:)==1));
            
            trueNegativeVector = guessLabels==0 & trueLabels==0;
            Ttn = (sum(trueNegativeVector(:)==1));
            
            falsePositiveVector = guessLabels==1 & trueLabels==0;
            Tfp = (sum(falsePositiveVector(:)==1));
            
            falseNegativeVector = guessLabels==0 & trueLabels==1;
            Tfn = (sum(falseNegativeVector(:)==1));
            
            stats(j).tp = Ttp; 
            stats(j).tn = Ttn; 
            stats(j).fp  = Tfp; 
            stats(j).fn = Tfn; 
            stats(j).acc = (Ttp+Ttn)/(Ttp+Ttn+Tfp+Tfn);
            stats(j).ppv = Ttp/(Ttp+Tfp);
            stats(j).sens = Ttp/(Ttp+Tfn);
            stats(j).spec = Ttn/(Tfp+Ttn);
            Pre = ((Ttp+Tfp)*(Ttp+Tfn) + (Ttn+Tfn)*(Ttn+Tfp)) / (Ttp+Ttn+Tfp+Tfn)^2;
            stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);        
            stats(j).Fscore = 2*Ttp/(2*Ttp+Tfp+Tfn);
        catch
        end
    end
end
            
fprintf('\nDONE.\n');