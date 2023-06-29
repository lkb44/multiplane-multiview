function stats = nFold_SVM(data_set,data_labels,shuffle,n,nIter,Subsets,SVM_kernel)
warning('may be outdated. check code and results');


% Using n-fold subset selection and SVM classifier
% Input:
%   data_set: data
%   data_labels: labels
%   classifier: string containing the classifier name (ie: 'SVM')
%   shuffle: 1 for random, 0 for non-random partition (Default: 0)
%   n: Number of folds to your cross-validation (Default: 3)
%   SVM_kernel: 'linear', 'quadratic','mlp','rbf', 'poly' (Default: 'linear')
%   nIter: Number of cross-validation iterations (Default: 1)
%   Subsets: pass your own training and testing subsets & labels (Default:
%   computer will generate using 'nFold')
%
% Output:
%   stats: struct containing true positive, accuracy
%
% The function was originally written by George Lee @2009
% Function has been last modified by Jacob Antunes @2-15-2018
%
% addpath('\\csehomes\csehome$\jta35\My Documents\MATLAB\Research\Madabhushi\REPO\repo_george\ClassificationPackage-2009\classifiers\svm_mex\svmmex_tools');
% addpath('\\csehomes\csehome$\jta35\My Documents\MATLAB\Research\Madabhushi\REPO\repo_george\ClassificationPackage-2009\classifiers\svm_mex\bin');
if nargin < 7
    SVM_kernel = 'linear'; %linear, rbf, polynomial, quadratic, mlp
end
if nargin < 6
    Subsets = {};
end
if nargin < 5
    nIter = 1;
end
if nargin < 4
    n = 3; % 3-fold cross-validation
end
if nargin < 3
    shuffle = 1; % randomized
end

if any(~xor(data_labels == 1, data_labels == -1)), error('Labels must be 1 and -1'); end
data_set = svm_scale(data_set);
params.kernel = SVM_kernel;

stats = struct; %cell(1,nIter);
for j=1:nIter
    fprintf('Iteration: %i\n',j);
    
    % reset total statistics
    Ttp = 0; Ttn = 0; Tfp = 0; Tfn = 0;
    
    if isempty(Subsets)
        [tra tes]=GenerateSubsets('nFold',data_set,data_labels,shuffle,n);
    else
        tra = Subsets(j).training;
        tes = Subsets(j).testing;
    end
    
    for i=1:n
        fprintf(['Fold # ' num2str(i) '\n']);
        
        training_set = data_set(tra{i},:);
        testing_set = data_set(tes{i},:);
        training_labels = data_labels(tra{i});
        testing_labels = data_labels(tes{i});
        
        [temp_stats,methodstring] = Classify_wrapper( 'SVM', training_set , testing_set, training_labels, testing_labels, params);
        Ttp = Ttp + temp_stats.tp;
        Ttn = Ttn + temp_stats.tn;
        Tfp = Tfp + temp_stats.fp;
        Tfn = Tfn + temp_stats.fn;
        decision(tes{i}) = temp_stats.decision;
        prediction(tes{i}) = temp_stats.prediction;

    end
    
    % output statistics
  
    % output statistics
    if numel(unique(data_labels))>1 %numel(unique(testing_labels))>1
        if n == 1
        [FPR,TPR,T,AUC,OPTROCPT,~,~] = perfcurve(data_labels(tes{i}),prediction(tes{i}),1);
        else
        [FPR,TPR,T,AUC,OPTROCPT,~,~] = perfcurve(data_labels,prediction,1);
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
    stats(j).decision = decision';
    stats(j).prediction = prediction';
    Pre = ((Ttp+Tfp)*(Ttp+Tfn) + (Ttn+Tfn)*(Ttn+Tfp)) / (Ttp+Ttn+Tfp+Tfn)^2;
    stats(j).kappa = (stats(j).acc - Pre) / (1 - Pre);
end