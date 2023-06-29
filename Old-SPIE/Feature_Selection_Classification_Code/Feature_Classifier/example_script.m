%% setup
clear;clc;
addpath(genpath('.'));
addpath('../feature_selection');

%% load data
load fisheriris

data_labels = zeros(size(species));
for i = 1:length(species)
    if strcmp(species{i},'versicolor')
        data_labels(i)=1; %data_labels must be -1 and 1;
    else
        data_labels(i)=-1;
    end
end
data = meas;

clear species meas

%% establish parameters
params.classifier='SVM';
params.fsname='wilcoxon';
params.num_top_feats=4;
params.shuffle = 1;
params.n = 3;
params.nIter = 25;
params.num_top_feats = 1;
params.subsets = {};
params.featnames = {'sepal length','spepal width','petal length','petal width'};

%% remove correlated features
num_features = ceil(1*size(data,2)); %what percent of the features
idx = 1:size(data,2); %but search through all available features
correlation_factor = 0.6;
correlation_metric = 'Spearman';
set_candiF=pick_best_uncorrelated_features(data,data_labels, idx, num_features,correlation_factor,correlation_metric);
params.feature_idxs = set_candiF; % pre-selected set of features

clear num_features idx correlation_factor correlation_metric set_candiF

%% evaluate performance of remaining features
stats = nFoldCV_withFS_v3(data,data_labels,params);

%% what were most frequently appearing features

%p-values
minLab = min(unique(data_labels));
maxLab = max(unique(data_labels));
for i = 1:size(data,2)
    if all(isnan(data(:,i)))
        p(i) = 1;
    else
        [p(i),~]=ranksum(data((data_labels==maxLab),i), data((data_labels==minLab),i)); 
    end
end

% count the number of times certain features were ranked. NOTE: only storing features that received at least 1 voted rank.
[ordered_count, ordered_featInds] = count_topfeatures(cat(1,stats.topfeatinds),1:size(data,2),'descend',1);
ordered_freq = ordered_count./(params.n*params.nIter); %proportion of times each feature was selected
ordered_names = params.featnames(ordered_featInds);
ordered_pvals = p(ordered_featInds);

%store as single struct
ordered_feat_results.count = ordered_count;
ordered_feat_results.freq = ordered_freq;
ordered_feat_results.featinds = ordered_featInds;
ordered_feat_results.featnames = ordered_names;
ordered_feat_results.pvals = ordered_pvals;
clear ordered_count ordered_freq ordered_featInds ordered_names ordered_pvals i minLab maxLab p;

%%
%\*THUS WE ARE LEFT WITH THREE VARIABLES OF INTEREST: 
% params
% stats
% ordered_feat_results
%*/
