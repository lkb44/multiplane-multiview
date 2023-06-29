addpath(genpath('.'));

clear;close all;clc;
load fisheriris

data = meas;
labels = zeros(size(species));
for i = 1:length(species)
    if strcmp(species{i},'versicolor')
        labels(i)=1; %data_labels must be -1 and 1;
    else
        labels(i)=-1;
    end
end

type = 'mrmr'; %wilcoxn, mrmr
K = 2; % return indices of top 2 features

[top_feat_inds] = top_feature_selection_ndim(data, labels, type, K)

