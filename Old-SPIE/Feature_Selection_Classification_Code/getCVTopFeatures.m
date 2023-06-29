function ordered_feat_results = getCVTopFeatures(features, CV_params, feat_stats, data_labels, output_path)
%GETCVTOPFEATURES Determine the top features from K-Fold Cross Validation.
%   After K-Fold Cross Validation, the top features are selected. For each
%   feature, the following information is computed:
%       - ordered_count = the number of times each feature was ranked.
%       - ordered_freq = the proportion of times each feature was selected.
%       - featinds = indices of the top ranked features in the input
%       feature matrix.
%       - pvals = p-values of each ranked features. Indicates whether the
%       selected value is statistically significant
%   The information for each ranked feature is saved into a .mat file.
%   Parameters:
%       features: matrix
%           Matrix of features used for K-Fold Cross Validation with
%           Feature Selection
%       CV_params: struct
%           Struct containing the parameters passed into K-Fold Cross
%           Validation with Feature Selection
%       feat_stats: struct
%           Struct containing the metrics from each iteration of K-Fold
%           Cross Validation with Feature Selection
%       data_labels: vector
%           mx1 vector, where m is the number of patients. Each row
%           corresponds to the data label of a patient.
%       output_path: str
%           Directory in which the information for the top ranked features
%           will be saved.
%   Returns:
%      ordered_feat_results: struct
%           Struct containing the top ranked features. The fields of the
%           struct include:
%                - ordered_count
%                - ordered_freq
%                - featinds
%                - pvals
%
%p-values
minLab = min(unique(data_labels));
maxLab = max(unique(data_labels));
for i = 1:size(features,2)
    if any(isnan(features(:,i)))
        p(i) = 1;
    else
        [p(i),~]=ranksum(features((data_labels==maxLab),i), features((data_labels==minLab),i)); 
    end
end

% count the number of times certain features were ranked. NOTE: only storing features that received at least 1 voted rank.
[ordered_count, ordered_featInds] = count_topfeatures(cat(1,feat_stats.topfeatinds),1:size(features,2),'descend',1);
ordered_freq = ordered_count./(CV_params.n*CV_params.nIter); %proportion of times each feature was selected
%ordered_names = CV_params.featnames(ordered_featInds);
ordered_pvals = p(ordered_featInds);

%store as single struct
ordered_feat_results.count = ordered_count;
ordered_feat_results.freq = ordered_freq;
ordered_feat_results.featinds = ordered_featInds;
%ordered_feat_results.featnames = ordered_names;
ordered_feat_results.pvals = ordered_pvals;

save(strcat(output_path,'Top_Features.mat'),'ordered_feat_results');
end

