function [CV_Results, meanValues] = analyze_CV_results(output_path, feature_ref, NumFeats)
%ANALYZE_CV_RESULTS Finds top feature information
%   This function takes the output of K-Fold cross validation (CV) as input.
%   featType dictates which feature K-Fold cross validation result is
%   loaded. NumFeats dictates how many top features are selected. This
%   function returns a struct that contains the top feature information.
%   Parameters:
%       output_path: str
%           Directory to save the following files:
%               - "CV_Results.xlsx" - quantitative metrics for CV
%               - "Top_Features_Information.mat" - information about the top N
%                  features saved as .mat file
%               - "Top_Features_Information.xlsx" - information about the
%                  top N features saved as .xlsx file
%       feature_ref: cell array of strings
%           1xM cell array of strings, where M is the number of features.
%           Each cell contains the name of a feature as a string.
%       NumFeats: int
%           Number of top features to select
%   Returns:
%       CV_Results: struct
%           Top feature results. These results include:
%               - Indices
%               - Names
%               - Counts
%               - Frequency
%               - p-values
%       meanValues: struct
%           The mean metrics for K-Fold CV. This struct includes:
%               - AUC
%               - Accuracy
%               - Positive Predictive Value (ppv)
%               - Sensitivity
%               - Specificity
%               - Kappa
%               - F-Score
%               - MCC
%               - Optimal Threshold

%% Load in stats files
feat_stats = load(strcat(output_path,'CV_Feature_Results.mat'));
feat_stats = feat_stats.feat_stats;

%% Calculate stats
total_itrs = length(feat_stats);
meanValues = getMeanFromStructField(feat_stats, total_itrs);

%% Save results to excel file
writetable(meanValues,strcat(output_path,'CV_Training_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Metrics');

%% Load in top feature indices
top_feats = load(strcat(output_path,'Top_Features.mat'));
top_feat_indices = top_feats.ordered_feat_results.featinds(1:NumFeats)'; % Get top feature indices
topResults.indices = top_feat_indices;
% topResults.names = feature_ref(top_feat_indices)'; % Get top  feature names
topResults.counts = top_feats.ordered_feat_results.count(1:NumFeats)'; % Get top  feature counts
topResults.freq = top_feats.ordered_feat_results.freq(1:NumFeats)'; % Get top feature frequencies
topResults.pvals = top_feats.ordered_feat_results.pvals(1:NumFeats)'; % Get top feature p-vals

CV_Results = topResults;

%% Save top feature information
save(strcat(output_path,'Top_Features_Information.mat'),'topResults');
topResults = struct2table(topResults);
writetable(topResults,strcat(output_path,'Top_Features_Information.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Top Features');



function meanValues = getMeanFromStructField(structVar, iters)
AUC_cell = cell(iters,1);
acc_cell = cell(iters,1);
ppv_cell = cell(iters,1);
sens_cell = cell(iters,1);
spec_cell = cell(iters,1);
kappa_cell = cell(iters,1);
fscore_cell = cell(iters,1);
mcc_cell = cell(iters,1);
thresh_cell = cell(iters,1);
for i=1:iters
    AUC_cell{i} = structVar(i).AUCs(1);
    acc_cell{i} = structVar(i).acc;
    ppv_cell{i} = structVar(i).ppv;
    sens_cell{i} = structVar(i).sens;
    spec_cell{i} = structVar(i).spec;
    kappa_cell{i} = structVar(i).kappa;
    fscore_cell{i} = structVar(i).Fscore;
    mcc_cell{i} = structVar(i).MCC;
    thresh_cell{i} = structVar(i).optThresh;
end
results.AUC = round(mean(cell2mat(AUC_cell)), 3);
results.AUC_std = round(std(cell2mat(AUC_cell)), 3);

results.Accuracy = round(mean(cell2mat(acc_cell)), 3);
results.Accuracy_std = round(std(cell2mat(acc_cell)), 3);

results.ppv_value = round(mean(cell2mat(ppv_cell)), 3);
results.ppv_value_std = round(std(cell2mat(ppv_cell)), 3);

results.sens = round(mean(cell2mat(sens_cell)), 3);
results.sens_std = round(std(cell2mat(sens_cell)), 3);

results.specif = round(mean(cell2mat(spec_cell)), 3);
results.specif_std = round(std(cell2mat(spec_cell)), 3);

results.kappa_value = round(mean(cell2mat(kappa_cell)), 3);
results.kappa_value_std = round(std(cell2mat(kappa_cell)), 3);

results.fscore = round(mean(cell2mat(fscore_cell)), 3);
results.fscore_std = round(std(cell2mat(fscore_cell)), 3);

results.MCC_value = round(mean(cell2mat(mcc_cell)), 3);
results.MCC_value_std = round(std(cell2mat(mcc_cell)), 3);

results.Opt_Thres = mean(cell2mat(thresh_cell));

results = struct2table(results);
meanValues = results;
end
end

