clear;close all;clc;

if isunix
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/general'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/nFold_cross_validation/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/training_and_testing_sets/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/classifiers/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/Scripts-DataPreProcessing/Util-preProcessing'));
    % addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/estpab.mexmac'));
end
region = {'Best_Fat'}; % Can be 'best_fat' or 'best_tumor'
split = {'missingCollage'};
scheme = {'mrmr_svm'};

matrix_root_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollage/';

axial_train_matrix_path = string(strcat(matrix_root_path, 'Axial_', region, '_Training.mat'));
axial_test1_matrix_path = string(strcat(matrix_root_path, 'Axial_', region, '_Testing1.mat'));
axial_test2_matrix_path = string(strcat(matrix_root_path, 'Axial_', region, '_Testing2.mat'));

coronal_train_matrix_path = string(strcat(matrix_root_path, 'Coronal_', region, '_Training.mat'));
coronal_test1_matrix_path = string(strcat(matrix_root_path, 'Coronal_', region, '_Testing1.mat'));
coronal_test2_matrix_path = string(strcat(matrix_root_path, 'Coronal_', region, '_Testing2.mat'));

axial_path_to_train = string(axial_train_matrix_path);
axial_path_to_test1 = string(axial_test1_matrix_path);
axial_path_to_test2 = string(axial_test2_matrix_path);

coronal_path_to_train = string(coronal_train_matrix_path);
coronal_path_to_test1 = string(coronal_test1_matrix_path);
coronal_path_to_test2 = string(coronal_test2_matrix_path);

label_path_root = '/Users/leobao/Documents/MultiPlanePipeline/Data/train_test_labels/MissingCollageLabels/';

training_label_path = string(strcat(label_path_root, region, '/train_labels.csv'));
testing1_label_path = string(strcat(label_path_root, region, '/test1_labels.csv'));
testing2_label_path = string(strcat(label_path_root, region, '/test2_labels.csv'));

feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

if strcmp(scheme, 'wilcoxon_qda')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_16_46_13_Axial_proxfat10_only_wilcoxon_qda/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_08_32_Coronal_proxfat10_only_wilcoxon_qda/Top_Features_Information.mat';
elseif strcmp(scheme, 'wilcoxon_lda')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_16_46_24_Axial_proxfat10_only_wilcoxon_lda/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_09_08_Coronal_proxfat10_only_wilcoxon_lda/Top_Features_Information.mat';
elseif strcmp(scheme, 'wilcoxon_rf')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_16_46_35_Axial_proxfat10_only_wilcoxon_rf/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_09_19_Coronal_proxfat10_only_wilcoxon_rf/Top_Features_Information.mat';
elseif strcmp(scheme, 'wilcoxon_svm')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_16_46_59_Axial_proxfat10_only_wilcoxon_svm/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_09_49_Coronal_proxfat10_only_wilcoxon_svm/Top_Features_Information.mat';
elseif strcmp(scheme, 'mrmr_qda')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_17_09_58_Axial_proxfat10_only_mrmr_qda/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_25_36_Coronal_proxfat10_only_mrmr_qda/Top_Features_Information.mat';
elseif strcmp(scheme, 'mrmr_lda')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_17_10_27_Axial_proxfat10_only_mrmr_lda/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_25_46_Coronal_proxfat10_only_mrmr_lda/Top_Features_Information.mat';
elseif strcmp(scheme, 'mrmr_rf')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_17_10_47_Axial_proxfat10_only_mrmr_rf/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_26_00_Coronal_proxfat10_only_mrmr_rf/Top_Features_Information.mat';
elseif strcmp(scheme, 'mrmr_svm')
    axial_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_17_11_13_Axial_proxfat10_only_mrmr_svm/Top_Features_Information.mat';
    coronal_top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_19_26_29_Coronal_proxfat10_only_mrmr_svm/Top_Features_Information.mat';
end

axial_top_feats = load(axial_top_feature_path).topResults;
coronal_top_feats = load(coronal_top_feature_path).topResults;

axial_top_indices = axial_top_feats.indices';
coronal_top_indices = coronal_top_feats.indices';

fprintf("Got directory paths to feature matrices! \n");

output_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/";
experiment_date = strrep(string(datetime("now")), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");
output_path = strcat(output_path, experiment_date, "_", region, "_", scheme, "/");
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

datasets = struct;

datasets.(strcat("axial_train_", region)) = load(axial_path_to_train);
datasets.(strcat("axial_test1_", region)) = load(axial_path_to_test1);
datasets.(strcat("axial_test2_", region)) = load(axial_path_to_test2);
datasets.(strcat("coronal_train_", region)) = load(coronal_path_to_train);
datasets.(strcat("coronal_test1_", region)) = load(coronal_path_to_test1);
datasets.(strcat("coronal_test2_", region)) = load(coronal_path_to_test2);

data_labels_training = readmatrix(training_label_path);
data_labels_holdout1 = readmatrix(testing1_label_path);
data_labels_holdout2 = readmatrix(testing2_label_path);

axial_training_features = datasets.axial_train_Best_Fat.feature_matrix;
axial_testing1_features = datasets.axial_test1_Best_Fat.feature_matrix;
axial_testing2_features = datasets.axial_test2_Best_Fat.feature_matrix;
coronal_training_features = datasets.coronal_train_Best_Fat.feature_matrix;
coronal_testing1_features = datasets.coronal_test1_Best_Fat.feature_matrix;
coronal_testing2_features = datasets.coronal_test2_Best_Fat.feature_matrix;

axial_top5_training_feats = zeros(size(axial_training_features, 1), 5);
axial_top5_testing1_feats = zeros(size(axial_testing1_features, 1), 5);
axial_top5_testing2_feats = zeros(size(axial_testing2_features, 1), 5);
coronal_top5_training_feats = zeros(size(coronal_training_features, 1), 5);
coronal_top5_testing1_feats = zeros(size(coronal_testing1_features, 1), 5);
coronal_top5_testing2_feats = zeros(size(coronal_testing2_features, 1), 5);

axial_top5_feature_names = cell(1, 5);
coronal_top5_feature_names = cell(1, 5);

for i = 1:5
    axial_index = axial_top_indices(i);
    coronal_index = coronal_top_indices(i);

    axial_top5_training_feats(:, i) = axial_training_features(:, axial_index);
    axial_top5_testing1_feats(:, i) = axial_testing1_features(:, axial_index);
    axial_top5_testing2_feats(:, i) = axial_testing2_features(:, axial_index);
    axial_top5_feature_names(:, i) = feature_column_names(:, axial_index);

    coronal_top5_training_feats(:, i) = coronal_training_features(:, coronal_index);
    coronal_top5_testing1_feats(:, i) = coronal_testing1_features(:, coronal_index);
    coronal_top5_testing2_feats(:, i) = coronal_testing2_features(:, coronal_index);
    coronal_top5_feature_names(:, i) = feature_column_names(:, coronal_index);
end

feature_column_names = [axial_top5_feature_names, coronal_top5_feature_names];
features_training = [axial_top5_training_feats, coronal_top5_training_feats];
features_holdout1 = [axial_top5_testing1_feats, coronal_top5_testing1_feats];
features_holdout2 = [axial_top5_testing2_feats, coronal_top5_testing2_feats];

feature_file_name = string(strcat(output_path, 'top10_features.xlsx'));
feature_table = cell2table(feature_column_names);
writetable(feature_table, feature_file_name);

holdout_test_size1 = size(features_holdout1);
holdout_test_size2 = size(features_holdout2);
assert(holdout_test_size1(1) == 14, "The size of the dataset is incorrect!");
assert(holdout_test_size2(1) == 27, "The size of the dataset is incorrect!");

fprintf("Loaded in train and holdout testing datasets! \n");
%% Whiten the data
features_training = whitenData(features_training, 'simple');
features_holdout1 = whitenData(features_holdout1, 'simple');
features_holdout2 = whitenData(features_holdout2, 'simple');

fprintf("Whitened the training and holdout testing datasets successfully! \n");

%% Initialize classifier params for classifier
saveParams = true;
if(isequal(saveParams,true))
    if strcmp(scheme, 'wilcoxon_qda')
        params.classifier='QDA';
        params.fsname='wilcoxon';
    elseif strcmp(scheme, 'wilcoxon_lda')
        params.classifier='LDA';
        params.fsname='wilcoxon';
    elseif strcmp(scheme, 'wilcoxon_rf')
        params.classifier='RANDOMFOREST';
        params.fsname='wilcoxon';
    elseif strcmp(scheme, 'wilcoxon_svm')
        params.classifier='SVM';
        params.fsname='wilcoxon';
    elseif strcmp(scheme, 'mrmr_qda')
        params.classifier='QDA';
        params.fsname='mrmr';
    elseif strcmp(scheme, 'mrmr_lda')
        params.classifier='LDA';
        params.fsname='mrmr';
    elseif strcmp(scheme, 'mrmr_rf')
        params.classifier='RANDOMFOREST';
        params.fsname='mrmr';
    elseif strcmp(scheme, 'mrmr_svm')
        params.classifier='SVM';
        params.fsname='mrmr';
    end
    params.shuffle = 1;
    params.n = 5;
    params.nIter = 50;
    params.num_top_feats = 5;
    params.threshmeth = 'euclidean';
    params.osname = 'SMOTE';
    params.fsprunecorr = 'false';
    params.featnames = feature_column_names;
    save(strcat(output_path, 'Classifier_params.mat'),'params');
else
    params = load('Classifier_params.mat');
    params = params.params;
end

%% Create experimental log
log_file_name = strcat(output_path, "Summary_of_Experiment_", experiment_date, ".txt");
fileid = fopen(log_file_name,'w');
fprintf(fileid, strcat("Experiment Date and Time: ", string(datetime("now")), '\n\r'));
fprintf(fileid, "FEATURE FAMILY: 3D Shape \n");
fprintf(fileid, strcat("Trained on ", int2str(length(features_training(:,1))), " patients \n\r"));
fprintf(fileid, strcat("Tested on ", int2str(length(features_holdout1(:,1)))," patients \n\r"));
fprintf(fileid, strcat("Selection Method: ", params.fsname, '\n\r'));
fprintf(fileid, strcat("Classifier: ", params.classifier, '\n\r'));
fprintf(fileid, strcat("Number of Top Features: ", num2str(params.num_top_feats), '\n\r'));
fprintf(fileid, strcat("Number of Folds: ", num2str(params.n), '\n\r'));
fclose(fileid);

%% Run K-fold Cross Validation
fprintf("Running K-Fold Cross Validation with Feature Selection... \n");
feat_stats = nFoldCV_withFS(features_training,data_labels_training,params);
save(strcat(output_path,'CV_Feature_Results.mat'),'feat_stats');
FSResults = getCVTopFeatures(features_training, params, feat_stats, data_labels_training, output_path);

%% Analyze the output of feature selection
[topResults,meanValues] = analyze_CV_results(output_path, feature_column_names, params.num_top_feats);
fprintf("Found the top ranked features! \n")

%% Set the optimal threshold from analyze_CV_results
classifier_thresh = meanValues.Opt_Thres;

%% Train and test the QDA classifier
train_feats = features_training(:,topResults.indices);
test1_feats = features_holdout1(:,topResults.indices);
test2_feats = features_holdout2(:,topResults.indices);
fprintf(strcat("Training on ",int2str(length(train_feats(1,:)))," features for ", int2str(length(features_training(:,1))), " patients \n"));
fprintf(strcat("Testing on ",int2str(length(test1_feats(1,:)))," features for ", int2str(length(features_holdout1(:,1)))," patients \n"));

[stats1, methodstring] = TrainTestModel(params.classifier, train_feats, test1_feats, data_labels_training, data_labels_holdout1,classifier_thresh);
save(strcat(output_path,'Holdout_Testing1_Dataset_Results.mat'),'stats1');

stats_for_export1 = struct;
stats_for_export1.ACC = round(stats1.acc, 3);
stats_for_export1.AUC = round(stats1.AUC, 3);
stats_for_export1.AUCPRC = round(stats1.AUPRC, 3);
stats_for_export1.SENS = round(stats1.sens, 3);
stats_for_export1.SPEC = round(stats1.spec, 3);
stats_for_export1.FSCORE = round(stats1.Fscore, 3);
stats_for_export1.TP = stats1.tp;
stats_for_export1.FP = stats1.fp;
stats_for_export1.FN = stats1.fn;
stats_for_export1.TN = stats1.tn;
stats_for_export1.KAPPA = round(stats1.kappa, 3);
stats_for_export1.MCC = round(stats1.MCC, 3);
    

QDAResults1 = struct2table(stats_for_export1);
writetable(QDAResults1,strcat(output_path,'Holdout_Testing1_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Entire_Test_Cohort_Results');

fprintf(strcat("Testing on ",int2str(length(test2_feats(1,:)))," features for ", int2str(length(features_holdout2(:,1)))," patients \n"));

[stats2, methodstring] = TrainTestModel(params.classifier, train_feats, test2_feats, data_labels_training, data_labels_holdout2,classifier_thresh);
save(strcat(output_path,'Holdout_Testing2_Dataset_Results.mat'),'stats2');

stats_for_export2 = struct;
stats_for_export2.ACC = round(stats2.acc, 3);
stats_for_export2.AUC = round(stats2.AUC, 3);
stats_for_export2.AUCPRC = round(stats2.AUPRC, 3);
stats_for_export2.SENS = round(stats2.sens, 3);
stats_for_export2.SPEC = round(stats2.spec, 3);
stats_for_export2.FSCORE = round(stats2.Fscore, 3);
stats_for_export2.TP = stats2.tp;
stats_for_export2.FP = stats2.fp;
stats_for_export2.FN = stats2.fn;
stats_for_export2.TN = stats2.tn;
stats_for_export2.KAPPA = round(stats2.kappa, 3);
stats_for_export2.MCC = round(stats2.MCC, 3);
    

QDAResults2 = struct2table(stats_for_export2);
writetable(QDAResults2,strcat(output_path,'Holdout_Testing2_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Entire_Test_Cohort_Results');
%% Separate and save UH and VA results

function results = compute_metrics(prediction, groundtruth)

    assert(size(prediction, 1) == size(groundtruth, 1), "Input prediction array and ground truth array are not the same size in the first axis!");
    assert(size(prediction, 2) == size(groundtruth, 2), "Input prediction array and ground truth array are not the same size in teh second axis!");
    
    tp = length(find(prediction(prediction == groundtruth) == 1));
    tn = length(find(prediction(prediction == groundtruth) == 0));
    fn = length(find(prediction(prediction ~= groundtruth) == 0));
    fp = length(find(prediction(prediction ~= groundtruth) == 1));
    
    [FPR,TPR,T,AUC,OPTROCPT,~,~] = perfcurve(groundtruth, prediction, 1, 'XVals', [0:0.02:1]);  % calculate AUC. 'perfcurve' can also calculate sens, spec etc. to plot the ROC curve.
    [RECALL,PRECISION,~,AUPRC,~,~,~] = perfcurve(groundtruth, prediction, 1, 'XVals', [0:0.02:1], 'xCrit', 'reca', 'yCrit', 'prec'); 
    
    ACC = (tp+tn)/(tp+tn+fp+fn);
    results.ACC = round(ACC, 3);
    
    results.AUC = round(AUC, 3);
    results.AUCPRC = round(AUPRC, 3);
    
    PPV = tp/(tp+fp);
    
    SENS = tp/(tp+fn);
    results.SENS = round(SENS, 3);
    
    SPEC = tn/(fp+tn);
    results.SPEC = round(SPEC, 3);
    
    FSCORE = 2*tp/(2*tp+fp+fn);
    results.FSCORE = round(FSCORE, 3);
    
    results.TP = tp;
    results.FP = fp;
    results.FN = fn;
    results.TN = tn;
    
    Pre = ((tp+fp)*(tp+fn) + (tn+fn)*(tn+fp)) / (tp+tn+fp+fn)^2;
    KAPPA = (results.ACC - Pre) / (1 - Pre);
    results.KAPPA = round(KAPPA, 3);
    
    MCC = (tp*tn - fp*fn) / (((tp+fp)*(fn+tn)*(fp+tn)*(tp+fn))^0.5);
    results.MCC = round(MCC, 3);

end