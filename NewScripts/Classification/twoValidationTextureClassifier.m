clear;close all;clc;

if isunix
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/general'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/nFold_cross_validation/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/training_and_testing_sets/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/classifiers/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/multiplane-multiview/Old-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
    % addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/estpab.mexmac'));
end
%% Specify dataset information
view = {'Coronal'}; % {Can only be 'Axial'}, {'Coronal'} or {'Axial','Coronal'}
rois = {'Proximal_Fat10'};
scheme = {'wilcoxon_rf'}; % {'wilcoxon_qda'}, {'wilcoxon_rf'}, {'mrmr_qda'}, or {'mrmr_rf'}

if strcmp(rois, "Proximal_Fat5")
    experiment_type = 'proxfat5_only';
    region = 'ProxFat5';
elseif strcmp(rois, "Proximal_Fat10")
    experiment_type = 'proxfat10_only';
    region = 'ProxFat10';
elseif strcmp(rois, "Proximal_Fat15")
    experiment_type = 'proxfat15_only';
    region = 'ProxFat15';
elseif strcmp(rois, "Fat")
    experiment_type = 'fat_only';
    region = 'Fat';
elseif strcmp(rois, "Tumor")
    experiment_type = 'tumor_only';
    region = 'Tumor';
end
% Print experiment types
fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf("EXPERIMENT TYPE: %s \n", experiment_type);
%% Specify paths to ground truth labels
label_path_root = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Labels/';
matrix_path_root = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/';

training_label_path = string(strcat(label_path_root, view, '_', region, '/train_labels.csv'));
testing1_label_path = string(strcat(label_path_root, view, '_', region, '/test1_labels.csv'));
testing2_label_path = string(strcat(label_path_root, view, '_', region, '/test2_labels.csv'));

train_matrix_path = strcat(matrix_path_root, 'TrainingTextureFeatures/', view, '_', region, '_Training.mat');
test1_matrix_path = strcat(matrix_path_root, 'Validation1TextureFeatures/', view, '_', region, '_Validation1.mat');
test2_matrix_path = strcat(matrix_path_root, 'Validation2TextureFeatures/', view, '_', region, '_Validation2.mat');

path_to_train = string(train_matrix_path);
path_to_test1 = string(test1_matrix_path);
path_to_test2 = string(test2_matrix_path);

feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

fprintf("Got directory paths to feature matrices! \n");

%% Create output directory


output_path = "/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Results/";
experiment_date = strrep(string(datetime("now")), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");
output_path = strcat(output_path, experiment_date, "_", view, "_", experiment_type, "_", scheme, "/");
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

%% Load in datasets
datasets = struct;

for i = 1:length(rois)
    datasets.(strcat("train_",rois{i})) = load(path_to_train);
    datasets.(strcat("test1_",rois{i})) = load(path_to_test1);
    datasets.(strcat("test2_",rois{i})) = load(path_to_test2);
end

data_labels_training = readmatrix(training_label_path);
% data_labels_training = data_labels_training(:, 2);

data_labels_holdout1 = readmatrix(testing1_label_path);
data_labels_holdout2 = readmatrix(testing2_label_path);
% data_labels_holdout = data_labels_holdout(:, 2);

if(length(view) == 1 && experiment_type == "tumor_only")
    features_training = datasets.train_Tumor.feature_matrix;
    features_holdout1 = datasets.test1_Tumor.feature_matrix;
    features_holdout2 = datasets.test2_Tumor.feature_matrix;
    feature_column_names = strcat(view{1}, "_tumor_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "fat_only")
    features_training = datasets.train_Fat.feature_matrix;
    features_holdout1 = datasets.test1_Fat.feature_matrix;
    features_holdout2 = datasets.test2_Fat.feature_matrix;
    feature_column_names = strcat(view{1}, "_fat_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "proxfat5_only")
    features_training = datasets.train_Proximal_Fat5.feature_matrix;
    features_holdout1 = datasets.test1_Proximal_Fat5.feature_matrix;
    features_holdout2 = datasets.test2_Proximal_Fat5.feature_matrix;
    feature_column_names = strcat(view{1}, "_proxfat5_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "proxfat10_only")
    features_training = datasets.train_Proximal_Fat10.feature_matrix;
    features_holdout1 = datasets.test1_Proximal_Fat10.feature_matrix;
    features_holdout2 = datasets.test2_Proximal_Fat10.feature_matrix;
    feature_column_names = strcat(view{1}, "_proxfat10_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "proxfat15_only")
    features_training = datasets.train_Proximal_Fat15.feature_matrix;
    features_holdout1 = datasets.test1_Proximal_Fat15.feature_matrix;
    features_holdout2 = datasets.test2_Proximal_Fat15.feature_matrix;
    feature_column_names = strcat(view{1}, "_proxfat15_", feature_column_names);
end

holdout_test_size1 = size(features_holdout1);
holdout_test_size2 = size(features_holdout2);
assert(holdout_test_size1(1) == 27, "The size of the dataset is incorrect!");
assert(holdout_test_size2(1) == 36, "The size of the dataset is incorrect!");

fprintf("Loaded in train and holdout testing datasets! \n");

if strcmp(rois(i), 'Proximal_Fat5')
    nan_vector_train = features_training(1, :);
    nan_indices_train = find(~isnan(nan_vector_train));
    features_training = features_training(:, nan_indices_train);
    features_holdout1 = features_holdout1(:, nan_indices_train);
    features_holdout2 = features_holdout2(:, nan_indices_train);
    feature_column_names = feature_column_names(:, nan_indices_train);
    nan_vector_test1 = features_holdout1(1, :);
    nan_indices_test1 = find(~isnan(nan_vector_test1));
    features_training = features_training(:, nan_indices_test1);
    features_holdout1 = features_holdout1(:, nan_indices_test1);
    features_holdout2 = features_holdout2(:, nan_indices_test1);
    feature_column_names = feature_column_names(:, nan_indices_test1);
    nan_vector_test2 = features_holdout2(1, :);
    nan_indices_test2 = find(~isnan(nan_vector_test2));
    features_training = features_training(:, nan_indices_test2);
    features_holdout1 = features_holdout1(:, nan_indices_test2);
    features_holdout2 = features_holdout2(:, nan_indices_test2);
    feature_column_names = feature_column_names(:, nan_indices_test2);
elseif strcmp(rois(i), 'Proximal_Fat10')
    nan_vector_train = features_training(1, :);
    nan_indices_train = find(~isnan(nan_vector_train));
    features_training = features_training(:, nan_indices_train);
    features_holdout1 = features_holdout1(:, nan_indices_train);
    features_holdout2 = features_holdout2(:, nan_indices_train);
    feature_column_names = feature_column_names(:, nan_indices_train);
    nan_vector_test1 = features_holdout1(1, :);
    nan_indices_test1 = find(~isnan(nan_vector_test1));
    features_training = features_training(:, nan_indices_test1);
    features_holdout1 = features_holdout1(:, nan_indices_test1);
    features_holdout2 = features_holdout2(:, nan_indices_test1);
    feature_column_names = feature_column_names(:, nan_indices_test1);
    nan_vector_test2 = features_holdout2(1, :);
    nan_indices_test2 = find(~isnan(nan_vector_test2));
    features_training = features_training(:, nan_indices_test2);
    features_holdout1 = features_holdout1(:, nan_indices_test2);
    features_holdout2 = features_holdout2(:, nan_indices_test2);
    feature_column_names = feature_column_names(:, nan_indices_test2);
elseif strcmp(rois(i), 'Proximal_Fat15')
    nan_vector_train = features_training(1, :);
    nan_indices_train = find(~isnan(nan_vector_train));
    features_training = features_training(:, nan_indices_train);
    features_holdout1 = features_holdout1(:, nan_indices_train);
    features_holdout2 = features_holdout2(:, nan_indices_train);
    feature_column_names = feature_column_names(:, nan_indices_train);
    nan_vector_test1 = features_holdout1(1, :);
    nan_indices_test1 = find(~isnan(nan_vector_test1));
    features_training = features_training(:, nan_indices_test1);
    features_holdout1 = features_holdout1(:, nan_indices_test1);
    features_holdout2 = features_holdout2(:, nan_indices_test1);
    feature_column_names = feature_column_names(:, nan_indices_test1);
    nan_vector_test2 = features_holdout2(1, :);
    nan_indices_test2 = find(~isnan(nan_vector_test2));
    features_training = features_training(:, nan_indices_test2);
    features_holdout1 = features_holdout1(:, nan_indices_test2);
    features_holdout2 = features_holdout2(:, nan_indices_test2);
    feature_column_names = feature_column_names(:, nan_indices_test2);
elseif strcmp(rois(i), 'Fat')
    nan_vector_train = features_training(1, :);
    nan_indices_train = find(~isnan(nan_vector_train));
    features_training = features_training(:, nan_indices_train);
    features_holdout1 = features_holdout1(:, nan_indices_train);
    features_holdout2 = features_holdout2(:, nan_indices_train);
    feature_column_names = feature_column_names(:, nan_indices_train);
    nan_vector_test1 = features_holdout1(1, :);
    nan_indices_test1 = find(~isnan(nan_vector_test1));
    features_training = features_training(:, nan_indices_test1);
    features_holdout1 = features_holdout1(:, nan_indices_test1);
    features_holdout2 = features_holdout2(:, nan_indices_test1);
    feature_column_names = feature_column_names(:, nan_indices_test1);
    nan_vector_test2 = features_holdout2(1, :);
    nan_indices_test2 = find(~isnan(nan_vector_test2));
    features_training = features_training(:, nan_indices_test2);
    features_holdout1 = features_holdout1(:, nan_indices_test2);
    features_holdout2 = features_holdout2(:, nan_indices_test2);
    feature_column_names = feature_column_names(:, nan_indices_test2);
elseif strcmp(rois(i), 'Tumor')
    nan_vector_train = features_training(1, :);
    nan_indices_train = find(~isnan(nan_vector_train));
    features_training = features_training(:, nan_indices_train);
    features_holdout1 = features_holdout1(:, nan_indices_train);
    features_holdout2 = features_holdout2(:, nan_indices_train);
    feature_column_names = feature_column_names(:, nan_indices_train);
    nan_vector_test1 = features_holdout1(1, :);
    nan_indices_test1 = find(~isnan(nan_vector_test1));
    features_training = features_training(:, nan_indices_test1);
    features_holdout1 = features_holdout1(:, nan_indices_test1);
    features_holdout2 = features_holdout2(:, nan_indices_test1);
    feature_column_names = feature_column_names(:, nan_indices_test1);
    nan_vector_test2 = features_holdout2(1, :);
    nan_indices_test2 = find(~isnan(nan_vector_test2));
    features_training = features_training(:, nan_indices_test2);
    features_holdout1 = features_holdout1(:, nan_indices_test2);
    features_holdout2 = features_holdout2(:, nan_indices_test2);
    feature_column_names = feature_column_names(:, nan_indices_test2);
end

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
    elseif strcmp(scheme, 'mrmr_qda')
        params.classifier='QDA';
        params.fsname='mrmr';
    elseif strcmp(scheme, 'mrmr_lda')
        params.classifier='LDA';
        params.fsname='mrmr';
    elseif strcmp(scheme, 'mrmr_rf')
        params.classifier='RANDOMFOREST';
        params.fsname='mrmr';
    end
    params.shuffle = 1;
    params.n = 3;
    params.nIter = 50;
    params.num_top_feats = 5;
    %params.num_top_feats = length(features_training(1,:));
    params.threshmeth = 'euclidean';
    params.osname = 'SMOTE';
    params.fsprunecorr = 'true';
    params.featnames = feature_column_names;

    %% Use same patients in each fold for every iteration of cross validation.
        % In every iteration of cross validation, the same XX patients are used
        % for training in each fold.
        % In every iteration of cross validation, the same XX patients are used
        % for testing in each fold.
    params.shuffle = 1;

    % for idx=1:params.n
    %     train_folds{idx} = train_folds{1};
    %     test_folds{idx} = test_folds{1};
    % end
%     for i=1:params.nIter
%         [train_folds, test_folds] = nFold([], data_labels_training, params.shuffle, params.n);
%         params.subsets(i).training = train_folds;
%         params.subsets(i).testing = test_folds;
%     end
    save(strcat(output_path, 'Classifier_params.mat'),'params');
else
    params = load('Classifier_params.mat');
    params = params.params;
end

%% Create experimental log
log_file_name = strcat(output_path, "Summary_of_Experiment_", experiment_date, ".txt");
fileid = fopen(log_file_name,'w');
fprintf(fileid, strcat("Experiment Date and Time: ", string(datetime("now")), '\n\r'));
fprintf(fileid, strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf(fileid, "EXPERIMENT TYPE: %s \n", experiment_type);
fprintf(fileid, "FEATURE FAMILY: 3D Shape \n");
fprintf(fileid, strcat("Trained on ", int2str(length(features_training(:,1))), " patients \n\r"));
fprintf(fileid, strcat("Tested on ", int2str(length(features_holdout1(:,1)))," patients \n\r"));
fprintf(fileid, strcat("View: ", view{:}, '\n\r'));
fprintf(fileid, strcat("ROIs Analyzed: ", experiment_type, '\n\r'));
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






