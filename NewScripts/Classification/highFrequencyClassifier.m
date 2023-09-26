clear
close all
clc

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
%% Specify dataset information
view = {'Axial'}; % {Can only be 'Axial'}, {'Coronal'} or {'Axial','Coronal'}
rois = {'HighFrequencyFeats'};
scheme = {'mrmr_qda'}; % {'wilcoxon_qda'}, {'wilcoxon_rf'}, {'mrmr_qda'}, or {'mrmr_rf'}
split = {'missingCollage'};
region = 'ProxFat10';

% Print experiment types
fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
%% Specify paths to ground truth labels
label_path_root = '/Users/leobao/Documents/MultiPlanePipeline/Data/train_test_labels/MissingCollageLabels/';
matrix_path_root = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/';

training_label_path = string(strcat(label_path_root, view, '_', region, '/train_labels.csv'));
testing1_label_path = string(strcat(label_path_root, view, '_', region, '/test1_labels.csv'));
testing2_label_path = string(strcat(label_path_root, view, '_', region, '/test2_labels.csv'));

t_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor_Training.mat';
t_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor/Testing1/HighFrequency/Axial_Tumor_Testing1.mat';
t_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor_Testing2.mat';

pf5_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat5_Training.mat';
pf5_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat5_Testing1.mat';
pf5_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat5_Testing2.mat';

pf10_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat10_Training.mat';
pf10_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat10_Testing1.mat';
pf10_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat10_Testing2.mat';

pf15_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat15_Training.mat';
pf15_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat15_Testing1.mat';
pf15_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat15_Testing2.mat';

feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

f_indices = {};
t_indices = {596};
pf5_indices = {79, 366, 588};
pf10_indices = {720};
pf15_indices = {863};

fprintf("Got directory paths to feature matrices! \n");

%% Create output directory

output_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/";
experiment_date = strrep(string(datetime("now")), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");
output_path = strcat(output_path, experiment_date, "_", view, "_", rois, "_", scheme, "/");
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

%% Load in datasets
datasets = struct;

datasets.tumor_train = load(t_path_to_train);
datasets.tumor_test1 = load(t_path_to_test1);
datasets.tumor_test2 = load(t_path_to_test2);
datasets.proxfat5_train = load(pf5_path_to_train);
datasets.proxfat5_test1 = load(pf5_path_to_test1);
datasets.proxfat5_test2 = load(pf5_path_to_test2);
datasets.proxfat10_train = load(pf10_path_to_train);
datasets.proxfat10_test1 = load(pf10_path_to_test1);
datasets.proxfat10_test2 = load(pf10_path_to_test2);
datasets.proxfat15_train = load(pf15_path_to_train);
datasets.proxfat15_test1 = load(pf15_path_to_test1);
datasets.proxfat15_test2 = load(pf15_path_to_test2);

data_labels_training = readmatrix(training_label_path);
data_labels_holdout1 = readmatrix(testing1_label_path);
data_labels_holdout2 = readmatrix(testing2_label_path);

t_training_features = datasets.tumor_train.feature_matrix;
t_testing1_features = datasets.tumor_test1.feature_matrix;
t_testing2_features = datasets.tumor_test2.feature_matrix;
pf5_training_features = datasets.proxfat5_train.feature_matrix;
pf5_testing1_features = datasets.proxfat5_test1.feature_matrix;
pf5_testing2_features = datasets.proxfat5_test2.feature_matrix;
pf10_training_features = datasets.proxfat10_train.feature_matrix;
pf10_testing1_features = datasets.proxfat10_test1.feature_matrix;
pf10_testing2_features = datasets.proxfat10_test2.feature_matrix;
pf15_training_features = datasets.proxfat15_train.feature_matrix;
pf15_testing1_features = datasets.proxfat15_test1.feature_matrix;
pf15_testing2_features = datasets.proxfat15_test2.feature_matrix;

t_training_feats = zeros(size(t_training_features, 1), length(t_indices));
t_testing1_feats = zeros(size(t_testing1_features, 1), length(t_indices));
t_testing2_feats = zeros(size(t_testing2_features, 1), length(t_indices));
pf5_training_feats = zeros(size(pf5_training_features, 1), length(pf5_indices));
pf5_testing1_feats = zeros(size(pf5_testing1_features, 1), length(pf5_indices));
pf5_testing2_feats = zeros(size(pf5_testing2_features, 1), length(pf5_indices));
pf10_training_feats = zeros(size(pf10_training_features, 1), length(pf10_indices));
pf10_testing1_feats = zeros(size(pf10_testing1_features, 1), length(pf10_indices));
pf10_testing2_feats = zeros(size(pf10_testing2_features, 1), length(pf10_indices));
pf15_training_feats = zeros(size(pf15_training_features, 1), length(pf15_indices));
pf15_testing1_feats = zeros(size(pf15_testing1_features, 1), length(pf15_indices));
pf15_testing2_feats = zeros(size(pf15_testing2_features, 1), length(pf15_indices));

t_feature_names = cell(1, length(t_indices));
pf5_feature_names = cell(1, length(pf5_indices));
pf10_feature_names = cell(1, length(pf10_indices));
pf15_feature_names = cell(1, length(pf15_indices));

for i = 1 : length(t_indices)
    index = t_indices{i};
    t_training_feats(:, i) = t_training_features(:, index);
    t_testing1_feats(:, i) = t_testing1_features(:, index);
    t_testing2_feats(:, i) = t_testing2_features(:, index);
    t_feature_names(:, i) = feature_column_names(:, index);
end

for i = 1 : length(pf5_indices)
    index = pf5_indices{i};
    pf5_training_feats(:, i) = pf5_training_features(:, index);
    pf5_testing1_feats(:, i) = pf5_testing1_features(:, index);
    pf5_testing2_feats(:, i) = pf5_testing2_features(:, index);
    pf5_feature_names(:, i) = feature_column_names(:, index);
end

for i = 1 : length(pf10_indices)
    index = pf10_indices{i};
    pf10_training_feats(:, i) = pf10_training_features(:, index);
    pf10_testing1_feats(:, i) = pf10_testing1_features(:, index);
    pf10_testing2_feats(:, i) = pf10_testing2_features(:, index);
    pf10_feature_names(:, i) = feature_column_names(:, index);
end

for i = 1 : length(pf15_indices)
    index = pf15_indices{i};
    pf15_training_feats(:, i) = pf15_training_features(:, index);
    pf15_testing1_feats(:, i) = pf15_testing1_features(:, index);
    pf15_testing2_feats(:, i) = pf15_testing2_features(:, index);
    pf15_feature_names(:, i) = feature_column_names(:, index);
end

feature_column_names = [t_feature_names, pf5_feature_names, pf10_feature_names, pf15_feature_names];
features_training = [t_training_feats, pf5_training_feats, pf10_training_feats, pf15_training_feats];
features_holdout1 = [t_testing1_feats, pf5_testing1_feats, pf10_testing1_feats, pf15_testing1_feats];
features_holdout2 = [t_testing2_feats, pf5_testing2_feats, pf10_testing2_feats, pf15_testing2_feats];

feature_file_name = string(strcat(output_path, 'highfrequency_features.xlsx'));
feature_table = cell2table(feature_column_names);
writetable(feature_table, feature_file_name);

holdout_test_size1 = size(features_holdout1);
holdout_test_size2 = size(features_holdout2);
assert(holdout_test_size1(1) == 14, "The size of the dataset is incorrect!");
assert(holdout_test_size2(1) == 28, "The size of the dataset is incorrect!");

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
    params.n = 3;
    params.nIter = 50;
    params.num_top_feats = 6;
    %params.num_top_feats = length(features_training(1,:));
    params.threshmeth = 'euclidean';
    params.osname = 'SMOTE';
    params.fsprunecorr = 'false';
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
fprintf(fileid, "FEATURE FAMILY: 3D Shape \n");
fprintf(fileid, strcat("Trained on ", int2str(length(features_training(:,1))), " patients \n\r"));
fprintf(fileid, strcat("Tested on ", int2str(length(features_holdout1(:,1)))," patients \n\r"));
fprintf(fileid, strcat("View: ", view{:}, '\n\r'));
fprintf(fileid, strcat("Selection Method: ", params.fsname, '\n\r'));
fprintf(fileid, strcat("Classifier: ", params.classifier, '\n\r'));
fprintf(fileid, strcat("Number of Top Features: ", num2str(params.num_top_feats), '\n\r'));
fprintf(fileid, strcat("Number of Folds: ", num2str(params.n), '\n\r'));
fclose(fileid);

%% Train and test the QDA classifier
train_feats = features_training;
test1_feats = features_holdout1;
test2_feats = features_holdout2;
fprintf(strcat("Training on ",int2str(length(train_feats(1,:)))," features for ", int2str(length(features_training(:,1))), " patients \n"));
fprintf(strcat("Testing on ",int2str(length(test1_feats(1,:)))," features for ", int2str(length(features_holdout1(:,1)))," patients \n"));

[stats1, methodstring] = TrainTestModel(params.classifier, train_feats, test1_feats, data_labels_training, data_labels_holdout1);
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

[stats2, methodstring] = TrainTestModel(params.classifier, train_feats, test2_feats, data_labels_training, data_labels_holdout2);
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

% clear
% close all
% clc
% 
% if isunix
%     addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/general'));
%     addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/nFold_cross_validation/'));
%     addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/training_and_testing_sets/'));
%     addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/classifiers/'));
%     addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/'));
%     addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
%     addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/Scripts-DataPreProcessing/Util-preProcessing'));
%     % addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/estpab.mexmac'));
% end
% %% Specify dataset information
% view = {'Axial'}; % {Can only be 'Axial'}, {'Coronal'} or {'Axial','Coronal'}
% rois = {'HighFrequencyFeats'};
% scheme = {'mrmr_qda'}; % {'wilcoxon_qda'}, {'wilcoxon_rf'}, {'mrmr_qda'}, or {'mrmr_rf'}
% split = {'missingCollage'};
% region = 'ProxFat10';
% 
% % Print experiment types
% fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
% %% Specify paths to ground truth labels
% label_path_root = '/Users/leobao/Documents/MultiPlanePipeline/Data/train_test_labels/MissingCollageLabels/';
% matrix_path_root = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/';
% 
% training_label_path = string(strcat(label_path_root, view, '_', region, '/train_labels.csv'));
% testing1_label_path = string(strcat(label_path_root, view, '_', region, '/test1_labels.csv'));
% testing2_label_path = string(strcat(label_path_root, view, '_', region, '/test2_labels.csv'));
% 
% t_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor_Training.mat';
% t_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor/Testing1/HighFrequency/Axial_Tumor_Testing1.mat';
% t_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor_Testing2.mat';
% 
% pf5_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat5_Training.mat';
% pf5_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat5_Testing1.mat';
% pf5_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat5_Testing2.mat';
% 
% pf10_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat10_Training.mat';
% pf10_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat10_Testing1.mat';
% pf10_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat10_Testing2.mat';
% 
% pf15_path_to_train = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat15_Training.mat';
% pf15_path_to_test1 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat15_Testing1.mat';
% pf15_path_to_test2 = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_ProxFat15_Testing2.mat';
% 
% feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
% feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
% feature_column_names = table2cell(feature_column_names);
% 
% f_indices = {};
% t_indices = {596};
% pf5_indices = {79, 366, 588};
% pf10_indices = {720};
% pf15_indices = {863};
% 
% fprintf("Got directory paths to feature matrices! \n");
% 
% %% Create output directory
% 
% output_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/";
% experiment_date = strrep(string(datetime("now")), " ", "_");
% experiment_date = strrep(experiment_date, ":", "_");
% output_path = strcat(output_path, experiment_date, "_", view, "_", rois, "_", scheme, "/");
% if(~exist(output_path, "dir"))
%     mkdir(output_path);
% end
% 
% fprintf("Created output directory for this experiment! \n");
% 
% %% Load in datasets
% datasets = struct;
% 
% datasets.tumor_train = load(t_path_to_train);
% datasets.tumor_test1 = load(t_path_to_test1);
% datasets.tumor_test2 = load(t_path_to_test2);
% datasets.proxfat5_train = load(pf5_path_to_train);
% datasets.proxfat5_test1 = load(pf5_path_to_test1);
% datasets.proxfat5_test2 = load(pf5_path_to_test2);
% datasets.proxfat10_train = load(pf10_path_to_train);
% datasets.proxfat10_test1 = load(pf10_path_to_test1);
% datasets.proxfat10_test2 = load(pf10_path_to_test2);
% datasets.proxfat15_train = load(pf15_path_to_train);
% datasets.proxfat15_test1 = load(pf15_path_to_test1);
% datasets.proxfat15_test2 = load(pf15_path_to_test2);
% 
% data_labels_training = readmatrix(training_label_path);
% data_labels_holdout1 = readmatrix(testing1_label_path);
% data_labels_holdout2 = readmatrix(testing2_label_path);
% 
% t_training_features = datasets.tumor_train.feature_matrix;
% t_testing1_features = datasets.tumor_test1.feature_matrix;
% t_testing2_features = datasets.tumor_test2.feature_matrix;
% pf5_training_features = datasets.proxfat5_train.feature_matrix;
% pf5_testing1_features = datasets.proxfat5_test1.feature_matrix;
% pf5_testing2_features = datasets.proxfat5_test2.feature_matrix;
% pf10_training_features = datasets.proxfat10_train.feature_matrix;
% pf10_testing1_features = datasets.proxfat10_test1.feature_matrix;
% pf10_testing2_features = datasets.proxfat10_test2.feature_matrix;
% pf15_training_features = datasets.proxfat15_train.feature_matrix;
% pf15_testing1_features = datasets.proxfat15_test1.feature_matrix;
% pf15_testing2_features = datasets.proxfat15_test2.feature_matrix;
% 
% t_training_feats = zeros(size(t_training_features, 1), length(t_indices));
% t_testing1_feats = zeros(size(t_testing1_features, 1), length(t_indices));
% t_testing2_feats = zeros(size(t_testing2_features, 1), length(t_indices));
% pf5_training_feats = zeros(size(pf5_training_features, 1), length(pf5_indices));
% pf5_testing1_feats = zeros(size(pf5_testing1_features, 1), length(pf5_indices));
% pf5_testing2_feats = zeros(size(pf5_testing2_features, 1), length(pf5_indices));
% pf10_training_feats = zeros(size(pf10_training_features, 1), length(pf10_indices));
% pf10_testing1_feats = zeros(size(pf10_testing1_features, 1), length(pf10_indices));
% pf10_testing2_feats = zeros(size(pf10_testing2_features, 1), length(pf10_indices));
% pf15_training_feats = zeros(size(pf15_training_features, 1), length(pf15_indices));
% pf15_testing1_feats = zeros(size(pf15_testing1_features, 1), length(pf15_indices));
% pf15_testing2_feats = zeros(size(pf15_testing2_features, 1), length(pf15_indices));
% 
% t_feature_names = cell(1, length(t_indices));
% pf5_feature_names = cell(1, length(pf5_indices));
% pf10_feature_names = cell(1, length(pf10_indices));
% pf15_feature_names = cell(1, length(pf15_indices));
% 
% for i = 1 : length(t_indices)
%     index = t_indices{i};
%     t_training_feats(:, i) = t_training_features(:, index);
%     t_testing1_feats(:, i) = t_testing1_features(:, index);
%     t_testing2_feats(:, i) = t_testing2_features(:, index);
%     t_feature_names(:, i) = feature_column_names(:, index);
% end
% 
% for i = 1 : length(pf5_indices)
%     index = pf5_indices{i};
%     pf5_training_feats(:, i) = pf5_training_features(:, index);
%     pf5_testing1_feats(:, i) = pf5_testing1_features(:, index);
%     pf5_testing2_feats(:, i) = pf5_testing2_features(:, index);
%     pf5_feature_names(:, i) = feature_column_names(:, index);
% end
% 
% for i = 1 : length(pf10_indices)
%     index = pf10_indices{i};
%     pf10_training_feats(:, i) = pf10_training_features(:, index);
%     pf10_testing1_feats(:, i) = pf10_testing1_features(:, index);
%     pf10_testing2_feats(:, i) = pf10_testing2_features(:, index);
%     pf10_feature_names(:, i) = feature_column_names(:, index);
% end
% 
% for i = 1 : length(pf15_indices)
%     index = pf15_indices{i};
%     pf15_training_feats(:, i) = pf15_training_features(:, index);
%     pf15_testing1_feats(:, i) = pf15_testing1_features(:, index);
%     pf15_testing2_feats(:, i) = pf15_testing2_features(:, index);
%     pf15_feature_names(:, i) = feature_column_names(:, index);
% end
% 
% feature_column_names = [t_feature_names, pf5_feature_names, pf10_feature_names, pf15_feature_names];
% features_training = [t_training_feats, pf5_training_feats, pf10_training_feats, pf15_training_feats];
% features_holdout1 = [t_testing1_feats, pf5_testing1_feats, pf10_testing1_feats, pf15_testing1_feats];
% features_holdout2 = [t_testing2_feats, pf5_testing2_feats, pf10_testing2_feats, pf15_testing2_feats];
% 
% feature_file_name = string(strcat(output_path, 'highfrequency_features.xlsx'));
% feature_table = cell2table(feature_column_names);
% writetable(feature_table, feature_file_name);
% 
% holdout_test_size1 = size(features_holdout1);
% holdout_test_size2 = size(features_holdout2);
% assert(holdout_test_size1(1) == 14, "The size of the dataset is incorrect!");
% assert(holdout_test_size2(1) == 28, "The size of the dataset is incorrect!");
% 
% fprintf("Loaded in train and holdout testing datasets! \n");
% %% Whiten the data
% features_training = whitenData(features_training, 'simple');
% features_holdout1 = whitenData(features_holdout1, 'simple');
% features_holdout2 = whitenData(features_holdout2, 'simple');
% 
% fprintf("Whitened the training and holdout testing datasets successfully! \n");
% 
% %% Initialize classifier params for classifier
% saveParams = true;
% if(isequal(saveParams,true))
%     if strcmp(scheme, 'wilcoxon_qda')
%         params.classifier='QDA';
%         params.fsname='wilcoxon';
%     elseif strcmp(scheme, 'wilcoxon_lda')
%         params.classifier='LDA';
%         params.fsname='wilcoxon';
%     elseif strcmp(scheme, 'wilcoxon_rf')
%         params.classifier='RANDOMFOREST';
%         params.fsname='wilcoxon';
%     elseif strcmp(scheme, 'wilcoxon_svm')
%         params.classifier='SVM';
%         params.fsname='wilcoxon';
%     elseif strcmp(scheme, 'mrmr_qda')
%         params.classifier='QDA';
%         params.fsname='mrmr';
%     elseif strcmp(scheme, 'mrmr_lda')
%         params.classifier='LDA';
%         params.fsname='mrmr';
%     elseif strcmp(scheme, 'mrmr_rf')
%         params.classifier='RANDOMFOREST';
%         params.fsname='mrmr';
%     elseif strcmp(scheme, 'mrmr_svm')
%         params.classifier='SVM';
%         params.fsname='mrmr';
%     end
%     params.shuffle = 1;
%     params.n = 3;
%     params.nIter = 50;
%     params.num_top_feats = 6;
%     %params.num_top_feats = length(features_training(1,:));
%     params.threshmeth = 'euclidean';
%     params.osname = 'SMOTE';
%     params.fsprunecorr = 'false';
%     params.featnames = feature_column_names;
% 
%     %% Use same patients in each fold for every iteration of cross validation.
%         % In every iteration of cross validation, the same XX patients are used
%         % for training in each fold.
%         % In every iteration of cross validation, the same XX patients are used
%         % for testing in each fold.
%     params.shuffle = 1;
% 
%     % for idx=1:params.n
%     %     train_folds{idx} = train_folds{1};
%     %     test_folds{idx} = test_folds{1};
%     % end
% %     for i=1:params.nIter
% %         [train_folds, test_folds] = nFold([], data_labels_training, params.shuffle, params.n);
% %         params.subsets(i).training = train_folds;
% %         params.subsets(i).testing = test_folds;
% %     end
%     save(strcat(output_path, 'Classifier_params.mat'),'params');
% else
%     params = load('Classifier_params.mat');
%     params = params.params;
% end
% 
% %% Create experimental log
% log_file_name = strcat(output_path, "Summary_of_Experiment_", experiment_date, ".txt");
% fileid = fopen(log_file_name,'w');
% fprintf(fileid, strcat("Experiment Date and Time: ", string(datetime("now")), '\n\r'));
% fprintf(fileid, strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
% fprintf(fileid, "FEATURE FAMILY: 3D Shape \n");
% fprintf(fileid, strcat("Trained on ", int2str(length(features_training(:,1))), " patients \n\r"));
% fprintf(fileid, strcat("Tested on ", int2str(length(features_holdout1(:,1)))," patients \n\r"));
% fprintf(fileid, strcat("View: ", view{:}, '\n\r'));
% fprintf(fileid, strcat("Selection Method: ", params.fsname, '\n\r'));
% fprintf(fileid, strcat("Classifier: ", params.classifier, '\n\r'));
% fprintf(fileid, strcat("Number of Top Features: ", num2str(params.num_top_feats), '\n\r'));
% fprintf(fileid, strcat("Number of Folds: ", num2str(params.n), '\n\r'));
% fclose(fileid);
% 
% %% Run K-fold Cross Validation
% fprintf("Running K-Fold Cross Validation with Feature Selection... \n");
% feat_stats = nFoldCV_withFS(features_training,data_labels_training,params);
% save(strcat(output_path,'CV_Feature_Results.mat'),'feat_stats');
% FSResults = getCVTopFeatures(features_training, params, feat_stats, data_labels_training, output_path);
% 
% %% Analyze the output of feature selection
% [topResults,meanValues] = analyze_CV_results(output_path, feature_column_names, params.num_top_feats);
% fprintf("Found the top ranked features! \n")
% 
% %% Set the optimal threshold from analyze_CV_results
% classifier_thresh = meanValues.Opt_Thres;
% 
% %% Train and test the QDA classifier
% train_feats = features_training(:,topResults.indices);
% test1_feats = features_holdout1(:,topResults.indices);
% test2_feats = features_holdout2(:,topResults.indices);
% fprintf(strcat("Training on ",int2str(length(train_feats(1,:)))," features for ", int2str(length(features_training(:,1))), " patients \n"));
% fprintf(strcat("Testing on ",int2str(length(test1_feats(1,:)))," features for ", int2str(length(features_holdout1(:,1)))," patients \n"));
% 
% [stats1, methodstring] = TrainTestModel(params.classifier, train_feats, test1_feats, data_labels_training, data_labels_holdout1,classifier_thresh);
% save(strcat(output_path,'Holdout_Testing1_Dataset_Results.mat'),'stats1');
% 
% stats_for_export1 = struct;
% stats_for_export1.ACC = round(stats1.acc, 3);
% stats_for_export1.AUC = round(stats1.AUC, 3);
% stats_for_export1.AUCPRC = round(stats1.AUPRC, 3);
% stats_for_export1.SENS = round(stats1.sens, 3);
% stats_for_export1.SPEC = round(stats1.spec, 3);
% stats_for_export1.FSCORE = round(stats1.Fscore, 3);
% stats_for_export1.TP = stats1.tp;
% stats_for_export1.FP = stats1.fp;
% stats_for_export1.FN = stats1.fn;
% stats_for_export1.TN = stats1.tn;
% stats_for_export1.KAPPA = round(stats1.kappa, 3);
% stats_for_export1.MCC = round(stats1.MCC, 3);
%     
% 
% QDAResults1 = struct2table(stats_for_export1);
% writetable(QDAResults1,strcat(output_path,'Holdout_Testing1_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Entire_Test_Cohort_Results');
% 
% fprintf(strcat("Testing on ",int2str(length(test2_feats(1,:)))," features for ", int2str(length(features_holdout2(:,1)))," patients \n"));
% 
% [stats2, methodstring] = TrainTestModel(params.classifier, train_feats, test2_feats, data_labels_training, data_labels_holdout2,classifier_thresh);
% save(strcat(output_path,'Holdout_Testing2_Dataset_Results.mat'),'stats2');
% 
% stats_for_export2 = struct;
% stats_for_export2.ACC = round(stats2.acc, 3);
% stats_for_export2.AUC = round(stats2.AUC, 3);
% stats_for_export2.AUCPRC = round(stats2.AUPRC, 3);
% stats_for_export2.SENS = round(stats2.sens, 3);
% stats_for_export2.SPEC = round(stats2.spec, 3);
% stats_for_export2.FSCORE = round(stats2.Fscore, 3);
% stats_for_export2.TP = stats2.tp;
% stats_for_export2.FP = stats2.fp;
% stats_for_export2.FN = stats2.fn;
% stats_for_export2.TN = stats2.tn;
% stats_for_export2.KAPPA = round(stats2.kappa, 3);
% stats_for_export2.MCC = round(stats2.MCC, 3);
%     
% 
% QDAResults2 = struct2table(stats_for_export2);
% writetable(QDAResults2,strcat(output_path,'Holdout_Testing2_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Entire_Test_Cohort_Results');
% %% Separate and save UH and VA results
% 
% function results = compute_metrics(prediction, groundtruth)
% 
%     assert(size(prediction, 1) == size(groundtruth, 1), "Input prediction array and ground truth array are not the same size in the first axis!");
%     assert(size(prediction, 2) == size(groundtruth, 2), "Input prediction array and ground truth array are not the same size in teh second axis!");
%     
%     tp = length(find(prediction(prediction == groundtruth) == 1));
%     tn = length(find(prediction(prediction == groundtruth) == 0));
%     fn = length(find(prediction(prediction ~= groundtruth) == 0));
%     fp = length(find(prediction(prediction ~= groundtruth) == 1));
%     
%     [FPR,TPR,T,AUC,OPTROCPT,~,~] = perfcurve(groundtruth, prediction, 1, 'XVals', [0:0.02:1]);  % calculate AUC. 'perfcurve' can also calculate sens, spec etc. to plot the ROC curve.
%     [RECALL,PRECISION,~,AUPRC,~,~,~] = perfcurve(groundtruth, prediction, 1, 'XVals', [0:0.02:1], 'xCrit', 'reca', 'yCrit', 'prec'); 
%     
%     ACC = (tp+tn)/(tp+tn+fp+fn);
%     results.ACC = round(ACC, 3);
%     
%     results.AUC = round(AUC, 3);
%     results.AUCPRC = round(AUPRC, 3);
%     
%     PPV = tp/(tp+fp);
%     
%     SENS = tp/(tp+fn);
%     results.SENS = round(SENS, 3);
%     
%     SPEC = tn/(fp+tn);
%     results.SPEC = round(SPEC, 3);
%     
%     FSCORE = 2*tp/(2*tp+fp+fn);
%     results.FSCORE = round(FSCORE, 3);
%     
%     results.TP = tp;
%     results.FP = fp;
%     results.FN = fn;
%     results.TN = tn;
%     
%     Pre = ((tp+fp)*(tp+fn) + (tn+fn)*(tn+fp)) / (tp+tn+fp+fn)^2;
%     KAPPA = (results.ACC - Pre) / (1 - Pre);
%     results.KAPPA = round(KAPPA, 3);
%     
%     MCC = (tp*tn - fp*fn) / (((tp+fp)*(fn+tn)*(fp+tn)*(tp+fn))^0.5);
%     results.MCC = round(MCC, 3);
% 
% end






