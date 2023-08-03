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
%% Specify dataset information
view = {'Axial'}; % {Can only be 'Axial'}, {'Coronal'} or {'Axial','Coronal'}
rois = {'Tumor'};
scheme = {'wilcoxon_qda'}; % {'wilcoxon_qda'}, {'wilcoxon_rf'}, {'mrmr_qda'}, or {'mrmr_rf'}
split = {'_TRG_T'};

if strcmp(rois, "Proximal_Fat5")
    experiment_type = 'proxfat5_only';
elseif strcmp(rois, "Proximal_Fat10")
    experiment_type = 'proxfat10_only';
elseif strcmp(rois, "Proximal_Fat15")
    experiment_type = 'proxfat15_only';
elseif strcmp(rois, "Fat")
    experiment_type = 'fat_only';
elseif strcmp(rois, "Tumor")
    experiment_type = 'tumor_only';
end
% Print experiment types
fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf("EXPERIMENT TYPE: %s \n", experiment_type);
%% Specify paths to ground truth labels
training_label_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/train_test_labels/train_labels_TRG_T.csv";
testing_label_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/train_test_labels/test_labels_TRG_T.csv';
train_matrix_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/AxialTrainingShapeFeatures_TRG_T/";
test_matrix_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/AxialTestingShapeFeatures_TRG_T/";


feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/2d_shape_feature_names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

fprintf("Got directory paths to feature matrices! \n");

%% Create output directory


output_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/";
experiment_date = strrep(string(datetime(now,'ConvertFrom','datenum')), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");
output_path = strcat(output_path, view, 'TrainingShapeFeatures', split, '/', experiment_date, "_", view, "_", split, "_", experiment_type, "_", scheme, "/");
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

%% Load in datasets

datasets = struct;

train_matrices = dir(fullfile(train_matrix_path));
train = {train_matrices.name};

fat_train = find(contains(train,'Axial_Fat'));
proxfat5_train = find(contains(train,'_ProxFat5'));
proxfat10_train = find(contains(train,'_ProxFat10'));
proxfat15_train = find(contains(train,'_ProxFat15'));
tumor_train = find(contains(train,'_Tumor'));

test_matrices = dir(fullfile(test_matrix_path));
test = {test_matrices.name};

fat_test = find(contains(test,'Axial_Fat'));
proxfat5_test = find(contains(test,'_ProxFat5'));
proxfat10_test = find(contains(test,'_ProxFat10'));
proxfat15_test = find(contains(test,'_ProxFat15'));
tumor_test = find(contains(test,'_Tumor'));

for i = 1:length(rois)
    if strcmp(rois(i), 'Fat')
        path_to_train = strcat(train_matrix_path, train{fat_train});
    end
    if strcmp(rois(i), 'Proximal_Fat5')
        path_to_train = strcat(train_matrix_path, train{proxfat5_train});
    end
    if strcmp(rois(i), 'Proximal_Fat10')
        path_to_train = strcat(train_matrix_path, train{proxfat10_train});
    end
    if strcmp(rois(i), 'Proximal_Fat15')
        path_to_train = strcat(train_matrix_path, train{proxfat15_train});
    end
    if strcmp(rois(i), 'Tumor')
        path_to_train = strcat(train_matrix_path, train{tumor_train});
    end

    if strcmp(rois(i), 'Fat')
        path_to_test = strcat(test_matrix_path, test{fat_test});
    end
    if strcmp(rois(i), 'Proximal_Fat5')
        path_to_test = strcat(test_matrix_path, test{proxfat5_test});
    end
    if strcmp(rois(i), 'Proximal_Fat10')
        path_to_test = strcat(test_matrix_path, test{proxfat10_test});
    end
    if strcmp(rois(i), 'Proximal_Fat15')
        path_to_test = strcat(test_matrix_path, test{proxfat15_test});
    end
    if strcmp(rois(i), 'Tumor')
        path_to_test = strcat(test_matrix_path, test{tumor_test});
    end
    datasets.(strcat("train_",rois{i})) = load(path_to_train);
    datasets.(strcat("test_",rois{i})) = load(path_to_test);
end

data_labels_training = readmatrix(training_label_path);
% data_labels_training = data_labels_training(:, 2);

data_labels_holdout = readmatrix(testing_label_path);
% data_labels_holdout = data_labels_holdout(:, 2);

if(length(view) == 1 && experiment_type == "tumor_only")
    features_training = datasets.train_Tumor.feature_matrix;
    features_holdout = datasets.test_Tumor.feature_matrix;
    feature_column_names = strcat(view{1}, "_tumor_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "fat_only")
    features_training = datasets.train_Fat.feature_matrix;
    features_holdout = datasets.test_Fat.feature_matrix;
    feature_column_names = strcat(view{1}, "_fat_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "proxfat5_only")
    features_training = datasets.train_Proximal_Fat5.feature_matrix;
    features_holdout = datasets.test_Proximal_Fat5.feature_matrix;
    feature_column_names = strcat(view{1}, "_proxfat5_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "proxfat10_only")
    features_training = datasets.train_Proximal_Fat10.feature_matrix;
    features_holdout = datasets.test_Proximal_Fat10.feature_matrix;
    feature_column_names = strcat(view{1}, "_proxfat10_", feature_column_names);

elseif(length(view) == 1 && experiment_type == "proxfat15_only")
    features_training = datasets.train_Proximal_Fat15.feature_matrix;
    features_holdout = datasets.test_Proximal_Fat15.feature_matrix;
    feature_column_names = strcat(view{1}, "_proxfat15_", feature_column_names);
end

holdout_test_size = size(features_holdout);
assert(holdout_test_size(1) == 20, "The size of the dataset is incorrect!");

fprintf("Loaded in train and holdout testing datasets! \n");
%% Whiten the data
features_training = whitenData(features_training, 'simple');
features_holdout = whitenData(features_holdout, 'simple');

if strcmp(rois(i), 'Proximal_Fat5')
    nan_vector = features_training(1, :);
    nan_indices = find(~isnan(nan_vector));
    features_training = features_training(:, nan_indices);
    features_holdout = features_holdout(:, nan_indices);
    feature_column_names = feature_column_names(:, nan_indices);

    inf_vector = features_training(1, :);
    inf_indices = find(~isinf(inf_vector));
    features_training = features_training(:, inf_indices);
    features_holdout = features_holdout(:, inf_indices);
    feature_column_names = feature_column_names(:, inf_indices);
elseif strcmp(rois(i), 'Proximal_Fat10')
    nan_vector = features_training(1, :);
    nan_indices = find(~isnan(nan_vector));
    features_training = features_training(:, nan_indices);
    features_holdout = features_holdout(:, nan_indices);
    feature_column_names = feature_column_names(:, nan_indices);

    inf_vector = features_training(1, :);
    inf_indices = find(~isinf(inf_vector));
    features_training = features_training(:, inf_indices);
    features_holdout = features_holdout(:, inf_indices);
    feature_column_names = feature_column_names(:, inf_indices);
elseif strcmp(rois(i), 'Proximal_Fat15')
    nan_vector = features_training(1, :);
    nan_indices = find(~isnan(nan_vector));
    features_training = features_training(:, nan_indices);
    features_holdout = features_holdout(:, nan_indices);
    feature_column_names = feature_column_names(:, nan_indices);

    inf_vector = features_training(1, :);
    inf_indices = find(~isinf(inf_vector));
    features_training = features_training(:, inf_indices);
    features_holdout = features_holdout(:, inf_indices);
    feature_column_names = feature_column_names(:, inf_indices);
elseif strcmp(rois(i), 'Fat')
    nan_vector = features_training(1, :);
    nan_indices = find(~isnan(nan_vector));
    features_training = features_training(:, nan_indices);
    features_holdout = features_holdout(:, nan_indices);
    feature_column_names = feature_column_names(:, nan_indices);

    inf_vector = features_training(1, :);
    inf_indices = find(~isinf(inf_vector));
    features_training = features_training(:, inf_indices);
    features_holdout = features_holdout(:, inf_indices);
    feature_column_names = feature_column_names(:, inf_indices);
elseif strcmp(rois(i), 'Tumor')
    nan_vector = features_training(1, :);
    nan_indices = find(~isnan(nan_vector));
    features_training = features_training(:, nan_indices);
    features_holdout = features_holdout(:, nan_indices);
    feature_column_names = feature_column_names(:, nan_indices);

    inf_vector = features_training(1, :);
    inf_indices = find(~isinf(inf_vector));
    features_training = features_training(:, inf_indices);
    features_holdout = features_holdout(:, inf_indices);
    feature_column_names = feature_column_names(:, inf_indices);
end

fprintf("Whitened the training and holdout testing datasets successfully! \n");

%% Initialize classifier params for classifier
saveParams = true;
if(isequal(saveParams,true))
    params.classifier='QDA';
    params.fsname='wilcoxon';
    params.shuffle = 1;
    params.n = 5;
    params.nIter = 50;
    params.num_top_feats = 5;
    %params.num_top_feats = length(features_training(1,:));
    params.threshmeth = 'euclidean';
%     params.osname = 'SMOTE';
%     params.fsprunecorr = 'true';
%     params.featnames = feature_column_names;

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
fprintf(fileid, strcat("Experiment Date and Time: ", string(datetime(now,'ConvertFrom','datenum')), '\n\r'));
fprintf(fileid, strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf(fileid, "EXPERIMENT TYPE: %s \n", experiment_type);
fprintf(fileid, "FEATURE FAMILY: 3D Shape \n");
fprintf(fileid, strcat("Trained on ", int2str(length(features_training(:,1))), " patients \n\r"));
fprintf(fileid, strcat("Tested on ", int2str(length(features_holdout(:,1)))," patients \n\r"));
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
test_feats = features_holdout(:,topResults.indices);
fprintf(strcat("Training on ",int2str(length(train_feats(1,:)))," features for ", int2str(length(features_training(:,1))), " patients \n"));
fprintf(strcat("Testing on ",int2str(length(test_feats(1,:)))," features for ", int2str(length(features_holdout(:,1)))," patients \n"));

[stats, methodstring] = TrainTestModel(params.classifier, train_feats, test_feats, data_labels_training, data_labels_holdout,classifier_thresh);
save(strcat(output_path,'Holdout_Testing_Dataset_Results.mat'),'stats');

stats_for_export = struct;
stats_for_export.ACC = round(stats.acc, 3);
stats_for_export.AUC = round(stats.AUC, 3);
stats_for_export.AUCPRC = round(stats.AUPRC, 3);
stats_for_export.SENS = round(stats.sens, 3);
stats_for_export.SPEC = round(stats.spec, 3);
stats_for_export.FSCORE = round(stats.Fscore, 3);
stats_for_export.TP = stats.tp;
stats_for_export.FP = stats.fp;
stats_for_export.FN = stats.fn;
stats_for_export.TN = stats.tn;
stats_for_export.KAPPA = round(stats.kappa, 3);
stats_for_export.MCC = round(stats.MCC, 3);
    

QDAResults = struct2table(stats_for_export);
writetable(QDAResults,strcat(output_path,'Holdout_Testing_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Entire_Test_Cohort_Results');

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






