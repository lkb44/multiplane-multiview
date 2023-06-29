clear;close all;clc;
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/general/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature Classifier/nFold_cross_validation/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature Classifier/training_and_testing_sets/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/feature_selection/mrmr_feature_select/'));
addpath(genpath("/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/mrmr_feature_select/"));
addpath(genpath("/Users/Tom/Documents/scripts/feature_selection/mrmr_feature_select/"));

%% Specify experiment information
feature_family = {'shape', 'texture'}; % shape or texture
dimension = 2;
view = {'Axial', 'Coronal'}; % Axial or coronal

% ROIs depends on feature family and view
% Shape - Axial: lumen, tumor, fat
% Shape - Coronal: lumen, rw, fat
% Texture - Axial: tumor, fat
% Texture - Coronal: rw, fat
shape_axial_rois = {'lumen', 'tumor', 'fat'};
shape_cor_rois = {'lumen', 'rw', 'fat'};
texture_axial_rois = {'tumor', 'fat'};
texture_cor_rois = {'rw', 'fat'};

% Experiment type depends on feature family and view
% Shape - Axial: lumen_only, tumor_only, or fat_only
% Shape - Coronal: lumen_only, rw_only, fat_only
% Texture - Axial: tumor_only, fat_only
% Texture - Coronal: rw_only, fat_only
experiment_type = 'all';

% Print experiment types
fprintf(strcat("USING THE FOLLOWING FEAUTRES: ", strjoin(feature_family), "\n"));
fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf("EXPERIMENT TYPE: %s \n", experiment_type);

%% Specify paths to Feature matrices
% Specify paths to ground truth labels
train_labels_path = "../Data/train_test_labels/train_labels.csv";
test_labels_path = "../Data/train_test_labels/test_labels.csv";
va_test_labels = "../Data/train_test_labels/va_test_labels.csv";

% Axial 2D Shape Features
ax_train_shape_path = "../Shape_Features/2D_Feature_Matrices/Axial/train/";
ax_test_shape_path = "../Shape_Features/2D_Feature_Matrices/Axial/test/";
ax_va_test_shape_path = "../Shape_Features/2D_Feature_Matrices/Axial/va_test/";

% Coronal 2D Shape Features
cor_train_shape_path = "../Shape_Features/2D_Feature_Matrices/Coronal/train/";
cor_test_shape_path = "../Shape_Features/2D_Feature_Matrices/Coronal/test/";
cor_va_test_shape_path = "../Shape_Features/2D_Feature_Matrices/Coronal/va_test/";

% Axial 2D Texture Features
ax_train_texture_path = "../Texture_Features/2D_Features/Axial/train/";
ax_test_texture_path = "../Texture_Features/2D_Features/Axial/test/";
ax_va_test_texture_path = "../Texture_Features/2D_Features/Axial/va_test/";

% Coronal 2D Texture Features
cor_train_texture_path = "../Texture_Features/2D_Features/Coronal/train/";
cor_test_texture_path = "../Texture_Features/2D_Features/Coronal/test/";
cor_va_test_texture_path = "../Texture_Features/2D_Features/Coronal/va_test/";

%% Read in Feature Names
shape_col_names_path = "../2D_Shape_Feature_Extraction_Code/2d_shape_feature_names.xlsx";
shape_col_names = readtable(shape_col_names_path, 'ReadVariableNames',false);
shape_col_names = table2cell(shape_col_names);

texture_col_names_path = "../Texture_Features/2D_Features/Axial/train/StatName.mat";
texture_col_names = load(texture_col_names_path);
texture_col_names = texture_col_names.nameStat;

% Combine cell arrays
%feature_column_names = cat(2, shape_feature_column_names, texture_feature_column_names);

%% Create output directory
% 
% output_path = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_and_texture_Feature_Results/";
% experiment_date = strrep(string(datetime(now,'ConvertFrom','datenum')), " ", "_");
% experiment_date = strrep(experiment_date, ":", "_");
% 
% % If working with only 1 view
% if(length(view) == 1)  
%     output_path = strcat(output_path, experiment_date, "_", strjoin(feature_family, "_"),"_", view{1}, "_", experiment_type, "/");
% else % Working with both views
%     output_path = strcat(output_path, experiment_date, "_", strjoin(feature_family, "_"), "_", strjoin(view, '_'), "_", experiment_type, "/");
% end
% 
% if(~exist(output_path, "dir"))
%     mkdir(output_path);
% end
% 
% fprintf("Created output directory for this experiment! \n");

%% Load in datasets

datasets = struct;

% Axial Shape
datasets = load_data({'shape'}, dimension, {'Axial'}, shape_axial_rois, ax_train_shape_path, datasets);
datasets = load_data({'shape'}, dimension, {'Axial'}, shape_axial_rois, ax_test_shape_path, datasets);
datasets = load_data({'shape'}, dimension, {'Axial'}, shape_axial_rois, ax_va_test_shape_path, datasets);

% Coronal Shape
datasets = load_data({'shape'}, dimension, {'Coronal'}, shape_cor_rois, cor_train_shape_path, datasets);
datasets = load_data({'shape'}, dimension, {'Coronal'}, shape_cor_rois, cor_test_shape_path, datasets);
datasets = load_data({'shape'}, dimension, {'Coronal'}, shape_cor_rois, cor_va_test_shape_path, datasets);

% Axial Texture
datasets = load_data({'texture'}, dimension, {'Axial'}, texture_axial_rois, ax_train_texture_path, datasets);
datasets = load_data({'texture'}, dimension, {'Axial'}, texture_axial_rois, ax_test_texture_path, datasets);
datasets = load_data({'texture'}, dimension, {'Axial'}, texture_axial_rois, ax_va_test_texture_path, datasets);    

% Coronal Texture
datasets = load_data({'texture'}, dimension, {'Coronal'}, texture_cor_rois, cor_train_texture_path, datasets);
datasets = load_data({'texture'}, dimension, {'Coronal'}, texture_cor_rois, cor_test_texture_path, datasets);
datasets = load_data({'texture'}, dimension, {'Coronal'}, texture_cor_rois, cor_va_test_texture_path, datasets);

%% Initialize training and holdout testing datasets
features_training = [];
features_holdout = [];
features_va = [];
feature_column_names = [];

%% Axial Shape Lumen
shape_ax_lumen_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_30_07_shape_Axial_lumen_only/Top_Features_Information.xlsx");
shape_ax_lumen_top_indices = shape_ax_lumen_top_feats.indices;
shape_ax_lumen_top_names = shape_ax_lumen_top_feats.names';

shape_ax_lumen_train = datasets.shape_Axial_train_lumen.feature_matrix;
shape_ax_lumen_train = shape_ax_lumen_train(:, shape_ax_lumen_top_indices);
features_training = shape_ax_lumen_train;

shape_ax_lumen_test = datasets.shape_Axial_test_lumen.feature_matrix;
shape_ax_lumen_test = shape_ax_lumen_test(:, shape_ax_lumen_top_indices);
features_holdout = shape_ax_lumen_test;

shape_ax_lumen_va = datasets.shape_Axial_va_test_lumen.feature_matrix;
shape_ax_lumen_va = shape_ax_lumen_va(:, shape_ax_lumen_top_indices);
features_va = shape_ax_lumen_va;

feature_column_names = shape_ax_lumen_top_names;

%% Axial Shape Tumor
shape_ax_tumor_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_31_07_shape_Axial_tumor_only/Top_Features_Information.xlsx");
shape_ax_tumor_top_indices = shape_ax_tumor_top_feats.indices;
shape_ax_tumor_top_names = shape_ax_tumor_top_feats.names';

shape_ax_tumor_train = datasets.shape_Axial_train_tumor.feature_matrix;
shape_ax_tumor_train = shape_ax_tumor_train(:, shape_ax_tumor_top_indices);
features_training = cat(2, features_training, shape_ax_tumor_train);

shape_ax_tumor_test = datasets.shape_Axial_test_tumor.feature_matrix;
shape_ax_tumor_test = shape_ax_tumor_test(:, shape_ax_tumor_top_indices);
features_holdout = cat(2, features_holdout, shape_ax_tumor_test);

shape_ax_tumor_va = datasets.shape_Axial_va_test_tumor.feature_matrix;
shape_ax_tumor_va = shape_ax_tumor_va(:, shape_ax_tumor_top_indices);
features_va = cat(2, features_va, shape_ax_tumor_va);

feature_column_names = cat(2, feature_column_names, shape_ax_tumor_top_names);

%% Axial Shape Fat
shape_ax_fat_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_32_11_shape_Axial_fat_only/Top_Features_Information.xlsx");
shape_ax_fat_top_indices = shape_ax_fat_top_feats.indices;
shape_ax_fat_top_names = shape_ax_fat_top_feats.names';

shape_ax_fat_train = datasets.shape_Axial_train_fat.feature_matrix;
shape_ax_fat_train = shape_ax_fat_train(:, shape_ax_fat_top_indices);
features_training = cat(2, features_training, shape_ax_fat_train);

shape_ax_fat_test = datasets.shape_Axial_test_fat.feature_matrix;
shape_ax_fat_test = shape_ax_fat_test(:, shape_ax_fat_top_indices);
features_holdout = cat(2, features_holdout, shape_ax_fat_test);

shape_ax_fat_va = datasets.shape_Axial_va_test_fat.feature_matrix;
shape_ax_fat_va = shape_ax_fat_va(:, shape_ax_fat_top_indices);
features_va = cat(2, features_va, shape_ax_fat_va);

feature_column_names = cat(2, feature_column_names, shape_ax_fat_top_names);

%% Coronal Shape Lumen
shape_cor_lumen_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_32_53_shape_Coronal_lumen_only/Top_Features_Information.xlsx");
shape_cor_lumen_top_indices = shape_cor_lumen_top_feats.indices;
shape_cor_lumen_top_names = shape_cor_lumen_top_feats.names';

shape_cor_lumen_train = datasets.shape_Coronal_train_lumen.feature_matrix;
shape_cor_lumen_train = shape_cor_lumen_train(:, shape_cor_lumen_top_indices);
features_training = cat(2, features_training, shape_cor_lumen_train);

shape_cor_lumen_test = datasets.shape_Coronal_test_lumen.feature_matrix;
shape_cor_lumen_test = shape_cor_lumen_test(:, shape_cor_lumen_top_indices);
features_holdout = cat(2, features_holdout, shape_cor_lumen_test);

shape_cor_lumen_va = datasets.shape_Coronal_va_test_lumen.feature_matrix;
shape_cor_lumen_va = shape_cor_lumen_va(:, shape_cor_lumen_top_indices);
features_va = cat(2, features_va, shape_cor_lumen_va);

feature_column_names = cat(2, feature_column_names, shape_cor_lumen_top_names);

%% Coronal Shape Rectal Wall
shape_cor_rw_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_33_23_shape_Coronal_rw_only/Top_Features_Information.xlsx");
shape_cor_rw_top_indices = shape_cor_rw_top_feats.indices;
shape_cor_rw_top_names = shape_cor_rw_top_feats.names';

shape_cor_rw_train = datasets.shape_Coronal_train_rw.feature_matrix;
shape_cor_rw_train = shape_cor_rw_train(:, shape_cor_rw_top_indices);
features_training = cat(2, features_training, shape_cor_rw_train);

shape_cor_rw_test = datasets.shape_Coronal_test_rw.feature_matrix;
shape_cor_rw_test = shape_cor_rw_test(:, shape_cor_rw_top_indices);
features_holdout = cat(2, features_holdout, shape_cor_rw_test);

shape_cor_rw_va = datasets.shape_Coronal_va_test_rw.feature_matrix;
shape_cor_rw_va = shape_cor_rw_va(:, shape_cor_rw_top_indices);
features_va = cat(2, features_va, shape_cor_rw_va);

feature_column_names = cat(2, feature_column_names, shape_cor_rw_top_names);

%% Coronal Shape Fat
shape_cor_fat_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_34_30_shape_Coronal_fat_only/Top_Features_Information.xlsx");
shape_cor_fat_top_indices = shape_cor_fat_top_feats.indices;
shape_cor_fat_top_names = shape_cor_fat_top_feats.names';

shape_cor_fat_train = datasets.shape_Coronal_train_fat.feature_matrix;
shape_cor_fat_train = shape_cor_fat_train(:, shape_cor_fat_top_indices);
features_training = cat(2, features_training, shape_cor_fat_train);

shape_cor_fat_test = datasets.shape_Coronal_test_fat.feature_matrix;
shape_cor_fat_test = shape_cor_fat_test(:, shape_cor_fat_top_indices);
features_holdout = cat(2, features_holdout, shape_cor_fat_test);

shape_cor_fat_va = datasets.shape_Coronal_va_test_fat.feature_matrix;
shape_cor_fat_va = shape_cor_fat_va(:, shape_cor_fat_top_indices);
features_va = cat(2, features_va, shape_cor_fat_va);

feature_column_names = cat(2, feature_column_names, shape_cor_fat_top_names);

%% Axial Texture Tumor
texture_ax_tumor_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_41_49_texture_Axial_tumor_only/Top_Features_Information.xlsx");
texture_ax_tumor_top_indices = texture_ax_tumor_top_feats.indices;
texture_ax_tumor_top_names = texture_ax_tumor_top_feats.names';

texture_ax_tumor_train = datasets.texture_Axial_train_tumor.feature_matrix;
texture_ax_tumor_train = texture_ax_tumor_train(:, texture_ax_tumor_top_indices);
features_training = cat(2, features_training, texture_ax_tumor_train);

texture_ax_tumor_test = datasets.texture_Axial_test_tumor.feature_matrix;
texture_ax_tumor_test = texture_ax_tumor_test(:, texture_ax_tumor_top_indices);
features_holdout = cat(2, features_holdout, texture_ax_tumor_test);

texture_ax_tumor_va = datasets.texture_Axial_va_test_tumor.feature_matrix;
texture_ax_tumor_va = texture_ax_tumor_va(:, texture_ax_tumor_top_indices);
features_va = cat(2, features_va, texture_ax_tumor_va);

feature_column_names = cat(2, feature_column_names, texture_ax_tumor_top_names);

%% Axial Texture Fat
texture_ax_fat_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_40_05_texture_Axial_fat_only/Top_Features_Information.xlsx");
texture_ax_fat_top_indices = texture_ax_fat_top_feats.indices;
texture_ax_fat_top_names = texture_ax_fat_top_feats.names';

texture_ax_fat_train = datasets.texture_Axial_train_fat.feature_matrix;
texture_ax_fat_train = texture_ax_fat_train(:, texture_ax_fat_top_indices);
features_training = cat(2, features_training, texture_ax_fat_train);

texture_ax_fat_test = datasets.texture_Axial_test_fat.feature_matrix;
texture_ax_fat_test = texture_ax_fat_test(:, texture_ax_fat_top_indices);
features_holdout = cat(2, features_holdout, texture_ax_fat_test);

texture_ax_fat_va = datasets.texture_Axial_va_test_fat.feature_matrix;
texture_ax_fat_va = texture_ax_fat_va(:, texture_ax_fat_top_indices);
features_va = cat(2, features_va, texture_ax_fat_va);

feature_column_names = cat(2, feature_column_names, texture_ax_fat_top_names);

%% Coronal Texture Rectal Wall
texture_cor_rw_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_42_46_texture_Coronal_rw_only/Top_Features_Information.xlsx");
texture_cor_rw_top_indices = texture_cor_rw_top_feats.indices;
texture_cor_rw_top_names = texture_cor_rw_top_feats.names';

texture_cor_rw_train = datasets.texture_Coronal_train_rw.feature_matrix;
texture_cor_rw_train = texture_cor_rw_train(:, texture_cor_rw_top_indices);
features_training = cat(2, features_training, texture_cor_rw_train);

texture_cor_rw_test = datasets.texture_Coronal_test_rw.feature_matrix;
texture_cor_rw_test = texture_cor_rw_test(:, texture_cor_rw_top_indices);
features_holdout = cat(2, features_holdout, texture_cor_rw_test);

texture_cor_rw_va = datasets.texture_Coronal_va_test_rw.feature_matrix;
texture_cor_rw_va = texture_cor_rw_va(:, texture_cor_rw_top_indices);
features_va = cat(2, features_va, texture_cor_rw_va);

feature_column_names = cat(2, feature_column_names, texture_cor_rw_top_names);

%% Coronal Texture Fat
texture_cor_fat_top_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_43_53_texture_Coronal_fat_only/Top_Features_Information.xlsx");
texture_cor_fat_top_indices = texture_cor_fat_top_feats.indices;
texture_cor_fat_top_names = texture_cor_fat_top_feats.names';


texture_cor_fat_train = datasets.texture_Coronal_train_fat.feature_matrix;
texture_cor_fat_train = texture_cor_fat_train(:, texture_cor_fat_top_indices);
features_training = cat(2, features_training, texture_cor_fat_train);

texture_cor_fat_test = datasets.texture_Coronal_test_fat.feature_matrix;
texture_cor_fat_test = texture_cor_fat_test(:, texture_cor_fat_top_indices);
features_holdout = cat(2, features_holdout, texture_cor_fat_test);

texture_cor_fat_va = datasets.texture_Coronal_va_test_fat.feature_matrix;
texture_cor_fat_va = texture_cor_fat_va(:, texture_cor_fat_top_indices);
features_va = cat(2, features_va, texture_cor_fat_va);

feature_column_names = cat(2, feature_column_names, texture_cor_fat_top_names);

%% Create holdout test set

features_holdout = cat(1, features_holdout, features_va);

%% Load in labels

data_labels_training = readmatrix(train_labels_path);
data_labels_training = data_labels_training(:, 2);

uh_data_labels_holdout = readmatrix(test_labels_path);
uh_data_labels_holdout = uh_data_labels_holdout(:, 2);

va_data_labels = readmatrix(va_test_labels);
va_data_labels = va_data_labels(:, 2);


%% Create holdout testing dataset
data_labels_holdout = vertcat(uh_data_labels_holdout, va_data_labels);

holdout_test_size = size(features_holdout);
assert(holdout_test_size(1) == 24, "The size of the dataset is incorrect!");

fprintf("Loaded in train and holdout testing datasets! \n");

%% Whiten the data
features_training = whitenData(features_training, 'modified');
features_holdout = whitenData(features_holdout, 'modified');
fprintf("Whitened the training and holdout testing datasets successfully! \n");

%% Initialize classifier params for classifier
saveParams = true;
if(isequal(saveParams,true))
    params.classifier='QDA';
    params.fsname='mrmr';
    params.shuffle = 1;
    params.n = 5;
    params.nIter = 50;
    params.num_top_feats = 5;
    %params.num_top_feats = length(features_training(1,:));
    params.threshmeth = 'euclidean';

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
fprintf(fileid, "FEATURE FAMILY: %dD %s \n", dimension, strjoin(feature_family));
fprintf(fileid, strcat("Trained on ", int2str(length(features_training(:,1))), " patients \n\r"));
fprintf(fileid, strcat("Tested on ", int2str(length(features_holdout(:,1)))," patients \n\r"));
fprintf(fileid, strcat("View: ", view{:}, '\n\r'));
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
num_va_pts = length(features_va(:, 1));

uh_predictions = stats.decisions(1:end-6);
va_predictions = stats.decisions(end-5:end); % VA patients were appended to the end of each feature matrix.


uh_metrics = compute_metrics(double(uh_predictions), uh_data_labels_holdout);
va_metrics = compute_metrics(double(va_predictions), va_data_labels);


% Save UH and VA results to separate files
uh_metrics_table = struct2table(uh_metrics);
writetable(uh_metrics_table, strcat(output_path,'Holdout_Testing_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'UH_Only');

va_metrics_table = struct2table(va_metrics);
writetable(va_metrics_table, strcat(output_path,'Holdout_Testing_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'VA_Only');


function results = compute_metrics(prediction, groundtruth)
%% COMPUTE_METRICS
%   Compute metrics obtained from the TrainTestClassifier.m
%   INPUTS:
%       prediction: boolean matrix
%           nx1 matrix of predictions generated by trained classifier. n is
%           the number of observations (i.e., rows/patients). Must have the
%           same dimension as prediction.
%       groundtruth: boolean matrix
%           nx1 matrix of ground truth classes. n is the number of
%           observations (i.e., rows/patients). Must have the same
%           dimension as prediction.
%   RETURNS:
%       results: struct
%           Struct containing values of metrics used to evaluate
%           classifier. Metrics include accuracy (ACC), AUC, area
%           under precision recall curve (AUPRC), positive predictive value 
%           (PPV), sensitivity (SENS), specificity (SPEC), F1-score
%           (FSCORE), kappa (KAPPA), MCC (MCC), true positives (TP), true
%           negatives (TN), false positives (FP), and false negatives (FN)

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

function dataset_struct = load_data(feature_family, dimension, view, rois, input_path, datasets)
%% LOAD_DATA
%   Loads in a dataset and stores dataset in a struct
%   INPUTS:
%       feature_family: cell array
%           Cell array of string with only one value: shape or texture.
%           Important for naming datasets.
%       dimension: int
%           Dimension of features. Currently supports only 2D features.
%           Important for naming datasets.
%      view: cell array
%           Cell array of string with only one value: Axial or Coronal.
%           Important for naming datasets.
%       rois: cell array
%           Cell array of strings describing regions from which features
%           were extracted. Important for naming datasets.
%       input_path: string
%           Base path to folder that contains dataset in .mat file
%       datasets struct
%           Struct into which datasets are loaded
%   RETURNS:
%       dataset_struct: struct
%           Struct that holds loaded datasets. Each struct field is a
%           dataset. The field names are based on feature family, view, and
%           roi.

    data_type = strsplit(input_path, '/');
    data_type = data_type{end-1};

    if(dimension == 2)
        for j=1:length(view)
            for i=1:length(rois)
                path_to_train = strcat(input_path, view{j}, "_", rois{i}, "_", data_type, "_2D_", feature_family{1}, "_features.mat");
                datasets.(strcat(feature_family{1}, "_", view{j}, "_", data_type, "_", rois{i})) = load(path_to_train);
            end
        end
        
    end
    
    dataset_struct = datasets;

end






