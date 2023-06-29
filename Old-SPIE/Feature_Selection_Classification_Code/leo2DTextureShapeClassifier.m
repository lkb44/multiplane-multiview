clear
close all
clc

addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/general'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/nFold_cross_validation/'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/training_and_testing_sets/'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/classification/classifiers/'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));

view = {'Coronal', 'Axial'};
rois = {'Proximal_Fat', 'Tumor', 'Fat'};
feature_family = {'Texture', 'Shape'};
scheme = {'wilcoxon_qda'};
numtop = {'noPruning'};
dimension = 2;

shape_axial_rois = {'Proximal_Fat', 'Tumor', 'Fat'};
shape_cor_rois = {'Proximal_Fat', 'Tumor', 'Fat'};
texture_axial_rois = {'Proximal_Fat', 'Tumor', 'Fat'};
texture_cor_rois = {'Proximal_Fat', 'Tumor', 'Fat'};

% if strcmp(rois, "Proximal_Fat")
%     experiment_type = 'proxfat_only';
% elseif strcmp(rois, "Fat")
%     experiment_type = 'fat_only';
% elseif strcmp(rois, "Tumor")
experiment_type = 'all';
% end

% Print experiment types
fprintf(strcat("USING THE FOLLOWING FEAUTRES: ", strjoin(feature_family), "\n"));
fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf("EXPERIMENT TYPE: %s \n", experiment_type);

%% Specify paths to Feature matrices
% Specify paths to ground truth labels
training_label_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/train_test_labels/train_labels.csv";
testing_label_path = '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/train_test_labels/test_labels.csv';

% Axial 2D Shape Features
ax_train_shape_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/";
ax_test_shape_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTestingShapeFeatures/";

% Coronal 2D Shape Features
cor_train_shape_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingShapeFeatures/";
cor_test_shape_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTestingShapeFeatures/";

% Axial 2D Texture Features
ax_train_texture_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/";
ax_test_texture_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTestingTextureFeatures/";

% Coronal 2D Texture Features
cor_train_texture_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/";
cor_test_texture_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTestingTextureFeatures/";

%% Read in Feature Names
shape_col_names_path = '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/2d_shape_feature_names.xlsx';
shape_col_names = readtable(shape_col_names_path, 'ReadVariableNames',false);
shape_col_names = table2cell(shape_col_names);

texture_col_names = readtable('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/column_feature_names.xlsx', 'ReadVariableNames',false);
% texture_col_names = load(texture_col_names_path);
% texture_col_names = texture_col_names.nameStat;

% Combine cell arrays
feature_column_names = cat(2, shape_col_names, texture_col_names);

%% Create output directory

output_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/";
experiment_date = strrep(string(datetime(now,'ConvertFrom','datenum')), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");
output_path = strcat(output_path, experiment_date, "_", experiment_type, "_", scheme, "_", numtop, "/");
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

%% Load in datasets

datasets = struct;

% Axial Shape
datasets = load_data({'shape'}, dimension, {'Axial'}, shape_axial_rois, ax_train_shape_path, datasets);
datasets = load_data({'shape'}, dimension, {'Axial'}, shape_axial_rois, ax_test_shape_path, datasets);

% Coronal Shape
datasets = load_data({'shape'}, dimension, {'Coronal'}, shape_cor_rois, cor_train_shape_path, datasets);
datasets = load_data({'shape'}, dimension, {'Coronal'}, shape_cor_rois, cor_test_shape_path, datasets);

% Axial Texture
datasets = load_data({'texture'}, dimension, {'Axial'}, texture_axial_rois, ax_train_texture_path, datasets);
datasets = load_data({'texture'}, dimension, {'Axial'}, texture_axial_rois, ax_test_texture_path, datasets);   

% Coronal Texture
datasets = load_data({'texture'}, dimension, {'Coronal'}, texture_cor_rois, cor_train_texture_path, datasets);
datasets = load_data({'texture'}, dimension, {'Coronal'}, texture_cor_rois, cor_test_texture_path, datasets);

% dataset_names = fieldnames(datasets);
% for k=1:length(dataset_names)
%     
%     if(feature_family{1} == "shape")
%         
%         matrix_size = size(datasets.(dataset_names{k}).feature_matrix);
%         assert( matrix_size(2) == 18)
%         
%     end
%     
%     if(feature_family{1} == "texture")
%         
%         matrix_size = size(datasets.(dataset_names{k}).feature_matrix);
%         assert( matrix_size(2) == 656)
%     end
%     
% end



%% Initialize training and holdout testing datasets
features_training = [];
features_holdout = [];
features_va = [];
feature_column_names = [];

%% Axial Shape Proximal Fat
shape_ax_proxfat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/16-Feb-2023_21_53_00_Axial_proxfat_only_wilcoxon_qda/Top_Features_Information.xlsx");
shape_ax_proxfat_top_indices = shape_ax_proxfat_top_feats.indices;
shape_ax_proxfat_top_names = shape_ax_proxfat_top_feats.names';

shape_ax_proxfat_train = datasets.shape_Axial_AxialTrainingShapeFeatures_Proximal_Fat.feature_matrix;
shape_ax_proxfat_train = shape_ax_proxfat_train(:, shape_ax_proxfat_top_indices);
features_training = shape_ax_proxfat_train;

shape_ax_proxfat_test = datasets.shape_Axial_AxialTestingShapeFeatures_Proximal_Fat.feature_matrix;
shape_ax_proxfat_test = shape_ax_proxfat_test(:, shape_ax_proxfat_top_indices);
features_holdout = shape_ax_proxfat_test;

feature_column_names = shape_ax_proxfat_top_names;
%% Axial Shape Tumor
shape_ax_tumor_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/16-Feb-2023_22_09_19_Axial_tumor_only_wilcoxon_qda/Top_Features_Information.xlsx");
shape_ax_tumor_top_indices = shape_ax_tumor_top_feats.indices;
shape_ax_tumor_top_names = shape_ax_tumor_top_feats.names';

shape_ax_tumor_train = datasets.shape_Axial_AxialTrainingShapeFeatures_Tumor.feature_matrix;
shape_ax_tumor_train = shape_ax_tumor_train(:, shape_ax_tumor_top_indices);
features_training = cat(2, features_training, shape_ax_tumor_train);

shape_ax_tumor_test = datasets.shape_Axial_AxialTestingShapeFeatures_Tumor.feature_matrix;
shape_ax_tumor_test = shape_ax_tumor_test(:, shape_ax_tumor_top_indices);
features_holdout = cat(2, features_holdout, shape_ax_tumor_test);

feature_column_names = cat(2, feature_column_names, shape_ax_tumor_top_names);

%% Axial Shape Fat
shape_ax_fat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/16-Feb-2023_21_53_12_Axial_fat_only_wilcoxon_qda/Top_Features_Information.xlsx");
shape_ax_fat_top_indices = shape_ax_fat_top_feats.indices;
shape_ax_fat_top_names = shape_ax_fat_top_feats.names';

shape_ax_fat_train = datasets.shape_Axial_AxialTrainingShapeFeatures_Fat.feature_matrix;
shape_ax_fat_train = shape_ax_fat_train(:, shape_ax_fat_top_indices);
features_training = cat(2, features_training, shape_ax_fat_train);

shape_ax_fat_test = datasets.shape_Axial_AxialTestingShapeFeatures_Fat.feature_matrix;
shape_ax_fat_test = shape_ax_fat_test(:, shape_ax_fat_top_indices);
features_holdout = cat(2, features_holdout, shape_ax_fat_test);

feature_column_names = cat(2, feature_column_names, shape_ax_fat_top_names);

%% Coronal Shape Proximal Fat
shape_cor_proxfat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingShapeFeatures/16-Feb-2023_22_02_25_Coronal_proxfat_only_wilcoxon_qda/Top_Features_Information.xlsx");
shape_cor_proxfat_top_indices = shape_cor_proxfat_top_feats.indices;
shape_cor_proxfat_top_names = shape_cor_proxfat_top_feats.names';

shape_cor_proxfat_train = datasets.shape_Coronal_CoronalTrainingShapeFeatures_Proximal_Fat.feature_matrix;
shape_cor_proxfat_train = shape_cor_proxfat_train(:, shape_cor_proxfat_top_indices);
features_training = cat(2, features_training, shape_cor_proxfat_train);

shape_cor_proxfat_test = datasets.shape_Coronal_CoronalTestingShapeFeatures_Proximal_Fat.feature_matrix;
shape_cor_proxfat_test = shape_cor_proxfat_test(:, shape_cor_proxfat_top_indices);
features_holdout = cat(2, features_holdout, shape_cor_proxfat_test);

feature_column_names = cat(2, feature_column_names, shape_cor_proxfat_top_names);
%% Coronal Shape Tumor
shape_cor_tumor_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingShapeFeatures/16-Feb-2023_22_05_54_Coronal_tumor_only_wilcoxon_qda/Top_Features_Information.xlsx");
shape_cor_tumor_top_indices = shape_cor_tumor_top_feats.indices;
shape_cor_tumor_top_names = shape_cor_tumor_top_feats.names';

shape_cor_tumor_train = datasets.shape_Coronal_CoronalTrainingShapeFeatures_Tumor.feature_matrix;
shape_cor_tumor_train = shape_cor_tumor_train(:, shape_cor_tumor_top_indices);
features_training = cat(2, features_training, shape_cor_tumor_train);

shape_cor_tumor_test = datasets.shape_Coronal_CoronalTestingShapeFeatures_Tumor.feature_matrix;
shape_cor_tumor_test = shape_cor_tumor_test(:, shape_cor_tumor_top_indices);
features_holdout = cat(2, features_holdout, shape_cor_tumor_test);

feature_column_names = cat(2, feature_column_names, shape_cor_tumor_top_names);

%% Coronal Shape Fat
shape_cor_fat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingShapeFeatures/16-Feb-2023_22_04_10_Coronal_fat_only_wilcoxon_qda/Top_Features_Information.xlsx");
shape_cor_fat_top_indices = shape_cor_fat_top_feats.indices;
shape_cor_fat_top_names = shape_cor_fat_top_feats.names';

shape_cor_fat_train = datasets.shape_Coronal_CoronalTrainingShapeFeatures_Fat.feature_matrix;
shape_cor_fat_train = shape_cor_fat_train(:, shape_cor_fat_top_indices);
features_training = cat(2, features_training, shape_cor_fat_train);

shape_cor_fat_test = datasets.shape_Coronal_CoronalTestingShapeFeatures_Fat.feature_matrix;
shape_cor_fat_test = shape_cor_fat_test(:, shape_cor_fat_top_indices);
features_holdout = cat(2, features_holdout, shape_cor_fat_test);

feature_column_names = cat(2, feature_column_names, shape_cor_fat_top_names);

%% Axial Texture Tumor
texture_ax_tumor_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/15-Feb-2023_22_50_32_Axial_tumor_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_ax_tumor_top_indices = texture_ax_tumor_top_feats.indices;
texture_ax_tumor_top_names = texture_ax_tumor_top_feats.names';

texture_ax_tumor_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Tumor.feature_matrix;
texture_ax_tumor_train = texture_ax_tumor_train(:, texture_ax_tumor_top_indices);
features_training = cat(2, features_training, texture_ax_tumor_train);

texture_ax_tumor_test = datasets.texture_Axial_AxialTestingTextureFeatures_Tumor.feature_matrix;
texture_ax_tumor_test = texture_ax_tumor_test(:, texture_ax_tumor_top_indices);
features_holdout = cat(2, features_holdout, texture_ax_tumor_test);

feature_column_names = cat(2, feature_column_names, texture_ax_tumor_top_names);

%% Axial Texture Fat
texture_ax_fat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/16-Feb-2023_22_53_12_Axial_fat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_ax_fat_top_indices = texture_ax_fat_top_feats.indices;
texture_ax_fat_top_names = texture_ax_fat_top_feats.names';

texture_ax_fat_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Fat.feature_matrix;
texture_ax_fat_train = texture_ax_fat_train(:, texture_ax_fat_top_indices);
features_training = cat(2, features_training, texture_ax_fat_train);

texture_ax_fat_test = datasets.texture_Axial_AxialTestingTextureFeatures_Fat.feature_matrix;
texture_ax_fat_test = texture_ax_fat_test(:, texture_ax_fat_top_indices);
features_holdout = cat(2, features_holdout, texture_ax_fat_test);

feature_column_names = cat(2, feature_column_names, texture_ax_fat_top_names);

%% Axial Texture Proximal Fat
texture_ax_proxfat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/16-Feb-2023_22_51_54_Axial_proxfat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_ax_proxfat_top_indices = texture_ax_proxfat_top_feats.indices;
texture_ax_proxfat_top_names = texture_ax_proxfat_top_feats.names';

texture_ax_proxfat_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Proximal_Fat.feature_matrix;
texture_ax_proxfat_train = texture_ax_proxfat_train(:, texture_ax_proxfat_top_indices);
features_training = cat(2, features_training, texture_ax_proxfat_train);

texture_ax_proxfat_test = datasets.texture_Axial_AxialTestingTextureFeatures_Proximal_Fat.feature_matrix;
texture_ax_proxfat_test = texture_ax_proxfat_test(:, texture_ax_proxfat_top_indices);
features_holdout = cat(2, features_holdout, texture_ax_proxfat_test);

feature_column_names = cat(2, feature_column_names, texture_ax_proxfat_top_names);

%% Coronal Texture Fat
texture_cor_fat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/15-Feb-2023_22_47_26_Coronal_fat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_cor_fat_top_indices = texture_cor_fat_top_feats.indices;
texture_cor_fat_top_names = texture_cor_fat_top_feats.names';


texture_cor_fat_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Fat.feature_matrix;
texture_cor_fat_train = texture_cor_fat_train(:, texture_cor_fat_top_indices);
features_training = cat(2, features_training, texture_cor_fat_train);

texture_cor_fat_test = datasets.texture_Coronal_CoronalTestingTextureFeatures_Fat.feature_matrix;
texture_cor_fat_test = texture_cor_fat_test(:, texture_cor_fat_top_indices);
features_holdout = cat(2, features_holdout, texture_cor_fat_test);

feature_column_names = cat(2, feature_column_names, texture_cor_fat_top_names);

%% Coronal Texture Proximal Fat
texture_cor_proxfat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/16-Feb-2023_22_41_46_Coronal_proxfat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_cor_proxfat_top_indices = texture_cor_proxfat_top_feats.indices;
texture_cor_proxfat_top_names = texture_cor_proxfat_top_feats.names';


texture_cor_proxfat_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Proximal_Fat.feature_matrix;
texture_cor_proxfat_train = texture_cor_proxfat_train(:, texture_cor_proxfat_top_indices);
features_training = cat(2, features_training, texture_cor_proxfat_train);

texture_cor_proxfat_test = datasets.texture_Coronal_CoronalTestingTextureFeatures_Proximal_Fat.feature_matrix;
texture_cor_proxfat_test = texture_cor_proxfat_test(:, texture_cor_proxfat_top_indices);
features_holdout = cat(2, features_holdout, texture_cor_proxfat_test);

feature_column_names = cat(2, feature_column_names, texture_cor_proxfat_top_names);

%% Coronal Texture Tumor
texture_cor_tumor_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/15-Feb-2023_22_48_48_Coronal_tumor_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_cor_tumor_top_indices = texture_cor_tumor_top_feats.indices;
texture_cor_tumor_top_names = texture_cor_tumor_top_feats.names';


texture_cor_tumor_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Tumor.feature_matrix;
texture_cor_tumor_train = texture_cor_tumor_train(:, texture_cor_tumor_top_indices);
features_training = cat(2, features_training, texture_cor_tumor_train);

texture_cor_tumor_test = datasets.texture_Coronal_CoronalTestingTextureFeatures_Tumor.feature_matrix;
texture_cor_tumor_test = texture_cor_tumor_test(:, texture_cor_tumor_top_indices);
features_holdout = cat(2, features_holdout, texture_cor_tumor_test);

feature_column_names = cat(2, feature_column_names, texture_cor_tumor_top_names);

%% Create holdout test set

% features_holdout = cat(1, features_holdout, features_va);

%% Load in labels

data_labels_training = readmatrix(training_label_path);
data_labels_holdout = readmatrix(testing_label_path);

%% Create holdout testing dataset
holdout_test_size = size(features_holdout);
assert(holdout_test_size(1) == 14, "The size of the dataset is incorrect!");

fprintf("Loaded in train and holdout testing datasets! \n");

%% Whiten the data
features_training = whitenData(features_training, 'simple');
features_holdout = whitenData(features_holdout, 'simple');

nan_vector = features_training(1, :);
nan_indices = find(~isnan(nan_vector));
features_training = features_training(:, nan_indices);
features_holdout = features_holdout(:, nan_indices);
feature_column_names = feature_column_names(:, nan_indices);

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
    params.osname = 'SMOTE';
    % params.fsprunecorr = 'true';
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
                path_to_train = strcat(input_path, view{j}, '_', rois{i}, "_2D_", feature_family{1}, "_features.mat");
                datasets.(strcat(feature_family{1}, "_", view{j}, "_", data_type, "_", rois{i})) = load(path_to_train);
            end
        end
        
    end
    
    dataset_struct = datasets;

end