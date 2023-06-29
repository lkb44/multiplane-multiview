clear
close all
clc

addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/evaluation'));

view = {'Coronal', 'Axial'};
rois = {'Proximal_Fat', 'Tumor', 'Fat'};
feature_family = {'Texture'};
scheme = {'wilcoxon_qda'};
extras = {'singleTopAUCFromEveryROI'};
dataset = "TrainingTextureFeatures/";
dimension = 2;
top30features = [];

texture_axial_rois = {'Proximal_Fat', 'Tumor', 'Fat'};
texture_cor_rois = {'Proximal_Fat', 'Tumor', 'Fat'};

experiment_type = 'all';

fprintf(strcat("USING THE FOLLOWING FEAUTRES: ", strjoin(feature_family), "\n"));
fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf("EXPERIMENT TYPE: %s \n", experiment_type);

%% Specify filepaths
% root_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/";
% 
% input_path = strcat(root_path, "Axial", dataset, roi, '/'); 
% output_path = strcat(root_path, "Axial", dataset, '/', "Axial", '_',roi, "_2D_texture_features.mat");
% 
% id_matrix_output_path = strcat(root_path, view, dataset);

%% Get patient IDs
% patients = dir(fullfile(input_path, strcat('*RectalCA_*')));
% patient_ids = {patients.name};
% % Axial 2D Texture Features
ax_train_texture_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/";

% Coronal 2D Texture Features
cor_train_texture_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/";

output_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/Stability/";
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

%% Specify paths to Feature matrices
% Specify paths to ground truth labels
training_label_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/train_test_labels/train_labels.csv";

datasets = struct;

% Axial Texture
datasets = load_data({'texture'}, dimension, {'Axial'}, texture_axial_rois, ax_train_texture_path, datasets);

% Coronal Texture
datasets = load_data({'texture'}, dimension, {'Coronal'}, texture_cor_rois, cor_train_texture_path, datasets);

%% Initialize training datasets
features_axial = [];
features_coronal = [];
feature_column_names_axial = [];
feature_column_names_coronal = [];

%% Axial Top Features
texture_ax_tumor_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/15-Feb-2023_22_48_48_Coronal_tumor_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_ax_tumor_top_indices = texture_ax_tumor_top_feats.indices;
texture_ax_tumor_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Tumor.feature_matrix;
texture_ax_tumor_train = texture_ax_tumor_train(:, texture_ax_tumor_top_indices);
features_axial = texture_ax_tumor_train;
top_indices = texture_ax_tumor_top_indices;

corresponding_coronal_tumor_indices = texture_ax_tumor_top_indices;
corresponding_texture_cor_tumor_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Tumor.feature_matrix;
corresponding_texture_cor_tumor_train = corresponding_texture_cor_tumor_train(:, corresponding_coronal_tumor_indices);
corresponding_features_axial = corresponding_texture_cor_tumor_train;

texture_ax_fat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/15-Feb-2023_22_47_26_Coronal_fat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_ax_fat_top_indices = texture_ax_fat_top_feats.indices;
texture_ax_fat_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Fat.feature_matrix;
texture_ax_fat_train = texture_ax_fat_train(:, texture_ax_fat_top_indices);
features_axial = cat(2, features_axial, texture_ax_fat_train);
top_indices = cat(2, top_indices, texture_ax_fat_top_indices);

corresponding_coronal_fat_indices = texture_ax_fat_top_indices;
corresponding_texture_cor_fat_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Fat.feature_matrix;
corresponding_texture_cor_fat_train = corresponding_texture_cor_fat_train(:, corresponding_coronal_fat_indices);
corresponding_features_axial = cat(2, corresponding_features_axial, corresponding_texture_cor_fat_train);

texture_ax_proxfat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/16-Feb-2023_23_05_52_Coronal_proxfat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
texture_ax_proxfat_top_indices = texture_ax_proxfat_top_feats.indices;
texture_ax_proxfat_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Proximal_Fat.feature_matrix;
texture_ax_proxfat_train = texture_ax_proxfat_train(:, texture_ax_proxfat_top_indices);
features_axial = cat(2, features_axial, texture_ax_proxfat_train);
top_indices = cat(2, top_indices, texture_ax_proxfat_top_indices);

corresponding_coronal_proxfat_indices = texture_ax_proxfat_top_indices;
corresponding_texture_cor_proxfat_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Proximal_Fat.feature_matrix;
corresponding_texture_cor_proxfat_train = corresponding_texture_cor_proxfat_train(:, corresponding_coronal_proxfat_indices);
corresponding_features_axial = cat(2, corresponding_features_axial, corresponding_texture_cor_proxfat_train);


% %% Axial Texture Tumor
% texture_ax_tumor_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/15-Feb-2023_22_50_32_Axial_tumor_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
% texture_ax_tumor_top_indices = texture_ax_tumor_top_feats.indices;
% % texture_ax_tumor_top_indices = texture_ax_tumor_top_indices(1);
% texture_ax_tumor_top_names = texture_ax_tumor_top_feats.names';
% % texture_ax_tumor_top_names = texture_ax_tumor_top_names(1);
% 
% texture_ax_tumor_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Tumor.feature_matrix;
% texture_ax_tumor_train = texture_ax_tumor_train(:, texture_ax_tumor_top_indices);
% features_axial = texture_ax_tumor_train;
% 
% feature_column_names_axial = texture_ax_tumor_top_names;

% %% Axial Texture Fat
% texture_ax_fat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/16-Feb-2023_22_53_12_Axial_fat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
% texture_ax_fat_top_indices = texture_ax_fat_top_feats.indices;
% % texture_ax_fat_top_indices = texture_ax_fat_top_indices(1);
% texture_ax_fat_top_names = texture_ax_fat_top_feats.names';
% % texture_ax_fat_top_names = texture_ax_fat_top_names(1);
% 
% texture_ax_fat_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Fat.feature_matrix;
% texture_ax_fat_train = texture_ax_fat_train(:, texture_ax_fat_top_indices);
% features_axial = cat(2, features_axial, texture_ax_fat_train);
% 
% feature_column_names_axial = cat(2, feature_column_names_axial, texture_ax_fat_top_names);

% %% Axial Texture Proximal Fat
% texture_ax_proxfat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/16-Feb-2023_22_51_54_Axial_proxfat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
% texture_ax_proxfat_top_indices = texture_ax_proxfat_top_feats.indices;
% % texture_ax_proxfat_top_indices = texture_ax_proxfat_top_indices(1);
% texture_ax_proxfat_top_names = texture_ax_proxfat_top_feats.names';
% % texture_ax_proxfat_top_names = texture_ax_proxfat_top_names(1);
% 
% texture_ax_proxfat_train = datasets.texture_Axial_AxialTrainingTextureFeatures_Proximal_Fat.feature_matrix;
% texture_ax_proxfat_train = texture_ax_proxfat_train(:, texture_ax_proxfat_top_indices);
% features_axial = cat(2, features_axial, texture_ax_proxfat_train);
% 
% feature_column_names_axial = cat(2, feature_column_names_axial, texture_ax_proxfat_top_names);

% %% Coronal Texture Tumor
% texture_cor_tumor_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/15-Feb-2023_22_48_48_Coronal_tumor_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
% texture_cor_tumor_top_indices = texture_cor_tumor_top_feats.indices;
% % texture_cor_tumor_top_indices = texture_cor_tumor_top_indices(1);
% texture_cor_tumor_top_names = texture_cor_tumor_top_feats.names';
% % texture_cor_tumor_top_names = texture_cor_tumor_top_names(1);
% 
% texture_cor_tumor_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Tumor.feature_matrix;
% texture_cor_tumor_train = texture_cor_tumor_train(:, texture_cor_tumor_top_indices);
% features_coronal = texture_cor_tumor_train;
% 
% feature_column_names_coronal = texture_cor_tumor_top_names;

% %% Coronal Texture Fat
% texture_cor_fat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/15-Feb-2023_22_47_26_Coronal_fat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
% texture_cor_fat_top_indices = texture_cor_fat_top_feats.indices;
% % texture_cor_fat_top_indices = texture_cor_fat_top_indices(1);
% texture_cor_fat_top_names = texture_cor_fat_top_feats.names';
% % texture_cor_fat_top_names = texture_cor_fat_top_names(1);
% 
% 
% texture_cor_fat_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Fat.feature_matrix;
% texture_cor_fat_train = texture_cor_fat_train(:, texture_cor_fat_top_indices);
% features_coronal = cat(2, features_coronal, texture_cor_fat_train);
% 
% feature_column_names_coronal = cat(2, features_coronal, texture_cor_fat_top_names);

% %% Coronal Texture Proximal Fat
% texture_cor_proxfat_top_feats = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/16-Feb-2023_23_05_52_Coronal_proxfat_only_wilcoxon_qda_withPruningAndOversampling/Top_Features_Information.xlsx");
% texture_cor_proxfat_top_indices = texture_cor_proxfat_top_feats.indices;
% % texture_cor_proxfat_top_indices = texture_cor_proxfat_top_indices(1);
% texture_cor_proxfat_top_names = texture_cor_proxfat_top_feats.names';
% % texture_cor_proxfat_top_names = texture_cor_proxfat_top_names(1);
% 
% 
% texture_cor_proxfat_train = datasets.texture_Coronal_CoronalTrainingTextureFeatures_Proximal_Fat.feature_matrix;
% texture_cor_proxfat_train = texture_cor_proxfat_train(:, texture_cor_proxfat_top_indices);
% features_coronal = cat(2, features_coronal, texture_cor_proxfat_train);
% 
% feature_column_names_coronal = cat(2, feature_column_names_coronal, texture_cor_proxfat_top_names);

%% Load in labels

feature_cell = cell(2, 1);
feature_cell{1} = features_axial;
feature_cell{2} = corresponding_features_axial;

data_labels_training = readmatrix(training_label_path);

[interDifScore intraDifScore, H] = measureStability(feature_cell, 100, output_path, 'AxialStability.mat', 0.05);

disp(interDifScore);
disp(intraDifScore);





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