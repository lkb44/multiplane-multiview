clear;close all;clc;
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/general/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature Classifier/nFold_cross_validation/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature Classifier/training_and_testing_sets/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/feature_selection/mrmr_feature_select/'));
addpath(genpath("/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/mrmr_feature_select/"));
addpath(genpath("/Users/Tom/Documents/scripts/feature_selection/mrmr_feature_select/"));
%% Specify dataset information
view = {'Axial', 'Coronal'}; % {Can only be 'Axial'}, {'Coronal'} or {'Axial','Coronal'}
rois = {'lumen','rw', 'fat'};

experiment_type = 'fat_only'; % https://www.mathworks.com/help/matlab/ref/nchoosek.html#d123e968841 nchoosek(rois, 2)

% Print experiment types
fprintf(strcat("USING THE FOLLOWING VIEW(S): ", strjoin(view), "\n"));
fprintf("EXPERIMENT TYPE: %s \n", experiment_type);

% Specify paths to ground truth labels
train_labels_path = "../Data/train_test_labels/train_labels.csv";
test_labels_path = "../Data/train_test_labels/test_labels.csv";
va_test_labels = "../Data/train_test_labels/va_test_labels.csv";

% Set directory path to the feature matrices depending on the view
if(length(view) == 1 && view{1} == "Axial")
    train_matrix_path = "../Shape_Features/Feature_Matrices/Axial/train/";
    test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/test/";
    va_test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/va_test/";
elseif(length(view) == 1 && view{1} == "Coronal")
    train_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/train/"; 
    test_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/test/";
    va_test_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/va_test/";
elseif(length(view) == 2)
    % Assert Axial is first in list and Coronal is second. This is
    % important for data loading and naming output logs. If view is not in
    % this format, the results will be confusing!!
    assert(strcmp(view{1}, 'Axial') && strcmp(view{2}, 'Coronal'), 'View Cell array must contain Axial first and Coronal second');
    ax_train_matrix_path = "../Shape_Features/Feature_Matrices/Axial/train/";
    ax_test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/test/";
    ax_va_test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/va_test/";
    cor_train_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/train/"; 
    cor_test_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/test/";
    cor_va_test_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/va_test/";
end

% Create cell array with names of each feature
feature_column_names = readtable("Shape_Feature_Col_Names.xlsx", 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

fprintf("Got directory paths to feature matrices! \n");

%% Create output directory

output_path = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_Feature_Results/";
experiment_date = strrep(string(datetime(now,'ConvertFrom','datenum')), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");

% If working with only 1 view
if(length(view) == 1)  
    output_path = strcat(output_path, experiment_date, "_", view{1}, "_", experiment_type, "/");
else % Working with both views
    output_path = strcat(output_path, experiment_date, "_", strjoin(view, '_'), "_", experiment_type, "/");
end

if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

%% Load in datasets

datasets = struct;

% If working with only 1 view
if(length(view) == 1)
    for i=1:length(rois)
        path_to_train = strcat(train_matrix_path, view{1}, "_", rois{i}, "_train_3D_shape_features.mat");
        datasets.(strcat("train_",rois{i})) = load(path_to_train);

        path_to_test = strcat(test_matrix_path, view{1}, "_", rois{i}, "_test_3D_shape_features.mat");
        datasets.(strcat("test_",rois{i})) = load(path_to_test);
        
        path_to_test = strcat(va_test_matrix_path, view{1}, "_", rois{i}, "_va_test_3D_shape_features.mat");
        datasets.(strcat("va_test_",rois{i})) = load(path_to_test);
    end
else % Working with both views
    for i=1:length(rois)
        
        % Axial datasets
        path_to_train = strcat(ax_train_matrix_path, view{1}, "_", rois{i}, "_train_3D_shape_features.mat");
        datasets.(strcat(view{1}, "_train_",rois{i})) = load(path_to_train);
        
        path_to_test = strcat(ax_test_matrix_path, view{1}, "_", rois{i}, "_test_3D_shape_features.mat");
        datasets.(strcat(view{1}, "_test_",rois{i})) = load(path_to_test);
        
        path_to_test = strcat(ax_va_test_matrix_path, view{1}, "_", rois{i}, "_va_test_3D_shape_features.mat");
        datasets.(strcat(view{1}, "_va_test_",rois{i})) = load(path_to_test);
        
        % Coronal Datasets
        path_to_train = strcat(cor_train_matrix_path, view{2}, "_", rois{i}, "_train_3D_shape_features.mat");
        datasets.(strcat(view{2}, "_train_",rois{i})) = load(path_to_train);
        
        path_to_test = strcat(cor_test_matrix_path, view{2}, "_", rois{i}, "_test_3D_shape_features.mat");
        datasets.(strcat(view{2}, "_test_",rois{i})) = load(path_to_test);
        
        path_to_test = strcat(cor_va_test_matrix_path, view{2}, "_", rois{i}, "_va_test_3D_shape_features.mat");
        datasets.(strcat(view{2}, "_va_test_",rois{i})) = load(path_to_test);
        
    end

end

%% Load in labels

data_labels_training = readmatrix(train_labels_path);
data_labels_training = data_labels_training(:, 2);

uh_data_labels_holdout = readmatrix(test_labels_path);
uh_data_labels_holdout = uh_data_labels_holdout(:, 2);

va_data_labels = readmatrix(va_test_labels);
va_data_labels = va_data_labels(:, 2);

data_labels_holdout = vertcat(uh_data_labels_holdout, va_data_labels);
%% Specify training and testing datasets

if(length(view) == 1 && experiment_type == "lumen_only") % Single view lumen features only
    
    features_training = datasets.train_lumen.feature_matrix;
    features_holdout = datasets.test_lumen.feature_matrix;
    features_va = datasets.va_test_lumen.feature_matrix;
    features_holdout = vertcat(features_holdout, features_va);
    feature_column_names = strcat(view{1}, "_lumen_", feature_column_names);
    
elseif(length(view) == 1 && experiment_type == "rw_only") % Single view rectal wall features only
    
    features_training = datasets.train_rw.feature_matrix;
    features_holdout = datasets.test_rw.feature_matrix;
    features_va = datasets.va_test_rw.feature_matrix;
    features_holdout = vertcat(features_holdout, features_va);
    feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
    
elseif(length(view) == 1 && experiment_type == "fat_only") % Single view fat features only
    
    features_training = datasets.train_fat.feature_matrix;
    features_holdout = datasets.test_fat.feature_matrix;
    features_va = datasets.va_test_fat.feature_matrix;
    features_holdout = vertcat(features_holdout, features_va);
    feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
    
elseif(length(view) == 1 && experiment_type == "lumen_rw") % Single view with lumen and rectal wall features
    
    lumen_train = datasets.train_lumen.feature_matrix;
    rw_train = datasets.train_rw.feature_matrix;
    features_training = [lumen_train rw_train];
    
    lumen_test = datasets.test_lumen.feature_matrix;
    rw_test = datasets.test_rw.feature_matrix;
    features_holdout = [lumen_test rw_test];
    
    va_lumen_test = datasets.va_test_lumen.feature_matrix;
    va_rw_test = datasets.va_test_rw.feature_matrix;
    features_va = [va_lumen_test va_rw_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    lumen_feature_column_names = strcat(view{1}, "_lumen_", feature_column_names);
    rw_feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
    feature_column_names = cat(2, lumen_feature_column_names, rw_feature_column_names);
    
elseif(length(view) == 1 && experiment_type == "lumen_fat") % Single view with lumen and fat features
    
    lumen_train = datasets.train_lumen.feature_matrix;
    fat_train = datasets.train_fat.feature_matrix;
    features_training = [lumen_train fat_train];
    
    lumen_test = datasets.test_lumen.feature_matrix;
    fat_test = datasets.test_fat.feature_matrix;
    features_holdout = [lumen_test fat_test];
    
    va_lumen_test = datasets.va_test_lumen.feature_matrix;
    va_fat_test = datasets.va_test_fat.feature_matrix;
    features_va = [va_lumen_test va_fat_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    lumen_feature_column_names = strcat(view{1}, "_lumen_", feature_column_names);    
    fat_feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
    feature_column_names = cat(2, lumen_feature_column_names, fat_feature_column_names);
    
elseif(length(view) == 1 && experiment_type == "rw_fat") % Single view with rectal wall and fat features
    
    rw_train = datasets.train_rw.feature_matrix;
    fat_train = datasets.train_fat.feature_matrix;
    features_training = [rw_train fat_train];
    
    rw_test = datasets.test_rw.feature_matrix;
    fat_test = datasets.test_fat.feature_matrix;
    features_holdout = [rw_test fat_test];
    
    va_rw_test = datasets.va_test_rw.feature_matrix;
    va_fat_test = datasets.va_test_fat.feature_matrix;
    features_va = [va_rw_test va_fat_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    rw_feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
    fat_feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
    feature_column_names = cat(2, rw_feature_column_names, fat_feature_column_names);
    
elseif(length(view) == 1 && experiment_type == "all") % Single view with all lumne, rectal wall, and fat features
    
   lumen_train = datasets.train_lumen.feature_matrix; 
   rw_train = datasets.train_rw.feature_matrix;
   fat_train = datasets.train_fat.feature_matrix;
   features_training = [lumen_train rw_train fat_train];
   
   lumen_test = datasets.test_lumen.feature_matrix;
   rw_test = datasets.test_rw.feature_matrix;
   fat_test = datasets.test_fat.feature_matrix;
   features_holdout = [lumen_test rw_test fat_test];
   
   va_lumen_test = datasets.va_test_lumen.feature_matrix;
   va_rw_test = datasets.va_test_rw.feature_matrix;
   va_fat_test = datasets.va_test_fat.feature_matrix;
   features_va = [va_lumen_test va_rw_test va_fat_test];
   features_holdout = vertcat(features_holdout, features_va);
   
   lumen_feature_column_names = strcat(view{1}, "lumen_", feature_column_names);  
   rw_feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
   fat_feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
   feature_column_names = cat(2, lumen_feature_column_names, rw_feature_column_names, fat_feature_column_names);
   
elseif(length(view) == 2 && experiment_type == "lumen_only") % Two views lumen features only
    
    ax_lumen_train = datasets.Axial_train_lumen.feature_matrix;
    cor_lumen_train = datasets.Coronal_train_lumen.feature_matrix;
    features_training = [ax_lumen_train cor_lumen_train];
    
    ax_lumen_test = datasets.Axial_test_lumen.feature_matrix;
    cor_lumen_test = datasets.Coronal_test_lumen.feature_matrix;
    features_holdout = [ax_lumen_test cor_lumen_test];
    
    ax_va_lumen_test = datasets.Axial_va_test_lumen.feature_matrix;
    cor_va_lumen_test = datasets.Coronal_va_test_lumen.feature_matrix;
    features_va = [ax_va_lumen_test cor_va_lumen_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    ax_lumen_feature_column_names = strcat(view{1}, "_lumen_", feature_column_names); 
    cor_lumen_feature_column_names = strcat(view{2}, "_lumen_", feature_column_names); 
    feature_column_names = cat(2, ax_lumen_feature_column_names, cor_lumen_feature_column_names);
    
elseif(length(view) == 2 && experiment_type == "rw_only") % Two views rectal wall features only
    
    ax_rw_train = datasets.Axial_train_rw.feature_matrix;
    cor_rw_train = datasets.Coronal_train_rw.feature_matrix;
    features_training = [ax_rw_train cor_rw_train];
    
    ax_rw_test = datasets.Axial_test_rw.feature_matrix;
    cor_rw_test = datasets.Axial_test_rw.feature_matrix;
    features_holdout = [ax_rw_test cor_rw_test];
    
    ax_va_rw_test = datasets.Axial_va_test_rw.feature_matrix;
    cor_va_rw_test = datasets.Coronal_va_test_rw.feature_matrix;
    features_va = [ax_va_rw_test cor_va_rw_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    ax_rw_feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
    cor_rw_feature_column_names = strcat(view{2}, "_rw_", feature_column_names);
    feature_column_names = cat(2, ax_rw_feature_column_names, cor_rw_feature_column_names);
    
elseif(length(view) == 2 && experiment_type == "fat_only") % Two views fat features only
    
    ax_fat_train = datasets.Axial_train_fat.feature_matrix;
    cor_fat_train = datasets.Coronal_train_fat.feature_matrix;
    features_training = [ax_fat_train cor_fat_train];
    
    ax_fat_test = datasets.Axial_test_fat.feature_matrix;
    cor_fat_test = datasets.Coronal_test_fat.feature_matrix;
    features_holdout = [ax_fat_test cor_fat_test];
    
    ax_va_fat_test = datasets.Axial_va_test_fat.feature_matrix;
    cor_va_fat_test = datasets.Coronal_va_test_fat.feature_matrix;
    features_va = [ax_va_fat_test cor_va_fat_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    ax_fat_feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
    cor_fat_feature_column_names = strcat(view{2}, "_fat_", feature_column_names);
    feature_column_names = cat(2, ax_fat_feature_column_names, cor_fat_feature_column_names);
    
elseif(length(view) == 2 && experiment_type == "lumen_rw") % Two views with lumen and rectal wall features
    
    ax_lumen_train = datasets.Axial_train_lumen.feature_matrix;
    cor_lumen_train = datasets.Coronal_train_lumen.feature_matrix;
    ax_rw_train = datasets.Axial_train_rw.feature_matrix;
    cor_rw_train = datasets.Coronal_train_rw.feature_matrix;
    features_training = [ax_lumen_train cor_lumen_train ax_rw_train cor_rw_train];
    
    ax_lumen_test = datasets.Axial_test_lumen.feature_matrix;
    cor_lumen_test = datasets.Coronal_test_lumen.feature_matrix;
    ax_rw_test = datasets.Axial_test_rw.feature_matrix;
    cor_rw_test = datasets.Axial_test_rw.feature_matrix;
    features_holdout = [ax_lumen_test cor_lumen_test ax_rw_test cor_rw_test];
    
    ax_va_lumen_test = datasets.Axial_va_test_lumen.feature_matrix;
    cor_va_lumen_test = datasets.Coronal_va_test_lumen.feature_matrix;
    ax_va_rw_test = datasets.Axial_va_test_rw.feature_matrix;
    cor_va_rw_test = datasets.Coronal_va_test_rw.feature_matrix;
    features_va = [ax_va_lumen_test cor_va_lumen_test ax_va_rw_test cor_va_rw_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    ax_lumen_feature_column_names = strcat(view{1}, "_lumen_", feature_column_names); 
    cor_lumen_feature_column_names = strcat(view{2}, "_lumen_", feature_column_names);
    ax_rw_feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
    cor_rw_feature_column_names = strcat(view{2}, "_rw_", feature_column_names);
    feature_column_names = cat(2, ax_lumen_feature_column_names, cor_lumen_feature_column_names, ax_rw_feature_column_names, cor_rw_feature_column_names);
    
elseif(length(view) == 2 && experiment_type == "lumen_fat") % Two views with lumen and fat features
    
    ax_lumen_train = datasets.Axial_train_lumen.feature_matrix;
    cor_lumen_train = datasets.Coronal_train_lumen.feature_matrix;
    ax_fat_train = datasets.Axial_train_fat.feature_matrix;
    cor_fat_train = datasets.Coronal_train_fat.feature_matrix;
    features_training = [ax_lumen_train cor_lumen_train ax_fat_train cor_fat_train];
    
    ax_lumen_test = datasets.Axial_test_lumen.feature_matrix;
    cor_lumen_test = datasets.Coronal_test_lumen.feature_matrix;
    ax_fat_test = datasets.Axial_test_fat.feature_matrix;
    cor_fat_test = datasets.Coronal_test_fat.feature_matrix;
    features_holdout = [ax_lumen_test cor_lumen_test ax_fat_test cor_fat_test];
    
    ax_va_lumen_test = datasets.Axial_va_test_lumen.feature_matrix;
    cor_va_lumen_test = datasets.Coronal_va_test_lumen.feature_matrix;
    ax_va_fat_test = datasets.Axial_va_test_fat.feature_matrix;
    cor_va_fat_test = datasets.Coronal_va_test_fat.feature_matrix;
    features_va = [ax_va_lumen_test cor_va_lumen_test ax_va_fat_test cor_va_fat_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    ax_lumen_feature_column_names = strcat(view{1}, "_lumen_", feature_column_names); 
    cor_lumen_feature_column_names = strcat(view{2}, "_lumen_", feature_column_names);
    ax_fat_feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
    cor_fat_feature_column_names = strcat(view{2}, "_fat_", feature_column_names);
    feature_column_names = cat(2, ax_lumen_feature_column_names,  cor_lumen_feature_column_names, ax_fat_feature_column_names, cor_fat_feature_column_names);
 
elseif(length(view) == 2 && experiment_type == "rw_fat") % Two views with rectal wall and fat features
    
    ax_rw_train = datasets.Axial_train_rw.feature_matrix;
    cor_rw_train = datasets.Coronal_train_rw.feature_matrix;
    ax_fat_train = datasets.Axial_train_fat.feature_matrix;
    cor_fat_train = datasets.Coronal_train_fat.feature_matrix;
    features_training = [ax_rw_train cor_rw_train ax_fat_train cor_fat_train];
    
    ax_rw_test = datasets.Axial_test_rw.feature_matrix;
    cor_rw_test = datasets.Axial_test_rw.feature_matrix;
    ax_fat_test = datasets.Axial_test_fat.feature_matrix;
    cor_fat_test = datasets.Coronal_test_fat.feature_matrix;
    features_holdout = [ax_rw_test cor_rw_test ax_fat_test cor_fat_test];

    ax_va_rw_test = datasets.Axial_va_test_rw.feature_matrix;
    cor_va_rw_test = datasets.Coronal_va_test_rw.feature_matrix;
    ax_va_fat_test = datasets.Axial_va_test_fat.feature_matrix;
    cor_va_fat_test = datasets.Coronal_va_test_fat.feature_matrix;
    features_va = [ax_va_rw_test cor_va_rw_test ax_va_fat_test cor_va_fat_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    ax_rw_feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
    cor_rw_feature_column_names = strcat(view{2}, "_rw_", feature_column_names);
    ax_fat_feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
    cor_fat_feature_column_names = strcat(view{2}, "_fat_", feature_column_names);
    feature_column_names = cat(2, ax_rw_feature_column_names, cor_rw_feature_column_names, ax_fat_feature_column_names, cor_fat_feature_column_names);
    
elseif(length(view) == 2 && experiment_type == "all") % Two views with all lumen, rectal wall, and fat features
    
    ax_lumen_train = datasets.Axial_train_lumen.feature_matrix;
    cor_lumen_train = datasets.Coronal_train_lumen.feature_matrix;
    ax_rw_train = datasets.Axial_train_rw.feature_matrix;
    cor_rw_train = datasets.Coronal_train_rw.feature_matrix;
    ax_fat_train = datasets.Axial_train_fat.feature_matrix;
    cor_fat_train = datasets.Coronal_train_fat.feature_matrix;
    features_training = [ax_lumen_train cor_lumen_train ax_rw_train cor_rw_train ax_fat_train cor_fat_train];
    
    ax_lumen_test = datasets.Axial_test_lumen.feature_matrix;
    cor_lumen_test = datasets.Coronal_test_lumen.feature_matrix;
    ax_rw_test = datasets.Axial_test_rw.feature_matrix;
    cor_rw_test = datasets.Axial_test_rw.feature_matrix;
    ax_fat_test = datasets.Axial_test_fat.feature_matrix;
    cor_fat_test = datasets.Coronal_test_fat.feature_matrix;
    features_holdout = [ax_lumen_test cor_lumen_test ax_rw_test cor_rw_test ax_fat_test cor_fat_test];
    
    ax_va_lumen_test = datasets.Axial_va_test_lumen.feature_matrix;
    cor_va_lumen_test = datasets.Coronal_va_test_lumen.feature_matrix;
    ax_va_rw_test = datasets.Axial_va_test_rw.feature_matrix;
    cor_va_rw_test = datasets.Coronal_va_test_rw.feature_matrix;
    ax_va_fat_test = datasets.Axial_va_test_fat.feature_matrix;
    cor_va_fat_test = datasets.Coronal_va_test_fat.feature_matrix;
    features_va = [ax_va_lumen_test cor_va_lumen_test ax_va_rw_test cor_va_rw_test ax_va_fat_test cor_va_fat_test];
    features_holdout = vertcat(features_holdout, features_va);
    
    ax_lumen_feature_column_names = strcat(view{1}, "_lumen_", feature_column_names); 
    cor_lumen_feature_column_names = strcat(view{2}, "_lumen_", feature_column_names);
    ax_rw_feature_column_names = strcat(view{1}, "_rw_", feature_column_names);
    cor_rw_feature_column_names = strcat(view{2}, "_rw_", feature_column_names);
    ax_fat_feature_column_names = strcat(view{1}, "_fat_", feature_column_names);
    cor_fat_feature_column_names = strcat(view{2}, "_fat_", feature_column_names);
    feature_column_names = cat(2, ax_lumen_feature_column_names, cor_lumen_feature_column_names, ax_rw_feature_column_names, cor_rw_feature_column_names,...
                               ax_fat_feature_column_names, cor_fat_feature_column_names);
    
end

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
    params.classifier='RANDOMFOREST';
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






