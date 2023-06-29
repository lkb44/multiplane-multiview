clear;close all;clc;
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/general/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature Classifier/nFold_cross_validation/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature Classifier/training_and_testing_sets/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/Feature_Classifier/feature_selection/mrmr_feature_select/'));
addpath(genpath("/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Feature_Selection_Classification_Code/mrmr_feature_select/"));
addpath(genpath("/Users/Tom/Documents/scripts/"));

%% Params
view = {'Axial'};
rois = {'tumor', 'fat'};

test_matrix_path = "../Texture_Features/2D_Features/Axial/test/";
va_test_matrix_path = "../Texture_Features/2D_Features/Axial/va_test/";

% Load datasets
for i=1:length(rois)

    path_to_test = strcat(test_matrix_path, view{1}, "_", rois{i}, "_test_2D_texture_features.mat");
    datasets.(strcat("test_",rois{i})) = load(path_to_test);

    path_to_test = strcat(va_test_matrix_path, view{1}, "_", rois{i}, "_va_test_2D_texture_features.mat");
    datasets.(strcat("va_test_",rois{i})) = load(path_to_test);
end
    
%% Load labels
test_labels_path = "../Data/train_test_labels/test_labels.csv";
va_test_labels = "../Data/train_test_labels/va_test_labels.csv";

uh_data_labels_holdout = readtable(test_labels_path);
uh_test_pts = uh_data_labels_holdout.Patient;
uh_data_labels_holdout = uh_data_labels_holdout.ypT;

va_data_labels = readtable(va_test_labels);
va_test_pts = va_data_labels.Patient;
va_data_labels = va_data_labels.ypT;

data_labels_holdout = vertcat(uh_data_labels_holdout, va_data_labels);
%test_pts = vertcat(uh_test_pts, va_test_pts);

%% Load Top Features
% UPDATE: need to adjust to load in Axial-texture-fat only
top_tumor_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_41_49_texture_Axial_tumor_only/Top_Features_Information.xlsx");
top_fat_feats = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_40_05_texture_Axial_fat_only/Top_Features_Information.xlsx");

%% Read in dataset
% UPDATE: need to create variables for fat
uh_tumor_top_feats = datasets.test_tumor.feature_matrix;
uh_tumor_top_feats = uh_tumor_top_feats(:, top_tumor_feats.indices);

va_tumor_top_feats = datasets.va_test_tumor.feature_matrix;
va_tumor_top_feats = va_tumor_top_feats(:, top_tumor_feats.indices);

uh_fat_top_feats = datasets.test_fat.feature_matrix;
uh_fat_top_feats = uh_fat_top_feats(:, top_fat_feats.indices);

va_fat_top_feats = datasets.va_test_fat.feature_matrix;
va_fat_top_feats = va_fat_top_feats(:, top_fat_feats.indices);

%all_top_feats = vertcat(uh_tumor_top_feats, va_tumor_top_feats);

%% Load texture features
% Load per-pixel texture features
% The texture feature array is a 1x91 cell array
% where each column represents a patient. Each patient cell array is a
% mx164 cell array, where m is the number of non-zero pixels in the mask
% (i.e. pixels corresponding to ground truth label of ROI) and 164 is the
% total number of 2D texture features extracted per slice.
% UPDATE: need to load in fat per-pixel texture features
% (feat_preUH_AxialFat.mat). 
root_path = '/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Data/Axial_resampled/';
tumor_texture_feats = load("../Texture_Features/2D_Features/Axial/All_Features/feat_preUH_AxialTumor.mat");
tumor_texture_feats = tumor_texture_feats.feat;

fat_texture_feats = load("../Texture_Features/2D_Features/Axial/All_Features/feat_preUH_AxialFat.mat");
fat_texture_feats = fat_texture_feats.feat;

% Load columns IDs. The matrix of per-pixel texture features is NOT split
% between training and testing. All patients are listed in this matrix as
% rows. This mat file indicates which rows correspond to testing patients
% in the texture_feats matrix
texture_feats_ids = load("../Texture_Features/2D_Features/Axial/All_Features/testId.mat");
texture_feats_ids = texture_feats_ids.i2;

% Load texture feature names. This is a 4x164 cell array, where 4 is the
% number of statistics computed for each texture feature and 164 is the
% number of texture features. The first row corresponds to the median
% statistic, second row is the variance statistic, third row is the
% kurtosis statistic, and fourth row is the skewness statistic.
texture_feat_names = load("../Texture_Features/2D_Features/Axial/All_Features/UHfeatStatsName.mat");
texture_feat_names = texture_feat_names.statnames;

% Filter texture feature array to the testing patients only
test_tumor_texture_fests = tumor_texture_feats(:, texture_feats_ids);
test_fat_texture_feats = fat_texture_feats(:, texture_feats_ids);

% Index the top_feats table for the feature that we want to plot on a
% heatmap.
haralick_info2_ws7 = top_tumor_feats(1,:);
Laws_R5L5 = top_fat_feats(1,:);

% Get the feature name. We need the feature name to identify which column
% it is in the texture_feat_names cell array. This will also be the column
% of the feature in the texture_feats array
haralick_info2_ws7_name = string(haralick_info2_ws7.names);
haralick_info2_ws7_name = strsplit(haralick_info2_ws7_name, "_"); % I added some extra stuff to the string. Remove it to get the name of the feature.
haralick_info2_ws7_name = haralick_info2_ws7_name(:, 3);

findIndex = cellfun(@(x)isequal(x,haralick_info2_ws7_name),texture_feat_names); % Returns a logical array, where 1 corresponds to the index of the feature name
[haralick_info2_ws7_row,haralick_info2_ws7_index_col] = find(findIndex); % Get the row and column of the feature name

Laws_R5L5_name = string(Laws_R5L5.names);
Laws_R5L5_name = strsplit(Laws_R5L5_name, "_");
Laws_R5L5_name = Laws_R5L5_name(:, 3);

findIndex = cellfun(@(x)isequal(x,Laws_R5L5_name),texture_feat_names);
[Laws_R5L5_row, Laws_R5L5_col] = find(findIndex);


%% Plot Heatmaps
% We loop through each patient in the holdout testing dataset. We are
% plotting the radiomic heatmap on the slice of largest tumor area. NOTE,
% UPDATE FOR FAT
% 
for i=1:length(uh_test_pts)
    
    % Specify paths to image and mask
    image_in = strcat(root_path, uh_test_pts{i}, '/vol/vol_largest_slice.mha');
    mask_in = strcat(root_path, uh_test_pts{i}, '/masks/mask_tumor_largest_slice.mha');
    %mask_in = strcat(root_path, uh_test_pts{i}, '/masks/mask_fat_largest_slice.mha');
    
    % Load in the image and mask
    vol = vv(image_in);
    mask = vv(mask_in);
    
    % MATLAB loads images in the wrong orientation. Fix it.
    vol = rot90(vol, 3);
    mask = rot90(mask, 3);
    
    pt_label = uh_data_labels_holdout(i,:); % Get the responder/non-responder label for this patient
    
    pt_feature_matrix = test_tumor_texture_fests{:,i};  % Get the column in the feature matrix that corresponds to this patient
    %pt_feature_matrix = test_fat_texture_feats{:,i};
    
    pt_feature_matrix = pt_feature_matrix{1,1}; % The texture feature matrix for this patient is in a cell array in the first row/first column.
    
    pt_feature_matrix = pt_feature_matrix(:, haralick_info2_ws7_index_col); % Extract the column that corresponds to the desired feature.
    %pt_feature_matrix = pt_feature_matrix(:, Laws_R5L5_col);
    
    % Plot the heatmap
    figure();
    ploth = feature_map_slice(vol,mask,pt_feature_matrix);
    colorbar;
    
    % Create title for the heatmap. Title is the patient name followed by
    % their label of responder/non-responder.
    pt_name = strrep(uh_test_pts{i}, "-"," ");
    if(pt_label == 0)
        label = "responder";
    elseif(pt_label == 1)
        label = "non-responder";
    else
        disp("label not recognized!");
    end
    
    tit = strcat(pt_name, " ", label);
    title(tit);
    
    % Save the figure
    % UPDATE: save to fat-specific folder
    out_path = "../Abstract/FIgures/Texture_Heatmaps/Tumor/";
    %out_path = "../Abstract/FIgures/Texture_Heatmaps/Fat/";
    output_name = strcat(out_path, tit, "_heatmap.png");
    saveas(gcf, output_name);
    
end
