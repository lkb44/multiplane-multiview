function [features_training,data_labels_training] = getTrainingDataset(trainingDatasetType, features_path, labels_path, stage_grouping)
%GETTRAININGDATASET Load in the training dataset
%   Loads in the training data and corresponding training labels
%
%   Paramters:
%      getTrainingDataset: str
%           Training set to load
%       features_path: str
%           Path to the features. Features should be saved as a .mat file.
%   Returns:
%      features_training: matrix
%           mxn matrix of training features, where m is the number of
%           patients and n is the number of features
%      data_labels_training: matrix
%           mx1 matrix of training labels, where m is the number of
%           patients. The labels are binary.

data_labels_training = load(strcat(labels_path,stage_grouping,'-Data_Labels_Discovery.mat'));
data_labels_training = data_labels_training.stage_labels;

if(isequal(trainingDatasetType, "curvature only"))
    %%% Load in curvature features only
    features_training = load(strcat(features_path,'Discovery_Features.mat')); % Load in the mat file
    features_training = features_training.feats;
    
    fprintf("Successfully loaded curvature TRAINING features \n");
    
elseif(isequal(trainingDatasetType, "fd only"))
    %%% Load in fractal dimension features only
    features_training = load(strcat(features_path,'Training_FD.mat')); % Load in the mat file
    features_training = features_training.FD_feats;
    
    fprintf("Successfully loaded fractal dimension TRAINING features \n");
    
elseif(isequal(trainingDatasetType, "combined"))
    %%% Load and combine curvature and fractal dimension features
    features_discovery_curv = load(strcat(features_path,'Discovery_Features.mat')); % Load in the mat file
    features_discovery_curv = features_discovery_curv.feats;

    features_discovery_FD = load(strcat(features_path,'Training_FD.mat')); % Load in the mat file
    features_discovery_FD = features_discovery_FD.FD_feats;
    
    features_training = [features_discovery_curv features_discovery_FD];
    
    fprintf("Successfully combined curvature and fractal dimension TRAINING features \n");

elseif(isequal(trainingDatasetType, "SI_stats"))
    %%% Figure out shape index feature statistics only
    features_training = load(strcat(features_path,'Discovery_Features.mat')); % Load in the mat file
    features_training = features_training.feats;
    
elseif(isequal(trainingDatasetType, "lumen_si_histogram"))
    %%% Load in Lumen Shape Index Histogram features
    features_training = load(strcat(features_path,'Training_SI_Histogram_features.mat'));
    features_training = features_training.features;
    
    features_training = features_training(:,2:end); % Remove the 1st percentile because it has so many NaNs
    
    nan_rows = any(isnan(features_training),2); % Remove patients with NaNs
    features_training(nan_rows,:) = [];
    
    data_labels_training(nan_rows,:) = []; % Remove labels for patients with NaNs
    
    fprintf("Successfully loaded Lumen Shape Index Histogram TRAINING features \n");
    
elseif(isequal(trainingDatasetType, "er_lumen"))
    %%% Load in histogram features extracted from lumen across entire
    %%% rectum
    features_training = load(strcat(features_path,'Training_ER_Lumen_Histogram_features.mat')); % Load in the mat file
    features_training = features_training.features;
    
    fprintf("Successfully loaded ENTIRE RECTUM-LUMEN Histogram TRAINING features \n");
    
elseif(isequal(trainingDatasetType, "er_orw"))
    %%% Load in histogram features extracted from outer rectal wall across
    %%% entire rectum
    features_training = load(strcat(features_path,'Training_ER_Outer Rectal Wall_Histogram_features.mat')); % Load in the mat file
    features_training = features_training.features;
    
    fprintf("Successfully loaded ENTIRE RECTUM-ORW Histogram TRAINING features \n");
  
elseif(isequal(trainingDatasetType, "ts_lumen"))
    %%% Load in histogram features extracted from lumen across the tumor
    %%% subvolume
    features_training = load(strcat(features_path,'Training_TS_Lumen_Histogram_features.mat')); % Load in the mat file
    features_training = features_training.features;
    
    clear data_labels_training;
    data_labels_training = load(strcat(labels_path,stage_grouping,'-TS_Discovery.mat'));
    data_labels_training = data_labels_training.stage_labels;

    fprintf("Successfully loaded TUMOR SUBVOLUME-LUMEN Histogram TRAINING features \n");
    
elseif(isequal(trainingDatasetType, "ts_orw"))
    %%% Load in histogram features extracted from outer rectal wall across
    %%% the tumor subvolume
    features_training = load(strcat(features_path,'Training_TS_Outer Rectal Wall_Histogram_features.mat')); % Load in the mat file
    features_training = features_training.features;
    
    clear data_labels_training;
    data_labels_training = load(strcat(labels_path,stage_grouping,'-TS_Discovery.mat'));
    data_labels_training = data_labels_training.stage_labels;
    
    fprintf("Successfully loaded TUMOR SUBVOLUME-ORW Histogram TRAINING features \n");
    
else
    error("The type of training dataset(s) to load is not recognized!!");
end

end

