function [features_holdout,data_labels_holdout] = getHoldoutDataset(testingDatasetType, features_path, labels_path, stage_grouping)
%GETHOLDOUTDATASET Load in the holdout testing dataset.
%   Loads in the testing data and corresponding testing labels
%
%   Parameters:
%       testingDatasetType: str
%           Holdout testing set to load
%       features_path: str
%           Path to the features. Features should be saved as a .mat file.
%   Returns:
%       features_holdout: matrix
%           mxn matrix of testing features, where m is the number of
%           patients and n is the number of features
%       data_labels_holdout:
%           mx1 matrix of testing labels, where m is the number of
%           patients. The labels are binary.

%% CURVATURE FEATURES
if(isequal(testingDatasetType, "curvature original"))
    %%% Load in the original holdout testing dataset for curvature features
    features_holdout = load(strcat(features_path,'Holdout_Features.mat')); % Load in the mat file
    features_holdout = features_holdout.feats;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded curvature ORIGINAL TESTING features \n");
    
elseif(isequal(testingDatasetType, "curvature excluded"))
    %%% Load in the excluded dataset for curvature features
    features_holdout = load(strcat(features_path,'Excluded_Features.mat'));
    features_holdout = features_holdout.feats;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Excluded.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded curvature EXCLUDED TESTING features \n");
    
elseif(isequal(testingDatasetType, "curvature combined"))
    %%% Combined the original holdout and excluded datasets for curvature features
    features_holdout_org = load(strcat(features_path,'Holdout_Features.mat')); % Load in the mat file
    features_holdout_org = features_holdout_org.feats;
    
    features_holdout = [features_holdout_org; excluded_holdout];
    
    data_labels_org = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_org = data_labels_org.stage_labels;
    
    data_labels_ex = load(strcat(labels_path,stage_grouping,'-Data_Labels_Excluded.mat'));
    data_labels_ex = data_labels_ex.stage_labels;
    
    data_labels_holdout = [data_labels_org; data_labels_ex];
    
    fprintf("Successfully loaded curvature COMBINED TESTING features \n");

%% FRACTAL DIMENSION FEATURES
elseif(isequal(testingDatasetType, "fd original"))
    %%% Load in the original holdout testing dataset for fractal dimension features
    features_holdout = load(strcat(features_path,'Testing_FD.mat')); % Load in the mat file
    features_holdout = features_holdout.FD_feats;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded fractal dimension ORIGINAL TESTING features \n");
    
elseif(isequal(testingDatasetType, "fd excluded"))
    %%% Load in the excluded dataset for fractal dimension features
    features_holdout = load(strcat(features_path,'Excluded_FD.mat'));
    features_holdout = features_holdout.FD_feats;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Excluded.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded fractal dimension EXCLUDED TESTING features \n");
    
elseif(isequal(testingDatasetType, "fd combined"))
    %%% Combine the original holdout and excluded datasets for fractal
    %%% dimension features
    features_holdout_org = load(strcat(features_path,'Testing_FD.mat')); % Load in the mat file
    features_holdout_org = features_holdout_org.FD_feats;
    
    features_holdout = [features_holdout_org; excluded_holdout];
    
    data_labels_org = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_org = data_labels_org.stage_labels;
    
    data_labels_ex = load(strcat(labels_path,stage_grouping,'-Data_Labels_Excluded.mat'));
    data_labels_ex = data_labels_ex.stage_labels;
    
    data_labels_holdout = [data_labels_org; data_labels_ex];
    
    fprintf("Successfully loaded fractal dimension COMBINED TESTING features \n");
    
%% CURVATURE AND FRACTAL DIMENSION FEATURES
elseif(isequal(testingDatasetType, "combined original"))
    %%% Load in the original holdout test datasets for curvature features
    %%% and fractal dimension features
    features_holdout_curv = load(strcat(features_path,'Holdout_Features.mat')); % Load in the mat file
    features_holdout_curv = features_holdout_curv.feats;
    
    features_holdout_FD = load(strcat(features_path,'Testing_FD.mat')); % Load in the mat file
    features_holdout_FD = features_holdout_FD.FD_feats;
    
    features_holdout = [features_holdout_curv features_holdout_FD];
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded curvature + fractal dimension ORIGINAL TESTING features \n");
elseif(isequal(testingDatasetType, "combined excluded"))
    %%% Load in the excluded datasets for curvature features and fractal
    %%% dimension features
    features_holdout_curv = load(strcat(features_path,'Excluded_Features.mat'));
    features_holdout_curv = features_holdout_curv.feats;
    
    features_holdout_FD = load(strcat(features_path,'Excluded_FD.mat'));
    features_holdout_FD = features_holdout_FD.FD_feats;
    
    features_holdout = [features_holdout_curv features_holdout_FD];
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Excluded.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded curvature + fractal dimension EXCLUDED TESTING features \n");
    
elseif(isequal(testingDatasetType, "combined combined"))
    %%% Load in the original and excluded holdout testing datasets for
    %%% curvature features and fractal dimension features
    features_holdout_curv = load(strcat(features_path,'Holdout_Features.mat')); % Load in the mat file
    features_holdout_curv = features_holdout_curv.feats;
    
    excluded_holdout_curv = load(strcat(features_path,'Excluded_Features.mat'));
    excluded_holdout_curv = excluded_holdout_curv.feats;
    
    features_holdout_curv = [features_holdout_curv; excluded_holdout_curv];
    
    features_holdout_FD = load(strcat(features_path,'Testing_FD.mat')); % Load in the mat file
    features_holdout_FD = features_holdout_FD.FD_feats;
    
    excluded_holdout_FD = load(strcat(features_path,'Excluded_FD.mat'));
    excluded_holdout_FD = excluded_holdout_FD.FD_feats;
    
    features_holdout_FD = [features_holdout_FD; excluded_holdout_FD];
    
    features_holdout = [features_holdout_curv features_holdout_FD];
    
    data_labels_org = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_org = data_labels_org.stage_labels;
    
    data_labels_ex = load(strcat(labels_path,stage_grouping,'-Data_Labels_Excluded.mat'));
    data_labels_ex = data_labels_ex.stage_labels;
    
    data_labels_holdout = [data_labels_org; data_labels_ex];
    
    fprintf("Successfully loaded curvature + fractal dimension COMBINED TESTING features \n");

%% Histogram Features
elseif(isequal(testingDatasetType, "lumen_si_histogram original"))
    features_holdout = load(strcat(features_path,'Testing_SI_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;
    
    features_holdout = features_holdout(:,2:end); % Remove the 1st percentile because it has so many NaNs in training
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded Lumen Shape Index Histogram ORIGINAL TESTING features \n");
    
elseif(isequal(testingDatasetType, "er_lumen_original"))
     %%% Load in the original holdout testing dataset for histogram
     %%% features extracted from lumen across entire rectum
    features_holdout = load(strcat(features_path,'Testing_ER_Lumen_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded ENTIRE RECTUM-LUMEN ORIGINAL TESTING features \n");
    
elseif(isequal(testingDatasetType, "er_orw_original"))
     %%% Load in the original holdout testing dataset for histogram
     %%% features extracted from outer rectal wall across entire rectum
    features_holdout = load(strcat(features_path,'Testing_ER_Outer Rectal Wall_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded ENTIRE RECTUM-ORW ORIGINAL TESTING features \n");
    
elseif(isequal(testingDatasetType, "er_lumen_excluded"))
     %%% Load in the excluded holdout testing dataset for histogram
     %%% features extracted from lumen across entire rectum
    features_holdout = load(strcat(features_path,'Excluded_ER_Lumen_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded ENTIRE RECTUM-LUMEN EXCLUDED TESTING features \n");
    
elseif(isequal(testingDatasetType, "er_orw_excluded"))
     %%% Load in the excluded holdout testing dataset for histogram
     %%% features extracted from outer rectal wall across entire rectum
    features_holdout = load(strcat(features_path,'Excluded_ER_Outer Rectal Wall_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;

    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-Data_Labels_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;

    fprintf("Successfully loaded ENTIRE RECTUM-ORW EXCLUDED TESTING features \n");
    
 elseif(isequal(testingDatasetType, "ts_lumen_original"))
    %%% Load in the original holdout testing dataset for histogram
    %%% features extracted from lumen across tumor subvolume
    features_holdout = load(strcat(features_path,'Testing_TS_Lumen_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-TS_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded TUMOR SUBVOLUME-LUMEN ORIGINAL TESTING features \n");
    
elseif(isequal(testingDatasetType, "ts_orw_original"))
     %%% Load in the original holdout testing dataset for histogram
     %%% features extracted from outer rectal wall across tumor subvolume
    features_holdout = load(strcat(features_path,'Testing_TS_Outer Rectal Wall_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-TS_Holdout.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded TUMOR SUBVOLUME-ORW ORIGINAL TESTING features \n");
    
elseif(isequal(testingDatasetType, "ts_lumen_excluded"))
     %%% Load in the excluded holdout testing dataset for histogram
     %%% features extracted from lumen across tumor volume
    features_holdout = load(strcat(features_path,'Excluded_TS_Lumen_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;
    
    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-TS_Excluded.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;
    
    fprintf("Successfully loaded TUMOR SUBVOLUME-LUMEN EXCLUDED TESTING features \n");
    
elseif(isequal(testingDatasetType, "ts_orw_excluded"))
     %%% Load in the excluded holdout testing dataset for histogram
     %%% features extracted from outer rectal wall across tumor subvolume
    features_holdout = load(strcat(features_path,'Excluded_TS_Outer Rectal Wall_Histogram_features.mat')); % Load in the mat file
    features_holdout = features_holdout.features;

    data_labels_holdout = load(strcat(labels_path,stage_grouping,'-TS_Excluded.mat'));
    data_labels_holdout = data_labels_holdout.stage_labels;

    fprintf("Successfully loaded TUMOR SUBVOLUME-ORW EXCLUDED TESTING features \n");
    
end
end

