function feature_ref = getReferenceTable(trainingDatasetType, ref_outpath)
%GETREFERENCETABLE Load in reference table for feature dataset.
%   For a feature dataset, there is a corresponding reference table. The
%   reference table contains the names of the columns (features). This is
%   used for getting the names of the top features after K-Fold Cross
%   Validation.
%
%   Parameters:
%      getTrainingDataset: str
%           Training set that is loaded. The name of the training set
%           corresponds to the name of the reference table.
%      ref_outpath: str
%           Path to the reference table.
%   Returns:
%       feature_ref: vector
%           Horizontal vector of feature names. If the feature matrix is
%           mxn, the horizontal vector of feature names should be 1xn,
%           where n is the number of features. m is the number of patients.

if(isequal(trainingDatasetType,"curvature only"))
    feature_ref = readtable(strcat(ref_outpath,'patients_col_ref.xlsx'));
    feature_ref = table2array(feature_ref);
    
elseif(isequal(trainingDatasetType,"fd only"))
    feature_ref = load(strcat(ref_outpath,'Training_FD_cols.mat'));
    feature_ref = feature_ref.FD_feat_cols;
    
elseif(isequal(trainingDatasetType,"combined"))
    feature_ref_curv = readtable(strcat(ref_outpath,'patients_col_ref.xlsx'));
    feature_ref_curv = table2array(feature_ref_curv);
    feature_ref_fd = load(strcat(ref_outpath,'Training_FD_cols.mat'));
    feature_ref_fd = feature_ref_fd.FD_feat_cols;
    feature_ref = [feature_ref_curv feature_ref_fd];
    
elseif(isequal(trainingDatasetType, "lumen_si_histogram"))
    feature_ref = readtable(strcat(ref_outpath,'Histogram_col_ref.xlsx'));
    feature_ref = table2array(feature_ref);
    feature_ref = feature_ref(:,2:end);
    
elseif(isequal(trainingDatasetType, "er_lumen") || isequal(trainingDatasetType, "er_orw") || isequal(trainingDatasetType, "ts_lumen") || isequal(trainingDatasetType, "ts_orw"))
    feature_ref = readtable(strcat(ref_outpath, 'Histogram_Descriptor_col_ref.xlsx'),'ReadVariableNames',false);
    feature_ref = table2array(feature_ref);
    
else
    error("Feature type is not recognized!");
end
end

