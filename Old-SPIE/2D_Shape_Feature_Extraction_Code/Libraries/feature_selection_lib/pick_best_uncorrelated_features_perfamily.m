function set_candiF = pick_best_uncorrelated_features_perfamily(data,classes,feat_family_idx,num_features_per_family,correlation_factor,correlation_metric)

% Eliminate correlated features from large feature set giving priority to
% significant features based on ranksum test
%
% Compare with: rankfeatures() in MATLAB using the CCWeighting parameter
%
% INPUTS:
%   data = MxN feature matrix with M observations and N variables
%   classes = Mx1 vector of class labels corresponding to observations in data
%   num_features_per_family = maximum size of uncorrelated feature set to be returned per feature family ([DEFAULT: the single best feature from each family]
%   feat_family_idx = indices to indicate start of each feature family (i.e. [1 30 65] would mean feature family A begins at IDX=1, feature family B begins at IDX=30, and feature family C begins at IDX=65))
%   correlation_factor = minimum correlation between features to trim [DEFAULT: 0.6]
%   correlation_metric = 'pearson' (parametric) or 'spearman' (nonparametric) [Default: 'spearman']
% 
% OUTPUTS:
%   set_candiF = vector of sortd feature indices of those features which were uncorrelated yet most discriminative
%
% Original Author: Jacob Antunes, 2018

%% CHECK INPUTS
if isa(classes,'cell'), error('CLASSES must be a non-cell array with labels 1 and -1.'); end
classes = double(classes);
classes(classes==0) = -1;
if any(~xor(classes == 1, classes == -1)), error('CLASSES must be 1 and -1'); end
if size(data,1)~=length(classes), error('CLASSES must contain same number of labels as observations in DATA'); end

%% DEFAULTS

if exist('feat_family_idx','var')~=1 || isempty(feat_family_idx)
    error('Must specify the feature_family_idx vector'); 
end
if exist('num_features_per_family','var')~=1 || isempty(num_features_per_family)
    num_features_per_family = 1; 
end
if exist('correlation_factor','var')~=1 || isempty(correlation_factor)
    correlation_factor=0.6;
end
if exist('correlation_metric','var')~=1 || isempty(correlation_metric)
    correlation_metric='Spearman';
end

%% Identify Best Features Per Family
set_candiF = [];

for i = 1:length(feat_family_idx)
    
    startz = feat_family_idx(i);
    if i<length(feat_family_idx)
        endz = feat_family_idx(i+1)-1;
    else
        endz = size(data,2);
    end

    F=pick_best_uncorrelated_features(data,classes,startz:endz,num_features_per_family,correlation_factor,correlation_metric);
    set_candiF = cat(2,set_candiF,F);
    
end

