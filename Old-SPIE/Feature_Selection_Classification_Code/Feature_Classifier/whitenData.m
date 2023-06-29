function [whitenedData,mean_vec,mad_vec] = whitenData(data, whiteType, varargin)
%WHITENDATA Whiten a dataset using simple whitening or modified whitening.
%   Data can be whitened independently or with respect to another dataset.
%   To whiten input data with respect to another dataset, you must input the
%   mean values and the MAD (mean absolute deviation) values.
%
%   Parameters:
%       data: matrix
%           Matrix of raw/unwhitened features
%       whiteType: str
%           Type of whitening to do. The following inputs are acceptable
%           inputs:
%               - simple
%               - modified
%       varargin: vectors
%           Additional arguments. Pass the mean vector and MAD vector from
%           another dataset to whiten the input dataset with respect to the
%           other dataset.
%   Returns:
%       whitenedData: matrix
%           Matrix of whitened features
%       mean_vec: vector
%           Vector of mean values. Each column corresponds to the mean of
%           a feature in whitenedData.
%       mad_vec: vector
%           Vector of MAD values. Each column corresponds to the MAD of 
%           a feature in whitenedData.
%% Initialize
addpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/pre-processing');

%% Whiten the data

% Whiten with respect to another dataset
if nargin>3
    input_mean_vec = varargin{1,1}; % Get input means
    input_mad_vec = varargin{1,2}; % Get input MADs
    fprintf("Received input mean values and MAD values\n")
    if(isequal(whiteType,'simple'))
        [whitenedData,mean_vec,mad_vec] = simplewhiten(data,input_mean_vec,input_mad_vec);
    elseif(isequal(whiteType,'modified'))
        [whitenedData,mean_vec,mad_vec] = simplewhitennooutlier(data,input_mean_vec,input_mad_vec);
    else
        error('Whitening type not recognized!')
    end
    fprintf("Data whitened successfully!\n");
    
% Uh-oh! Can't whiten data independently or with respect to another dataset
elseif nargin>2
    error('Error Missing mean_vec or mad_vec. Must path both');

% Whiten data independently
else
    if(isequal(whiteType,'simple'))
        [whitenedData,mean_vec,mad_vec] = simplewhiten(data);
    elseif(isequal(whiteType,'modified'))
        [whitenedData,mean_vec,mad_vec] = simplewhitennooutlier(data);
    else
        error('Whitening type not recognized!')
    end
    fprintf("Data whitened successfully!\n");
end


end

