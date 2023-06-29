function [top_feat_inds] = top_feature_selection_ndim(data, labels, type, K, varargin)
% INPUTS:
%           data: M x N matrix with M observations and N variables
%         labels: N x 1 column vector with N class labels corresponding to each observation
%           type: (optional) string specifying which feature selection method to use ('mrmr','wilcoxn' DEFAULT='wilcoxn')
%              K: (optional) specifies the number of top feature indices to be returned (DEFAULT = 1);
%       varargin: (optional) any additional feature selection parameters
%
% OUTPUTS:
%  top_feat_inds: indices of the top features selected
% 
% Code written by Jacob Antunes @02-15-2018

if nargin < 4
    K = 1;
end

if nargin < 3
    type = 'wilcxon';
end

if isa(labels,'cell'), error('LABELS must be a non-cell array with labels 1 and -1.'); end

if ismember(type,{'mrmr','wilcoxn'})
        
    if strcmp(type,'wilcoxn')
        method = 'ranksum_ndim';
    end
    if strcmp(type,'mrmr')
        method = 'mrmr_mid_d';
    end
    
    top_feat_inds = feval(method,data, labels, K, varargin{:});
    
else
    error([type ' is an invalid feature selection method name']);
end
