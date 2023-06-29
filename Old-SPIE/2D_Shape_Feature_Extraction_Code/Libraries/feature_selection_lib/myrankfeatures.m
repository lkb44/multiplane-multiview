function [fea, p, h] = myrankfeatures(data,labels,K,varargin)
% INPUTS
%      data: a N*M matrix, indicating N samples, each having M dimensions.
%    labels: a N*1 matrix (vector), indicating the class/category of the N samples. Must be categorical.
%         K: (optional) the number of top features need to be returned
%  varargin: any additional parameters to wilcoxn
%
% OUTPUTS
%     fea: (if K specified): indices of the top K features based on rank (in descending order); otherwise, just sequential indices
%       p: p-values for significance testing for each feature
%       h: whether or not the null hypothesis was rejected (1) or failed to reject (0) for each feature

% Code written by Jacob Antunes @02-15-2018

if nargin<3
    K = [];
end

f = size(data,2); %# of features
    p = zeros(1,f);
    h = zeros(1,f);

for i = 1:f
    idx2rem = isnan(data(:,i));
    tempdata = data(~idx2rem,i);
    templabels = logical(labels(~idx2rem));
    [p(i),h(i),~]=ranksum(tempdata(templabels), tempdata(~templabels), varargin{:}); 
end

fea = 1:f; %sequential indices of all features considered
if ~isempty(K) %SORT!
    [p_sorted, fea_sorted] = sort(p,'ascend');
    h_sorted = h(fea_sorted);
    
    fea = fea_sorted(1:K);
    p = p_sorted(1:K);
    h = h_sorted(1:K);
end

