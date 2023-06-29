function [methodstring,outputs] = SVM(training_set, testing_set, training_labels, testing_labels, varargin)

%function [methodstring,outputs.tp,outputs.tn,outputs.fp,outputs.fn,outputs.prob_estimates] = SVM( training_set , testing_set, training_labels, testing_labels,params)

%%sample data
% data_set = [ 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 5 4 3 2 1; 5 4 3 2 1; 5 4 3 2 1];
% data_labels = [ -1 -1 -1 1 1 1]';

% Runs a demo on a given data set.
% This is an example on what should be done to train and predict on a data set. The information
% is not stored externally but it prints important information on screen. It gives a quick idea
% on how the SVM is doing.
% The data set is divided internally in the training set and testing set. Options are provided to
% select the training set randomly or just select a specific number of postive and negative examples
% for training.
%
% Parameters
%   data_set    list of eamples arranged by rows.
%       - N observations x P features
%   data_labels list of labels corresponding to the examples.
%       - labels [-1 -1 -1  1  1  1]'

% The function was originally written by George Lee @2009
% Function has been last modified by Jacob Antunes @2-15-2018
%
training_labels = training_labels(:);
testing_labels = testing_labels(:);

% if size(data_labels,1) == 1, data_labels = data_labels'; end
if nargin > 4
    kernel = varargin{1};
    params.kernel = lower(kernel);
    switch params.kernel
        case {'linear','quadratic','rbf','poly','mlp'}
        otherwise 
            fprintf('Invalid kernel selected.  Using linear kernel by default.\n');
            params.kernel = 'linear';
    end
else
    fprintf('No kernel selected.  Using linear kernel by default.\n');
    params.kernel = 'linear';
end

methodstring = ['SVM classifier (' params.kernel ' kernel) '];

outputs = svm_main(training_set, training_labels, testing_set, testing_labels,params);