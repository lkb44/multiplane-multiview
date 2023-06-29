function [training, testing] = nFold_multiclass(data_set, data_labels, shuffle, n)

% INPUTS
%      data_set: (unused, so just make empty [])
%   data_labels: M x 1 cell array of class labels
%       shuffle: 0 (unrandomized) or 1 (randomized) selection of training/testing groups per fold
%             n: number of folds
%
% OUTPUTS
%      training: 1 x n cell array of training observation indices used in each fold
%       testing: 1 x n cell array of testing observation indices used in each fold
    
    u = unique(data_labels);
    for i = 1:length(data_labels)
        dl(i) = find(contains(u,data_labels{i}));
    end
    
data_labels = dl';

% 1. First, we acquire the training set, training labels, testing set and testing labels.
%    For this, we will divide our data set in two. We will find positive (a)
%    and negative (b) examples to create a balanced training set.


%***HARDCODED FOR 3 CLASSES***%
a = find( data_labels == 1);
b = find( data_labels == 2);
c = find( data_labels == 3);
    
testing=cell(1,n);
training=cell(1,n);

if n==length(data_labels)
    x = 1:n;
    
    for i = 1:n
    d = setdiff(x,i);
    testing{i} = i;
    training{i} = d;
    end 
    
else

    if shuffle
        % commit a random portion of the dataset for training
        a_shuffle = randperm(length(a));    %randomize index
        b_shuffle = randperm(length(b));    %randomize index
        c_shuffle = randperm(length(c));    %randomize index
    else
        % or don't shuffle
        a_shuffle = 1:length(a);    %same index
        b_shuffle = 1:length(b);    %same index
        c_shuffle = 1:length(c);    %same index
    end

    % define n sets
    a_cuts = [0 round((1:n-1)/n*length(a)) size(a,1)];
    b_cuts = [0 round((1:n-1)/n*length(b)) size(b,1)];
    c_cuts = [0 round((1:n-1)/n*length(b)) size(c,1)];

    for i=1:n
        a_ind = a_shuffle(a_cuts(i)+1:a_cuts(i+1));
        b_ind = b_shuffle(b_cuts(i)+1:b_cuts(i+1));
        c_ind = c_shuffle(c_cuts(i)+1:c_cuts(i+1));
        
        a_values = a(a_ind);
        b_values = b(b_ind);
        c_values = c(c_ind);

        a_notvalues = a;
        b_notvalues = b;
        c_notvalues = c;
        
        a_notvalues(a_ind) = [];
        b_notvalues(b_ind) = [];
        c_notvalues(c_ind) = [];

        values = [a_values ; b_values; c_values]; %indices for testing set
        notvalues = [a_notvalues ; b_notvalues; c_notvalues]; %indices for training set


        % set testing set and labels
    %     testing_set{i} = data_set(values,:);
        testing{i} = values;
    %     testing_labels{i} = data_labels(values);

        % set training set as all samples not included in testing set
    %     training_set{i} = data_set;
        training{i} = notvalues;
    %     training_labels{i} = data_labels;

    %     temp = training_set{i};
    %     temp(values,:) = [];
    %     training_set{i} = temp;
    %     
    %     temp = training_labels{i};
    %     temp(values) = [];
    %     training_labels{i} = temp;
    end
    
end