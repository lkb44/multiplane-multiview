function [out_features, out_labels] = UNDERSAMPLE(in_features, in_labels)

% This function undersamples the majority class to imrpove class balance
%
% INPUTS:
%----------
%in_features:
%(N \times P) matrix of numerical features. each row is one example, each
%column is one feature, hence there are N examples with P features each.
%
%in_labels:
%boolean N-vector of labels, defining the classes to which the examples in
%in_features belong.
%

% OUTPUTS:
%----------
%out_features, out_labels:
%features and labels of ONLY the undersampled data.
%note that each entry of out_labels is the label of the majority class
%since only examples of the majority class are sampled.
%concatenating [in_features out_features] and [in_labels out_labels]
%gives a new example set with the desired class balance.

mkNeg = 0; %handle -1s as 0s
if ~all(in_labels==0 | in_labels==1)
    if any(in_labels==-1)
        in_labels(in_labels==-1) = 0;
        mkNeg = 1;
    else
        error('UNDERSAMPLE: in_labels may contain only the values 0 and 1.');
    end
end

numZeros = sum(in_labels==0);
numOnes  = sum(in_labels==1);

if numOnes == numZeros
    %nothing needs to be done because if classes are already balanced, then
    %for any in_beta, this is already the desired result.
    out_features = [];
    out_labels   = [];
    return;
else
    if numZeros > numOnes
        majLabel = logical(0);
        minLabel = logical(1);
    else
        majLabel = logical(1);
        minLabel = logical(0);
    end
end

%rename:
S = in_features;
clear in_features;

%feature sets by class:
Smin = S(in_labels==minLabel,:);
Smaj = S(in_labels==majLabel,:);

%handle boundary cases:
if size(Smin,1)==0
    warning('UNDERSAMPLE: there were no examples of the minority class in the data. hence balancing is not possible. Returning empty matrices.');
    out_features = [];
    out_labels   = [];
    return;
end

if size(Smin,1)==1
    warning('UNDERSAMPLE: there was only one example of the minority class in the data. Hence returning 1 majority class sample for balancing.');
    out_features = datasample(Smaj, 1);
    out_labels   = logical(majLabel*1);
    return;
end

%Number of minority examples to determine how many samples to draw from majority class
k = size(Smin,1);

out_features = datasample(Smaj,k,'Replace',false);
out_labels = (majLabel * ones([size(out_features,1) 1]));

if mkNeg
    out_labels(out_labels==0) = -1;
else
    out_labels = logical(out_labels);
end

end