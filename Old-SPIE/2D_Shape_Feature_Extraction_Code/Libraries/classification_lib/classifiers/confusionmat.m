function [cmat,colclasses] = confusionmat(truelabels,estimatedlabels)

%setup to handle multiclass as well
% tp = []; tn = []; fp = []; fn = []; %doesnt really hold as well for multiclass
colclasses = unique(truelabels);
cmat = zeros(length(colclasses));
if isa(truelabels,'cell')
    for i = 1:length(truelabels)
        u_true = find(contains(colclasses,truelabels{i}));
        u_est = find(contains(colclasses,estimatedlabels{i}));
        cmat(u_est,u_true) = cmat(u_est,u_true)+1;
    end
else
    error('not implemented for anything other than cell array of strings');
end

