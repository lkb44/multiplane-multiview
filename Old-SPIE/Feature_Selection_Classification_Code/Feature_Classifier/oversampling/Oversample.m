function [out_data, out_data_labels] = Oversample(method, in_data, in_data_labels, varargin)

% Makes function calls to oversample/undersample the minority class to establish class balance.
% Returns the new data (original data concatenated with newly oversampled minority class data).
% Class balance may not be perfectly 1:1 due to rounding in certain algorithms.

% make results reproducible by resetting the random number generator:
% rng('default');
%***warning*** this may make subsets the same across
% nIterations!!! If you want reproducibility, rng('default') should be
% called before running classification altogether!

if ismember(upper(method),{'SMOTE','ADASYN','UNDERSAMPLE'})
    if strcmp(upper(method),'UNDERSAMPLE')
        [out_majordata,out_majorlabel]=feval(method, in_data, in_data_labels, varargin{:});
        if isempty(out_majorlabel) %there was already class balance!
            out_data = in_data; out_data_labels = in_data_labels;
            return;
        end
        out_minordata = in_data(in_data_labels~=out_majorlabel(1),:);
        out_minorlabel = in_data_labels(in_data_labels~=out_majorlabel(1));
        out_data = [out_majordata; out_minordata];
        out_data_labels = [out_majorlabel; out_minorlabel];
    else
        [out_minordata,out_minorlabel]=feval(method, in_data, in_data_labels, varargin{:});
        if isempty(out_minorlabel) %there was already class balance!
            out_data = in_data; out_data_labels = in_data_labels;
            return;
        end
        out_data = [in_data; out_minordata];
        out_data_labels = [in_data_labels; out_minorlabel];
    end
else
    error('OVERSAMPLE: invalid method option.');
end
