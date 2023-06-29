function slice_idxs = findMaxNSlices(volmask,n)
% Returns n slice indices corresponding to slices surrounding (and including) the slice with the largest ROI.

% INPUTS
% volmask: 3D label volume (label > 0)
% n: number of desired slices to identify

% OUTPUTS
% slice_idxs: the indices of the slices identified

%% Check inputs
if nargin ~=2
    error('Incorrect number of input arguments.');
end
volmask = double(volmask);
if isempty(volmask) || isempty(n) || nnz(volmask)==0 || numel(size(volmask))~=3 || ~isinteger(uint8(n)) || n <= 0
    error('Inputs must be a (1) non-empty 3D matrix and (2) non-empty integer')
end

if n == 1
    slice_idxs = findMaxSlice(volmask);
    return;
end

%% Identify slice with maximum ROI
non_empties = []; %non empty slices
nz_entries = []; %# non-zero pixels on non empty slices
for i = 1:size(volmask,3)
    if nnz(volmask(:,:,i))>0
        non_empties = cat(1,non_empties,i);
        nz_entries = cat(1,nz_entries,nnz(volmask(:,:,i)));
    end
end
        
if length(non_empties)<n 
    error('VOLMASK does not contain the desired number of non-zero slices.');
end

%% Get slice w/ maxROI & surroundng slices
[~, max_idx] = max(nz_entries);
if max_idx == 1 %if the largest ROI is the beginning of the ROI... 
    slice_idxs = non_empties(max_idx:max_idx+n-1);
elseif max_idx == length(non_empties) %if the largest ROI is the end of the ROI... 
    slice_idxs = non_empties(max_idx-n+1:max_idx);
else
    slice_idxs = non_empties(max_idx)-1:non_empties(max_idx)+1;
end

        