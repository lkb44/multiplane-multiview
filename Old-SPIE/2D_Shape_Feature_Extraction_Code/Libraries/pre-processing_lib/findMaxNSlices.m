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
idxs = max_idx;
if max_idx == 1 %if the largest ROI is the end of the ROI... (this should be updated later to check for this when assigning j
    idxs = 1:n;
elseif max_idx == length(nz_entries)
    idxs = length(nz_entries)-n+1:length(nz_entries);
else
    for j = 1:floor(n/2)
        nz_1 = nz_entries(max_idx-j); nz_2 = nz_entries(max_idx+j); 
        [~, next_max_idx] = max([nz_1,nz_2]); %determine which j adjacent slice to add first
        if next_max_idx == 1
            idxs = cat(1,idxs,max_idx-j);
        else
            idxs = cat(1,idxs,max_idx+j);
        end  
        if length(idxs)<n && next_max_idx == 1, idxs = cat(1,idxs,max_idx+j); end %in case of n being odd, add the other j adjacent slice
        if length(idxs)<n && next_max_idx == 2, idxs = cat(1,idxs,max_idx-j); end
    end
end
idxs = sort(idxs,'ascend');
slice_idxs = non_empties(idxs);

        