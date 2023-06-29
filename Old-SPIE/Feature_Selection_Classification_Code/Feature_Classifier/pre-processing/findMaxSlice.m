function maxslice_idx = findMaxSlice(vol)
% Returns the index of the maximum-data-containing slice

% INPUTS
% vol = 3D matrix, with background = 0; typically this would be a binary mask

if numel(size(vol))~=3
%     error('findMaxSlice.m is only meant for 3D inputs.');
    maxslice_idx = 1;
end

vol = double(vol);
vol = double(vol>0);

maxslice_idx = 1;
maxvals = sum(sum(vol(:,:,1)));
for i = 2:size(vol,3)
    vals = sum(sum(vol(:,:,i)));
    if vals>maxvals
        maxvals = vals;
        maxslice_idx = i;
    end
end

end




