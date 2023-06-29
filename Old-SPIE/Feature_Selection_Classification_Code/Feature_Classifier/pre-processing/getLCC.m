function lcc_mask = getLCC(mask)
%function LCC_MASK = getLCC(MASK);
%
% Goal: Get the Largest Connected Component (LCC) from a mask
%
% INPUTS: 
%   + MASK: a 2D or 3D mask of 1s and 0s (1 = ROI)
% 
% OUTPUTS:
%   + LCC_MASK: a mask (same size as MASK) containing only the largest
%               connected component;
%

BW = (mask ==1);
CC = bwconncomp(BW);
l = [];
for i = 1:length(CC.PixelIdxList)
    l = [l length(CC.PixelIdxList{i})];
end
[maxPix, maxIdx] = max(l);
lcc_mask = zeros(size(mask));
lcc_mask(CC.PixelIdxList{maxIdx}) = 1;



