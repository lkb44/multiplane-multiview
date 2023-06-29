function paddedmask = padbox(mask,n)
% Create a new padded mask using a bounding box of n-pixels outside the current mask boundary

% example call:
% pm = padmask(mask,10); 

%{
INPUTS
mask = original mask mask
n = number of pixels away from outer edges of mask (width of bounding box) in x-y plane
%}

%{
OUTPUTS
paddedmask = an image (same size as mask) but now "dilated" n pixels in x-y
%}

%{
 ____________________________________________
|(1,1)          ...           (1,size(vol,2)|
|                                           |
|                                           |
|                                           |
|                                           |
|(size(vol,1),1)...(size(vol,1),size(vol,2))|
_____________________________________________

%}

% @Jacob Antunes, 2015
% Last Updated: 04-23-2018


%% make all zero values in the original volume nonzero

if ~isempty(mask) %use the mask to create a bounding box
    [row, col, slice] = ind2sub(size(mask),find(mask ~= 0)); 
else
    error('mask is empty.');
end

%% find dimensions of padded box
bottom = min(max(row)+n,size(mask,1));
top = max(min(row)-n,1);
left = max(min(col)-n,1);
right = min(max(col)+n,size(mask,2));
front = min(slice);
back = max(slice);


%% apply padding

paddedmask = mask;
paddedmask(top:bottom,left:right,front:back) = 10;
paddedmask(mask~=0) = mask(mask~=0);

%testing purposes
% subplot(2,1,1);
% imagesc(mask(:,:,10));
% subplot(2,1,2);
% imagesc(paddedmask(:,:,10));

%view results
% v(:,:,:,1) = vol;
% v(:,:,:,2) = croppedV;
% mosaic4D(v);

end


