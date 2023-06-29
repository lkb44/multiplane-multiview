function ploth = feature_map_slice_v2(img,featmap)
% Function for visualizing feature heatmaps overlaid onto imaging data
% See also: feature_map_slice.m for another version

% INPUTS
%       1- img: 2D slice or 3D imgume of interest
%       2- featmap: 2D heatmap of feature values extracted within mask
% OUTPUTS
%       1- ploth: handle to current plot

% example call:
%   feature_map(img,featmap)


%% check inputs
if nargin ~= 2
    error('Incorrect number of arguments');
end

if ~all(ismember(size(img),size(featmap)))
    error('Size of FEATMAP must equal size of IMG');
end
if numel(size(img))~=2 
    error('FEATURE_MAP_SLICE_V2.M handles 2D images only. See FEATURE_MAP.M for 3D visualization.');
end

%% Preapre Data
img = double(img);
featslice = double(featmap);

%% ----overlay of feature map onto img------%
cla;colorbar off;

imgslice = img/max(img(:));
rgbslice = imgslice(:,:,[1 1 1]);

bwROIlocations = ~isnan(featslice);
g = imagesc(featslice);
colormap(gca,'jet'); c = colorbar('east');
c.Color = 'w'; c.FontSize = 12;

alpha(g,1); 
hold on; h = imagesc(rgbslice);
set(h,'AlphaData',~bwROIlocations);
axis image
axis off

if nargout == 1
    ploth = h;
end

end



