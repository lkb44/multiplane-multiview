%% resampleScalarVolume
function resampleScalarVolume(newXres, newYres, newZres, img_file_in, img_path_out, mask_file_in, mask_path_out)

% Goal: Batch-automated resampling of 3-D image data to a desired, fixed
%       resolution across a cohort. Currently, parameters have been set to
%       1by1byz resolution, where z is the original resolution in
%       z-direction. Call's 3D Slicer's "ResampleScalarVolume" module

% INPUTS
%   newXres: desired new resolution in X
%   newYres: desired new resolution in X
%   newZres: desired new resolution in X
%   img_file_in: full file path (with extension) to a 3D image file
%   img_path_out: directory to where resampled 3D image file should be saved to
%   mask_file_in: (optional) full file path (with extension) to a 3D mask file
%   mask_path_out: (optional) directory to where resampled 3D mask file should be saved to

% **** WARNING #1****
% If mask is to be resampled, please pass in the corresponding image!
% Otherwise nearest neighbor resampling may yield slightly different size
% *****************

% **** WARNING #2****
% Directories probably can't be on shared / external drives. i.e. need to be on local or C:/ drive!
% *****************

% @author Jacob Antunes
% @date 09-17-2019

%% FOR MORE INFO
% ResampleScalarVolume --help

%***IN ORDER TO WORK***%
%{
1- identify where your Slicer program folder is:
    'C:\Program Files\Slicer 4.8.1';
2- on WINDOWS, need to add the following subdirectories to PATH environment variable:
    'C:\Program Files\Slicer 4.8.1'
    'C:\Program Files\Slicer 4.8.1\lib\Slicer-4.8'
    'C:\Program Files\Slicer 4.8.1\lib\Slicer-4.8\cli-modules'
3-restart MATLAB after applying these changes
%}

%% Resample Image

out_dir = img_path_out;
myinterp = 'linear'; %linear, for T2w MRI sequences, anyways
[~,name,ext] = fileparts(img_file_in); 
module = 'ResampleScalarVolume ';
interpolator = ['-i ' myinterp ' '];
spacing = ['-s ' num2str(newXres) ',' num2str(newYres) ',' num2str(newZres) ' '];  %new X, Y, Z resolutiosn
img_out_file = fullfile(out_dir,[name '_resampled' ext]); %img_file_out is same name as in img_file_in, but with "_resampled" tag
command = [module,interpolator,spacing,img_file_in,' ',img_out_file]; 
system(command);

%% Resample Mask
if nargin > 5 && ~isempty(mask_file_in)
    out_dir = mask_path_out;
    myinterp = 'nearestNeighbor'; %nearestNeighbor, for masks, anyways
    [~,name,ext] = fileparts(mask_file_in); 
    module = 'ResampleScalarVolume ';
    interpolator = ['-i ' myinterp ' '];
    spacing = ['-s ' num2str(newXres) ',' num2str(newYres) ',' num2str(newZres) ' '];  %new X, Y, Z resolutiosn
    mask_out_file = fullfile(out_dir,[name '_resampled' ext]); %img_file_out is same name as in img_file_in, but with "_resampled" tag
    command = [module,interpolator,spacing,mask_file_in,' ',mask_out_file]; 
    system(command);
    
    %***CHECK IMG, MASK DIMENSIONS***
    adjustMaskDims(img_out_file,mask_out_file)
end

end

%% Subfunction
function adjustMaskDims(img_file,mask_file)
funcname = 'resampleScalarVolume.m';
funcpath = which(funcname);
codepath = funcpath(1:end-length(funcname));

addpath([codepath '../images']);
addpath([codepath '../mha']);

vol = vv(img_file);
[mask,~,maskinfo] = vv(mask_file);

if all(size(vol)==size(mask)), return; end %dont need to edit anything

% error4: mask is 1 slice too short & 1 row/1column short
if (size(mask,1) == size(vol,1)-1 || size(mask,2) == size(vol,2)-1) && size(mask,3) == size(vol,3)-1
    maskcopy = zeros(size(vol));
    maskcopy(2:end,2:end,1:end-1) = mask(:,:,:);

%error 1: mask is 1 row / 1 column short
elseif size(mask,1) == size(vol,1)-1 || size(mask,2) == size(vol,2)-1
    maskcopy = zeros(size(vol));
    maskcopy(2:end,2:end,:) = mask(:,:,:);

% error 2: mask is 1 slice too many
elseif size(mask,3) == size(vol,3)+1
    maskcopy = mask(:,:,1:end-1);
    
% error3: mask is 1 slice too short
elseif size(mask,3) == size(vol,3)-1
    maskcopy = zeros(size(vol));
    maskcopy(:,:,1:end-1) = mask;
else
    error('Unhandled error. Are you sure mask and img files came from same patient?');
end

% testing: visualize to be sure we are happy with our correction
figure;vv(vol,maskcopy)
mask = maskcopy;

% save it back out
mha_write_volume(mask_file,mask,maskinfo.PixelDimensions,maskinfo.Offset,maskinfo.TransformMatrix);

rmpath([codepath '../images']);
rmpath([codepath '../mha']);
end