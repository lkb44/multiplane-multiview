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