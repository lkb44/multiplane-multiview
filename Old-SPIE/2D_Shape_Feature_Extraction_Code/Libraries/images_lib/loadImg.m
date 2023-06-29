function img = loadImg(filename)
% Generic function for loading filename into MATLAB workspace

if strcmp(filename(end-3:end),'.mha')
    img = cast(mha_read_volume(mha_read_header(filename)),'double');
elseif strcmp(filename(end-6:end),'.nii.gz') || strcmp(filename(end-3:end),'.nii')
    try
        nii_vol = load_nii(filename);
    catch
        nii_vol = load_nii_untouched(filename);
    end
    img = double(nii_vol.img);
else
    error('Currently only implemented for mha and nifti extensions.');
end

