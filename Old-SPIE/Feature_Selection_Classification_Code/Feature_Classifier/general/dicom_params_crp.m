function [params_c,params_t]=dicom_params_crp(params)

params_t{1,1}={'ScanningSequence'};
params_t{1,2}={'TR'};
params_t{1,3}={'TE'};
params_t{1,4}={'ImagingFreq'};
params_t{1,5}={'Field'};
params_t{1,6}={'PixelSpacing'};
params_t{1,7}={'Resolution_[Y,X]'};
params_t{1,8}={'SliceSpacing'};

params_c{1,1}=params.ScanningSequence;
params_c{1,2}=params.RepetitionTime;
params_c{1,3}=params.EchoTime;
params_c{1,4}=params.ImagingFrequency;
params_c{1,5}=params.MagneticFieldStrength;
params_c{1,6}=[params.PixelSpacing(1),params.PixelSpacing(2)];
params_c{1,7}=double([params.Rows,params.Columns]);
params_c{1,8}=params.SliceThickness;


