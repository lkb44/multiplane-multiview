import os
import SimpleITK as sitk
import nibabel as nib

def convert_mha_to_nifti(input_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    file_list = os.listdir(input_path)

    for file_name in file_list:
        if file_name.endswith('.mha'):
            mha_path = os.path.join(input_path, file_name)
            output_name = os.path.splitext(file_name)[0] + '.nii.gz'
            nifti_path = os.path.join(output_path, output_name)

            image = sitk.ReadImage(mha_path)
            sitk.WriteImage(image, nifti_path)

if __name__ == "__main__":
    input_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Coronal_Image"
    output_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Coronal_Image_NIFTY"
    convert_mha_to_nifti(input_path, output_path)