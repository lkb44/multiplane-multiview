import os
import SimpleITK as sitk
import numpy as np

def convert_mha_to_npy(input_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    file_list = os.listdir(input_path)

    for file_name in file_list:
        if file_name.endswith('.mha'):
            mha_path = os.path.join(input_path, file_name)
            output_name = os.path.splitext(file_name)[0] + '.npy'
            npy_path = os.path.join(output_path, output_name)

            image = sitk.ReadImage(mha_path)
            image_array = sitk.GetArrayFromImage(image)

            np.save(npy_path, image_array)

if __name__ == "__main__":
    input_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Axial_ProxFat10"
    output_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Axial_ProxFat10_NPY"
    convert_mha_to_npy(input_path, output_path)
