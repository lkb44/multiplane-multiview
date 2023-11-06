import medviz as viz
from pathlib import Path
from scipy import io
import numpy as np
import os
import re


ROI = 'Fat'
roi = 'fat'
PLANE = 'Axial'
plane = 'ax'

dir_path = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/ResumeNPY/' + PLANE + '/'
out_dir_path = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/CollageFeatures/' + PLANE + '/'

image_dir_path = dir_path + 'Image/'
mask_dir_path = dir_path + ROI + '/'

patient_data = {}
window_sizes = ['3', '5', '7', '9', '11']

for filename in os.listdir(image_dir_path):
    if not filename.startswith('.'):  # Exclude hidden files
        patient_id = re.search(r'Patient-(\d+)', filename).group(1)
        image_path = image_dir_path + 'Patient-' + patient_id + '_' + plane + '_ls.npy'
        mask_path = mask_dir_path + 'Patient-' + patient_id + '_' + plane + '_label_' + roi + '_ls.npy'
        out_path = out_dir_path + ROI
        file_name = 'Patient-' + patient_id + '_' + plane

        # Load the two .npy arrays
        mask = np.load(mask_path)
        image = np.load(image_path)

        # Get the dimensions of each array
        mask_dim = mask.shape
        image_dim = image.shape


        # Check if the dimensions match
        if mask_dim == image_dim:
            print("The dimensions of both arrays match.")
            print("Extracting collage features from Patient-" + patient_id) 
            image, mask = viz.read_image_mask(image = image_path, mask = mask_path)
            viz.feats.collage2d(image, mask, window_sizes = [3, 5, 7, 9, 11], save_path = out_path, out_name = file_name)

        else:
            print("The dimensions of the arrays do not match. Attempting to convert 3D to 2D.")

            array_3d = np.load(mask_path)
            # Extract the 2D array from the 3D array
            array_2d = array_3d[0]

            # Save the 2D array back over the original file
            np.save(mask_path, array_2d)

            print("Conversion complete and saved to the original file. Now attempting to extract features")
            image, mask = viz.read_image_mask(image = image_path, mask = mask_path)
            viz.feats.collage2d(image, mask, window_sizes = [3, 5, 7, 9, 11], save_path = out_path, out_name = file_name)

print("Features extracted into .npy and .pkl format. Combining window sizes and converting to .mat now.")

npy_path = out_dir_path + ROI
mat_path = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/CollageFeaturesMAT/' + PLANE + '/' + ROI + '/'
if not os.path.exists(mat_path):
    os.makedirs(mat_path)

for file_name in os.listdir(npy_path):
    if not file_name.startswith('.'):  # Exclude hidden files
        if file_name.endswith('ws_3.npy'):
            # Extract patient ID and from the file name
            parts = file_name.split('_')
            patient_id = parts[2]

            print("Attempting to combine window sizes/convert: " + patient_id)

            # Specify window file paths
            ws3_file_name = f'Feats_Col_{patient_id}_{plane}_ws_3.npy'
            ws5_file_name = f'Feats_Col_{patient_id}_{plane}_ws_5.npy'
            ws7_file_name = f'Feats_Col_{patient_id}_{plane}_ws_7.npy'
            ws9_file_name = f'Feats_Col_{patient_id}_{plane}_ws_9.npy'
            ws11_file_name = f'Feats_Col_{patient_id}_{plane}_ws_11.npy'

            # Load in npy arrays for each window size
            ws3_npy_data = np.load(os.path.join(npy_path, ws3_file_name), allow_pickle=True).item()
            ws5_npy_data = np.load(os.path.join(npy_path, ws5_file_name), allow_pickle=True).item()
            ws7_npy_data = np.load(os.path.join(npy_path, ws7_file_name), allow_pickle=True).item()
            ws9_npy_data = np.load(os.path.join(npy_path, ws9_file_name), allow_pickle=True).item()
            ws11_npy_data = np.load(os.path.join(npy_path, ws11_file_name), allow_pickle=True).item()

            # Create array lists for each feature in each window size array
            ws3_array_list = [ws3_npy_data[key] for key in ws3_npy_data]
            ws5_array_list = [ws5_npy_data[key] for key in ws5_npy_data]
            ws7_array_list = [ws7_npy_data[key] for key in ws7_npy_data]
            ws9_array_list = [ws9_npy_data[key] for key in ws9_npy_data]
            ws11_array_list = [ws11_npy_data[key] for key in ws11_npy_data]

            concatenated_arrays = []

            array_lists = [ws3_array_list, ws5_array_list, ws7_array_list, ws9_array_list, ws11_array_list]

            # Concatenate array lists from each window size and reshape into a vector
            for array_list in array_lists:
                concatenated_array = np.vstack(array_list).reshape(1, -1)
                concatenated_arrays.append(concatenated_array)

            # Reshape final array into a vector
            concatenated_array = np.concatenate(concatenated_arrays).reshape(1, -1)

            # Create the output file name
            output_file_name = f'Feats_Col_{patient_id}_{plane}.mat'
            output_file_path = os.path.join(mat_path, output_file_name)

            # Save as a .mat file
            io.savemat(output_file_path, {'array_data': concatenated_array})
            print(f"Saved transformed array as {output_file_path}")