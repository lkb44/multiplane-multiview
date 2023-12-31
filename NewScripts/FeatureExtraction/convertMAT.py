import os
import numpy as np
from scipy import io

ROI = 'ProxFat10'
roi = 'proxfat10'
PLANE = 'Axial'
plane = 'ax'

out_dir_path = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/CollageFeatures/' + PLANE + '/'

patient_data = {}
window_sizes = ['3', '5', '7', '9', '11']

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