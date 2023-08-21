import os
import numpy as np
from scipy import io

# Define the plane
plane = 'cor'

# Define the input directory and window sizes
#input_directory = '/Volumes/Crucial X6/MP_COLLAGE/Axial_Tumor_NPY_PKL'
input_directory = '/Volumes/Crucial X6/Resume2'
window_sizes = ['3', '5', '7', '9', '11']

# Define the output directory
output_directory = '/Volumes/Crucial X6/CollageFeatures/Coronal_ProxFat15/'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Loop through the input directory
for file_name in os.listdir(input_directory):
    if not file_name.startswith('.'):  # Exclude hidden files
        if file_name.endswith('ws_3.npy'):
            # Extract patient ID and from the file name
            parts = file_name.split('_')
            patient_id = parts[2]

            # Specify window file paths
            ws3_file_name = f'Feats_Col_{patient_id}_{plane}_ls_ws_3.npy'
            ws5_file_name = f'Feats_Col_{patient_id}_{plane}_ls_ws_5.npy'
            ws7_file_name = f'Feats_Col_{patient_id}_{plane}_ls_ws_7.npy'
            ws9_file_name = f'Feats_Col_{patient_id}_{plane}_ls_ws_9.npy'
            ws11_file_name = f'Feats_Col_{patient_id}_{plane}_ls_ws_11.npy'

            # Load in npy arrays for each window size
            ws3_npy_data = np.load(os.path.join(input_directory, ws3_file_name), allow_pickle=True).item()
            ws5_npy_data = np.load(os.path.join(input_directory, ws5_file_name), allow_pickle=True).item()
            ws7_npy_data = np.load(os.path.join(input_directory, ws7_file_name), allow_pickle=True).item()
            ws9_npy_data = np.load(os.path.join(input_directory, ws9_file_name), allow_pickle=True).item()
            ws11_npy_data = np.load(os.path.join(input_directory, ws11_file_name), allow_pickle=True).item()

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
            output_file_path = os.path.join(output_directory, output_file_name)

            # Save as a .mat file
            io.savemat(output_file_path, {'array_data': concatenated_array})
            print(f"Saved transformed array as {output_file_path}")