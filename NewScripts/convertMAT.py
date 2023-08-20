import os
import numpy as np
from scipy import io

# Define the input directory and window sizes
input_directory = '/Volumes/Crucial X6/Resume'
window_sizes = ['3', '5', '7', '9', '11']

# Define the output directory
output_directory = '/Volumes/Crucial X6/CollageFeatures/Axial_Fat/'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Loop through the input directory
for file_name in os.listdir(input_directory):
    if not file_name.startswith('.'):  # Exclude hidden files
        if file_name.endswith('.npy'):
            # Extract patient ID and window size from the file name
            parts = file_name.split('_')
            patient_id = parts[2]
            window_size = parts[-1].replace('.npy', '')

            if window_size in window_sizes:
                # Load the .npy array
                npy_data = np.load(os.path.join(input_directory, file_name), allow_pickle=True).item()

                # Extract and concatenate the arrays
                array_list = [npy_data[key] for key in npy_data]
                concatenated_array = np.vstack(array_list)

                # Transform to have a shape of (1, 52)
                transformed_array = concatenated_array.reshape(1, -1)

                # Create the output file name
                output_file_name = f'Feats_Col_Patient-{patient_id}_ax_ls_{window_size}.mat'
                output_file_path = os.path.join(output_directory, output_file_name)

                # Save as a .mat file
                io.savemat(output_file_path, {'array_data': transformed_array})
                print(f"Saved transformed array as {output_file_path}")