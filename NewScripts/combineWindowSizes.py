import os
import numpy as np
import re

input_directory = "/Volumes/Crucial X6/CollageFeatures/Axial_Fat/"
output_directory = "/Users/leobao/Documents/MultiPlanePipeline/CollageReshaped/Axial_Fat/"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

patient_data = {}

for filename in os.listdir(input_directory):
    if not filename.startswith('.'):  # Exclude hidden files
        if filename.endswith(".mat"):
            patient_id = re.search(r'Patient-(\d+)', filename).group(1)
            window_size = int(re.search(r'ws_(\d+)', filename).group(1))
            
            if patient_id not in patient_data:
                patient_data[patient_id] = {}
            
            patient_data[patient_id][window_size] = np.loadtxt(os.path.join(input_directory, filename))

final_patient_data = {}

for patient_id, windows in patient_data.items():
    combined_array = np.concatenate([windows[3], windows[5], windows[7], windows[9], windows[11]], axis=1)
    final_patient_data[patient_id] = combined_array

for patient_id, combined_array in final_patient_data.items():
    output_filename = os.path.join(output_directory, f"Patient_{patient_id}_combined.npy")
    np.save(output_filename, combined_array)

