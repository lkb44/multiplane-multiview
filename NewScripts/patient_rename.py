import os
import csv

parent_path = '/Users/leobao/Documents/MultiPlanePipeline/Renaming'
spreadsheet_path = '/Users/leobao/Documents/MultiPlanePipeline/multiplane-multiview/T0-2vT3-4cohort_split.csv'

# Read the spreadsheet and store the mapping of original folder names to new folder names
folder_mapping = {}
with open(spreadsheet_path, 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row
    for row in reader:
        original_folder = row[0].strip()
        new_folder = row[1].strip()
        folder_mapping[original_folder] = new_folder

# Iterate over the patient folders and rename them accordingly
for folder_name in os.listdir(parent_path):
    folder_path = os.path.join(parent_path, folder_name)
    if os.path.isdir(folder_path) and folder_name.startswith(('UH', 'VA')) and '-' in folder_name:
        if folder_name in folder_mapping:
            new_folder_name = folder_mapping[folder_name]
            
            # Rename the patient folder
            new_folder_path = os.path.join(parent_path, new_folder_name)
            os.rename(folder_path, new_folder_path)
            print(f'Renamed folder "{folder_name}" to "{new_folder_name}"')
        else:
            print(f'No mapping found for original folder name "{folder_name}" in the spreadsheet.')