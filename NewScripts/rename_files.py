import os
import sys

data_path = "/Users/leobao/Documents/MultiPlanePipeline/multiplane-multiview/Coronal_Raw"

dir_list = os.listdir(data_path)
dir_list = sorted([ d for d in dir_list if "." not in d])
    
for folder in dir_list:
    files = os.listdir(os.path.join(data_path, folder))
    for old_name in files:
        if "label" not in old_name and "DS" not in old_name:
            new_name = old_name.replace(old_name, folder + "_raw_cor.mha")
        if "segmentation" in old_name or "Segmentation" in old_name or "mask" in old_name or "label" in old_name:
            new_name = old_name.replace(old_name, folder + "_raw_cor_label.mha")
        os.rename(os.path.join(data_path, folder, old_name), os.path.join(data_path, folder, new_name))