import os

data_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Axial_Image/"

dir_list = os.listdir(data_path)
dir_list = sorted([d for d in dir_list if "." not in d])

for folder in dir_list:
    files = os.listdir(os.path.join(data_path, folder))
    for old_name in files:
        if "Patient" in old_name:
            new_name = old_name.replace("_", "-")
            os.rename(os.path.join(data_path, folder, old_name), os.path.join(data_path, folder, new_name))
#        if "label" not in old_name and "DS" not in old_name and "Label" not in old_name:
#            new_name = folder + "_ax.mha"
#            os.rename(os.path.join(data_path, folder, old_name), os.path.join(data_path, folder, new_name))
#        elif "Label" in old_name or "label" in old_name:
#            new_name = folder + "_ax_label.mha"
#            os.rename(os.path.join(data_path, folder, old_name), os.path.join(data_path, folder, new_name))