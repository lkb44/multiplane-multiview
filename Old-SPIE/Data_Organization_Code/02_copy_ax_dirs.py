import os
import shutil
import pandas as pd


data = pd.read_excel("../Dataset_Information.xlsx", sheet_name = "Data_to_copy")
patients = data["All Patient ID"].to_list()


#source_directory = "/Volumes/GoogleDrive/Shared drives/INVent/sxs2658/Rectal_Data/UH_Resamp/Pre/"
source_directory = "/Volumes/GoogleDrive/Shared drives/INVent/sxs2658/Rectal_Data/UH_Resamp/New_Patients2020/Pre/"
dir_list = os.listdir(source_directory)


dir_list = sorted([d for d in dir_list if "." not in d])
dir_list_formatted = dir_list
dir_list_formatted = [ele.replace("RectalCA_", "UH-RectalCA-") for ele in dir_list_formatted]
print(dir_list_formatted)

destination_directory = "/Volumes/GoogleDrive/My Drive/RSNA_2022/Morphology/Data/Axial/"

i = 0
for index, directory in enumerate(dir_list_formatted):
    if directory in patients:
        src = os.path.join(source_directory, dir_list[index])
        dest = os.path.join(destination_directory, directory)
        shutil.copytree(src, dest)
        print("Moved " + src + " to " + dest + "/n")
        i += 1
