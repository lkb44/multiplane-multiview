import os
import shutil

source_directory = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialOut"
dir_list = os.listdir(source_directory)


dir_list = sorted([d for d in dir_list if "." not in d])
print(dir_list)

destination_directory = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialMasks"

for directory in dir_list:
    src = os.path.join(source_directory, directory)
    dest = os.path.join(destination_directory, directory)
    shutil.copytree(src, dest)
    print("Moved " + src + " to " + dest + "/n")
