import os
import sys
path = os.getcwd()
import numpy as np
from inspect import getmembers, isfunction

import medviz as viz

np.load(r"/Volumes/Crucial X6/MacMiniData882023/Data/Axial_Image_NPY/Patient-002_ax_ls.npy")
#load a test image to make sure orientation is working the way you need it

images_path="/Volumes/Crucial X6/MacMiniData882023/Data/Coronal_Image_NPY"
masks_path="/Volumes/Crucial X6/MacMiniData882023/Data/Coronal_Fat_NPY"
save_path="dataset_schema.csv"

#viz.preprocess.match_image_masks(images_path, masks_path,image_extension=".nii",mask_extension=".nii",image_id_func= lambda x: x.split("_")[0],mask_id_func= lambda x: x.split("_")[0],save_path="dataset_schema.csv")  
viz.preprocess.match_image_masks(images_path, 
                                masks_path,
                                image_extension=".npy",
                                mask_extension=".npy",
                                image_id_func= lambda x: x.split("-")[-1][:3],
                                mask_id_func= lambda x: x.split("-")[-1][:3],
                                save_path="dataset_schema.csv")

viz.feats.compute_stats_collage(
    csv_path="dataset_schema.csv",
    stats_save_path="test_stats/",
)
print("Finish Data Scheme")