#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 14:22:59 2022

@author: Tom
"""

import os
import SimpleITK as sitk
import numpy as np
import shutil
import glob

directory = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Data/Axial/"
# directory = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Data/Coronal_resampled/" # UNCOMMENT FOR CORONAL

directories = os.listdir(directory)
directories = [d for d in directories if os.path.isdir(os.path.join(directory, d))]

mask_titles = ["mask","label", "Segmentation"]

for folder in directories:
    
    # Find volume file
    files = glob.glob(os.path.join(directory, folder) + "/*.mha")
    
    vol_file = [f for f in files if "mask" not in f and "label" not in f and "Segmentation" not in f]
     
    # # DEBUGGING ONLY - REMOVE VOL FOLDERS
    # if os.path.isdir(os.path.join(directory, folder, "vol")):
    #     shutil. rmtree(os.path.join(directory, folder, "vol"))
    
    # Create vol directory
    if not os.path.isdir(os.path.join(directory, folder, "vol")):
        os.mkdir(os.path.join(directory, folder, "vol"))
    
    # Move volume file
    if(len(vol_file) != 1):
        raise Exception ("Too many volumes found!")
    else:
        vol_file = vol_file[0]
        
    src = os.path.join(directory, folder, vol_file)
    
    output = os.path.join(directory, folder, 'vol')
    
    shutil.copy2(src, output)
    
    print("Moving " + src + " to " + output)

