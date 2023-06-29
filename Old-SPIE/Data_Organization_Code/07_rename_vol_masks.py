#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 17:15:55 2022

@author: Tom
"""

import os

directory = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Data/Axial/"

directories = os.listdir(directory)
directories = [d for d in directories if os.path.isdir(os.path.join(directory, d))]

for folder in directories:
    vol_files = os.listdir(os.path.join(directory, folder, "vol"))
    mask_files = os.listdir(os.path.join(directory, folder, "masks"))
    
    os.chdir(os.path.join(directory, folder, "vol"))
    os.rename(vol_files[0], "vol.mha")
    print("Renamed " + vol_files[0])
    
    # mask_org = [i for i in mask_files if "rw" in i]
    # os.chdir(os.path.join(directory, folder, "masks"))
    # os.rename(mask_org[0], "mask_rw.mha")
    # print("Renamed " + mask_org[0])
    
    # mask_org = [i for i in mask_files if "lumen" in i]
    # os.chdir(os.path.join(directory, folder, "masks"))
    # os.rename(mask_org[0], "mask_lumen.mha")
    # print("Renamed " + mask_org[0])
    
    # mask_org = [i for i in mask_files if "fat" in i]
    # os.chdir(os.path.join(directory, folder, "masks"))
    # os.rename(mask_org[0], "mask_fat.mha")
    # print("Renamed " + mask_org[0])