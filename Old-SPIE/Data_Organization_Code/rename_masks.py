#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:47:20 2022

@author: Tom
"""

import os
import sys

data_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/Axial"

dir_list = os.listdir(data_path)
dir_list = sorted([ d for d in dir_list if "." not in d])
    
for folder in dir_list:
    files = os.listdir(os.path.join(data_path, folder))
    for old_name in files:
        if "label" in old_name or "Label" in old_name:
            new_name = old_name.replace(old_name, folder + "_pre_ax_label_resampled.mha")
        if "label" not in old_name and "DS" not in old_name:
            new_name = old_name.replace(old_name, folder + "_pre_ax_resampled.mha")
        #new_name = new_name.replace("-", "_")
        os.rename(os.path.join(data_path, folder, old_name), os.path.join(data_path, folder, new_name))
        
        
dir_list = os.listdir(data_path)
dir_list = sorted([ d for d in dir_list if "." not in d])    
    
for folder in dir_list:
    files = os.listdir(os.path.join(data_path, folder))
    for old_name in files:
        new_name = old_name.replace("raw_", "")
        os.rename(os.path.join(data_path, folder, old_name), os.path.join(data_path, folder, new_name))
        
    
        
        