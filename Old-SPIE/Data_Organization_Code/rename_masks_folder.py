#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:47:20 2022

@author: Tom
"""

import os
import sys

data_path = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialOut"

dir_list = os.listdir(data_path)
dir_list = sorted([ d for d in dir_list if "." not in d])
    
for folder in dir_list:
    files = os.listdir(os.path.join(data_path, folder))
    for old_name in files:
        if "masks" in old_name:
            new_name = old_name.replace(old_name, folder + "_masks")
            os.rename(os.path.join(data_path, folder, old_name), os.path.join(data_path, folder, new_name))
        
        
dir_list = os.listdir(data_path)
dir_list = sorted([ d for d in dir_list if "." not in d])    