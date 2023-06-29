#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 18:47:03 2022

@author: Tom
"""

import os

directory = "/Volumes/GoogleDrive/Shared drives/INVent/Dhruv/Multi_View_Shape_Experiments/Data/Coronal/"

directories = os.listdir(directory)
directories = [d for d in directories if os.path.isdir(os.path.join(directory, d))]

for folder in directories:
    if not os.path.isdir(os.path.join(directory, folder, "features")):
        os.mkdir(os.path.join(directory, folder, "features"))
    