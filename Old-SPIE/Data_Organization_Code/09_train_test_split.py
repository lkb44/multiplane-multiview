#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on thu May 12 09:21:57 2022

@author: tom
"""

import pandas as pd
import shutil
import glob
import os

# Specify input
dataset = "train" # train, test, or va_test
view = "Coronal" # Axial or Coronal
dataset_file = "../Data/train_test_labels/" + dataset + "_labels.csv"
src = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_Features/" + view + "/All/"

# Read csv file with data labels
data = pd.read_csv(dataset_file)

# Get a list of patients beloning to dataset
pt_list = data["Patient"].tolist()

# Copy patient directories with feature text files
for pt in pt_list:
    dest = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_Features/" + view + "/" + dataset + "/" + pt + "/"
    if os.path.exists(dest):
        shutil.rmtree(dest)
    os.mkdir(dest)
    for file in glob.glob(os.path.join(src, pt, "*.txt")):
        shutil.copy2(file,dest)
    print("Finished copying text files for " + pt)






