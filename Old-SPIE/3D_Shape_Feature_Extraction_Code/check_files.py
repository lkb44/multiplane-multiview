#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 18:21:43 2022

@author: Tom
"""

import os
import glob

root = "/scratch/pbsjobs/job.16967448.hpc/Axial/"
dirs = sorted(os.listdir(root))

issues = []
for thedir in dirs:
    thepath = os.path.join(root, thedir)
    thepath += "/*.txt"
    files = glob.glob(thepath)
    if not files:
        issues.append(thepath)

print("The following directories have no text files: " + str(issues))

with open('issues.txt', 'a') as f:
    f.write("The following directories have no text files: \n")
    for i in issues:
        f.write(i + "\n")
