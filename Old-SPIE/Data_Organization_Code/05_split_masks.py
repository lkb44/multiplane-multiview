#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 14:22:59 2022

@author: Tom
"""

import os
import SimpleITK as sitk
import numpy as np


def write_mha(image, spacing, origin, direction, output_filename):
    
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    image.SetDirection(direction)
    
    sitk.WriteImage(image, output_filename)
    
def check_mask(input_mask):
    
    result = np.unique(input_mask)
    result = sorted(result.tolist())
    if(result != [0, 1]):
        print("Mask is not binary!")
        return False
    else:
        print("Mask is binary.")
        return True
        
    # for i in range(input_mask.shape[0]):
    #     check_slice = input_mask[i, :, :]              
                
    #     for y in range(check_slice.shape[1]):
    #         for x in range(check_slice.shape[0]):
    #             if(check_slice[x,y] > 1):
    #                 raise ValueError("Mask is not binary!")
    #             else:
    #                 binary = True
    # if(binary is True):
        
                    
                
         

#directory = "../Data/Coronal_resampled/"
directory = "/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialOut"

directories = os.listdir(directory)
directories = [d for d in directories if os.path.isdir(os.path.join(directory, d))]

# directories = [d for d in directories if "VA" in d] # UNCOMMENT FOR VA

# List of patients without binary masks
issues = []

for folder in directories:
    
    # Find mask file
    files = os.listdir(os.path.join(directory, folder))
    mask_file = [f for f in files if "mask" in f or "label" in f or "Segmentation" in f]
    volume_name = folder.replace("-","_") + "_pre_ax_label_resampled_aniresampled.mha" # Replace with _pre_ax_resampled.mha for UH axial, _pre_cor_resampled for UH cor, and or _cor_resampled.mha for VA
    vol_file = os.path.join(directory, folder, volume_name)
    
    # Read in mask file
    src_mask = sitk.ReadImage(os.path.join(directory, folder, mask_file[0]))
    try:
        vol = sitk.ReadImage(vol_file)
        print("Splitting " + mask_file[0])
        
        # Get Spacing and Origin
        spacing = vol.GetSpacing()
        origin = vol.GetOrigin()
        direction = vol.GetDirection()
        #direction = [round(i, 6) for i in direction]
        
        # Convert mask to numpy array
        src_mask = sitk.GetArrayFromImage(src_mask)
        vol = sitk.GetArrayFromImage(vol)
        
        # create empty masks for tumor
        src_mask_shape = vol.shape
        num_slices = src_mask_shape[0]
        width = src_mask_shape[1]
        height = src_mask_shape[2]
        tumor_mask = np.zeros((len(src_mask), width, height))
        rectal_wall_mask = np.zeros((len(src_mask), width, height))
        fat_mask = np.zeros((len(src_mask), width, height))
        slice_area_tumor = []

        
        # Extract tumor from each slice
        for idx in range(len(src_mask)):
            
            img = src_mask[idx, :, :]
            
            tumor_img = np.copy(img)
            rw_img = np.copy(img)
            lumen_img = np.copy(img)
            fat_img = np.copy(img)

            tumor_img[tumor_img != 1] = 0
            tumor_img[tumor_img == 1] = 1
            rw_img[rw_img == 7] = 1
            rw_img[rw_img == 8] = 1
            rw_img[rw_img > 1] = 0
            fat_img[fat_img != 4] = 0
            fat_img[fat_img == 4] = 1

            tumor_mask[idx, :, :] = tumor_img
            rectal_wall_mask[idx, :, :] = rw_img
            fat_mask[idx, :, :] = fat_img
            
           

        
        # Check masks
        binary_tumor = check_mask(tumor_mask)
        if(binary_tumor is False):
            issues.append(mask_file[0] + "_tumor")
            raise ValueError("Tumor mask for " + mask_file[0] + " is NOT binary!" )
        binary_wall = check_mask(rectal_wall_mask)
        if(binary_wall is False):
            issues.append(mask_file[0] + "_rw")
            raise ValueError("Rectal wall mask for " + mask_file[0] + " is NOT binary!" )
        binary_fat = check_mask(fat_mask)
        if(binary_fat is False):
            issues.append(mask_file[0] + "_fat")
            raise ValueError("Fat mask for " + mask_file[0] + " is NOT binary!" )
        
        
        # convert masks into sitk images
        tumor_mask = sitk.GetImageFromArray(tumor_mask)
        rectal_wall_mask = sitk.GetImageFromArray(rectal_wall_mask)
        fat_mask = sitk.GetImageFromArray(fat_mask)
        
        # Create masks directory
        if not os.path.isdir(os.path.join(directory, folder, "masks")):
            os.mkdir(os.path.join(directory, folder, "masks"))
        
        # Save tumor mask
        out_name = mask_file[0].replace(".mha","_tumor.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(tumor_mask, spacing, origin ,direction, output)

        # Save rectal wall mask
        out_name = mask_file[0].replace(".mha","_rw.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(rectal_wall_mask, spacing, origin ,direction, output)
        
        # Save fat mask
        out_name = mask_file[0].replace(".mha","_fat.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(fat_mask, spacing, origin, direction, output)
    
    except Exception as e:
        print(e)
    
    
    


# FOR DEBUGGING ONLY: DELETE  "MASKS" DIRECTORIES IN ALL PATIENT FOLDERS
# import shutil
# import os
# directory = "../Data/Axial/"

# directories = os.listdir(directory)
# directories = [d for d in directories if os.path.isdir(os.path.join(directory, d))]
# for folder in directories:
#     print("Deleting masks folder in " + folder)
#     dir_to_rm = os.path.join(directory, folder, "masks")
#     shutil.rmtree(dir_to_rm)
    
# FOR DEBUGGING ONLY: UNCOMMENT TO CHECK THE DIRECTION OF A SPECIFIC MASK
# import SimpleITK as sitk
# volume = "/Volumes/GoogleDrive/My Drive/RSNA_2022/Morphology/Data/Axial/UH-RectalCA-002/vol/vol.mha"
# mask = "/Volumes/GoogleDrive/My Drive/RSNA_2022/Morphology/Data/Axial/UH-RectalCA-070/masks/UH_RectalCA_002_pre_ax_label_resampled_lumen.mha"

# vol =sitk.ReadImage(volume)
# mask = sitk.ReadImage(mask)

# mask_direction = mask.GetDirection()
# image_direction = mask.GetDirection()

# mask_floats = []
# image_floats = []

# import decimal 

# for i in mask_direction:
#     mask_floats.append(decimal.Decimal.from_float(i))
    
# for i in image_direction:
#     image_floats.append(decimal.Decimal.from_float(i))

# for idx in range(len(mask_floats)):
#     if(mask_floats[idx] == image_floats[idx]):
#         print("Directions are equal")
#     else:
#         print("Directions are NOT equal")