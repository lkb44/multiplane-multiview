"""
Created on Thu Jun 29 8:58:58 2023

@author: Leo
"""

import SimpleITK as sitk
import pandas as pd
import numpy as np
import sys
import os
import cv2

def write_mha(image, spacing, origin, output_filename):
    
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    #image.SetDirection(direction)
    
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

directory = "/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/lkb44/SPIE/SPIE Data/AxialSingle/"
dir_list = os.listdir(directory)
dir_list = sorted([ d for d in dir_list if "." not in d])


slice_tracker = pd.DataFrame(columns = ['Patient', 'Slice_Number', 'Tumor_Size', 'Dilated_Mask_Size', 'Fat_Mask_Size'])

# List of patients without binary masks
issues = []

for folder in dir_list:
    
    # Find mask file
    files = os.listdir(os.path.join(directory, folder))
    mask_file = [f for f in files if "label" in f]
    volume_name = folder.replace("-","_") + "_pre_ax_resampled_resampled.mha" # Replace with _pre_ax_resampled.mha for UH axial, _pre_cor_resampled for UH cor, and or _cor_resampled.mha for VA
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
        
        # create empty masks for rectal wall and lumen
        src_mask_shape = vol.shape
        num_slices = src_mask_shape[0]
        width = src_mask_shape[1]
        height = src_mask_shape[2]
        
        tumor_mask = np.zeros((width, height))
        lumen_mask = np.zeros((width, height))
        fat_mask = np.zeros((width, height))
        rectal_wall_mask = np.zeros((len(src_mask), width, height))
        dilated5_mask = np.zeros((width, height))
        dilated10_mask = np.zeros((width, height))
        dilated15_mask = np.zeros((width, height))
        new_vol = np.zeros((width, height))
        
        largest_slice_mask = np.zeros((1, width, height))
        
        tumor_area = 0
        dilate5_area = 0
        dilate10_area = 0
        dilate15_area = 0
        fat_area = 0
        tumor_area_slice = 0
        dilate5_kernel = np.ones((5, 5), np.uint8)
        dilate10_kernel = np.ones((10, 10), np.uint8)
        dilate15_kernel = np.ones((15, 15), np.uint8)

        for idx in range(len(src_mask)):
            
            img = src_mask[idx, :, :]
            
            slice_candidate = sum(np.in1d([1,2,4], img)) # Slice must have tumor and fat
            
            if(slice_candidate == 3):
                current_tumor_area = np.count_nonzero( img == 1 )
                if(current_tumor_area > tumor_area):
                    tumor_area = current_tumor_area
                    tumor_area_slice = idx
                    
                    tumor_img = np.copy(img)
                    lumen_img = np.copy(img)
                    fat_img = np.copy(img)
                    rw_img = np.copy(img)
                    zero_img = np.copy(img)
                    
                    
                    tumor_img[tumor_img == 1] = 1
                    tumor_img[tumor_img != 1] = 0
                    tumor_mask[:, :] = tumor_img
                    dilated5_img = cv2.dilate(tumor_img, dilate5_kernel, iterations = 1)
                    dilated5_img[dilated5_img == 1] = 11
                    dilated10_img = cv2.dilate(tumor_img, dilate10_kernel, iterations = 1)
                    dilated10_img[dilated5_img == 1] = 11
                    dilated15_img = cv2.dilate(tumor_img, dilate15_kernel, iterations = 1)
                    dilated15_img[dilated15_img == 1] = 11
                    
                    lumen_img[lumen_img != 2] = 0
                    lumen_img[lumen_img == 2] = 1
                    lumen_mask[:, :] = lumen_img
                    
                    fat_img[fat_img != 4] = 0
                    fat_img[fat_img == 4] = 1
                    fat_mask[:, :] = fat_img
                    fat_area = np.count_nonzero( fat_img == 1 )

                    rw_img[rw_img == 7] = 1
                    rw_img[rw_img == 8] = 1
                    rw_img[rw_img > 1] = 0
                    rectal_wall_mask[idx, :, :] = rw_img

                    #zero_img[zero_img != 0] = 1

                    dilated5_img = img + dilated5_img
                    dilated5_img[dilated5_img != 15] = 0
                    dilated5_img[dilated5_img == 15] = 1
                    dilated5_mask[:, :] = dilated5_img

                    dilated10_img = img + dilated10_img
                    dilated10_img[dilated10_img != 15] = 0
                    dilated10_img[dilated10_img == 15] = 1
                    dilated10_mask[:, :] = dilated10_img

                    dilated15_img = img + dilated15_img
                    dilated15_img[dilated15_img != 15] = 0
                    dilated15_img[dilated15_img == 15] = 1
                    dilated15_mask[:, :] = dilated15_img

                    new_vol[:, :] = vol[idx, :, :]
                    dilate5_area = np.count_nonzero( dilated5_img == 1 )
                    dilate10_area = np.count_nonzero( dilated10_img == 1 )
                    dilate15_area = np.count_nonzero( dilated15_img == 1 )
                    
                    print(folder, ": Slice " + str(idx) + " has the largest tumor area of " + str(tumor_area))
        
        slice_tracker = slice_tracker.append({"Patient" : folder, "Slice_Number" : tumor_area_slice, "Tumor_Size" : tumor_area, "Dilated_Mask_Size" : dilate_area, "Fat_Mask_Size" : fat_area}, ignore_index=True) 
        
        # Check masks
        #binary_lumen = check_mask(lumen_mask)
        #if(binary_lumen is False):
        #    issues.append(mask_file[0] + "_lumen")
        #    raise ValueError("Lumen mask for " + mask_file[0] + " is NOT binary!" )
        #binary_tumor = check_mask(tumor_mask)
        #if(binary_tumor is False):
        #    issues.append(mask_file[0] + "_tumor")
        #    raise ValueError("Rectal wall mask for " + mask_file[0] + " is NOT binary!" )
        #binary_fat = check_mask(fat_mask)
        #if(binary_fat is False):
        #    issues.append(mask_file[0] + "_fat")
        #    raise ValueError("Fat mask for " + mask_file[0] + " is NOT binary!" )
            
        # convert masks into sitk images
        tumor_mask = sitk.GetImageFromArray(tumor_mask)
        lumen_mask = sitk.GetImageFromArray(lumen_mask)
        fat_mask = sitk.GetImageFromArray(fat_mask)
        rectal_wall_mask = sitk.GetImageFromArray(rectal_wall_mask)
        dilated5_mask = sitk.GetImageFromArray(dilated5_mask)
        dilated10_mask = sitk.GetImageFromArray(dilated10_mask)
        dilated15_mask = sitk.GetImageFromArray(dilated15_mask)
        new_vol = sitk.GetImageFromArray(new_vol)
        
        # Create masks directory
        if not os.path.isdir(os.path.join(directory, folder, "masks")):
            os.mkdir(os.path.join(directory, folder, "masks"))
        
        # Save tumor mask
        out_name = mask_file[0].replace(".mha","_tumor_ls.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(tumor_mask, spacing, origin, output)
        
        # Save lumen mask
        out_name = mask_file[0].replace(".mha","_lumen_ls.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(lumen_mask, spacing, origin, output)
        
        # Save fat mask
        out_name = mask_file[0].replace(".mha","_fat_ls.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(fat_mask, spacing, origin, output)

        # Save rectal wall mask
        out_name = mask_file[0].replace(".mha","_rw_ls.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(rectal_wall_mask, spacing, origin, output)

        # Save proximal fat mask 5 pixels
        out_name = mask_file[0].replace(".mha","_proxfat5_ls.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(dilated5_mask, spacing, origin, output)

        # Save proximal fat mask 10 pixels
        out_name = mask_file[0].replace(".mha","_proxfat10_ls.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(dilated10_mask, spacing, origin, output)

        # Save proximal fat mask 15 pixels
        out_name = mask_file[0].replace(".mha","_proxfat15_ls.mha")
        output = os.path.join(directory, folder, "masks", out_name)
        write_mha(dilated15_mask, spacing, origin, output)
        
        # Save scan
        out_name = volume_name.replace("_pre_ax_resampled_resampled", "_pre_ax_resampled_ls")
        output = os.path.join(directory, folder, out_name)
        write_mha(new_vol, spacing, origin, output)
        
    except Exception as e:
        print(e)
        sys.exit(1)
        
    
        # Extract rectal wall, lumen, and fat from each slice
    
slice_tracker.to_excel("LTS_Problems_Tracker.xlsx", index = False)