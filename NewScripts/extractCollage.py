import medviz as viz
from pathlib import Path
import numpy as np


image_path = '/Users/leobao/Documents/MultiPlanePipeline/CollageMAT/Patient-012/Patient-012_cor_ls.npy'
mask_path = '/Users/leobao/Documents/MultiPlanePipeline/CollageMAT/Patient-012/Patient-012_cor_label_proxfat10_ls.npy'

# Load the two .npy arrays
mask = np.load(mask_path)
image = np.load(image_path)

# Get the dimensions of each array
mask_dim = mask.shape
image_dim = image.shape

# Check if the dimensions match
if mask_dim == image_dim:
    print("The dimensions of both arrays match.")
    image, mask = viz.read_image_mask(image = image_path, mask = mask_path)
    viz.feats.collage2d(image, mask, window_sizes = [3, 5, 7, 9, 11], save_path = '/Users/leobao/Documents/MultiPlanePipeline/CollageReadytoReshape', out_name = 'Patient-012_cor')
else:
    print("The dimensions of the arrays do not match. Attempting to convert 3D to 2D.")

    array_3d = np.load(mask_path)
    # Extract the 2D array from the 3D array
    array_2d = array_3d[0]

    # Save the 2D array back over the original file
    np.save(mask_path, array_2d)

    print("Conversion complete and saved to the original file. Now attempting to extract features")
    image, mask = viz.read_image_mask(image = image_path, mask = mask_path)
    viz.feats.collage2d(image, mask, window_sizes = [3, 5, 7, 9, 11], save_path = '/Users/leobao/Documents/MultiPlanePipeline/CollageReadytoReshape', out_name = 'Patient-012_cor')