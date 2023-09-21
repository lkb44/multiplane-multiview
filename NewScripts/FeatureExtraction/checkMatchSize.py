import medviz as viz
from pathlib import Path
import numpy as np

# Load the two .npy arrays
mask = np.load('/Users/leobao/Documents/MultiPlanePipeline/CollageMAT/Patient-049/Patient-049_cor_label_fat_ls.npy')
image = np.load('/Users/leobao/Documents/MultiPlanePipeline/CollageMAT/Patient-049/Patient-049_cor_ls.npy')

# Get the dimensions of each array
mask_dim = mask.shape
image_dim = image.shape

# Print the dimensions of each array
print(f"Dimensions of mask: {mask_dim}")
print(f"Dimensions of image: {image_dim}")

# Check if the dimensions match
if mask_dim == image_dim:
    print("The dimensions of both arrays match.")
else:
    print("The dimensions of the arrays do not match.")