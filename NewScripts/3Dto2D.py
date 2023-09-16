import numpy as np

# Load the 3D numpy array from the .npy file
filename = '/Users/leobao/Documents/MultiPlanePipeline/CollageMAT/Patient-049/Patient-049_cor_label_fat_ls.npy'
array_3d = np.load(filename)

# Check the shape of the 3D array
shape = array_3d.shape

# Ensure the input is 1 x n x n
if len(shape) != 3 or shape[0] != 1:
    raise ValueError("Input array should have shape (1, n, n)")

# Extract the 2D array from the 3D array
array_2d = array_3d[0]

# Save the 2D array back over the original file
np.save(filename, array_2d)

print("Conversion complete and saved to the original file:", filename)
