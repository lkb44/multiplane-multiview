import SimpleITK as sitk
import cv2
import numpy as np


# Read .mha file using SimpleITK
image = sitk.ReadImage("/Users/leobao/Downloads/RectalCA_011_pre_ax_label_resampled_aniresampled_fat_largest_slice.mha")
image = sitk.GetArrayFromImage(image)

# Convert the image data to a numpy array
image = np.array(image, dtype=np.uint8)

# Define the structuring element for dilation
kernel = np.ones((5, 5), np.uint8)

# Apply dilation to the image
dilated_image = cv2.dilate(image, kernel, iterations=1)

# Show the original and dilated images
cv2.imshow("Original", image)
cv2.imshow("Dilated", dilated_image)