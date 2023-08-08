import numpy as np
import matplotlib.pyplot as plt

def load_and_display_volume(image_path, mask_path):
    image_array = np.load(image_path)
    mask_array = np.load(mask_path)

    # Display the dimensions of the image and mask arrays
    print("Image Dimensions:", image_array.shape)
    print("Mask Dimensions:", mask_array.shape)

    # Display the image slice
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.imshow(image_array, cmap='gray')
    ax1.set_title("Image Slice")
    ax1.axis('off')

    # Display the mask slice
    ax2.imshow(mask_array, cmap='jet')  # Using 'jet' colormap for better mask visualization
    ax2.set_title("Mask Slice")
    ax2.axis('off')

    plt.tight_layout()
    plt.show()
    
if __name__ == "__main__":
    image_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Coronal_Image_NPY/Patient-002_cor_ls.npy"
    mask_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Coronal_Fat_NPY/Patient-002_cor_label_fat_ls.npy"

    load_and_display_volume(image_path, mask_path)