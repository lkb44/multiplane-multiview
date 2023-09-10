import medviz as viz
from pathlib import Path

image_path = '/Users/leobao/Documents/MultiPlanePipeline/CollageMAT/Patient-011/Patient-011_cor_ls.npy'
mask_path = '/Users/leobao/Documents/MultiPlanePipeline/CollageMAT/Patient-011/Patient-011_cor_label_fat_ls.npy'

image, mask = viz.read_image_mask(image = image_path, mask = mask_path)

viz.feats.collage2d(image, mask, window_sizes = [3, 5, 7, 9, 11], save_path = '/Users/leobao/Documents/MultiPlanePipeline/CollageReadytoReshape', out_name = 'Patient-011_cor')