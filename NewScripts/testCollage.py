import medviz as viz
from pathlib import Path

viz.preprocess.match_image_masks(
        images_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Coronal_Image_NPY",
        masks_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Coronal_Fat_NPY",
        image_id_func = lambda x: x.split("-")[-1][:3],
        mask_id_func = lambda x: x.split("-")[-1][:3],
        save_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/CoronalCollageFeatures/Coronal_Fat_mask_match.csv",

    )

viz.compute_stats_collage(
    csv_path="/Users/leobao/Documents/MultiPlanePipeline/Data/CoronalCollageFeatures/Coronal_Fat_mask_match.csv",
    stats_save_path="/Users/leobao/Documents/MultiPlanePipeline/Data/CoronalCollageFeatures",
)
