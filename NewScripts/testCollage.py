import medviz as viz
from pathlib import Path

viz.preprocess.match_image_masks(
        images_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Axial_Image_NIFTY",
        masks_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Axial_Fat_NIFTY",
        image_id_func = lambda x: x.split("-")[-1][:3],
        mask_id_func = lambda x: x.split("-")[-1][:3],
        save_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/AxialCollageFeatures/Axial_Fat_mask_match.csv",

    )

viz.compute_stats_collage(
    csv_path="/Users/leobao/Documents/MultiPlanePipeline/Data/AxialCollageFeatures/Axial_Fat_mask_match.csv",
    stats_save_path="/Users/leobao/Documents/MultiPlanePipeline/Data/AxialCollageFeatures",
)
