import medviz as viz

viz.match_image_mask(
        images_path = "/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/lkb44/Multiplane-Multiview/Multiplane Data/Test/RectalCA-001/RectalCA_001_pre_ax_resampled_largest_slice.mha",
        masks_path = "/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/lkb44/Multiplane-Multiview/Multiplane Data/Test/RectalCA-001/masks/RectalCA_001_pre_ax_label_resampled_aniresampled_tumor_largest_slice.mha",
        image_func = lambda x: x.split("_")[1],
        mask_func = lambda x: x.split("_")[1],
        save_path = "/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/lkb44/Multiplane-Multiview/Multiplane Data/Test/RectalCA-001/tumor_ax_collage.csv",

    )


