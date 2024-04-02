# """
# inference pipeline for cell segmentation
# """

# import hydra
# import numpy as np
# from omegaconf import OmegaConf

# from src.data_prep.postprocessing import (
#     match_input_shape,
#     stitching_instance_segmentation,
# )
# from src.data_prep.preprocessing import Preprocessing
# from src.utils.data_loader import load_images, load_pos_list
# from src.utils.data_saver import save_data
# from src.utils.model_loader import cellpose_model, predict


# def predict_pipeline() -> None:
#     """Predict pipeline (Pipeline to load an H&E image and output a segmented mask)"""
#     hydra.core.global_hydra.GlobalHydra.instance().clear()
#     hydra.initialize(config_path="../../../conf")
#     params = OmegaConf.to_container(hydra.compose(config_name="config"))

#     # Loading configurations for image_seg
#     image_seg = params["image_segmentation"]
#     # Load Image
#     # img_arr = data_loader.load_images(image_seg["paths"]["image_path"])
#     img_arr = load_images(image_seg["paths"]["image_path"])
#     # Preprocess Image
#     # tissue_pos_csv = data_loader.load_pos_list(image_seg["paths"]["tissue_pos_path"])
#     tissue_pos_csv = load_pos_list(image_seg["paths"]["tissue_pos_path"])
#     image_processor = Preprocessing(
#         img_arr, image_seg["params"]["patch_size"], image_seg["params"]["overlap"]
#     )
#     detected_tissue_region = image_processor.detect_tissue_bbox(tissue_pos_csv)
#     cropped_region = image_processor.crop_tissue_region(img_arr, detected_tissue_region)
#     normalized_image = image_processor.norm_stain_remove(cropped_region)
#     downsampled_image = image_processor.downsample(
#         normalized_image, image_seg["params"]["downsample_factor"]
#     )
#     resized_image_dim = image_processor.resize_img_dim(downsampled_image.shape)
#     resized_image = image_processor.resize_img(
#         img_arr, downsampled_image.shape, resized_image_dim
#     )
#     patched_arr = image_processor.generate_patches(
#         resized_image, image_seg["params"]["n_channels"]
#     )
#     ## Model Inference
#     model = cellpose_model(image_seg["params"]["model_type"])
#     masks, flows, styles, diams = predict(
#         model,
#         patched_arr,
#         image_seg["params"]["batch_size"],
#         image_seg["params"]["diameter"],
#         image_seg["params"]["channels"],
#     )
#     ## PostProcess
#     full_mask_arr = stitching_instance_segmentation(
#         np.asarray(masks),
#         (image_seg["params"]["patch_size"], image_seg["params"]["patch_size"]),
#         image_seg["params"]["overlap"],
#         (resized_image.shape[0], resized_image.shape[1]),
#         minimum_pixels=6,
#     )
#     final_mask, upsampled_mask = match_input_shape(
#         full_mask_arr,
#         image_processor.img_shape,
#         detected_tissue_region,
#     )
#     save_data("numpy", image_seg["paths"]["home_path"], full_mask_arr, "wsi_seg")
#     save_data("numpy", image_seg["paths"]["home_path"], final_mask, "full_size_wsi_seg")


# if __name__ == "__main__":
#     predict_pipeline()
