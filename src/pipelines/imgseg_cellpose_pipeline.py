"""
inference pipeline for cell segmentation
"""
# pylint:disable=no-member

import os
from datetime import datetime

import hydra
import numpy as np
from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.data_prep.postprocessing import (
    match_input_shape,
    stitching_instance_segmentation,
)
from src.data_prep.preprocessing import Preprocessing
from src.pipelines.dataio_pipeline import DataIO
from src.utils.data_loader import load_images, load_pos_list, load_segmentation
from src.utils.data_saver import save_data
from src.utils.logger import Logger
from src.utils.model_loader import cellpose_model, predict

load_dotenv()
logger = Logger()


class ImageSeg:
    """
    Python class for image segmentation
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "imgseg"

    def load_model_configs_from_flow(self, pipeline_configs_from_flow):
        """load_model_configs_from_flow(pipeline_configs)

        Args:
            model_configs (_type_): _description_
        """
        self.model_name = pipeline_configs_from_flow["model"]["name"]
        self.model_version = pipeline_configs_from_flow["model"]["version"]
        self.config_key = self.pipeline_name + "." + self.model_name
        self.configs = OmegaConf.to_container(hydra.compose(config_name="config_tree"))[
            self.pipeline_name
        ][self.config_key]
        if "params" in pipeline_configs_from_flow["model"]:
            self.configs["params"] = pipeline_configs_from_flow["model"]["params"]

    def load_unittest_data(self, cache_path):
        """load_unittest_data
        Load unittest data from S3
        """
        logger.info("<load_data> load image array")
        self.resized_image = load_segmentation(cache_path + "unit_test_image.npy")
        logger.info("<load_data> load patched image array")
        self.patched_arr = load_segmentation(cache_path + "unit_test_image_patched.npy")
        logger.info("<load_data> load segmentation array")
        self.seg_arr = load_segmentation(cache_path + "unit_test_seg.npy")

    def load_data(self, prep_dir):
        """Load data as needed for img seg"""
        # TODO: need parameterize this
        # for file_name in self.configs["inputs"].keys():
        #     file_ext = pipeline_configs["inputs"][input]["extension"]
        #     file_path = cache_path + file_name + "." + file_ext

        logger.info("<load_data> load image array")
        self.img_arr = load_images(prep_dir + "wsi.tif")

        logger.info("<load_data> load tissue position table")
        self.tissue_pos_csv = load_pos_list(prep_dir + "tissue_positions_list.csv")

    def pre_process(self):
        """Preporcessing of image"""
        self.image_processor = Preprocessing(
            self.img_arr,
            self.configs["params"]["patch_size"],
            self.configs["params"]["overlap"],
        )
        logger.info("<pre_process> detected tissue")
        self.detected_tissue_region = self.image_processor.detect_tissue_bbox(
            self.tissue_pos_csv
        )
        # Padding the image all around with patch_size pixels
        self.padded_img = self.image_processor.pad_all_around(self.img_arr)
        # Shifting the tissue boundary to account for the padding
        self.adjusted_tissue_boundary = self.image_processor.adjust_tissue_boundary(
            self.detected_tissue_region
        )

        logger.info("<pre_process> cropping region")
        cropped_region = self.image_processor.crop_tissue_region(
            self.padded_img,
            self.adjusted_tissue_boundary,
            self.configs["params"]["patch_size"],
        )
        logger.info("<pre_process> remove background artifact")
        self.cleaned_tissue = self.image_processor.tissue_detector(cropped_region)
        logger.info("<pre_process> normalize image")

        normalized_image = self.image_processor.norm_stain_remove(self.cleaned_tissue)
        logger.info("<pre_process> downsample image")
        downsampled_image = self.image_processor.downsample(
            normalized_image, self.configs["params"]["downsample_factor"]
        )
        logger.info("<pre_process> resize image")
        resized_image_dim = self.image_processor.resize_img_dim(downsampled_image.shape)
        self.resized_image = self.image_processor.resize_img(
            downsampled_image, resized_image_dim, resized_image_dim
        )
        logger.info("<pre_process> get patched array")
        self.patched_arr = self.image_processor.generate_patches(
            self.resized_image, self.configs["params"]["n_channels"]
        )
        logger.info("<pre_process> DONE---")

    def predict(self, gpu_enabled=True):
        """Run inference for cell segmentation"""

        logger.info("<predict> get model")
        self.model = cellpose_model(
            self.configs["params"]["model_type"], gpu_enabled=gpu_enabled
        )
        logger.info("<predict> predict")
        self.masks, self.flows, self.styles, self.diams = predict(
            self.model,
            self.patched_arr,
            self.configs["params"]["batch_size"],
            self.configs["params"]["diameter"],
            self.configs["params"]["channels"],
            self.configs["params"]["flow_threshold"],
        )
        logger.info("<predict> DONE---")

    def post_process(self, exp_id):
        """Post processing"""

        logger.info("<post-process> stitch segmentation")
        self.full_mask_arr = stitching_instance_segmentation(
            np.asarray(self.masks),
            (
                self.configs["params"]["patch_size"],
                self.configs["params"]["patch_size"],
            ),
            self.configs["params"]["overlap"],
            (self.resized_image.shape[0], self.resized_image.shape[1]),
            minimum_pixels=6,
        )
        if exp_id != "unittest":
            logger.info("<post-process> match input shape")
            self.final_mask, self.upsampled_mask = match_input_shape(
                self.full_mask_arr,
                self.padded_img.shape,
                self.adjusted_tissue_boundary,
                self.cleaned_tissue.shape,
                self.configs["params"]["downsample_factor"],
            )
            # Undo padding
            logger.info("<post-process> Undo padding")
            self.final_mask = self.image_processor.undo_padding(self.final_mask)
            logger.info(
                f"<post-process> Number of detected cells: {self.final_mask.max()}"
            )

        # im.save("cells_layer.png")
        logger.info("<post-process> DONE---")

    # do we need to save data? can we pass directly to next flow?
    def save_data(self, results_dir):
        """Save data outputs"""
        logger.info("<save-data> save (wsi_seg, full_size_wsi_seg)")
        save_data(
            "numpy",
            results_dir,
            self.full_mask_arr,
            "cellpose_cell_segmentation_crop.npy",
        )
        save_data(
            "numpy",
            results_dir,
            self.final_mask,
            "cellpose_cell_segmentation_full.npy",
        )
        logger.info("<save-data> DONE---")


if __name__ == "__main__":
    exp_id = "ZZ268953_V11T09-086D"
    user_id = "sunaal"
    # flow_id = "12345"
    flow_id = str(datetime.timestamp(datetime.now())).replace(".", "")
    hydra.core.global_hydra.GlobalHydra.instance().clear()
    hydra.initialize(config_path="../../data/conf")
    config_flow = OmegaConf.to_container(
        hydra.compose(config_name="visium_config_flow")
    )
    logger.info(f"exp_id: {exp_id}, user_id: {user_id}")
    data_io = DataIO(
        exp_id=exp_id,
        user_id=user_id,
        flow_id=flow_id,
    )
    image_seg = ImageSeg()
    image_seg.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[image_seg.pipeline_name]
    )

    # logger.info(f"[{image_seg.pipeline_name}] Fetch relevant files")
    # data_io.fetch_input_files(pipeline_configs=image_seg.configs)

    logger.info(f"[{image_seg.pipeline_name}] Start pipeline")
    image_seg.load_data(prep_dir=data_io.prep_dir)
    image_seg.pre_process()
    image_seg.predict()
    image_seg.post_process(exp_id=exp_id)
    image_seg.save_data(results_dir=data_io.results_dir)

    # data_io.push_output_files(pipeline_configs=image_seg.configs)
