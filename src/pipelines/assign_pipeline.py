"""
Inference pipeline for morphological clustering
"""
# pylint:disable=no-member
from datetime import datetime

import hydra
from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.image_analysis.clustering.feature_extraction import extract_image_features
from src.image_analysis.clustering.spot_clustering import (
    assign_celltype,
    assign_probabilistic_celltype,
    assign_random_celltype,
    create_barcode_dict,
)
from src.image_analysis.clustering.transform import dimensionality_reduction, scaling
from src.pipelines.dataio_pipeline import DataIO
from src.utils.data_loader import load_images, load_pos_list, load_segmentation
from src.utils.data_saver import save_data
from src.utils.logger import Logger

logger = Logger()
load_dotenv()


class SpotCluster:
    """
    Python class for morphological clustering
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "assign"

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
        """Load data as needed for img seg"""
        logger.info("<load_data> load truncated master table")
        self.unit_test_master_table = load_pos_list(
            cache_path + "unit_test_master_table.csv"
        )
        self.master_table = self.unit_test_master_table
        # self.master_table = unit_test_master_table.iloc[1:, :]
        logger.info("<load segmentation> load cropped segmentation")
        self.img_arr = np.zeros
        self.mask = load_segmentation(cache_path + "unit_test_seg.npy")
        self.img_arr = np.zeros(self.mask.shape)
        # logger.info("<load_data> load tissue position table")
        # self.tissue_pos_csv = load_pos_list(self.configs["paths"]["unittest_pos_path"])

    # TODO: can improve by dynamically creating variables
    def load_data(self, prep_dir, results_dir, config_flow):
        """Load mask"""
        logger.info("<load_data> load the wsi image")
        imgseg_method = config_flow["imgseg"]["model"]["name"]
        self.img_arr = load_images(prep_dir + "wsi.tif")
        logger.info("<load_data> load the mask")
        self.mask = load_segmentation(
            results_dir + f"{imgseg_method}_cell_segmentation_full.npy"
        )

    def extract_features(self):
        """Extract features for clustering"""
        logger.info("<fit_cluster> extract image features")
        self.master_table = extract_image_features(self.mask, self.img_arr)

    def preprocess_data(self):
        """Perform scaling and dimensionality reduction"""
        logger.info("<fit_cluster> scale training data")
        scaled_training_data = scaling(self.master_table)
        logger.info("<fit_cluster> reduce dimensions -> pca, reduced_df")
        pca, self.reduced_df = dimensionality_reduction(scaled_training_data)

    def spot_cluster(self, results_dir, deconv_method="cell2location", unittest=False):
        "Perform Clustering by spot and assign cell types to corresponding cells"
        if unittest:
            logger.info("<Barcode dictionary> creating unit test barcode dictionary")
            barcode_dict = create_barcode_dict(
                results_dir + "unit_test_cell2location_results.csv"
            )
            if self.model_name == "local":
                logger.info(
                    "<Spot clustering> cell type  assigning with local clustering"
                )
                self.barcodes = assign_celltype(
                    results_dir + "unit_test_cells_df.csv",
                    self.reduced_df,
                    barcode_dict,
                )
            elif self.model_name == "global":
                logger.info(
                    "<Spot clustering> cell type assigning with probabilistic algorithm"
                )
                self.barcodes = assign_probabilistic_celltype(
                    results_dir + "unit_test_cells_df.csv",
                    self.reduced_df,
                    barcode_dict,
                    self.configs["params"]["n_clusters"],
                    self.configs["params"]["clustering_batch_size"],
                    self.configs["params"]["random_state"],
                )
            elif self.model_name == "naive":
                logger.info(
                    "<Spot clustering> cell type assigning with naive algorithm"
                )
                self.barcodes = assign_random_celltype(
                    results_dir + "unit_test_cells_df.csv",
                    self.reduced_df,
                    barcode_dict,
                    self.configs["params"]["random_state"],
                )
            else:
                logger.error(
                    f"<Spot clustering>  - Algorithm {self.model_name} not supported"
                )
        else:
            logger.info("<Barcode dictionary> creating barcode dictionary")
            barcode_dict = create_barcode_dict(
                results_dir + f"{deconv_method}_results.csv"
            )
            if self.model_name == "local":
                logger.info(
                    "<Spot clustering> cell type  assigning with local clustering"
                )
                self.barcodes = assign_celltype(
                    results_dir + "cells_df.csv",
                    self.reduced_df,
                    barcode_dict,
                )
            elif self.model_name == "global":
                logger.info(
                    "<Spot clustering> cell type assigning with probabilistic algorithm"
                )
                self.barcodes = assign_probabilistic_celltype(
                    results_dir + "cells_df.csv",
                    self.reduced_df,
                    barcode_dict,
                    self.configs["params"]["n_clusters"],
                    self.configs["params"]["clustering_batch_size"],
                    self.configs["params"]["random_state"],
                )
            elif self.model_name == "naive":
                logger.info(
                    "<Spot clustering> cell type assigning with naive algorithm"
                )
                self.barcodes = assign_random_celltype(
                    results_dir + "cells_df.csv",
                    self.reduced_df,
                    barcode_dict,
                    self.configs["params"]["random_state"],
                )
            else:
                logger.error(
                    f"<Spot clustering>  - Algorithm {self.model_name} not supported"
                )

    def save_data(self, results_dir):
        """Save data from clustering"""
        logger.info("<save-data> save assigned cell types")
        save_data(
            "csv",
            results_dir,
            self.barcodes,
            "spot_clusters.csv",
        )
        logger.info("<save-data> DONE---")


if __name__ == "__main__":
    exp_id = "V10S24-023A"
    user_id = "Ana"
    flow_id = "12345"
    flow_id = str(datetime.timestamp(datetime.now())).replace(".", "")
    hydra.core.global_hydra.GlobalHydra.instance().clear()
    hydra.initialize(config_path="../../conf")
    config_flow = OmegaConf.to_container(
        hydra.compose(config_name="visium_config_flow")
    )
    logger.info(f"exp_id: {exp_id}, user_id: {user_id}")
    data_io = DataIO(
        exp_id=exp_id,
        user_id=user_id,
        flow_id=flow_id,
    )
    assign_celltype = SpotCluster()
    assign_celltype.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[assign_celltype.pipeline_name]
    )

    # logger.info(f"[{morph_cluster.pipeline_name}] Fetch relevant files")
    # data_io.fetch_input_files(pipeline_configs=morph_cluster.configs)

    assign_celltype.load_data(
        prep_dir=data_io.prep_dir,
        results_dir=data_io.results_dir,
        config_flow=config_flow,
    )
    assign_celltype.extract_features()
    assign_celltype.preprocess_data()
    assign_celltype.spot_cluster(
        results_dir=data_io.results_dir,
        deconv_method=config_flow["celldeconv"]["model"]["name"],
    )
    assign_celltype.save_data(results_dir=data_io.results_dir)

    # # Push all output files from local cache to s3
    # data_io.push_output_files(pipeline_configs=assign_celltype.configs)
