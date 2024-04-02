"""
Image segmentation Metaflow pipeline
"""
import json
import os
from datetime import datetime

# pylint:disable=no-member, no-pylint
import hydra
import pytz
import scanpy as sc
from omegaconf import OmegaConf

from src.pipelines.assign_pipeline import SpotCluster
from src.pipelines.cell2spot_pipeline import Cell2Spot
from src.pipelines.celldeconv_card_pipeline import CARD
from src.pipelines.celldeconv_cell2location_pipeline import Cell2Location
from src.pipelines.clustering_pipeline import MorphCluster
from src.pipelines.dataio_pipeline import DataIO
from src.pipelines.datamerge_pipeline import DataMerge
from src.pipelines.imgseg_cellpose_pipeline import ImageSeg
from src.pipelines.imgseg_hovernet_pipeline import HoverNet

from src.pipelines.param_pipeline import ParamPipeline
from src.pipelines.qc_pipeline import QC
from src.pipelines.spatial_one_flow import SpatialOneFlow
from src.pipelines.spatialanalysis_pipeline import SpatialAnalysis
from src.utils.env_variables import Env
from src.utils.logger import Logger

# Note: Always load the environment variables before importing MetaFlow
# since some of these variable are required by MetaFlow for proper
# initialization
e = Env()
from metaflow import (  # pylint:disable=no-name-in-module, wrong-import-position, wrong-import-order, line-too-long; kubernetes,; schedule,
    Parameter,
    Run,
    current,
    project,
    step,
)

logger = Logger()


@project(name=e.project_name)
# @schedule(cron=e.schedule)  # schedule pipeline (import schedule)
class VisiumFlow(SpatialOneFlow):
    """
    Run the spatial one pipeline in "local" mode in your workbench with the following command (configs can be modified in conf/config_flow.yaml)
    e.g. PYTHONPATH=${PYTHONPATH}:. python3 -m src.pipelines.visium_flow run --max-workers 1
    """

    @step
    def start(self):
        """Start
        :param None(rtype):
        :return None(rtype):
        """
        self.spatial_one_flow_start()
        self.next(self.get_params, foreach="exp_ids")

    @step
    def get_params(self):
        """Upload config_flow to analysis
        :param None(rtype):
        :return None(rtype):
        """
        self.spatial_one_flow_get_params()
        self.next(self.image_seg)

    @step
    def image_seg(self):
        """Image segmentation
        :param None(rtype):
        :return None(rtype):
        """
        image_seg = ImageSeg()
        if (self.input is not None) and (
            self.pipelines_enabled.get(image_seg.pipeline_name, False)
        ):
            exp_id = self.input
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )
            # Check for model type
            if self.configs["imgseg"]["model"]["name"] == "cellpose":
                logger.info("Setting up cellpose -----------")
                image_seg = ImageSeg()
                # get pipeline_confgis from config_flow
                image_seg.load_model_configs_from_flow(
                    pipeline_configs_from_flow=self.configs[image_seg.pipeline_name]
                )

                logger.info(f"[{image_seg.pipeline_name}] Start pipeline")
                image_seg.load_data(prep_dir=data_io.prep_dir)
                image_seg.pre_process()
                image_seg.predict()
                image_seg.post_process(exp_id=exp_id)
                image_seg.save_data(results_dir=data_io.results_dir)

            elif self.configs["imgseg"]["model"]["name"] == "hovernet":
                logger.info("Setting up data_hovernet -----------")
                image_seg = HoverNet()
                image_seg.load_model_configs_from_flow(
                    pipeline_configs_from_flow=self.configs[image_seg.pipeline_name]
                )
                logger.info(
                    f"[{image_seg.pipeline_name}] POST predict API request Hovernet"
                )
                image_seg.predict(data_io)
        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")

        self.next(self.cell_to_spot)

    @step
    def cell_to_spot(self):
        """Cell to spot
        :param None(rtype):
        :return None(rtype):
        """
        cell2spot = Cell2Spot()
        if (self.input is not None) and (
            self.pipelines_enabled.get(cell2spot.pipeline_name, False)
        ):
            exp_id = self.input
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )

            cell2spot.load_model_configs_from_flow(
                pipeline_configs_from_flow=self.configs[cell2spot.pipeline_name]
            )

            cell2spot.load_data(
                prep_dir=data_io.prep_dir,
                results_dir=data_io.results_dir,
                config_flow=self.configs,
            )
            cell2spot.convert_mask_to_polygons()
            cell2spot.generate_visium_spot_polygons()
            cell2spot.create_cell_search_tree()
            cell2spot.get_cells_per_spot()
            cell2spot.save_data(results_dir=data_io.results_dir)

        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")

        self.next(self.cell_deconv)

    @step
    def cell_deconv(self):
        """Cell deconv
        :param None(rtype):
        :return None(rtype):
        """
        celldeconv = (
            Cell2Location()
        )  # create any celldeconv class to get the pipeline_name dynamically
        if (self.input is not None) and (
            self.pipelines_enabled.get(celldeconv.pipeline_name, False)
        ):
            exp_id = self.input
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )
            # Check model type
            if self.configs["celldeconv"]["model"]["name"] == "card":
                card_celldeconv = CARD()
                logger.info("CARD")
                # Load configs from flow
                card_celldeconv.load_model_configs_from_flow(
                    pipeline_configs_from_flow=self.configs[
                        card_celldeconv.pipeline_name
                    ]
                )
                logger.info(
                    f"[{card_celldeconv.pipeline_name}] POST predict API request CARD"
                )
                respone = False
                response = card_celldeconv.predict(data_io=data_io)

                # data_io.poll_for_update(file_path="card_results.csv")
            # Cell2location deconv
            elif self.configs["celldeconv"]["model"]["name"] == "cell2location":
                cell2loc = Cell2Location()
                # Load configs from flow
                cell2loc.load_model_configs_from_flow(
                    pipeline_configs_from_flow=self.configs[cell2loc.pipeline_name]
                )

                # If it cell_sig exists and there is no retraining of cell_sig requried, then load cell_sig from S3
                if cell2loc.configs["params"]["retrain_cell_sig"] == False:
                    logger.info(
                        f"[{cell2loc.pipeline_name}] Fetch trained cell signature data from reference directory"
                    )
                    cell_sig = sc.read_h5ad(
                        data_io.reference_dir
                        + cell2loc.configs["params"]["atlas_type"]
                        + "_cell_sig.h5ad",
                    )
                else:
                    logger.info(
                        f"[{cell2loc.pipeline_name}] Fetch single cell atlas data from s3 for training cell signature"
                    )
                    cell_sig = cell2loc.save_train_sig(
                        reference_dir=data_io.reference_dir
                    )  ### train sc data and save to reference directory

                logger.info(f"[{cell2loc.pipeline_name}] Start pipeline")
                cell2loc.load_data(
                    prep_dir=data_io.prep_dir,
                    results_dir=data_io.results_dir,
                    reference_dir=data_io.reference_dir,
                )

                cell2loc.train_location_model(cell_sig)
                adata_st = cell2loc.get_cell_type_proportion()
                cell2loc.save_data(adata_st, results_dir=data_io.results_dir)

        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")

        self.next(self.morph_cluster)

    @step
    def morph_cluster(self):
        """Morphological clustering
        :param None(rtype):
        :return None(rtype):
        """
        morph_cluster = MorphCluster()
        if (self.input is not None) and (
            self.pipelines_enabled.get(morph_cluster.pipeline_name, False)
        ):
            exp_id = self.input
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )

            morph_cluster.load_model_configs_from_flow(
                pipeline_configs_from_flow=self.configs[morph_cluster.pipeline_name]
            )

            morph_cluster.load_data(
                prep_dir=data_io.prep_dir,
                results_dir=data_io.results_dir,
                config_flow=self.configs,
            )
            morph_cluster.extract_features()
            morph_cluster.preprocess_data()
            morph_cluster.fit_cluster()
            morph_cluster.save_data(results_dir=data_io.results_dir)
        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")

        self.next(self.assign_celltype)

    @step
    def assign_celltype(self):
        """Cell type assignation
        :param None(rtype):
        :return None(rtype):
        """
        assign_celltype = SpotCluster()
        if (self.input is not None) and (
            self.pipelines_enabled.get(assign_celltype.pipeline_name, False)
        ):
            exp_id = self.input
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )

            assign_celltype.load_model_configs_from_flow(
                pipeline_configs_from_flow=self.configs[assign_celltype.pipeline_name]
            )

            assign_celltype.load_data(
                prep_dir=data_io.prep_dir,
                results_dir=data_io.results_dir,
                config_flow=self.configs,
            )

            assign_celltype.extract_features()
            assign_celltype.preprocess_data()
            assign_celltype.spot_cluster(
                results_dir=data_io.results_dir,
                deconv_method=self.configs["celldeconv"]["model"]["name"],
            )

            assign_celltype.save_data(results_dir=data_io.results_dir)

        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")

        self.next(self.spot_qc)

    @step
    def spot_qc(self):
        """Spot-Level QC
        :param None(rtype):
        :return None(rtype):
        """
        qc = QC()
        if (self.input is not None) and (
            self.pipelines_enabled.get(qc.pipeline_name, False)
        ):
            exp_id = self.input
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )
            qc.load_model_configs_from_flow(
                pipeline_configs_from_flow=self.configs[qc.pipeline_name]
            )

            qc.load_data(prep_dir=data_io.prep_dir, results_dir=data_io.results_dir)
            qc.get_qc_metrics(exp_id=exp_id)
            qc.save_data(results_dir=data_io.results_dir)
        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")

        self.next(self.data_merge)

    @step
    def data_merge(self):
        """Data merge
        :param None(rtype):
        :return None(rtype):
        """
        data_merger = DataMerge()
        if (self.input is not None) and (
            self.pipelines_enabled.get(data_merger.pipeline_name, False)
        ):
            exp_id = self.input
            logger.info(f"Flow-id for {exp_id}: {self.flow_id}")
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )

            data_merger.load_model_configs_from_flow(
                pipeline_configs_from_flow=self.configs[data_merger.pipeline_name],
                tissue_type=self.experimental_params["tissue_type"],
                species_type=self.experimental_params["species"],
            )

            logger.info(f"[{data_merger.pipeline_name}] Start pipeline")
            data_merger.load_data(
                prep_dir=data_io.prep_dir,
                results_dir=data_io.results_dir,
                config_flow=self.configs,
                pipelines_enabled=self.pipelines_enabled,
            )
            data_merger.create_cell_outlines_layer(self.pipelines_enabled)
            data_merger.create_visium_spots_layer(self.pipelines_enabled)
            data_merger.filter_gene_data()
            data_merger.create_cells_adata()
            data_merger.create_spots_adata()
            data_merger.summarize_spots_adata_to_dataframe()
            data_merger.summarize_cells_adata_to_dataframe()
            data_merger.create_tmap(
                data_io.results_dir,
                exp_id,
                data_merger.configs["params"]["annotation_file"],
            )
            data_merger.save_data(
                results_dir=data_io.results_dir,
                pipelines_enabled=self.pipelines_enabled,
            )

            data_io.save_run_configs(configs=self.configs)
        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")
        self.next(self.spatial_analysis)

    @step
    def spatial_analysis(self):
        """Spatial analysis
        :param None(rtype):
        :return None(rtype):
        """
        spatanalysis = SpatialAnalysis()
        if (self.input is not None) and (
            self.pipelines_enabled.get(spatanalysis.pipeline_name, False)
        ):
            exp_id = self.input
            flow_id = self.flow_id
            logger.info(f"Run-id for {exp_id}: {self.flow_id}")
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )

            spatanalysis.load_model_configs_from_flow(
                pipeline_configs_from_flow=self.configs[spatanalysis.pipeline_name]
            )
            logger.info(f"[{spatanalysis.pipeline_name}] Start pipeline")
            spatanalysis.load_data(
                prep_dir=data_io.prep_dir, results_dir=data_io.results_dir
            )
            # Region-level analysis
            logger.info(f"[{spatanalysis.pipeline_name}] Region level analysis")
            region_names = spatanalysis.run_region_spatial_analysis(
                exp_id=exp_id, flow_id=flow_id
            )
            # Whole-tissue analysis
            logger.info(
                f"[{spatanalysis.pipeline_name}] Running spatial structure analysis for whole tissue"
            )
            report, cells_adata, spots_adata = spatanalysis.run_spatial_analysis(
                cells_df=spatanalysis.cells_df,
                spots_df=spatanalysis.spots_df,
                level="tissue",
                exp_id=exp_id,
                flow_id=flow_id,
            )
            spatanalysis.save_data(
                report=report,
                cells_adata=cells_adata,
                spots_adata=spots_adata,
                report_name="",
            )
            # Pushing output files to S3
            # spatanalysis.add_files_to_configs(
            #     file_names=[f"{region_name}_report" for region_name in region_names],
            #     bucket="results",
            #     extension="html",
            #     file_type="outputs",
            # )

            data_io.save_run_configs(configs=self.configs)
        else:
            logger.info(f"<{self.input}> Pipeline not requested, SKIP >>>>>")

        self.next(self.join)

    @step
    def join(self, inputs):
        """Join
        :param None(rtype):
        :return None(rtype):
        """
        self.next(self.end)

    @step
    def end(self):
        """End Flow"""
        logger.info("Flow END ---")


if __name__ == "__main__":
    VisiumFlow()
