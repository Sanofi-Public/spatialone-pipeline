"""
Image segmentation Metaflow pipeline
"""
import json
import os
from datetime import datetime

# pylint:disable=no-member
import hydra
import pytz
import scanpy as sc
from omegaconf import OmegaConf

from src.pipelines.assign_pipeline import SpotCluster
from src.pipelines.cell2spot_pipeline import Cell2Spot
from src.pipelines.celldeconv_cell2location_pipeline import Cell2Location
from src.pipelines.clustering_pipeline import MorphCluster
from src.pipelines.dataio_pipeline import DataIO
from src.pipelines.datamerge_pipeline import DataMerge
from src.pipelines.imgseg_cellpose_pipeline import ImageSeg

from src.pipelines.param_pipeline import ParamPipeline
from src.pipelines.qc_pipeline import QC
from src.pipelines.spatialanalysis_pipeline import SpatialAnalysis
from src.utils.env_variables import Env
from src.utils.logger import Logger

# Note: Always load the environment variables before importing MetaFlow
# since some of these variable are required by MetaFlow for proper
# initialization
e = Env()
from metaflow import (  # pylint:disable=no-name-in-module, wrong-import-position, wrong-import-order, line-too-long; kubernetes,; schedule,
    FlowSpec,
    Parameter,
    Run,
    current,
    project,
    step,
)

logger = Logger()


class SpatialOneFlow(FlowSpec):
    """Class defining common functions for visium and phenocycler flows"""

    # Specify Metaflow parameters
    config_flow_filename = Parameter(
        "config_flow_filename",
        help="File name for configuation and parameters",
        default="visium_config_flow",
    )
    env = Parameter(
        "env",
        help="environment requested for spatial one run",
        default="dev",
    )

    def spatial_one_flow_start(self):
        """Start
        :param None(rtype):
        :return None(rtype):
        """

        # log start_datetime
        self.start_datetime = datetime.now(pytz.timezone("UTC"))
        logger.info(f"Flow start datetime [UTC]: {self.start_datetime}")
        logger.info(f"Requested environment: {self.env}")
        # logger.info(f"project name:{current.project_name}")
        # logger.info(f"project branch:{current.branch_name}")
        # logger.info(f"Is this a production run?: {current.is_production}")

        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../data/conf")

        self.configs = OmegaConf.to_container(
            hydra.compose(config_name=self.config_flow_filename)
        )

        if self.configs is not None:
            # load variables from configs
            self.exp_ids = self.configs["experiment.ids"]
            self.user_id = self.configs["user.id"]
            self.pipelines_enabled = self.configs["pipelines.enabled"]
            self.version = os.environ.get("DOCKER_IMAGE_VERSION")
            # update configs to have version
            self.configs["version"] = self.version
            # TODO: determine if wee need to change to only single exp.id and remove "foreach" in flow steps
            if type(self.exp_ids) != list:
                self.exp_ids = [self.exp_ids]

            # Spatial specific parameters
            logger.info(f"experiment_ids: {self.exp_ids}")
            logger.info(f"pipeline list: {self.pipelines_enabled}")
            logger.info(f"version : {self.version}")

            # add tags to current run
            # Note: metaflow's run_id is considered to be the "flow_id" for spatial_one
            self.flow_id = current.run_id
            current_flow = Run(f"{current.flow_name}/{self.flow_id}")
            logger.info(f"current_flow : {current_flow}")
            logger.info(f"flow_id : {self.flow_id}")

            current_flow.add_tags(
                [
                    f"user:{self.user_id}",
                    f"expid:{self.exp_ids}",
                    f"v.{self.version}",
                ]
            )

    def spatial_one_flow_get_params(self):
        """Upload config_flow to analysis
        :param None(rtype):
        :return None(rtype):
        """
        if self.input is not None:
            # upload the config_flow.yaml or {self.configs} to s3
            exp_id = self.input
            data_io = DataIO(
                exp_id=exp_id,
                user_id=self.user_id,
                flow_id=self.flow_id,
            )
            param_pipeline = ParamPipeline()
            param_pipeline.load_configs()
            logger.info(f"[{param_pipeline.pipeline_name}] Fetch relevant files")
            # data_io.fetch_input_files(
            #     pipeline_configs=param_pipeline.configs,
            # )
            logger.info(f"[{param_pipeline.pipeline_name}] Load params")
            param_pipeline.load_data(prep_dir=data_io.prep_dir)

            # Save the parameters from prep/"parameters.csv" in self.configs
            self.experimental_params = param_pipeline.params
            self.configs["experimental.params"] = param_pipeline.params

    def spatial_one_log_end(self):
        """Update flow info if triggered from rabbitmq"""
        self.end_datetime = datetime.now(pytz.timezone("UTC"))
        logger.info(f"Flow({self.flow_id}) end datetime: {self.end_datetime}")
        logger.info(
            f"Flow({self.flow_id}) duration: {self.end_datetime - self.start_datetime}"
        )


if __name__ == "__main__":
    SpatialOneFlow()
