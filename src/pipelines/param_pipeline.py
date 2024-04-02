"""
Pipeline for managing parameters/configurations from experiments
"""
import csv
from datetime import datetime

# pylint:disable=no-member
import hydra
import pandas as pd
from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.pipelines.dataio_pipeline import DataIO
from src.utils.logger import Logger

load_dotenv()
logger = Logger()


class ParamPipeline:
    """
    Python class for image segmentation
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "param"

    def load_configs(self):
        """load_configs"""
        self.model_name = "default"
        self.config_key = self.pipeline_name + "." + self.model_name
        self.configs = OmegaConf.to_container(hydra.compose(config_name="config_tree"))[
            self.pipeline_name
        ][self.config_key]

    def load_data(self, prep_dir):
        """load_data"""
        # Create params dict and store the parameter_name (column=0) and parameter_value (column=1)
        self.params = {}
        with open(prep_dir + "parameters.csv", "r") as file:
            csv_reader = csv.reader(file)
            next(csv_reader)  # skip the first row
            for row in csv_reader:
                self.params[row[0]] = row[1]
        logger.info(f"experimental.params: {self.params}")


if __name__ == "__main__":
    exp_id = "test"
    user_id = "sunaal"
    flow_id = str(datetime.timestamp(datetime.now())).replace(".", "")
    hydra.initialize(config_path="../../data/conf")
    # config_flow = OmegaConf.to_container(hydra.compose(config_name="visium_config_flow"))
    logger.info(f"exp_id: {exp_id}, user_id: {user_id}")
    data_io = DataIO(
        exp_id=exp_id,
        user_id=user_id,
        flow_id=flow_id,
    )

    param_pipeline = ParamPipeline()
    param_pipeline.load_configs()

    logger.info(f"[{param_pipeline.pipeline_name}] Fetch relevant files")
    # data_io.fetch_input_files(
    #     pipeline_configs=param_pipeline.configs,
    # )

    # logger.info(f"[{param_pipeline.pipeline_name}] Load params")
    param_pipeline.load_data(prep_dir=data_io.prep_dir)
