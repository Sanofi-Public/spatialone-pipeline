"""
inference pipeline for cell deconvolution
model: cell2location
"""
# pylint:disable=no-member

import json

import hydra
import requests
from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.pipelines.dataio_pipeline import DataIO
from src.utils.logger import Logger

load_dotenv()
logger = Logger()


class CARD:
    """
    Python class for cell deconvolution using CARD
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "celldeconv"

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

        # Setup atlas property
        key_cell_atlas_file = list(self.configs["cell_atlas"].keys())[0]
        self.atlas_file_name = (
            self.configs["params"]["atlas_type"]
            + key_cell_atlas_file
            + "."
            + self.configs["cell_atlas"][key_cell_atlas_file]["extension"]
        )
        logger.info(
            f"<{self.pipeline_name}/{self.model_name}/{self.model_version}> cell atlas file name: {self.atlas_file_name}"
        )

    def predict(self, data_io) -> bool:
        """_summary_

        Args:
            data_io (_type_): _description_

        Returns:
            bool: _description_
        """

        endpoint = self.configs["deployment"]["predict_url"]
        logger.info(f"Predict URL: {endpoint}")
        configs_json = json.dumps(self.configs)
        logger.info(f"JSON configs: {configs_json}")

        headers = {"Keep-Alive": "timeout=2200, max=2200"}
        post_params = {
            "configs": configs_json,
            "atlas_filename_ext": self.atlas_file_name,
            "prep_path": data_io.prep_dir,
            "results_path": data_io.results_dir,
            "reference_path": data_io.reference_dir,
            "sample_id": data_io.exp_id,
        }
        logger.info(f"Params for POST predict request: {post_params}")

        try:
            response = requests.post(
                endpoint,
                headers=headers,
                data=post_params,
                proxies={
                    "http": None,
                    "https": None,
                },  # no proxy for cross container communication
            )
            response.raise_for_status()
        except requests.exceptions.ConnectionError as e:
            logger.info("Connection error occured")
            return False
        except requests.exceptions.HTTPError as e:
            logger.info("HTTP error occured")
            return False
        else:
            logger.info(f"Success: {response.text}")
            return True


if __name__ == "__main__":
    # exp_id = os.environ["TEST_exp_id"]
    exp_id = "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma"
    user_id = "sunaal"
    flow_id = "12345"
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
    card_celldeconv = CARD()
    card_celldeconv.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[card_celldeconv.pipeline_name]
    )
    card_celldeconv.predict(data_io=data_io)
