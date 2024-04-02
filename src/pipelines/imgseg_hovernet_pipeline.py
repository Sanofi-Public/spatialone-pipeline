import hydra
import requests
from omegaconf import OmegaConf

from src.pipelines.dataio_pipeline import DataIO
from src.utils.logger import Logger

logger = Logger()


class HoverNet:
    """
    Python class for HoverNet implementation
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
        print(self.configs)

    def predict(self, data_io):
        """_summary_

        Args:
            data_io (_type_): _description_
        """
        endpoint = self.configs["deployment"]["predict_url"]
        logger.info(f"Predict URL: {endpoint}")

        data = {
            "configs": self.configs,
            "exp_id": data_io.exp_id,
            "input_path": data_io.prep_dir,
            "output_path": data_io.results_dir,
            "reference_path": data_io.reference_dir,
            "model_filename": self.configs["params"]["model_name"],
            "type_filename": self.configs["params"]["type_filename"],
        }
        logger.info(f"The post parameter payload is: \n {data}")

        try:
            response = requests.post(
                endpoint, json=data, proxies={"http": None, "https": None}
            )
            response.raise_for_status()
        except requests.exceptions.ConnectionError as e:
            logger.info("Connection error occured")
        except requests.exceptions.HTTPError as e:
            logger.info("HTTP error occured")
        else:
            logger.info("Success")


if __name__ == "__main__":
    exp_id = "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma"
    user_id = "sunaal"
    user_group = "demo_user"
    storage_env = "dev"
    flow_id = "12345"
    hydra.core.global_hydra.GlobalHydra.instance().clear()
    hydra.initialize(config_path="../../data/conf")
    config_flow = OmegaConf.to_container(
        hydra.compose(config_name="visium_config_flow")
    )
    logger.info(
        f"exp_id: {exp_id}, user_id: {user_id}, user_group: {user_group}, storage_env: {storage_env}"
    )
    data_io = DataIO(
        exp_id=exp_id,
        user_id=user_id,
        user_group=user_group,
        storage_env=storage_env,
        flow_id=flow_id,
    )
    hovernet_imgseg = HoverNet()
    hovernet_imgseg.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[hovernet_imgseg.pipeline_name]
    )
    logger.info(f"[{hovernet_imgseg.pipeline_name}] Start pipeline")
    logger.info(f"Running HoverNet Inference...")
    hovernet_imgseg.predict(data_io)
