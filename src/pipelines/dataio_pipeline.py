"""
Created By  : Mena Kamel
Created Date: 2022/11/23
Python class to handle push and pull from storage
"""
import datetime

# pylint:disable=no-member
import json
import os
import shutil
import time

import hydra
from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.utils.logger import Logger

logger = Logger()
load_dotenv()


class DataIO:
    """
    Python class to compile data from all analysis streams and build a .tmap file
    """

    def __init__(self, exp_id, user_id, flow_id, pipeline_name=""):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.configs = OmegaConf.to_container(hydra.compose(config_name="config_tree"))[
            "dataio"
        ]

        self.exp_id = exp_id
        self.user_id = user_id
        self.flow_id = flow_id
        # self.time_stamp = time.strftime("%m/%d/%Y, %H:%M:%S", time.localtime())
        self.time_stamp = "0000000_00:00:00"
        dataio_paths_keyname = "dataio.paths"  # defining literal once
        for key, val in self.configs[dataio_paths_keyname].items():
            val = (
                val.replace("<exp_id>", self.exp_id)
                .replace("<user_id>", self.user_id)
                .replace("<flow_id>", flow_id)
            )
            self.configs[dataio_paths_keyname][key] = val

        # for local
        self.data_dir = self.configs[dataio_paths_keyname]["data_dir"]
        self.prep_dir = (
            self.data_dir + self.configs[dataio_paths_keyname]["local_prep_path"]
        )
        self.results_dir = (
            self.data_dir + self.configs[dataio_paths_keyname]["local_results_path"]
        )
        # Create results dir if doesn't exist
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

        self.reference_dir = (
            self.data_dir + self.configs[dataio_paths_keyname]["local_reference_path"]
        )

        self.config_flow_dir = self.configs[dataio_paths_keyname]["config_flow_dir"]
        # self.unittest_path = self.configs[dataio_paths_keyname]["unittest_path"]

    def validate_param(self, param_name, param_value):
        """_summary_

        Args:
            param_name (string): name of parameter
            param_value (string): parameter value

        Raises:
            ValueError: _description_
        """
        if param_value not in self.configs["dataio.params"][param_name]:
            print(param_value)
            raise ValueError(
                f"Parameter value of {param_value} is not a supported parameter for {param_name}. "
                f"Valid parameters are {self.configs['dataio.params']}"
            )
        else:
            # print(f"validated {param_name}: {param_value}")
            return param_value

    # def fetch_unittest_files(self, pipeline_configs):
    #     """_fetch_unittest_files_

    #     Args:
    #         pipeline_configs (dict): dictionary of pipeline configurations
    #     """
    #     unittest_file_names = pipeline_configs["unittest"].keys()
    #     for file in unittest_file_names:
    #         # Grab metadata of file from pipeline configs
    #         # Create input file path to load from cloud provider
    #         file_name = f"{file}.{file_ext}"
    #         target_file_path = self.cache_dir + file_name
    #         source_file_path = self.unittest_path + file_name
    #         logger.info(f"source_path: {source_file_path}")
    #         logger.info(f"target_path: {target_file_path}")

    def move_file(self, source_file_path, target_file_path):  # pragma: no cover
        """moving filesystem

        Args:
            source_file_path (string): source file path
            target_file_path (string): target file path
        """
        shutil.copy(source_file_path, target_file_path)

    def save_run_configs(self, configs):
        """update configs"""

        logger.info("<save_run_configs> save self.configs to cache dir")
        # Write JSON data to the file
        file_path = self.results_dir + "run_configs.json"
        logger.info(file_path)
        with open(file_path, "w") as json_file:
            json.dump(configs, json_file, indent=4)

    def poll_for_update(self, file_path, polling_interval=60, poll_timeout_min=300):

        # TODO check if we can use requests without timeout
        if os.path.exists(file_path):
            initial_modification_time = os.path.getmtime(file_path)
        else:
            initial_modification_time = datetime.datetime.now()
        start_time = time.time()
        while True:
            logger.info(f"({datetime.datetime.now()} polling for {file_path}")

            current_modification_time = os.path.getmtime(file_path)

            if current_modification_time > initial_modification_time:
                print(f"Finished polling for {file_path}")
                return

            if time.time() - start_time > poll_timeout_min:
                print(f"Timeout reached. File '{file_path}' was not found.")
                return

            time.sleep(polling_interval)


if __name__ == "__main__":
    exp_id = "ZZ268953_V11T09-086D"
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

    # data_io.save_run_configs(config_flow)

    # data_io.poll_for_update(
    #     bucket_name="results",
    #     source_base_path=data_io.analysis_path,
    #     parent_folder="",
    #     file_name_ext="card_results.csv"
    # )
