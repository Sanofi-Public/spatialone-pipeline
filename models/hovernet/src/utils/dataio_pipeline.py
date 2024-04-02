"""
Python class to handle push and pull from storage
"""
# pylint:disable=no-member
import json
import os
import shutil
import time
from pathlib import Path

import hydra

# from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.utils.logger import Logger

logger = Logger()


class DataIO:
    """
    DataIO
    """

    def __init__(self, exp_id):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.configs = OmegaConf.to_container(hydra.compose(config_name="config"))[
            "storage"
        ]["storage.data_io"]
        self.exp_id = exp_id
        self.input_path = self.configs["paths"]["input_path"].replace(
            "<exp_id>", str(exp_id)
        )
        self.output_path = self.configs["paths"]["output_path"].replace(
            "<exp_id>", str(exp_id)
        )
        self.cache_dir = self.configs["paths"]["cache_dir"]
        self.reference_path = self.configs["paths"]["reference_path"]
