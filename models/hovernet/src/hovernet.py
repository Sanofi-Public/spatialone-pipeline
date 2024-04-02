"""
Hovernet wrapper
"""
import glob
import logging
import os

import hydra
import torch

from src.misc.utils import log_info
from src.utils.dataio_pipeline import DataIO
from src.utils.logger import Logger

logger = Logger()


class HoverNet:
    """
    Python class for HoverNet implementation
    """

    def __init__(self, configs):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../conf")
        self.pipeline_name = "imgseg"
        self.configs = configs

    def predict(
        self,
        exp_id,
        input_path,
        output_path,
        cache_dir,
        reference_path,
        model_name,
        type_filename,
    ):
        """predict"""
        model_path = reference_path + model_name
        tile_cli = """
        Arguments for processing tiles.

        usage:
            tile (--input_dir=<path>) (--output_dir=<path>) \
                [--draw_dot] [--save_qupath] [--save_raw_map] [--mem_usage=<n>]

        options:
        --input_dir=<path>     Path to input data directory. Assumes the files are not nested within directory.
        --output_dir=<path>    Path to output directory..

        --mem_usage=<n>        Declare how much memory (physical + swap) should be used for caching.
                                By default it will load as many tiles as possible till reaching the
                                declared limit. [default: 0.2]
        --draw_dot             To draw nuclei centroid on overlay. [default: False]
        --save_qupath          To optionally output QuPath v0.2.3 compatible format. [default: False]
        --save_raw_map         To save raw prediction or not. [default: False]
        """

        wsi_cli = """
        Arguments for processing wsi

        usage:
            wsi (--input_dir=<path>) (--output_dir=<path>) [--proc_mag=<n>]\
                [--cache_path=<path>] [--input_mask_dir=<path>] \
                [--ambiguous_size=<n>] [--chunk_shape=<n>] [--tile_shape=<n>] \
                [--save_thumb] [--save_mask]

        options:
            --input_dir=<path>      Path to input data directory. Assumes the files are not nested within directory.
            --output_dir=<path>     Path to output directory.
            --cache_path=<path>     Path for cache. Should be placed on SSD with at least 100GB. [default: cache]
            --mask_dir=<path>       Path to directory containing tissue masks.
                                    Should have the same name as corresponding WSIs. [default: '']

            --proc_mag=<n>          Magnification level (objective power) used for WSI processing. [default: 40]
            --ambiguous_size=<int>  Define ambiguous region along tiling grid to perform re-post processing. [default: 128]
            --chunk_shape=<n>       Shape of chunk for processing. [default: 10000]
            --tile_shape=<n>        Shape of tiles for processing. [default: 2048]
            --save_thumb            To save thumb. [default: False]
            --save_mask             To save mask. [default: False]
        """
        sub_cli_dict = {"tile": tile_cli, "wsi": wsi_cli}
        # hovernet = docopt(__doc__, help=False, options_first=True,
        #                 version='HoVer-Net Pytorch Inference v1.0')
        if self.configs["params"]["wsi"]:
            sub_cmd = "wsi"
        else:
            sub_cmd = "tile"
        # sub_cmd = args.pop('<command>')
        # sub_cmd_args = args.pop('<hovernet>')

        # ! TODO: where to save logging
        logging.basicConfig(
            level=logging.INFO,
            format="|%(asctime)s.%(msecs)03d| [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d|%H:%M:%S",
            handlers=[logging.FileHandler("debug.log"), logging.StreamHandler()],
        )

        if self.configs["params"]["help"]:
            if sub_cmd in sub_cli_dict:
                print(sub_cli_dict[sub_cmd])
            else:
                print(__doc__)
            exit()
        # Commenting out (29.05.2023, since it is redundant above)
        # if self.configs["params"]["help"]:
        #     print(__doc__)
        #     exit()

        # hovernet = docopt(sub_cli_dict[sub_cmd], argv=sub_cmd_args, help=True)

        # hovernet.pop('--version')
        gpu_list = self.configs["params"]["gpu"]
        os.environ["CUDA_VISIBLE_DEVICES"] = gpu_list

        nr_gpus = torch.cuda.device_count()
        log_info("Detect #GPUS: %d" % nr_gpus)

        # hovernet = {k.replace('--', '') : v for k, v in hovernet.items()}
        # hovernet = {k.replace('--', '') : v for k, v in hovernet.items()}

        if not glob.glob(model_path):
            raise Exception(
                "A model path must be supplied as an argument with --model_path."
            )

        nr_types = (
            int(self.configs["params"]["nr_types"])
            if int(self.configs["params"]["nr_types"]) > 0
            else None
        )
        method_args = {
            "method": {
                "model_args": {
                    "nr_types": nr_types,
                    "mode": self.configs["params"]["model_mode"],
                },
                "model_path": model_path,
            },
            "type_info_path": None
            if reference_path + type_filename == ""
            else reference_path + type_filename,
        }

        # ***
        run_args = {
            "batch_size": int(self.configs["params"]["batch_size"]) * nr_gpus,
            "nr_inference_workers": int(self.configs["params"]["nr_inference_workers"]),
            "nr_post_proc_workers": int(self.configs["params"]["nr_post_proc_workers"]),
        }

        if self.configs["params"]["model_mode"] == "fast":
            run_args["patch_input_shape"] = 256
            run_args["patch_output_shape"] = 164
        else:
            run_args["patch_input_shape"] = 270
            run_args["patch_output_shape"] = 80

        if self.configs["params"]["tile"]:
            run_args.update(
                {
                    "exp_id": exp_id,
                    "input_dir": input_path,
                    "output_dir": output_path,
                    "mem_usage": float(self.configs["tile"]["mem_usage"]),
                    "draw_dot": self.configs["tile"]["draw_dot"],
                    "save_qupath": self.configs["tile"]["save_qupath"],
                    "save_raw_map": self.configs["tile"]["save_raw_map"],
                }
            )

        if self.configs["params"]["wsi"]:
            run_args.update(
                {
                    "exp_id": exp_id,
                    "input_dir": input_path,
                    "output_dir": output_path,
                    "input_mask_dir": input_path,
                    "cache_path": cache_dir,
                    "proc_mag": int(self.configs["wsi"]["proc_mag"]),
                    "ambiguous_size": int(self.configs["wsi"]["ambiguous_size"]),
                    "chunk_shape": int(self.configs["wsi"]["chunk_shape"]),
                    "tile_shape": int(self.configs["wsi"]["tile_shape"]),
                    "save_thumb": self.configs["wsi"]["save_thumb"],
                    "save_mask": self.configs["wsi"]["save_mask"],
                }
            )
        # ***

        if self.configs["params"]["tile"]:
            from src.infer.tile import InferManager

            infer = InferManager(**method_args)
            infer.process_file_list(run_args)
        else:
            from src.infer.wsi import InferManager

            infer = InferManager(**method_args)
            infer.process_wsi_list(run_args)


if __name__ == "__main__":
    hovernet = HoverNet(configs)
    data_io = DataIO("CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma")
    print(data_io.output_path)
    print(data_io.input_path)

    model_name = configs["model_filename"]
    type_filename = configs["type_filename"]

    if not model_name:
        model_name = "hovernet_original_consep_type_tf2pytorch.tar"
    if not type_filename:
        type_filename = "type_info.json"

    hovernet.predict(
        "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma",
        data_io.input_path,
        data_io.output_path,
        data_io.cache_dir,
        data_io.reference_path,
        model_name,
        type_filename,
    )
