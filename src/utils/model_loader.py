"""
Model loader for segmentation
"""

import numpy as np
import torch
from cellpose import models

from src.utils.logger import Logger

logger = Logger()


# PUT THIS FUNC IN TRAIN FOLDER
def cellpose_model(model_type: str, gpu_enabled=True):
    """Predict on patched array with model of choice
    Arguments:
        model_type: define type of model from nuclei, cyto, cyto2
    Returns:
        model
    """
    if gpu_enabled and torch.cuda.is_available():
        device = torch.device("cuda:0")
        print("Detected GPU device, using GPU")
    else:
        device = torch.device("cpu")
        logger.info("No GPU device found, using CPU")
        gpu_enabled = False

    model = models.Cellpose(gpu=gpu_enabled, model_type=model_type, device=device)

    return model


# LOAD IN ELTON'S MODEL WEIGHTS
def predict(
    model,
    patches: np.ndarray,
    batch_size: int,
    diameter,
    channels: list,
    flow_threshold: float,
):
    """Predict on patched array with model of choice
    Arguments:
        model -- Segmentation model
        patches {np.ndarray} -- array of patched original image
        batch_size {int} -- batch_size for processing segmentation through model on GPU
        diameter {int or str} -- size of cell diameter, if set to 0 (same as setting to None), then diameter is automatically estimated if size model is loaded
        flow_threshold {float} - the maximum allowed error of the flows for each mask
        channels {list} --  list of channels(ex: [0,0,0])
    Returns:
        masks {np.ndarray} -- labelled image
        flows {list} -- XY flows at each pixel
        styles {list} -- style vector summarizing each image, also used to estimate size of objects in image
        diams {list} -- list of diameters
    """
    masks, flows, styles, diams = model.eval(
        patches,
        batch_size=batch_size,
        diameter=diameter,
        channels=channels,
        normalize=True,
        flow_threshold=flow_threshold,
    )
    return masks, flows, styles, diams
