"""
Feature extraction
"""

import numpy as np
import pandas as pd
from skimage.measure import label, regionprops_table


def extract_image_features(mask: np.ndarray, wsi: np.ndarray) -> pd.DataFrame:
    """Extract features from segmented mask with labeled cells.
    Arguments:
        mask {np.ndarray} -- Segmented Mask
    Returns:
        pd.DataFrame -- Dataframe with collected features
    """
    labeled_mask = label(mask, return_num=False)
    features = regionprops_table(
        labeled_mask,
        wsi,
        properties=[
            "label",
            "area",
            "area_bbox",
            "area_convex",
            "area_filled",
            "axis_major_length",
            "axis_minor_length",
            "eccentricity",
            "euler_number",
            "extent",
            "inertia_tensor",
            "inertia_tensor_eigvals",
            "moments_hu",
            "orientation",
            "perimeter",
            "perimeter_crofton",
            "solidity",
            "intensity_max",
            "intensity_mean",
            "intensity_min",
        ],
    )
    master_table = pd.DataFrame(features)
    return master_table
