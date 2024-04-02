""" metrics.py
"""
import numpy as np


def dice_coeff(img1, img2):
    """Compute Dice coefficient between two binary mask segmentations
    Arguments:
        img1 {np.array} : masked segmentation array (binary)
        img2 {np.array} : masked segmentation array (binary)
    Returns:
        float -- dice score
    """
    intersection = np.logical_and(img1, img2)
    union = np.logical_or(img1, img2)
    dice = (2 * np.sum(intersection)) / (np.sum(union) + np.sum(intersection))
    return dice
