"""
Saves data to various types
"""
import os

import numpy as np
import skimage as skio


def save_data(file_type, file_path, item, file_name):
    """Save data to various data types
    Arguments:
        file_type {string} : string which denotes the saved file's type
        file_path {string} : path to where the file will be saved
        item : data that will be saved. Can be a numpy array, a pandas dataframe.
        file_name {string} : name of the file
    """
    if file_type == "numpy":
        np.save(os.path.join(file_path, file_name), item)
    if file_type == "csv":
        item.to_csv(os.path.join(file_path, file_name))
    if file_type == "png":
        skio.imsave(os.path.join(file_path, file_name), item)
