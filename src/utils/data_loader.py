"""
Loads whole slide images, segmentation masks, and Visium spot coordinates
"""
import csv

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import tables
import tifffile as tifi


def load_images(file_path):
    """Load image from given file path
    Arguments:
        file_path {string} : path to the file that we are trying to load
    Returns:
        np.array -- loaded image as numpy array
    """
    img_arr = tifi.imread(file_path)
    return img_arr


def load_segmentation(file_path):
    """Load segmented image from given file path
    Arguments:
        file_path {string} : path to the file that we are trying to load
    Returns:
        np.array -- loaded image as numpy array
    """
    img_arr = np.load(file_path)
    return img_arr


def load_pos_list(tissue_pos_path, header=0):
    """Load tissue position list from given file path
    Arguments:
        file_path {string} : path to the file that we are trying to load
    Returns:
        pd.DataFrame() -- loaded csv containing the positions of tissue in the slide image
    """
    tissue_pos_csv = pd.read_csv(tissue_pos_path, header=header)

    # commenting this hovernet fix out due to issues with reading csv files in cluster, trying header=0 above instead
    # with open(tissue_pos_path, 'r') as csvfile:
    #     dialect = csv.Sniffer().sniff(csvfile.read(1024))
    #     print(dialect)
    #     csvfile.seek(0)
    #     reader = csv.reader(csvfile, dialect)

    # if reader:
    #     tissue_pos_csv = pd.read_csv(tissue_pos_path)
    # else:
    #     tissue_pos_csv = pd.read_csv(tissue_pos_path, header= None)
    return tissue_pos_csv


def get_matrix_from_h5(filename):
    """Reads the gene-by-bead matrix given by the path filename

    Args:
        filename (str): Path to gene matrix (.h5).

    Returns:
        pandas.DataFrame: gene matrix
    """
    with tables.open_file(filename, "r") as gene_file:
        mat_group = gene_file.get_node(gene_file.root, "matrix")
        barcodes = gene_file.get_node(mat_group, "barcodes").read()
        data = getattr(mat_group, "data").read()
        indices = getattr(mat_group, "indices").read()
        indptr = getattr(mat_group, "indptr").read()
        shape = getattr(mat_group, "shape").read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        feature_ref = {}
        feature_group = gene_file.get_node(mat_group, "features")
        feature_ids = getattr(feature_group, "id").read()
        feature_names = getattr(feature_group, "name").read()
        feature_types = getattr(feature_group, "feature_type").read()
        feature_ref["id"] = feature_ids
        feature_ref["name"] = feature_names
        feature_ref["feature_type"] = feature_types
        tag_keys = getattr(feature_group, "_all_tag_keys").read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()
        gene_df = pd.DataFrame(
            matrix.toarray(), columns=np.char.decode(barcodes, "utf-8")
        ).T
        gene_df.columns = np.char.decode(feature_ref["name"], "utf-8")
        gene_df = gene_df[
            gene_df.columns.drop(list(gene_df.filter(regex="DEPRECATED")))
        ]
        gene_matrix = gene_df.loc[:, gene_df.sum(axis=0) > 0]
        return gene_matrix
