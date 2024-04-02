"""
Python class to run QC analysis
"""
# pylint:disable=no-member

import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime

import hydra
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.pipelines.dataio_pipeline import DataIO
from src.utils.data_loader import get_matrix_from_h5, load_pos_list
from src.utils.logger import Logger


load_dotenv()

logger = Logger()


class QC:
    """
    Class to compute QC metrics
    """

    def __init__(self):
        """Initializes QC class"""
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "qc"

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

    def load_data(self, prep_dir, results_dir):
        """Load data as needed for QC computation"""
        logger.info("<load_data> load gene matrix")
        self.gene_matrix = get_matrix_from_h5(
            os.path.join(prep_dir, "raw_feature_bc_matrix.h5")
        ).T
        self.spot_ids = self.gene_matrix.columns.tolist()
        self.gene_ids = self.gene_matrix.index.tolist()

        logger.info("<load_data> load Visium spot info")
        with open(
            prep_dir + "scalefactors_json.json", "r", encoding="utf-8"
        ) as spot_geometry_file:
            logger.info("<load_data> load scale factors")
            self.scale_factors = json.loads(spot_geometry_file.read())

        logger.info("<load_data> load spot coordinates")
        pos_filename = "spots_df.csv"
        self.spots_df = pd.read_csv(results_dir + pos_filename)

    def get_qc_metrics(self, np_array=None, gene_ids=None, spot_ids=None, exp_id=None):
        """Generates QC metrics at Spot, Gene and overall level
        In order to run this function mitochondrial genes filter parameters  and metrics labels
        need to be defined at the config file.

        Args:
            np_array (pandas.DataFrame or numpy array): Gene matrix (geneXspot). Defaults to None.
            spot_ids (list, optional): List of barcodes. Defaults to None.
            gene_ids (list, optional): List of genes. Defaults to None.
            exp_id (list, optional): Experiment id to be used as id in overall_qc outputs.

        Returns:
            pandas.DataFrame: Spot, Gene and overall level dataframes with QC metrics.
        """
        if np_array is None or type(np_array) == pd.DataFrame:
            np_array = np.array(self.gene_matrix)
        if spot_ids is None:
            spot_ids = self.spot_ids
        if gene_ids is None:
            gene_ids = self.gene_ids

        logger.info("<get_qc_metrics_by_axis> Generating spot level QC metrics")
        # getting mitohondrial metrics
        mit_n_df, mit_cts_df, mit_mask_df = self.get_gene_artifacts(
            np_array,
            gene_ids,
            start_with=self.configs["params"]["mit_genes_filter"]["start_with"],
            contains=self.configs["params"]["mit_genes_filter"]["list_of_gene_ids"],
            label=self.configs["params"]["mit_genes_filter"]["label"],
        )

        spot_qc = self.get_metrics_by_axis(np_array, spot_ids, 0)
        spot_qc = _get_renamed_metrics(
            spot_qc, self.configs["params"]["metrics_labels"], "spot_qc_labels"
        )
        spot_qc_df = pd.DataFrame(spot_qc).reset_index(drop=True)
        spot_qc_df = pd.concat([spot_qc_df, mit_n_df, mit_cts_df], axis=1)
        spot_qc_df = spot_qc_df.merge(
            self.spots_df[
                [
                    "barcode",
                    "in_tissue",
                    "pxl_row_in_fullres",
                    "pxl_col_in_fullres",
                    "num_contained_cells",
                ]
            ],
            how="inner",
            on="barcode",
        )
        spot_qc_df = spot_qc_df[spot_qc_df.in_tissue.astype(int) == 1]
        spot_qc_df = spot_qc_df.rename(
            columns={"num_contained_cells": "qc_num_contained_cells"}
        )
        self.spot_qc = spot_qc_df.set_index("barcode")

        logger.info("<get_qc_metrics_by_axis> Generating gene level QC metrics")
        gene_qc = self.get_metrics_by_axis(np_array, gene_ids, 1)
        gene_qc = _get_renamed_metrics(
            gene_qc, self.configs["params"]["metrics_labels"], "gene_qc_labels"
        )
        gene_qc_df = pd.DataFrame(gene_qc).reset_index(drop=True)
        gene_qc_df = pd.concat([gene_qc_df, mit_mask_df], axis=1)
        self.gene_qc = gene_qc_df.set_index("gene")

        logger.info("<get_overall_metrics> Generating overall level QC metrics")
        overall_qc = self.get_overall_metrics(np_array, exp_id=exp_id)
        overall_qc = _get_renamed_metrics(
            overall_qc, self.configs["params"]["metrics_labels"], "overall_qc_labels"
        )
        overall_qc = clean_null_keys(overall_qc)
        self.overall_qc = pd.DataFrame(overall_qc, index=[0])

    @staticmethod
    def get_metrics_by_axis(np_array: np.array, ids: list, axis_n: int):
        """This method can be used to generate QC metrics at Spot or Gene level accoring to the specified axis.

        Args:
            np_array (np.array): np_array (pandas.DataFrame or numpy array): Gene matrix (geneXspot).
            ids (list): List of genes or barcodes.
            axis_n (int): axis indicating which axis will be used to generate QC metrics.

        Returns:
            dict: a dictionary with QC metrics at the specified axis.
        """
        if type(np_array) == pd.DataFrame:
            np_array = np.array(np_array)

        # generate a binary array for UMI computation
        np_array_bin = np_array.copy()
        np_array_bin[np_array_bin > 0] = 1
        UMI = np_array_bin.sum(axis=axis_n)
        UMI_cts = np_array.sum(axis=axis_n)

        metrics = {
            "id": ids,
            # raw metrics
            "UMI": UMI,
            "UMI_prop": list(np.round(UMI / np_array_bin.shape[axis_n], 3)),
            "UMI_cts": UMI_cts,
            "UMI_cts_prop": UMI_cts / np_array.sum(),
            # scaled metrics
            "cts_by_UMI": list(np.round(UMI_cts / UMI, 3)),
            "Saturation": list(
                np.round((1 - UMI / UMI_cts) * (UMI / np_array_bin.shape[axis_n]), 3)
            ),
        }
        return metrics

    @staticmethod
    def get_overall_metrics(np_array: np.array, exp_id=None):
        """This method can be used to generate overall QC metrics using both axis from gene matrix table

        Args:
            np_array (np.array): np_array (pandas.DataFrame or numpy array): Gene matrix (geneXspot).
            exp_id (str, optional): String indicating the experiment ID. Defaults to None.

        Returns:
            dictionary: Dictionary with overall QC metrics.
        """
        if type(np_array) == pd.DataFrame:
            np_array = np.array(np_array)

        UMI_cts = np_array.sum()
        UMI = np.sum(np_array > 0)
        metrics = {
            "id": exp_id,
            "UMI": UMI,
            "UMI_cts": UMI_cts,
            "cts_by_UMI": np.round(UMI_cts / UMI, 3),
            "Saturation": np.round(1 - UMI / UMI_cts, 3),
        }
        return metrics

    @staticmethod
    def get_gene_artifacts(
        np_array: np.array, gene_ids: list, start_with=None, contains=None, label=None
    ):
        """This function generates  a report for specific gene artifacts. For instance can be used to filter out mitochondrial genes.
        A starting-with string (start_with) or a list of gene names (contains) can be used to mask and get metrics from the specific gene artifact.
        This method use two additional auxiliary methods to get gene ids and counts of genes at each spots (_get_list_of_reached_genes) and for
        masking (_get_gene_group_mask).

        Args:
            np_array (np.array): numpy array / pd.Dataframe gene matrix containing with barcodes as columns and genes as rows.
            gene_ids (list): list of gene ids.
            start_with (str, optional): String with a prefix asociated to specific genes. Defaults to None.
            contains ( list, optional): List of gene ids used for filtering. Defaults to None.
            label (str, optional): String. Defaults to None.

        Returns:
            pd.dataFrame : It returns dataframes with the number of genes,  the total number of counts for each gene and the used mask.
        """
        if type(np_array) == pd.DataFrame:
            np_array = np.array(np_array)
        np_array = np_array.transpose()
        reached_genes = _get_list_of_reached_genes(np_array)
        gene_mask = _get_gene_group_mask(gene_ids, start_with, contains, label)

        logger.info(
            f" <get_gene_artifacts> getting counts from {label} genes for each spots"
        )
        art_n = [
            np.sum(
                np.array(gene_mask.get("mask"))[
                    reached_genes["reached_genes_idx"].get(spot)
                ]
            )
            for spot in range(len(reached_genes["reached_genes_idx"]))
        ]
        art_cts = [
            np.sum(
                np.array(reached_genes["reached_genes_cts"].get(spot))[
                    np.array(gene_mask.get("mask"))[
                        reached_genes["reached_genes_idx"].get(spot)
                    ]
                ]
            )
            if reached_genes["reached_genes_idx"].get(spot)
            else 0
            for spot in range(len(reached_genes["reached_genes_idx"]))
        ]
        mask = gene_mask.get("mask")
        return (
            pd.DataFrame({gene_mask.get("mask_label") + "_n": art_n}),
            pd.DataFrame({gene_mask.get("mask_label") + "_cts": art_cts}),
            pd.DataFrame({gene_mask.get("mask_label"): mask}),
        )

    def save_data(self, results_dir):
        """Save output file

        Args:
            cache_path (str): Cache path.
        """
        logger.info("<save-data> save (spot_qc)")
        self.spot_qc.to_csv(os.path.join(results_dir, "spot_qc.csv"))
        self.gene_qc.to_csv(os.path.join(results_dir, "gene_qc.csv"))
        self.overall_qc.to_csv(os.path.join(results_dir, "overall_qc.csv"))


def _get_list_of_reached_genes(np_array: np.array):
    """This method gets a dictionary  gene indexes and counts detected at each spots.
    Args:
        np_array (np.array): numpy array / pd.Dataframe gene matrix containing with barcodes as columns and genes as rows.

    Returns:
        dict : dictionaty with gene indexes ("reached_genes_idx")  and their respective counts ("reached_genes_cts") reachec at each spot index.
    """
    if type(np_array) == pd.DataFrame:
        np_array = np.array(np_array)
    logger.info(
        f" <get_list_of_reached_genes> getting list of reached genes for each spots"
    )
    coord_rows = np.where(np_array > 0)[0]
    coord_cols = np.where(np_array > 0)[1]
    list_of_tuples = list(zip(coord_rows, coord_cols))
    genes_reached_by_spot_idx = defaultdict(list)
    genes_reached_by_spot_cts = defaultdict(list)
    for k, v in list_of_tuples:
        genes_reached_by_spot_idx[k].append(v)
        genes_reached_by_spot_cts[k].append(np_array[k, v])
    return {
        "reached_genes_idx": dict(genes_reached_by_spot_idx),
        "reached_genes_cts": dict(genes_reached_by_spot_cts),
    }


def _get_gene_group_mask(
    gene_list: list, start_with: str, contains: list, group_name: str
):
    """Thi method gets a mask for the specific group of genes.
    Args:
        gene_list (list): input list of genes.
        start_with (str): string with a prefix asociated to specific genes.
        contains (list): a list of gene ids used for filtering.
        group_name (str): label for the mask.

    Returns:
        dict: It returns a dictionary with mask name and the list with the masked specific group of genes.
    """

    logger.info(
        f" <get_gene_group_mask> getting a mask for the specific group of genes: {group_name}"
    )
    starting_list = np.zeros(len(gene_list), dtype=bool)
    containing_list = np.zeros(len(gene_list), dtype=bool)
    if start_with is not None:
        starting_list = np.array(
            [
                gene.upper().startswith(start.upper())
                for gene in gene_list
                for start in start_with
            ]
        )
    if contains is not None:
        containing_list = np.array([gene in contains for gene in gene_list])
    return {"mask_label": group_name, "mask": list(starting_list | containing_list)}


def _get_renamed_metrics(metrics: dict, metrics_label: dict, label: str):
    """Auxiliary method used to rename metrics labels from get_metrics_by_axis method outputs.
    Args:
        metrics (dict): dictionary with metrics obtained using  get_metrics_by_axis method.
        metrics_label (dict): a dictionary with metrics labels to rename metrics dict.
        label (str): label indicating level of metrics (gene, spot or overall)

    Returns:
        dict: a dictionary with renamed metrics.
    """
    return {
        metrics_label[label].get("id"): metrics.get("id"),
        metrics_label[label].get("UMI"): metrics.get("UMI"),
        metrics_label[label].get("UMI_prop"): metrics.get("UMI_prop"),
        metrics_label[label].get("UMI_cts"): metrics.get("UMI_cts"),
        metrics_label[label].get("UMI_cts_prop"): metrics.get("UMI_cts_prop"),
        metrics_label[label].get("cts_by_UMI"): metrics.get("cts_by_UMI"),
        metrics_label[label].get("Saturation"): metrics.get("Saturation"),
    }


def clean_null_keys(dt: dict):
    """This auxiliar function cleans the empty keys of input dictionary.

    Args:
        dt (dict): input dictionary.

    Returns:
        dict: cleaned dictionary.
    """
    for k, v in list(dt.items()):
        if k is None:
            del dt[k]
    return dt


if __name__ == "__main__":
    exp_id = "ZZ268953_V11T09-086D"
    user_id = "sunaal"
    flow_id = str(datetime.timestamp(datetime.now())).replace(".", "")
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
    qc = QC()
    qc.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[qc.pipeline_name]
    )

    logger.info(f"[{qc.pipeline_name}] Fetch relevant files")
    # data_io.fetch_input_files(pipeline_configs=qc.configs)
    qc.load_data(prep_dir=data_io.prep_dir, results_dir=data_io.results_dir)
    # print(qc.configs["params"]["mit_genes_filter"]["start_with"] )
    # print(qc.configs["params"]["metrics_labels"]["spot_qc_labels"].get("id") )
    # print(qc.configs["params"]["metrics_labels"]["gene_qc_labels"].get("id") )
    qc.get_qc_metrics(exp_id=exp_id)
    print(qc.spot_qc.iloc[0, :5].transpose())
    print(qc.gene_qc.iloc[0, :5].transpose())
    print(qc.overall_qc.iloc[0, :5].transpose())

    qc.save_data(results_dir=data_io.results_dir)
    # # Push all output files from local cache to s3
    # # data_io.push_output_files(pipeline_configs=qc.configs)
