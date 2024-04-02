"""
Python class to compile data from all analysis pipelines and build a .tmap file
"""
# pylint:disable=no-member

import os
import random
import warnings
from datetime import datetime

import anndata
import cv2
import hydra
import numpy as np
import pandas as pd
import rasterio
import scanpy
import tifffile
from dotenv import load_dotenv
from omegaconf import OmegaConf
from PIL import Image
from shapely import wkt

from src.pipelines.dataio_pipeline import DataIO
from src.utils.data_loader import get_matrix_from_h5
from src.utils.env_variables import Env
from src.utils.logger import Logger
from src.utils.tmap_generator import tmapGenerator

warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)
load_dotenv()
logger = Logger()


class DataMerge:
    """
    Class to compile data from all pipelines and build a .tmap file for TissUUmaps
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "datamerge"

    def load_model_configs_from_flow(
        self, pipeline_configs_from_flow, tissue_type, species_type
    ):
        """loading configs from config_flow to pipeline.configs
        Args:
            pipeline_configs_from_flow (dict): pipeline configs from config_flow
            tissue_type (string): e.g. "lung"
            species_type (string): e.g. "human"
        """
        if "version" in pipeline_configs_from_flow["model"].keys():
            self.model_version = pipeline_configs_from_flow["model"]["version"]

        if "run_data_merge_only" in pipeline_configs_from_flow.keys():
            self.run_data_merge_only = pipeline_configs_from_flow["run_data_merge_only"]

        # Access paths from default.yaml confg
        self.config_key = f"{self.pipeline_name}.default"
        self.configs = OmegaConf.to_container(hydra.compose(config_name="config_tree"))[
            self.pipeline_name
        ][self.config_key]

        if (
            "name" in pipeline_configs_from_flow["model"]
            and pipeline_configs_from_flow["model"]["name"]
        ):
            self.merge_setting = pipeline_configs_from_flow["model"]["name"]
        else:
            logger.warning(
                "Model name for data merge has not been properly defined. The generic configuration will be used."
            )
            self.merge_setting = "generic_config"

        self.config_key = f"{self.pipeline_name}.{self.merge_setting}"
        param_configs = OmegaConf.to_container(
            hydra.compose(config_name="config_tree")
        )[self.pipeline_name][self.config_key]

        # Check if user wants to override params
        if "params" in pipeline_configs_from_flow["model"]:
            logger.info("loading params from visium_config_flow")
            self.configs["params"] = pipeline_configs_from_flow["model"]["params"]
        else:
            logger.info(f"loading params from {self.pipeline_name}.default")
            self.configs["params"].update(param_configs["params"])

    def load_data(
        self, prep_dir, results_dir, config_flow, pipelines_enabled, unittest=False
    ):
        """Loads the data from all the upstream pipelines

        Args:
            prep_dir (_type_): _description_
            results_dir (_type_): _description_
            config_flow (_type_): _description_
            unittest (bool, optional): _description_. Defaults to False.
        """

        logger.info(self.configs["params"])
        self.target_genes = self.configs["params"]["target_genes"]
        self.spot_index = self.configs["params"]["spot_index"]
        self.cell_index = self.configs["params"]["cell_index"]

        # Loading all upstream pipeline outputs
        self.load_cell2spot_results(results_dir, pipelines_enabled)
        self.load_deconvolution_results(
            results_dir,
            config_flow,
            pipelines_enabled,
        )
        self.load_morph_clustering_results(results_dir, pipelines_enabled)
        self.load_qc_results(results_dir, pipelines_enabled)
        self.load_tsne_results(prep_dir)
        self.load_umap_results(prep_dir)
        self.load_graphical_clustering_results(prep_dir)
        self.load_kmeans_clustering_results(prep_dir)
        self.load_tissue_positions_list(prep_dir, unittest)
        self.load_gene_matrix(prep_dir)

        # Get WSI size
        if not unittest and ("wsi_shape" not in self.__dict__):
            logger.info("<load_data> Get WSI image size")
            with rasterio.open(prep_dir + "wsi.tif") as wsi:
                self.wsi_shape = wsi.shape

    def load_cell2spot_results(self, results_dir, pipelines_enabled):
        """
        Function to load the cell2spot results

        Args:
            results_dir (_type_): _description_
            pipelines_enabled (_type_): _description_
        """
        logger.info(
            "<load_data> Load cell2spot pipeline outputs: cell-level and "
            "spot-level dataframes"
        )
        spatial_cols = ["cell_x", "cell_y"]
        if pipelines_enabled.get("cell2spot", False) or self.run_data_merge_only:
            unnamed_cols = ["Unnamed: 0"]
            cells_df = pd.read_csv(results_dir + "cells_df.csv")
            cells_df = cells_df.drop(unnamed_cols, axis=1)
            cells_df["barcode"] = cells_df["barcode"].fillna("no barcode")
            self.cells_df = cells_df.sort_values(self.cell_index).set_index(
                self.cell_index
            )
            self.cell_spatial_df = self.cells_df[spatial_cols].astype(int)
            self.spots_df = pd.read_csv(results_dir + "spots_df.csv")
            self.spots_df = self.spots_df.drop(unnamed_cols, axis=1)
            self.spots_df = self.spots_df.sort_values(self.spot_index).set_index(
                self.spot_index
            )
        else:
            self.cells_df = None
            self.spots_df = None
            self.cell_spatial_df = None

    def load_deconvolution_results(self, results_dir, config_flow, pipelines_enabled):
        """Load cell deconvolution results

        Args:
            results_dir (_type_): _description_
            config_flow (_type_): _description_
            pipelines_enabled (_type_): _description_
            unittest (_type_): _description_
        """
        logger.info("<load_data> Load cell deconvolution results")
        if pipelines_enabled.get("celldeconv", False) or self.run_data_merge_only:
            deconv_method = config_flow["celldeconv"]["model"]["name"]
            deconv_results_df = pd.read_csv(
                results_dir + f"{deconv_method}_results.csv"
            )
            deconv_results_df = deconv_results_df.sort_values(
                self.spot_index
            ).set_index(self.spot_index)
            deconv_results_df.columns = [
                col.lower() for col in deconv_results_df.columns
            ]
            deconv_results_df["deconv_cell_counts"] = deconv_results_df.sum(axis=1)
            deconv_results_df = deconv_results_df.div(
                deconv_results_df["deconv_cell_counts"], axis=0
            )
            deconv_results_df = deconv_results_df.fillna(0.0)
            self.deconv_results_df = deconv_results_df
        else:
            self.deconv_results_df = None

    def load_morph_clustering_results(self, results_dir, pipelines_enabled):
        """Load clustering pipeline outputs

        Args:
            results_dir (_type_): _description_
            pipelines_enabled (_type_): _description_
        """
        logger.info("<load_data> Load morphological & spot-level clustering results")
        if pipelines_enabled.get("cluster", False) or self.run_data_merge_only:
            morph_clustering = pd.read_csv(results_dir + "morphological_clusters.csv")
            # Renaming the columns
            morph_clustering = morph_clustering.rename(
                columns={
                    "label": self.cell_index,
                    "clusters": "morphological_cluster",
                }
            ).drop(["Unnamed: 0"], axis=1)
            # Setting the index
            self.morph_clustering = (
                morph_clustering.sort_values(self.cell_index)
                .set_index(self.cell_index)
                .astype(str)
            )
            self.morph_clustering = self.morph_clustering[["morphological_cluster"]]
            self.spot_clustering = (
                pd.read_csv(results_dir + "spot_clusters.csv")
                .sort_values(self.cell_index)
                .set_index(self.cell_index)
            )
        else:
            self.morph_clustering = None
            self.spot_clustering = None

    def load_qc_results(self, results_dir, pipelines_enabled):
        """Load QC pipeline outputs

        Args:
            results_dir (_type_): _description_
            pipelines_enabled (_type_): _description_
        """
        logger.info("<load_data> Load QC results")
        if pipelines_enabled.get("qc", False) or self.run_data_merge_only:
            spot_qc = pd.read_csv(results_dir + "spot_qc.csv")
            qc_columns = [col for col in spot_qc.columns if "qc_" in col]
            spot_qc = spot_qc[[self.spot_index] + qc_columns]
            self.spot_qc = spot_qc.sort_values(self.spot_index).set_index(
                self.spot_index
            )
            self.gene_qc = (
                pd.read_csv(results_dir + "gene_qc.csv")
                .sort_values("gene")
                .drop_duplicates(subset="gene")
                .set_index("gene")
            )
            self.overall_qc = pd.read_csv(results_dir + "overall_qc.csv")
        else:
            self.spot_qc = None
            self.gene_qc = None
            self.overall_qc = None

    def load_tsne_results(self, prep_dir):
        """Load tSNE projection (SpaceRanger)

        Args:
            prep_dir (_type_): _description_
        """
        logger.info("<load_data> Load tSNE Projection")
        try:
            self.tsne_df = (
                pd.read_csv(prep_dir + "tsne_projection.csv")
                .rename(columns={"Barcode": self.spot_index})
                .sort_values(self.spot_index)
                .set_index(self.spot_index)
            )
        except FileNotFoundError:
            logger.warning("File not found: " + prep_dir + "tsne_projection.csv")
            self.tsne_df = None

    def load_umap_results(self, prep_dir):
        """Load UMAP projection (SpaceRanger)

        Args:
            prep_dir (_type_): _description_
        """
        logger.info("<load_data> Load UMAP Projection")
        try:
            self.umap_df = (
                pd.read_csv(prep_dir + "umap_projection.csv")
                .rename(columns={"Barcode": self.spot_index})
                .sort_values(self.spot_index)
                .set_index(self.spot_index)
            )
        except FileNotFoundError:
            logger.warning("File not found: " + prep_dir + "umap_projection.csv")
            self.umap_df = None

    def load_graphical_clustering_results(self, prep_dir):
        """Load graphical clustering results (SpaceRanger)


        Args:
            prep_dir (_type_): _description_
        """
        logger.info("<load_data> Load graphical clustering results")
        graph_clust_path = prep_dir + "gene_graph_clusters.csv"
        try:
            graph_clusters = pd.read_csv(graph_clust_path)
            graph_clusters.columns = [
                self.spot_index,
                "graphical_clustering",
            ]
            self.graph_clusters = graph_clusters.sort_values(self.spot_index).set_index(
                self.spot_index
            )

        except FileNotFoundError:
            logger.warning("File not found: " + graph_clust_path)
            self.graph_clusters = None

    def load_tissue_positions_list(self, prep_dir, unittest):
        """Load tissue position list


        Args:
            prep_dir (_type_): _description_
            unittest (_type_): _description_
        """
        logger.info("<load_data> Load tissue position table")
        vs_headings = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]
        pos_filename = "tissue_positions_list.csv"

        tissue_pos_csv = pd.read_csv(
            prep_dir + pos_filename,
            header=0,
            names=vs_headings,
        )
        self.tissue_pos_csv = (
            tissue_pos_csv[tissue_pos_csv["in_tissue"].astype(str) == "1"]
            .sort_values(self.spot_index)
            .set_index(self.spot_index)
        )
        self.in_tissue_barcodes = self.tissue_pos_csv.index.tolist()

    def load_gene_matrix(self, prep_dir):
        """Load spot-by-gene matrix (SpaceRanger)

        Args:
            prep_dir (_type_): _description_
        """
        logger.info("<load_data> Load spot-by-gene matrix")
        gene_matrix_path = prep_dir + "raw_feature_bc_matrix.h5"
        self.gene_matrix_adata = scanpy.read_10x_h5(gene_matrix_path)
        gene_matrix_df = get_matrix_from_h5(gene_matrix_path)
        self.gene_matrix_df = gene_matrix_df.sort_index()
        self.gene_matrix_df.index.name = self.spot_index

    def identify_top_n_genes(self, n=1000, mode="expression"):
        """
        Identifies and returns the top N genes based on either total expression ('expression' mode)
        or variability ('variability' mode) from the gene matrix.

        Args:
        - n (int): Number of top genes to return. Default is 1000.
        - mode (str): Selection mode - 'expression' for total expression or 'variability' for variance. Default is 'expression'.

        Returns:
        - List[str]: Names of the top N genes.
        """
        gene_df = self.gene_matrix_df

        if mode.lower() == "expression":
            col_sums = gene_df.sum()
        elif mode.lower() == "variability":
            col_sums = gene_df.var()
        else:
            logger.warning(
                f"Unrecognized mode '{mode}' when filtering genes. 'expression' mode will be used"
            )
            col_sums = gene_df.sum()

        return col_sums.nlargest(n).index.tolist()

    def filter_gene_data(self):
        """
        Formats the gene-by-bead matrix into a regular dataframe only
        containing the target genes
        """

        ## Logic to complement list of genes in case the user wants to
        ## a) analyze all genes
        ## b) analyze top n most epxressed genes
        ## c) analyze top n genes with more expression variability

        if "use_all_genes" in self.configs["params"]:
            if self.configs["params"]["use_all_genes"] is True:
                logger.warn(
                    "Spatial analysis will include all genes. Downstream analysis will take longer than usual"
                )
                self.target_genes.extend(self.gene_matrix_df.columns)

        if "n_top_expressed_genes" in self.configs["params"]:
            top_n = self.configs["params"]["n_top_expressed_genes"]
            # Logic that identifies the top_n most expressed genes
            logger.info(
                f"Spatial analysis will include the {top_n} most expressed genes."
            )
            top_n_genes = self.identify_top_n_genes(top_n, mode="expression")
            self.target_genes.extend(top_n_genes)

        if "n_top_variability_genes" in self.configs["params"]:
            top_n_var = self.configs["params"]["n_top_variability_genes"]
            # Logic that identifies the top_n_var genes with the highest variability in their expression
            logger.info(
                f"Spatial analysis will include the {top_n_var} genes with highest variability accross spots."
            )
            top_var_genes = self.identify_top_n_genes(top_n_var, mode="variability")
            self.target_genes.extend(top_var_genes)

        # Remove duplicates and sort alphabetically
        self.target_genes = list(set(self.target_genes))
        self.target_genes.sort()

        logger.info(f"Target list length: {len(self.target_genes)}")

        logger.info("<filter_gene_data> Filtering gene matrix for target genes")
        gene_df = self.gene_matrix_df
        gene_df.columns = [s.upper() for s in gene_df.columns]
        self.target_genes = [s.upper() for s in self.target_genes]
        matching = []
        for target_gene in self.target_genes:
            matching += [s for s in gene_df.columns if target_gene in s]
        if len(matching) == 0:
            matching = random.choices(gene_df.columns, k=5)

        # Remove duplicates by converting to a set and back to a list
        matching = sorted(list(set(matching)))

        logger.info(f"Matching list length: {len(matching)}")

        filtered_gene_matrix = gene_df[matching]

        # dropping duplicated columns
        filtered_gene_matrix = filtered_gene_matrix.loc[
            :, ~filtered_gene_matrix.columns.duplicated()
        ].copy()
        self.filtered_gene_matrix = filtered_gene_matrix

    def create_cell_outlines_layer(self, pipelines_enabled):
        """
        Creates a cells outlines layer
        """
        if (pipelines_enabled.get("cell2spot", False)) or self.run_data_merge_only:
            logger.info("<create_cell_outlines_layer> Creating cells layer")
            cells_mask = np.zeros(shape=self.wsi_shape)
            cells_df = self.cells_df.copy()
            # Converting the cell_polygons column from strings to proper Shapely format
            cells_df["cell_polygons"] = cells_df["cell_polygons"].apply(wkt.loads)
            cells_polys = (
                cells_df["cell_polygons"]
                .apply(lambda row: np.int32(row.exterior.coords))
                .tolist()
            )
            # Draw the cell outlines on the mask
            cv2.polylines(cells_mask, cells_polys, 1, color=(255, 0, 0))
            self.cells_layer = Image.fromarray(cells_mask.astype(np.uint8))
        else:
            logger.warning(
                "<create_cell_outlines_layer> Cells layer not created since "
                "cells_df is not available."
            )

    def create_visium_spots_layer(self, pipelines_enabled):
        """
        Creates a Visium spot outlines layer
        """
        if (pipelines_enabled.get("cell2spot", False)) or self.run_data_merge_only:
            logger.info("<create_visium_spots_layer> Creating spots layer")
            spots_mask = np.zeros(shape=self.wsi_shape)
            spots_df = self.spots_df.copy()
            # Converting the spot_polygons column from strings to proper Shapely format
            spots_df["spot_polygons"] = spots_df["spot_polygons"].apply(wkt.loads)
            spot_polys = (
                spots_df["spot_polygons"]
                .apply(lambda row: np.int32(row.exterior.coords))
                .tolist()
            )
            # Draw the spot outlines on the mask
            cv2.fillPoly(spots_mask, spot_polys, color=(255, 0, 0))
            self.spots_layer = Image.fromarray(spots_mask.astype(np.uint8))
        else:
            logger.warning(
                "<create_visium_spots_layer> Spots layer not created since "
                "spots_df is not available."
            )

    def load_kmeans_clustering_results(self, prep_dir):
        """
        Function to load kmeans files
        """
        kmeans_files = [fname for fname in os.listdir(prep_dir) if "kmeans" in fname]
        # Load kmeans dataframes into one dataframe
        kmeans_df = pd.DataFrame()
        # Loop through all the dataframes and concatenate them
        for kmeans_file in kmeans_files:
            f_contents = pd.read_csv(prep_dir + kmeans_file)
            cluster_name = kmeans_file.replace(".csv", "")
            # Renaming columns to avoid conflicts
            f_contents.columns = [self.spot_index, cluster_name]
            f_contents = f_contents.set_index(self.spot_index).astype(str)
            f_contents[cluster_name] = (
                f_contents[cluster_name].astype(str) + " (Gene cluster)"
            )
            # Concatenating the dataframes
            kmeans_df = pd.concat([kmeans_df, f_contents], axis=1)
        if self.validate_df(kmeans_df):
            self.kmeans_df = kmeans_df.sort_index()
        else:
            self.kmeans_df = None

    def validate_df(self, data_frame):
        """Function to validate the dataframe

        Args:
            data_frame (_type_): _description_

        Returns:
            _type_: _description_
        """
        return data_frame is not None and not data_frame.empty

    def create_spots_adata(self):
        """Function to create a spots anndata object from the spots_df dataframe

        cells_adata:
            uns:
                overal_qc              # overall tissue QC metrics
            varm:
                qc                     # gene-level QC metrics
            obsm:
                spatial                # spatial coordinates
                cell2spot              # cell to spot pipeline results
                deconv_results         # deconvolution results
                qc                     # spot-level QC metrics
                filtered_genes         # filtered gene matrix based on n_top_<expressed/variability>_genes
                UMAP                   # UMAP coordinates
                TSNE                   # tSNE coordinates
                kmeans_clustering      # kmeans clustering results
                graphical_clustering   # graph clustering results
        """
        # Filter the gene_matrix anndata to only include the in-tissue barcodes
        if self.validate_df(self.tissue_pos_csv):
            spatial_coords_df = self.tissue_pos_csv
            in_tissue_barcodes = self.in_tissue_barcodes
            adata = self.gene_matrix_adata[in_tissue_barcodes]
        else:
            spatial_coords_df = None
            adata = self.gene_matrix_adata
            in_tissue_barcodes = adata.obs.index

        adata.var_names_make_unique()
        adata.X = adata.X.A

        # Adding data to spots AnnData object
        if self.validate_df(spatial_coords_df):
            # Adding the spot coordinates
            adata.obsm["spatial"] = spatial_coords_df[
                ["pxl_row_in_fullres", "pxl_col_in_fullres"]
            ].astype(np.int16)
        if self.validate_df(self.spots_df):
            # Adding the cell2spot assignment pipeline results
            adata.obsm["cell2spot"] = self.spots_df[["num_contained_cells"]].astype(
                np.int16
            )
        if self.validate_df(self.deconv_results_df):
            # Adding the deconvolution pipeline results
            # Check if indices are identical
            if not self.deconv_results_df.index.equals(pd.Index(adata.obs_names)):
                obs_names_index = pd.Index(adata.obs_names)
                # Find indices in deconv_results_df not in adata.obs_names
                diff_df_not_in_obs = self.deconv_results_df.index.difference(
                    obs_names_index
                )
                # Find indices in adata.obs_names not in deconv_results_df
                diff_obs_not_in_df = obs_names_index.difference(
                    self.deconv_results_df.index
                )

                logger.warning(
                    f"{len(diff_df_not_in_obs)} barcodes do not match between the spatial "
                    "coordinates and the cell deconvolution outputs"
                )
                logger.warning(
                    "Only the barcodes present in the spatial coordinates list will be kept."
                )
                logger.warning("Any missing barcode will be filled with 0s.")

                # Reindex self.deconv_results_df based on adata.obs_names, this automatically aligns,
                # adds missing indices and drops extra indices.
                # Fill with 0s the missing data
                self.deconv_results_df = self.deconv_results_df.reindex(adata.obs_names)
                self.deconv_results_df.fillna(0, inplace=True)
                logger.warning("Barcodes have been aligned")

            adata.obsm["deconv_results"] = self.deconv_results_df.astype(np.float16)
        if self.validate_df(self.spot_qc):
            # Adding gene QC results
            gene_qc = self.gene_qc
            gene_qc["qc_mitochondrial_genes"] = gene_qc[
                "qc_mitochondrial_genes"
            ].astype(float)
            gene_qc_df = pd.DataFrame(index=adata.var.index).merge(
                gene_qc, how="left", left_index=True, right_index=True
            )
            adata.varm["qc"] = gene_qc_df.loc[adata.var.index]
            # Adding spot QC results
            adata.obsm["qc"] = self.spot_qc
            # Adding overall QC results
            adata.uns["overall_qc"] = self.overall_qc

        if self.validate_df(self.filtered_gene_matrix):
            # Adding filtered gene matrix (for spatial analysis)
            adata.obsm["filtered_genes"] = self.filtered_gene_matrix.loc[
                in_tissue_barcodes
            ]
        if self.validate_df(self.umap_df):
            # Adding UMAP projection results
            adata.obsm["UMAP"] = self.umap_df.astype(np.float16)
        if self.validate_df(self.tsne_df):
            # Adding t-SNE projection results
            adata.obsm["TSNE"] = self.tsne_df.astype(np.float16)
        if self.validate_df(self.kmeans_df):
            # Adding kmeans clustering results
            adata.obsm["kmeans_clustering"] = self.kmeans_df.astype("category")
        if self.validate_df(self.graph_clusters):
            # Adding graph clustering results
            adata.obsm["graphical_clustering"] = self.graph_clusters.astype("category")
        self.spots_adata = adata

    def create_cell_assignment_df(self, adata):
        """Function to create a cell assignment dataframe

        Args:
            adata (_type_): _description_
        """
        cell_assign_cols = ["barcode", "spot_clusters", "cell_type"]
        # Create an empty dataframe with the same index as the adata
        cell_assign_df = pd.DataFrame(adata.obs.index).set_index(self.cell_index)
        # Adding the spot clustering results to the dataframe
        cell_assign_df = cell_assign_df.merge(
            self.spot_clustering, how="left", on=self.cell_index
        )
        # Computing the feature columns
        feature_cols = set(cell_assign_df.columns) - set(cell_assign_cols)

        # Breaking the dataframe into a results and a features dataframe
        cell_assignment_results = cell_assign_df[["cell_type"]].fillna("n/a")
        cell_assignment_results["cell_type"] = cell_assignment_results[
            "cell_type"
        ].str.lower()
        cell_assignment_features = cell_assign_df[list(feature_cols)]
        return cell_assignment_results.astype(str), cell_assignment_features

    def create_cells_adata(self):
        """Function to create a cells anndata object from the cells_df dataframe
        morphological cluster

        cells_adata:
            obsm:
                spatial                         # spatial coordinates
                cell2spot                       # cell to spot pipeline results
                morphological_clustering        # morphological clustering results
                cell_assignment                 # cell assignment results
                cell_assignment_features        # cell assignment features matrix
        """
        adata = anndata.AnnData()
        if self.validate_df(self.cells_df):
            # Adding cell2spot pipeline results
            adata = anndata.AnnData(self.cell_spatial_df, dtype="int16")
            adata.obs.index = self.cell_spatial_df.index.astype(int)
            adata.obsm["spatial"] = self.cell_spatial_df.astype(np.int16)
            adata.obsm["cell2spot"] = self.cells_df[["barcode"]].astype("category")
            if self.validate_df(self.morph_clustering):
                # Adding morphological clustering results
                adata.obsm["morphological_clustering"] = self.morph_clustering.loc[
                    adata.obs.index
                ].astype("category")
            if self.validate_df(self.spot_clustering):
                # Adding spot clustering results
                results_df, _ = self.create_cell_assignment_df(adata)
                adata.obsm["cell_assignment"] = results_df.astype("category")
        adata.obs_names = adata.obs_names.astype(str)
        self.cells_adata = adata

    def summarize_spots_adata_to_dataframe(self):
        """
        Summarizes the spots anndata object to a dataframe
        """
        empty_df = pd.DataFrame()
        # Concatenating the obsm entries in the AnnData object into a dataframe
        merged_spots_df = pd.concat(
            [
                self.spots_adata.obsm.get("spatial", empty_df),
                self.spots_adata.obsm.get("cell2spot", empty_df).get(
                    "num_contained_cells", empty_df
                ),
                self.spots_adata.obsm.get("TSNE", empty_df),
                self.spots_adata.obsm.get("UMAP", empty_df),
                self.spots_adata.obsm.get("graphical_clustering", empty_df),
                self.spots_adata.obsm.get("kmeans_clustering", empty_df),
                self.spots_adata.obsm.get("qc", empty_df),
                self.spots_adata.obsm.get("deconv_results", empty_df),
                self.spots_adata.obsm.get("filtered_genes", empty_df),
            ],
            axis=1,
        )
        self.merged_spots_df = merged_spots_df.reset_index(names=self.spot_index)

    def summarize_cells_adata_to_dataframe(self):
        """
        Summarizes the spots anndata object to a dataframe
        """
        empty_df = pd.DataFrame()
        # Concatenating the obsm entries in the AnnData object into a dataframe
        merged_cells_df = pd.concat(
            [
                self.cells_adata.obsm.get("spatial", empty_df),
                self.cells_adata.obsm.get("cell2spot", empty_df).get(
                    "barcode", empty_df
                ),
                self.cells_adata.obsm.get("morphological_clustering", empty_df).get(
                    "morphological_cluster", empty_df
                ),
                self.cells_adata.obsm.get("cell_assignment", empty_df).get(
                    "cell_type", empty_df
                ),
            ],
            axis=1,
        )
        self.merged_cells_df = merged_cells_df.reset_index(names=self.cell_index)

    def create_tmap(self, results_dir, exp_id, annotation_file):
        """Function to create a TMAP object from the results directory

        Args:
            results_dir (str): Path to the results directory
            exp_id (str): Experiment ID
        """
        tmap_gen = tmapGenerator()
        tmap_gen.create_tmap_file(
            run_ids=[exp_id],
            spots_adata=self.spots_adata,
            cells_adata=self.cells_adata,
            results_dir=results_dir,
            annotation_file=annotation_file,
        )

    def save_data(self, results_dir, pipelines_enabled):
        """Save data outputs"""
        self.spots_adata.write(results_dir + "spots_adata.h5", compression="gzip")
        self.cells_adata.write(results_dir + "cells_adata.h5", compression="gzip")
        self.merged_spots_df.to_csv(results_dir + "merged_spots_df.csv")
        self.merged_cells_df.to_csv(results_dir + "merged_cells_df.csv")
        if (pipelines_enabled.get("cell2spot", False)) or self.run_data_merge_only:
            logger.info(
                "<save-data> save (spots_adata, cells_adata, merged_spots_df, "
                "merged_cells_df, spots_layer, cells_layer)"
            )
            self.spots_layer.save(results_dir + "spots_layer.png", format="png")
            self.cells_layer.save(results_dir + "cells_layer.png", format="png")
        else:
            logger.warning(
                "<save_data> Did not save cell and spot layers since cell2spot pipeline "
                "is disabled."
            )


if __name__ == "__main__":
    exp_id = "ZZ268953_V11T09-086D"
    user_id = "sunaal"
    tissue_type = "lung"
    species = "human"
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
    data_merger = DataMerge()
    data_merger.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[data_merger.pipeline_name],
        tissue_type=tissue_type,
        species_type=species,
    )

    logger.info(f"[{data_merger.pipeline_name}] Start pipeline")
    pipelines_enabled = config_flow["pipelines.enabled"]
    # data_io.fetch_input_files(pipeline_configs=data_merger.configs)

    data_merger.load_data(
        prep_dir=data_io.prep_dir,
        results_dir=data_io.results_dir,
        config_flow=config_flow,
        pipelines_enabled=pipelines_enabled,
    )
    data_merger.create_cell_outlines_layer(pipelines_enabled)
    data_merger.create_visium_spots_layer(pipelines_enabled)
    data_merger.filter_gene_data()
    data_merger.create_cells_adata()
    data_merger.create_spots_adata()
    data_merger.summarize_spots_adata_to_dataframe()
    data_merger.summarize_cells_adata_to_dataframe()
    data_merger.create_tmap(
        data_io.results_dir,
        exp_id,
        data_merger.configs["params"]["annotation_file"],
    )
    data_merger.save_data(
        results_dir=data_io.results_dir, pipelines_enabled=pipelines_enabled
    )

    import code

    code.interact(local=dict(globals(), **locals()))

    # data_io.save_run_configs(config_flow)
    # data_io.push_output_files(
    #     pipeline_configs=data_merger.configs, push_to_output=False
    # )
    # data_io.backup_experiment_files()
