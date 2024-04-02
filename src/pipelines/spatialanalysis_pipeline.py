"""
Python class to run spatial structure analysis.
"""
# pylint:disable=no-member

import json
import os
from collections import defaultdict
from datetime import datetime

import anndata as ad
import hydra
import numpy as np
import pandas as pd
import squidpy as sq
from dotenv import load_dotenv
from omegaconf import OmegaConf
from PIL import Image
from shapely.geometry import Point, Polygon

from src.pipelines.dataio_pipeline import DataIO
from src.spatial_report_generator.comparative_analysis import ComparativeAnalysis
from src.spatial_report_generator.plot_generator import PlotGenerator
from src.utils.logger import Logger

load_dotenv()

logger = Logger()


class SpatialAnalysis:
    """
    Python class to run spatial structure analysis.
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "spatialanalysis"

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

    def load_data(self, prep_dir, results_dir, unittest=False):
        """
        Loading required files and saving them as class variables.
        """
        self.prep_dir = prep_dir
        self.results_dir = results_dir
        logger.info("<load_data> Loading cells and spots data frames")
        self.cells_df = pd.read_csv(results_dir + "merged_cells_df.csv")
        self.cells_df["cell_type"] = self.cells_df.cell_type.fillna("not deconvolved")
        self.spots_df = pd.read_csv(results_dir + "merged_spots_df.csv")
        self.clustering_column = (
            "graphical_clustering"
            if "graphical_clustering" in self.spots_df.columns
            else "visium_cluster"
        )

        logger.info("<load_data> Loading low and high resolution tissue images")
        with Image.open(prep_dir + "tissue_hires_image.png") as hires_img:
            self.hires = np.asarray(hires_img, dtype="uint8")
        with Image.open(prep_dir + "tissue_lowres_image.png") as lowres_img:
            self.lowres = np.asarray(lowres_img, dtype="uint8")

        logger.info("<load_data> load Visium spot info")
        sf_path = prep_dir + "scalefactors_json.json"
        with open(sf_path, "r", encoding="utf-8") as spot_geometry_file:
            self.scale_factors = json.loads(spot_geometry_file.read())

        logger.info("<load_data> load expert annotations file")
        annotations_path = prep_dir + "annotations.geojson"
        if os.path.exists(annotations_path):
            with open(annotations_path) as annotations:
                file_contents = annotations.read()
                self.expert_annotations = json.loads(file_contents)
                if type(self.expert_annotations) == dict:
                    self.expert_annotations = self.expert_annotations["features"]
        else:
            self.expert_annotations = None

        logger.info("<load_data> load region-level data")
        self.load_region_level_data()
        logger.info("<load_data> load cluster-level data")
        self.load_cluster_level_data()
        self.report_template_path = "templates/report_template.html"
        self.report_directory = os.path.join(self.results_dir, "reports")
        os.makedirs(os.path.join(self.results_dir, "reports"), exist_ok=True)

    def get_subset_tables(self, region_shapely, cells_df=None, spots_df=None):
        """
        Getting subset of cells and spots that are within the given polygon.

        Args:
            region_shapely (shapely.geometry.Polygon): Polygon to be used for subset.

        Returns:
            cells_subs (pandas.DataFrame): Subset of cells.
            spots_subs (pandas.DataFrame): Subset of spots.
        """
        # Getting the subset of cells and spots that are within the region
        min_x, min_y, max_x, max_y = region_shapely.bounds
        cells_df = self.cells_df if cells_df is None else cells_df
        spots_df = self.spots_df if spots_df is None else spots_df
        cells_subs = cells_df[
            cells_df.cell_x.between(min_x, max_x)
            & cells_df.cell_y.between(min_y, max_y)
        ]
        spots_subs = spots_df[
            spots_df.pxl_col_in_fullres.between(min_x, max_x)
            & spots_df.pxl_row_in_fullres.between(min_y, max_y)
        ]
        contained_cells_ix = []
        for i, row in cells_subs.iterrows():
            if region_shapely.contains(Point(row.cell_x, row.cell_y)):
                contained_cells_ix.append(i)
        contained_spots_ix = []
        for i, row in spots_subs.iterrows():
            if region_shapely.contains(
                Point(row.pxl_col_in_fullres, row.pxl_row_in_fullres)
            ):
                contained_spots_ix.append(i)
        cells_df = cells_df.loc[contained_cells_ix].reset_index(drop=True)
        spots_df = spots_df.loc[contained_spots_ix].reset_index(drop=True)
        return cells_df, spots_df

    def load_region_level_data(self):
        """Iterates through the regions and extracts the cell and spot level data.

        Returns:
            regions (collections.defaultdict): Dictionary of regions with the following structure:
                region_class
                `- region_nameq
                    |- num_cells
                    |- num_spots
                    |- cells_df
                    `- spots_df
        """
        regions = defaultdict(dict)
        if self.expert_annotations is None:
            self.region_level_data = regions
            return

        for i, region in enumerate(self.expert_annotations):
            region_class = region["properties"]["classification"]["name"]
            region_poly = np.array(region["geometry"]["coordinates"])
            # Handling region names
            region_name = region["properties"]["name"]
            if region_poly.dtype is not np.dtype("object"):
                region_shapely = Polygon(region_poly.reshape((-1, 2)))
            else:
                # If region has invalid structure, we skip
                # (region has two sets of coordinates)
                continue
            # Getting the subset of cells and spots that are within the region
            cells_df, spots_df = self.get_subset_tables(region_shapely)
            # Adding the cells and spots subsets to the regions dictionary
            if len(spots_df) > 1:
                regions[region_class][region_name] = {
                    "num_cells": len(cells_df),
                    "num_spots": len(spots_df),
                    "cells_df": cells_df,
                    "spots_df": spots_df,
                    "region_shapely": region_shapely,
                }
            else:
                continue
        self.region_level_data = regions

    def load_cluster_level_data(self):
        """Iterates through the visium clusters and extracts the cell and spot level data."""
        clusters = sorted(self.spots_df[self.clustering_column].unique())
        regions = defaultdict(dict)
        for cluster in clusters:
            spots_df = self.spots_df[
                self.spots_df[self.clustering_column] == cluster
            ].reset_index(drop=True)
            barcodes = spots_df.barcode
            cells_df = self.cells_df[self.cells_df.barcode.isin(barcodes)].reset_index(
                drop=True
            )
            region_name = cluster
            regions[cluster][region_name] = {
                "num_cells": len(cells_df),
                "num_spots": len(spots_df),
                "cells_df": cells_df,
                "spots_df": spots_df,
            }
        self.cluster_level_data = regions

    def cells_df_to_adata(self, exp_id, df=None, drop_not_deconvolved=True):
        """
        Converts a cells dataframe to an AnnData object.

        Args:
            exp_id(str): experiment ID of sample.
            drop_not_deconvolved(bool, optional): whether to drop cells that have not been deconvolved. Defaults to True.
        """
        logger.info("<cells_df_to_adata> Converting cells dataframe to AnnData")
        if df is None:
            df = self.cells_df.copy()
        # Cleaning up the cells dataframe before proceeding
        df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
        # Removing not deconvolved cells

        if drop_not_deconvolved:
            non_deconv_label = self.configs["params"]["non_deconv_label"]
            df = df[df.cell_type != non_deconv_label].reset_index(drop=True)

        # Generating distinct colors for the different cell types
        colors = PlotGenerator.generate_distinct_colors(len(df.cell_type.unique()))

        # Converting the cell types column to categorical variables
        cell_types = pd.Categorical(
            df.cell_type.str.replace(" ", "-").str.replace(",", "")
        )
        if len(cell_types) > 0:
            cell_type_one_hot = pd.get_dummies(cell_types)
            df = df.join(cell_type_one_hot)

            # Create Anndata from the cells dataframe
            adata = ad.AnnData(df[list(cell_type_one_hot.columns)], dtype=np.float32)

            # Defining the observation names (barcodes in this case)
            adata.obs_names = [f"Cell_{i:d}" for i in df.cell_ids]
            # Defining the variable names (genes in this case)
            adata.var_names = cell_types.categories

            # Adding cell-level observations
            cell_cluster_key = self.configs["params"]["cell_cluster_key"]
            adata.obs[cell_cluster_key] = cell_types

            # Adding cell type colors to Anndata unstructured data
            adata.uns[f"{cell_cluster_key}_colors"] = colors

            # Adding cell coordinates to the observation matrix
            adata.obsm["spatial"] = df[["cell_x", "cell_y"]].to_numpy()
            adata.uns["spatial"] = {
                exp_id: {
                    "metadata": {},
                    "images": {"hires": self.hires, "lowres": self.lowres},
                    "scalefactors": self.scale_factors,
                }
            }
            return adata
        else:
            return None

    def spots_df_to_adata(self, exp_id, df=None):
        """
        Converts a spots dataframe to an AnnData object.

        Args:
            exp_id(str): experiment ID of sample
        """
        logger.info("<spots_df_to_adata> Converting spots dataframe to AnnData")
        if df is None:
            df = self.spots_df.copy()
        # Converting each visium cluster into a categorical variable
        cluster_categories = pd.Categorical(
            df[self.clustering_column].astype(str).str.replace(" ", "-")
        )
        df[self.clustering_column] = cluster_categories
        # Getting the gene columns - Look for better ways of doing this. Currently looks for capitalized columns
        cell_columns = list(self.cells_df.cell_type.unique())
        cell_columns.remove(self.configs["params"]["non_deconv_label"])
        cell_columns = list(set(cell_columns) & set(df.columns))

        gene_columns = df.columns[
            [
                (col.isupper() and "-" not in col and "." not in col)
                for col in df.columns
            ]
        ].tolist()
        gene_columns = list(dict.fromkeys(gene_columns))
        gene_columns = list(set(gene_columns) - set(cell_columns))
        gene_columns = list(set(gene_columns) & set(df.columns))

        columns = gene_columns + cell_columns

        # Defining distinct colors for the different Visium clusters
        colors = PlotGenerator.generate_distinct_colors(
            len(cluster_categories.categories)
        )

        # Create Anndata from the gene matrix
        spot_data = df[columns]
        spot_data.index = df["barcode"]
        adata = ad.AnnData(spot_data, dtype=np.float32)
        # Defining the observation names (barcodes in this case)
        adata.obs_names = df.barcode.tolist()
        # Defining the variable names (genes in this case)
        adata.var_names = columns

        # Adding spot-level observations - coordinates, number of contained cells, t-SNE, UMAP, cluster info
        obs_columns = ["num_contained_cells", "visium_cluster"]

        # Check if 'TSNE-1', 'TSNE-2', 'UMAP-1', and 'UMAP-2' columns exist in the DataFrame
        if "TSNE-1" in df.columns and "TSNE-2" in df.columns:
            obs_columns.extend(["TSNE-1", "TSNE-2"])
        if "UMAP-1" in df.columns and "UMAP-2" in df.columns:
            obs_columns.extend(["UMAP-1", "UMAP-2"])

        # Adding gene variables - total number of times gene is counted
        gene_sums = df[columns].sum(axis=0)
        adata.var["sum_counts"] = pd.DataFrame(gene_sums)

        # Adding cluster colors to Anndata unstructured data
        adata.uns["visum_cluster_colors"] = colors

        # Adding Visium spot coordinates to the observation matrix
        adata.obsm["spatial"] = df[
            ["pxl_col_in_fullres", "pxl_row_in_fullres"]
        ].to_numpy()

        # Storing the cells and gene columns for later reference
        adata.uns["gene_columns"] = gene_columns
        adata.uns["cell_columns"] = cell_columns

        # Adding images and scale factors
        adata.uns["spatial"] = {
            exp_id: {
                "metadata": {},
                "images": {"hires": self.hires, "lowres": self.lowres},
                "scalefactors": self.scale_factors,
            }
        }
        return adata

    def build_neighbours_graph(self, adata, marker_type="cells", level="tissue"):
        """
        Builds neighbours graph for cells and spot tables.

        Args:
            marker_type (str, optional): type of markers system. Options: ["cells", "spots"]. Defaults to "cells".
        """
        logger.info("<build_neighbours_graph> Building spatial neighbours graph")
        params = self.configs["params"][level]
        if marker_type == "cells":
            # Need to check if we have enough cells to construct a neighborhood graph
            if len(adata.obs) < params["n_cells_neighbors"]:
                return adata
            sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
        else:
            n_rings = params["n_rings"]
            n_neighs = params["n_neighs"]
            # Need to check if we have enough spots to construct a neighborhood graph
            if len(adata.obs) < n_neighs + 1:
                return adata
            sq.gr.spatial_neighbors(
                adata, n_rings=n_rings, coord_type="grid", n_neighs=n_neighs
            )
        return adata

    def neighborhood_enrichment_analysis(self, adata, level="tissue"):
        """
        Runs neighborhood enrichment analysis on the cell-level AnnData object.
        Prerequisite: run SpatialStructAnalysis.build_neighbours_graph().
        """
        logger.info(
            "<neighborhood_enrichment_analysis> Running neighbourhood enrichment analysis on the cells AnnData"
        )
        params = self.configs["params"][level]
        if not (params["cell_net_matrix"] or params["cell_net_chord"]) or (
            "spatial_neighbors" not in adata.uns.keys()
        ):
            return adata

        # Running neighborhood enrichment analysis
        cell_cluster_key = self.configs["params"]["cell_cluster_key"]
        sq.gr.nhood_enrichment(adata, cluster_key=cell_cluster_key, seed=1)
        # Post-processing output Z scores to replace -inf and +inf with the min and max, respectively
        zscore_df = pd.DataFrame(
            adata.uns[f"{cell_cluster_key}_nhood_enrichment"]["zscore"]
        )
        min_zscore = zscore_df.min().min()
        zscore_df = zscore_df.fillna(min_zscore).replace([np.inf, -np.inf], min_zscore)
        adata.uns[f"{cell_cluster_key}_nhood_enrichment"][
            "zscore"
        ] = zscore_df.to_numpy()
        return adata

    def cooccurence_analysis(self, adata, level="tissue"):
        """
        Runs co-occurence analysis on the cells AnnData object.
        """
        logger.info(
            "<cooccurence_analysis> Running co-occurence analysis on the cells AnnData"
        )
        params = self.configs["params"][level]
        cell_cluster_key = self.configs["params"]["cell_cluster_key"]
        # PUT ELSE STATEMENTS
        if not params["cell_cooccur"]:
            return adata
        else:
            n_splits = params["n_splits"]
            n_intervals = params["n_intervals"]
            intervals = [int((0.5 * i) ** 3) for i in range(1, n_intervals)]
            sq.gr.co_occurrence(
                adata,
                cluster_key=cell_cluster_key,
                n_splits=n_splits,
                interval=intervals,
            )
            return adata

    def morans_i_analysis(self, adata, level="tissue", mode="gene"):
        """
        Runs Moran's I analysis on the spots AnnData object.
        Prerequisite: run SpatialStructAnalysis.build_neighbours_graph().
        """
        logger.info(
            "<morans_i_analysis> Running Moran's I analysis on the spots AnnData object"
        )
        params = self.configs["params"][level]
        run_flag = params[f"morans_{mode}_heatmap"] or params[f"morans_{mode}_bar"]
        if not run_flag or ("spatial_neighbors" not in adata.uns.keys()):
            return adata
        markers = adata.uns[f"{mode}_columns"]
        sq.gr.spatial_autocorr(
            adata,
            mode="moran",
            genes=markers,
            n_perms=params["n_perms"],
            n_jobs=8,
            seed=1,
        )
        adata.uns[f"{mode}_moranI"] = adata.uns["moranI"]
        return adata

    def summarize_cells_within_spots(self, cells_df=None, drop_not_deconvolved=True):
        """
        Summarize cells_df at barcode and cell_type level.
        Generate statistics for each cell column in df_barcode (mean, median, max, min, std).

        Args:
            drop_not_deconvolved(bool, optional): whether to drop cells that have not been deconvolved. Defaults to True.
        """
        logger.info("<summarize_cells_within_spots> Summarizing cells within spots")
        if cells_df is None:
            cells_df = self.cells_df
        df = cells_df.copy()
        if drop_not_deconvolved:
            # df[~df.datecolumn.isin(a)]
            non_deconv_label = self.configs["params"]["non_deconv_label"]
            df = df[df.cell_type != non_deconv_label].reset_index(drop=True)
            df = df[df.cell_type != "undefined"].reset_index(drop=True)
            df = df[df.cell_type != "Undefined"].reset_index(drop=True)
        # summarize df at barcode and cell_type level
        df_barcode = pd.DataFrame()
        df_barcode["cell_type"] = df.cell_type
        df_barcode["cell_type"] = df_barcode["cell_type"].astype("category")
        df_barcode["barcode"] = df.barcode
        df_barcode = (
            df_barcode.groupby(["barcode", "cell_type"])
            .size()
            .reset_index(name="counts")
        )
        df_barcode = df_barcode.pivot(
            index="barcode", columns="cell_type", values="counts"
        )
        df_barcode = df_barcode.fillna(0).astype(int)

        # generate statistics for each cell column in df_barcode (mean, median, max, min, std)
        df_barcode_stats = pd.DataFrame()
        df_barcode_stats["cell_count"] = df_barcode.sum(axis=0)
        df_barcode_stats["avg_spot"] = df_barcode.mean(axis=0)
        df_barcode_stats["median_spot"] = df_barcode.median(axis=0)
        df_barcode_stats["max_spot"] = df_barcode.max(axis=0)
        df_barcode_stats["min_spot"] = df_barcode.min(axis=0)
        df_barcode_stats["std_spot"] = df_barcode.std(axis=0)
        df_barcode_stats["p90_spot"] = df_barcode.quantile(0.9, axis=0)
        df_barcode_stats["p10_spot"] = df_barcode.quantile(0.1, axis=0)
        df_barcode_stats["cell_type"] = df_barcode_stats.index.values
        df_barcode_stats = df_barcode_stats.reindex(
            columns=[
                "cell_type",
                "cell_count",
                "avg_spot",
                "median_spot",
                "max_spot",
                "min_spot",
                "std_spot",
                "p90_spot",
                "p10_spot",
            ]
        )
        df_barcode_stats = df_barcode_stats.round(2)
        df_barcode_stats = df_barcode_stats.sort_values("cell_count", ascending=False)
        return df_barcode_stats, df_barcode

    def expand_shrink_region(self, region_shapely):
        """Expands and shrinks the given region to create distance intervals for study.

        Args:
            region_shapely (shapely.Polygon): Region polygon.

        Returns:
            list: list of dictionaries storing the polygons that have been expanded or
            shrunken, and their corresponding distance from the original polygon.
        """
        logger.info("<expand_shrink_region> Creating distance intervals for study")
        params = self.configs["params"]["region"]
        num_levels = params["infilt_levels"]
        step_size = int(self.scale_factors["spot_diameter_fullres"])
        level = num_levels * step_size
        shrink = False
        region = region_shapely.buffer(level, join_style=3)
        regions = [
            {
                "region": region,
                "level": level,
                "level_range": (level - step_size, level),
            }
        ]
        while True:
            shrink = True if level == 0 else False
            if shrink:
                level = -np.abs(level + step_size)
                region = region_shapely.buffer(level, join_style=3)
            else:
                level = level - step_size
                region = region_shapely.buffer(level, join_style=3)
            if region.is_empty:
                break
            regions += [{"region": region, "level": level}]
        return regions

    def infiltration_analysis(self, region_shapely):
        """Run infiltration analysis for the given region.

        Args:
            region_shapely (shapely.Polygon): Region polygon.

        Returns:
            pandas.DataFrame: dataframe with the results of infiltration analysis.
        """
        logger.info("<infiltration_analysis> Running infiltration analysis on region")
        if region_shapely is None:
            return None
        regions = self.expand_shrink_region(region_shapely)
        cells_df = self.cells_df.copy()
        spots_df = self.spots_df.copy()
        non_deconv_label = self.configs["params"]["non_deconv_label"]
        cells_df = cells_df[cells_df["cell_type"] != non_deconv_label].reset_index(
            drop=True
        )
        region_cells_df, region_spots_df = self.get_subset_tables(
            regions[0]["region"], cells_df, spots_df
        )
        level_cells = region_cells_df.copy()
        level_spots = region_spots_df.copy()
        region_cells_df = region_cells_df.set_index("cell_ids")
        region_spots_df = region_spots_df.set_index("barcode")
        region_cells_df["infilt_level"] = None
        region_spots_df["infilt_level"] = None
        for region in regions:
            level_cells, level_spots = self.get_subset_tables(
                region["region"], level_cells, level_spots
            )
            region_cells_df.loc[level_cells["cell_ids"], "infilt_level"] = region[
                "level"
            ]
            region_spots_df.loc[level_spots["barcode"], "infilt_level"] = region[
                "level"
            ]
        return region_cells_df

    def run_spatial_analysis(
        self,
        exp_id,
        flow_id,
        cells_df=None,
        spots_df=None,
        region_name=None,
        region=None,
        level="tissue",
    ):
        """Runs spatial analysis on the provides cells and spots dataframes

        Args:
            cells_df (pandas.DataFrame, optional): Cell-level dataframe. Defaults to None.
            spots_df (pandas.DataFrame, optional): Spot=level dataframe. Defaults to None.
            level (str, optional): "tissue" or "region" level analysis. Defaults to "tissue".

        Returns:
            report (str): HTML report for the spatial analysis.
            cells_adata (AnnData): Cell-level AnnData object.
            spots_adata (AnnData): Spot-level AnnData object.
        """
        logger.info("<run_spatial_analysis> Creating spatial analysis report")
        if cells_df is None and spots_df is None:
            cells_df = self.cells_df.copy()
            spots_df = self.spots_df.copy()
        # Cell-level analysis
        logger.info("<run_spatial_analysis> Cell-level spatial structure analysis")
        cells_adata = self.cells_df_to_adata(exp_id=exp_id, df=cells_df)
        if cells_adata is None:
            return None, None, None
        cells_adata = self.build_neighbours_graph(
            adata=cells_adata, marker_type="cells", level=level
        )
        cells_adata = self.neighborhood_enrichment_analysis(
            adata=cells_adata, level=level
        )
        cells_adata = self.cooccurence_analysis(adata=cells_adata, level=level)
        spot_summary = self.summarize_cells_within_spots(cells_df=cells_df)

        # Spot-level analysis
        logger.info("<run_spatial_analysis> Spot-level spatial structure analysis")
        spots_adata = self.spots_df_to_adata(exp_id=exp_id, df=spots_df)
        if spots_adata is None:
            return None, None, None
        spots_adata = self.build_neighbours_graph(
            adata=spots_adata, marker_type="spots", level=level
        )
        spots_adata = self.morans_i_analysis(
            adata=spots_adata, level=level, mode="gene"
        )
        spots_adata = self.morans_i_analysis(
            adata=spots_adata, level=level, mode="cell"
        )

        # Comparative analysis
        params = self.configs["params"][level]
        if level == "tissue":
            genes_list = spots_adata.uns["gene_columns"]
            comp_analysis = ComparativeAnalysis(configs=self.configs)
            if params["diff_exp_annotations"]:
                spots_adata.uns[
                    "diff_exp_annotations"
                ] = comp_analysis.compare_gene_exp(
                    regions_dict=self.region_level_data,
                    genes_list=genes_list,
                    level="tissue",
                )
            if params["diff_exp_clusters"]:
                spots_adata.uns["diff_exp_clusters"] = comp_analysis.compare_gene_exp(
                    regions_dict=self.cluster_level_data,
                    genes_list=genes_list,
                    level="tissue",
                )
        # Infiltration Analysis
        if params["infilt_comparing_cell_types"] or params["infilt_comparing_levels"]:
            cells_adata.uns["infilt_analysis"] = self.infiltration_analysis(
                region_shapely=region
            )
        # Generating plots
        report = self.generate_spatial_analysis_report(
            prep_dir=self.prep_dir,
            exp_id=exp_id,
            flow_id=flow_id,
            cells_adata=cells_adata,
            spots_adata=spots_adata,
            spot_summary=spot_summary,
            cells_df=cells_df,
            spots_df=spots_df,
            region_name=region_name,
            level=level,
        )
        return report, cells_adata, spots_adata

    def run_region_spatial_analysis(self, exp_id, flow_id):
        """Runs spatial analysis on a list of regions.

        Args:
            regions (collections.defaultdict): Dictionary of regions with the following structure:
                region_class
                `- region_name
                    |- num_cells
                    |- num_spots
                    |- cells_df
                    `- spots_df

        Returns:
            region_names (str): Names of of the regions.
        """
        logger.info("<run_region_spatial_analysis> Creating reports for regions")
        region_names = []
        regions = self.region_level_data
        for region_class in regions:
            for region_name in regions[region_class]:
                logger.info(
                    f"[{self.pipeline_name}] Running spatial structure analysis for Region: {region_name}"
                )
                cells_df = regions[region_class][region_name]["cells_df"]
                spots_df = regions[region_class][region_name]["spots_df"]
                region_shapely = regions[region_class][region_name]["region_shapely"]
                report, cells_adata, spots_adata = self.run_spatial_analysis(
                    cells_df=cells_df,
                    spots_df=spots_df,
                    level="region",
                    region_name=region_name,
                    region=region_shapely,
                    exp_id=exp_id,
                    flow_id=flow_id,
                )
                if cells_adata is None or spots_adata is None:
                    continue
                self.save_data(
                    report=report,
                    cells_adata=cells_adata,
                    spots_adata=spots_adata,
                    report_name=region_name,
                    level="region",
                )
                region_names += [region_name]
        return region_names

    def get_report_name(self, region_name):
        """Gets the report name

        Args:
            region_name (str): Name of the region
        """
        if region_name is None:
            return "Whole Tissue Report"
        else:
            return f"Region: {region_name}"

    def generate_spatial_analysis_report(
        self,
        prep_dir,
        exp_id,
        flow_id,
        cells_adata,
        spots_adata,
        spot_summary,
        cells_df=None,
        spots_df=None,
        region_name=None,
        level="tissue",
    ):
        """Calls the methods in report_generator.py to create a spatial structure analysis report.

        Args:
            exp_id (str): experiment ID
            cells_adata (anndata.AnnData): cells anndata
            spots_adata (anndata.AnnData): spots anndata
            spot_summary (pandas.DataFrame): spot summary statistics table. (Results from summarize_cells_within_spots())
            cells_df (pandas.DataFrame, optional): cells dafaframe. Defaults to None.
            spots_df (pandas.DataFrame, optional): spots dafaframe. Defaults to None.
            region_name (str, optional): region name. Defaults to None.
            level (str, optional): Options ["tissue", "region"]. Defaults to "tissue".
        """

        logger.info(
            "<generate_spatial_analysis_report> Generating spatial analysis report"
        )
        spot_stats, _ = spot_summary
        if cells_df is None:
            cells_df = self.cells_df
        if spots_df is None:
            spots_df = self.spots_df

        plot_gen = PlotGenerator(
            cells_adata,
            spots_adata,
            configs=self.configs,
            level=level,
            cache_path=self.report_directory,
            report_template_path=self.report_template_path,
            scale_factors=self.scale_factors,
        )
        plot_gen.add_qc_report_button(exp_id=exp_id, flow_id=flow_id)
        plot_gen.cell_summary_table(df=spot_stats)
        plot_gen.cell_count_plot(cells_df=cells_df)
        plot_gen.cells_within_spot_plot(spot_summary=spot_summary)
        plot_gen.avg_gene_exp_plot(spots_df=spots_df)
        plot_gen.neighborhood_enrichment_plots()
        plot_gen.chord_plot(spot_stats=spot_stats)
        plot_gen.morans_i_exp_plots(task_name="morans_gene_heatmap", mode="gene")
        plot_gen.morans_i_bar_plot(task_name="morans_gene_bar", mode="gene")
        plot_gen.morans_i_exp_plots(task_name="morans_cell_heatmap", mode="cell")
        plot_gen.morans_i_bar_plot(task_name="morans_cell_bar", mode="cell")
        plot_gen.cooccurrence_plots()
        plot_gen.volcano_plots(task_name="diff_exp_annotations")
        plot_gen.volcano_plots(task_name="diff_exp_clusters")
        _, infiltrating_cells = plot_gen.abundance_vs_level()
        plot_gen.abundance_vs_cells_plot(infiltrating_cells=infiltrating_cells)
        report_name = self.get_report_name(region_name)
        report = plot_gen.get_report(exp_id=exp_id, report_name=report_name)
        if level == "tissue":
            # Exporting figure content to .csv files
            self.export_downstream_analysis_to_csv(plot_gen.inline_figs)
        return report
        return report

    def add_files_to_configs(self, bucket, file_names, extension, file_type="outputs"):
        """Adds extra files to config flow

        Args:
            file_names (list): List of file names
            extension (str): file extensions
            bucket (str): target bucket
            file_type (str, optional): One of ["outputs", "inputs"]. Defaults to "outputs".
        """
        logger.info("<add_files_to_configs> Adding files to configs")
        for file_name in file_names:
            self.configs[file_type][file_name] = {
                "bucket": bucket,
                "extension": extension,
            }

    def export_downstream_analysis_to_csv(self, fig_configs):
        """
        Exports the downstream analysis results into .csv

        Args:
            fig_configs (dict): Inline figures configurations dictionary
        """
        # Looping through all inline figures that have a "data" object and
        # saving them to a .csv file
        for fig_config in fig_configs:
            if "data" in fig_config.keys():
                data_df = fig_config["data"]
                f_name = fig_config["name"]
                data_df.to_csv(os.path.join(self.report_directory, f"{f_name}.csv"))

    def save_data(
        self, report, cells_adata, spots_adata, report_name="", level="tissue"
    ):
        """Save data outputs"""
        logger.info("<save-data> save report")
        if report_name == "":
            report_path = os.path.join(self.report_directory, "report.html")
        else:
            report_path = (
                os.path.join(
                    self.report_directory, self.configs["params"]["region_report_name"]
                )
            ).format(report_name=report_name)
        with open(report_path, mode="w", encoding="utf-8") as report_file:
            report_file.write(report)
        if level == "tissue":
            spots_adata.write(
                self.report_directory + "/spots_downstream_results.h5",
                compression="gzip",
            )
            cells_adata.write(
                self.report_directory + "/cells_downstream_results.h5",
                compression="gzip",
            )


if __name__ == "__main__":
    exp_id = "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma"
    user_id = "blanka"
    tissue_type = "skin"
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
    spatanalysis = SpatialAnalysis()
    spatanalysis.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[spatanalysis.pipeline_name]
    )
    # logger.info(f"[{spatanalysis.pipeline_name}] Fetch relevant files")
    # data_io.fetch_input_files(pipeline_configs=spatanalysis.configs)

    logger.info(f"[{spatanalysis.pipeline_name}] Start pipeline")
    spatanalysis.load_data(prep_dir=data_io.prep_dir, results_dir=data_io.results_dir)
    # Region-level analysis
    # Note: This runs on regions that have been pre-defined by the user. If no regions are
    # defined, it does nothing
    logger.info(f"[{spatanalysis.pipeline_name}] Region level analysis")
    region_names = spatanalysis.run_region_spatial_analysis(
        exp_id=exp_id, flow_id=flow_id
    )
    # Whole-tissue analysis
    logger.info(
        f"[{spatanalysis.pipeline_name}] Running spatial structure analysis for whole tissue"
    )
    report, cells_adata, spots_adata = spatanalysis.run_spatial_analysis(
        cells_df=spatanalysis.cells_df,
        spots_df=spatanalysis.spots_df,
        level="tissue",
        exp_id=exp_id,
        flow_id=flow_id,
    )
    spatanalysis.save_data(
        report=report, cells_adata=cells_adata, spots_adata=spots_adata, report_name=""
    )
