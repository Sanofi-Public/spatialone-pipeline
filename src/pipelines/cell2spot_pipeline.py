"""
Python class to run cell to spot assignment
"""
# pylint:disable=no-member

import json
import os
from datetime import datetime

import geopandas
import hydra
import numpy as np
import pandas as pd
import shapely.geometry
import skimage
from dotenv import load_dotenv
from omegaconf import OmegaConf
from rasterio import features
from rtree import index
from scipy import stats

from src.pipelines.dataio_pipeline import DataIO
from src.utils.data_saver import save_data
from src.utils.env_variables import Env
from src.utils.logger import Logger

load_dotenv()

logger = Logger()


class Cell2Spot:
    """
    Class to compute cell to spot assignment
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "cell2spot"

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

    def load_data(self, prep_dir, results_dir, config_flow, unittest=False):
        """Load data as needed for cell to spot assignment"""
        logger.info("<load_data> load cell segmentation results")
        self.imgseg_method = config_flow["imgseg"]["model"]["name"]
        if unittest:
            self.cell_seg_npy = np.load(prep_dir + "cell_segmentation_patch.npy")
        else:
            self.cell_seg_npy = np.load(
                results_dir + f"{self.imgseg_method}_cell_segmentation_full.npy"
            )

        # Load tissue position csv
        logger.info("<load_data> load tissue position table")
        vs_headings = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]
        if unittest:
            pos_filename = "tissue_positions_list_patch.csv"
        else:
            pos_filename = "tissue_positions_list.csv"

        vs_spots = pd.read_csv(
            prep_dir + pos_filename,
            header=0,
            names=vs_headings,
        )
        self.vs_spots = vs_spots[vs_spots["in_tissue"].astype(int) == 1]
        # Load Visium spot info
        logger.info("<load_data> load Visium spot info")
        with open(
            prep_dir + "scalefactors_json.json", "r", encoding="utf-8"
        ) as spot_geometry_file:
            self.spot_geometry_info = json.loads(spot_geometry_file.read())

    def convert_mask_to_polygons(self):
        """
        Method that converts the a given cell instance segmentation mask to a GeoPandas
        dataframe with the following columns: Cell shapely polygon (polygons), and cell labels (cell_ids)
        :param : None
        :return (GeoPandas dataframe): gdf - Geopandas dataframe containing the cells, the morphological cluster
        they belong to, and the cell labels
        """
        logger.info(
            "<convert_mask_to_polygons> converting cell segmentation masks to Shapely polygons"
        )
        label, _ = skimage.measure.label(self.cell_seg_npy, return_num=True)
        polygons = []
        cell_x = []
        cell_y = []
        cell_ids = []
        cell_areas = []
        for geom, val in features.shapes(label.astype("int32")):
            if val == 0:
                continue
            poly = shapely.geometry.shape(geom)
            cell_ids.append(val)
            polygons.append(poly)
            cell_x.append(poly.centroid.xy[0].tolist()[0])
            cell_y.append(poly.centroid.xy[1].tolist()[0])
            cell_areas.append(poly.area)
        gdf = geopandas.GeoDataFrame(polygons)
        gdf["cell_ids"] = np.array(cell_ids).astype("int32")
        gdf["cell_x"] = np.array(cell_x)
        gdf["cell_y"] = np.array(cell_y)
        gdf["area"] = np.array(cell_areas)
        gdf.columns = ["cell_polygons", "cell_ids", "cell_x", "cell_y", "area"]
        gdf = gdf.set_geometry("cell_polygons")
        self.cells_df = gdf.drop_duplicates(
            subset=["cell_ids"], keep="first", ignore_index=True
        )
        
    def filter_by_area(self, mode="quantile", lower_q=0.01, upper_q=0.99, std=3):
        """
        Filteres the cells dataframe based on area using two methods: quantile and std
        """
        logger.info("<filter_by_area> Filtering cells by area")
        cells_df = self.cells_df.copy()
        if mode == "quantile":
            # Only keeps cells that have areas between the 0.01 and 0.99 quantile
            q_low = cells_df["area"].quantile(lower_q)
            q_hi = cells_df["area"].quantile(upper_q)
            df_filtered = cells_df[
                (cells_df["area"] < q_hi) & (cells_df["area"] > q_low)
            ]
        elif mode == "std":
            # Computes the z score and only keeps the entries within 3 standard
            # deviations of the mean
            df_filtered = cells_df[np.abs(stats.zscore(cells_df["area"])) < std]
        else:
            df_filtered = cells_df.copy()
        df_filtered = df_filtered.reset_index(drop=True)
        logger.info(f"Number of filtered cells = {len(cells_df) - len(df_filtered)}")
        return df_filtered

    def generate_visium_spot_polygons(self):
        """
        Creates polygons for each of the Visium spot
        """
        logger.info(
            "<generate_visium_spot_polygons> generating Shapely polygons for each Visium spot"
        )
        spot_radius = self.spot_geometry_info["spot_diameter_fullres"] // 2
        spot_polygons = []
        for _, row in self.vs_spots.iterrows():
            cell_x = row["pxl_col_in_fullres"]
            cell_y = row["pxl_row_in_fullres"]
            spot_polygon = shapely.geometry.Point(cell_x, cell_y).buffer(spot_radius)
            spot_polygons.append(spot_polygon)
        self.vs_spots["spot_polygons"] = spot_polygons
        self.spots_df = self.vs_spots.reset_index(drop=True)

    def create_cell_search_tree(self):
        """
        Creates a search tree (Rtree) to be used for cell to spot assignment
        """
        print("<create_cell_search_tree> building cell search tree")
        tree = index.Index()
        for i, line in enumerate(self.cells_df["cell_polygons"]):
            tree.insert(i, line.bounds)
        self.cell_search_tree = tree

    def get_cells_per_spot(self):
        """
        Function that computes the cell to spot assignment. It updates self.spots_df
        by adding an additional column "num_contained_cells", storing the number of
        cells within each spots. Also updates self.cells_df by adding a "barcode"
        specifying the barcode that the cells belong to.
        """
        logger.info("<get_cells_per_spot> assigning cells to Visium spots")
        num_contained_cells_list = []
        self.cells_df["barcode"] = None
        for _, row in self.spots_df.iterrows():
            spot_polygon = row["spot_polygons"]
            contained_cells_ix = []
            for i in self.cell_search_tree.intersection(spot_polygon.bounds):
                if spot_polygon.contains(self.cells_df["cell_polygons"][i].centroid):
                    contained_cells_ix.append(i)
            num_contained_cells_list.append(len(contained_cells_ix))
            self.cells_df.loc[contained_cells_ix, "barcode"] = row["barcode"]
        self.spots_df["num_contained_cells"] = num_contained_cells_list

    def save_data(self, results_dir):
        """Save data outputs"""
        logger.info("<save-data> save (cells_df, spots_df)")
        self.spots_df.to_csv(results_dir + "spots_df.csv")
        self.cells_df.to_csv(results_dir + "cells_df.csv")


if __name__ == "__main__":
    # sr_id = os.environ["TEST_SR_ID"]
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
    cell2spot = Cell2Spot()
    cell2spot.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[cell2spot.pipeline_name]
    )

    # logger.info(f"[{cell2spot.pipeline_name}] Fetch relevant files")
    # data_io.fetch_input_files(pipeline_configs=cell2spot.configs)

    cell2spot.load_data(
        prep_dir=data_io.prep_dir,
        results_dir=data_io.results_dir,
        config_flow=config_flow,
    )
    cell2spot.convert_mask_to_polygons()
    cell2spot.generate_visium_spot_polygons()
    cell2spot.create_cell_search_tree()
    cell2spot.get_cells_per_spot()
    cell2spot.save_data(results_dir=data_io.results_dir)
    # Push all output files from local cache to s3
    # data_io.push_output_files(pipeline_configs=cell2spot.configs)
