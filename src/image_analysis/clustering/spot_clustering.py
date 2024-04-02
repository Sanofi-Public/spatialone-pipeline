"""
Bead clustering functions
"""

import random

import numpy as np
import pandas as pd
from numpy import unique
from sklearn.cluster import MiniBatchKMeans

from src.image_analysis.clustering.constrained_clustering import DeterministicAnnealing
from src.image_analysis.clustering.global_assignment import (
    estimate_cluster_probability,
    estimate_spot_probability,
    get_best_allocation_simanel,
    get_combined_probabilities,
    load_cells_in_spots,
)
from src.utils.logger import Logger

logger = Logger()
# Create rng
# TODO: to be moved as a parameter
seed = 0
rng = np.random.default_rng(seed)


def create_barcode_dict(cell2location_results_path) -> dict:
    """Create a formatted dictionary with associated cell types per barcode
    Arguments:
        cell2location_results - path to cell_type csv file
    Returns:
        dict -- dictionary with barcodes and associated cell types
    """
    # add path for cell_type table
    cell_type_table = pd.read_csv(cell2location_results_path)
    cell_type_table["cell_n"] = cell_type_table.sum(axis=1)
    # )  # pylint: disable=no-member
    # Drop columns where cell_n value is 0
    cell_type = cell_type_table[cell_type_table["cell_n"] != 0]
    barcode_dict = cell_type.set_index("barcode").to_dict("index")
    for key, value in barcode_dict.copy().items():
        for cell_k, cell_n in value.copy().items():
            if cell_n == 0:
                value.pop(cell_k)
        value["distinct_cell_types"] = len(value.keys()) - 1
        if len(value.keys()) - 1 == 0:
            barcode_dict.pop(key)
    return barcode_dict


def get_barcode(group):
    """Get barcode for pandas groupby object instance
    Arguments:
        group {pd.DataFrame} -- instance of groupby object
    Returns:
        string -- barcode
    """
    tag = group.barcode.unique().tolist()
    return tag[0]


# APPLY CLUSTERING ON EACH BEAD AND APPEND RESULTS TO DATA
def bead_clustering(group, barcode_dict: dict):
    """Get barcode for pandas groupby object instance
    Arguments:
        group {pd.DataFrame} -- instance of groupby object
        barcode_dict {dict} -- dictionary with barcodes and associated cell types
    Returns:
        group {pd.DataFrame} -- modified group dataframe with clustered results
    """
    row_count = group.shape[0]
    barcode = get_barcode(group)
    barcode_dist = get_distribution(barcode_dict)
    try:
        nclusters = barcode_dict[barcode]["distinct_cell_types"]
        if nclusters == 0:
            nclusters = 1
    except KeyError:
        nclusters = 1
    # Less data than cell type count check
    if nclusters > row_count:
        nclusters = row_count
    input_data = group.drop(["cell_ids", "barcode"], axis=1)
    bead_training_data = input_data.to_numpy()
    distribution = barcode_dist.get(barcode)
    if distribution is None:
        distribution = [1]
    model = DeterministicAnnealing(
        n_clusters=nclusters, distribution=distribution, np_seed=rng
    )
    model.fit(bead_training_data, enforce_cluster_distribution=True)
    # centers = model.cluster_centers_
    labels = model.labels_
    clusters = unique(labels)
    group["spot_clusters"] = labels
    cluster_to_ctype = celltype_map(group, barcode_dict, barcode, clusters)
    group["cell_type"] = group["spot_clusters"].map(cluster_to_ctype)
    return group


def get_distribution(barcode_dict: dict):
    """Calculate distribution of cell types
    Arguments:
        barcode_dict {dict} -- dictionary with barcodes and associated cell types
    Returns:
        barcodes_dist {dict} -- dictionary with barcodes and associated cell type distribution lists
    """
    barcode_dist = {}
    for key, value in barcode_dict.items():
        # cell_types = value.get("distinct_cell_types")
        cell_n = value.get("cell_n")
        if cell_n is None:
            cell_n = 1
        spot_distribution = [
            count / cell_n
            for cell_name, count in value.items()
            if cell_name not in ["cell_n", "distinct_cell_types"]
        ]
        barcode_dist[key] = spot_distribution
    return barcode_dist


def assign_celltype(cell_df_path, reduced_df, barcode_dict: dict) -> pd.DataFrame:
    """Assign cell_type to individual cells
    Arguments:
        cell_df_path -- path to visium spot information table
        reduced_df {np.ndarray} -- transformed extracted feature table
        barcode_dict {dict} -- dictionary with barcodes and associated cell types
    Returns:
        barcodes {pd.DataFrame} -- modified dataframe with clustered results
    """
    spot_table = pd.read_csv(cell_df_path)
    feature_table = pd.DataFrame(reduced_df)
    feature_table.reset_index(inplace=True)
    feature_table.rename(columns={"index": "cell_ids"}, inplace=True)
    spot_table = spot_table.drop(["cell_x", "cell_y", "cell_polygons"], axis=1)

    if np.min(spot_table.cell_ids) > np.min(
        feature_table.cell_ids
    ):  # TODO: There seems to be a bug here. Need to better see how to handle it
        # TODO: Need to revise if it's just an indexing start issue (0 vs 1) or if it has further implciations
        logger.warn("spot_table and spot cell ids are disaligned.")
        dif = np.min(spot_table.cell_ids) - np.min(feature_table.cell_ids)
        feature_table.cell_ids = feature_table.cell_ids + dif
        feature_table.index = feature_table.index + dif
        logger.info("spot_table, spot cell and feature table have been aligned")

    spot = pd.merge(
        feature_table, spot_table[["cell_ids", "barcode"]], on="cell_ids", how="left"
    )

    logger.info(
        "spot_tablecells_id values - Max:"
        + str(np.max(spot_table.cell_ids))
        + " Min: "
        + str(np.min(spot_table.cell_ids))
    )
    logger.info(
        "spot cells_id values (result from merge) - Max:"
        + str(np.max(spot.cell_ids))
        + " Min: "
        + str(np.min(spot.cell_ids))
    )
    logger.info(
        "Feature table cells_id values- Max:"
        + str(np.max(feature_table.cell_ids))
        + " Min: "
        + str(np.min(feature_table.cell_ids))
    )

    barcodes = spot.groupby(by=["barcode"], as_index=False).apply(
        bead_clustering, (barcode_dict)
    )
    return barcodes


def celltype_map(group, barcode_dict: dict, barcode: str, clusters: np.ndarray) -> dict:
    """Mapping cell type to grouped cells in each spot
    Arguments:
        group {pd.DataFrame} -- instance of groupby object
        barcode_dict {dict} -- dictionary with barcodes and associated cell types
        barcode {str} -- barcode tag for the spot
        clusters {np.ndarray} -- unique clusters
    Returns:
        cluster_to_ctype {dict} -- spot clustered dictionary with key as cluster and value as cell type
    """
    cluster_to_ctype = {}
    deconvolved_spots = list(barcode_dict.keys())
    ctype_deconv = barcode_dict.get(barcode)
    if barcode not in deconvolved_spots:
        cell_group = "not deconvolved"
        cluster_to_ctype[0] = cell_group
        return cluster_to_ctype
    if ctype_deconv is None:
        cell_group = "not deconvolved"
        cluster_to_ctype[0] = cell_group
        return cluster_to_ctype
    cluster_counts = group.spot_clusters.value_counts().to_dict()
    ctype_deconv.pop("cell_n")
    ctype_deconv.pop("distinct_cell_types")
    ctypes = list(ctype_deconv.values())
    if len(ctypes) != len(cluster_counts):
        ctype_deconv.copy().popitem()
        ctypes = list(ctype_deconv.values())
    for cell_count in ctypes:
        cell_group = [
            cell_k for cell_k, cell_n in ctype_deconv.items() if cell_n == cell_count
        ]
        if len(cell_group) > 1:
            cell_group = str(rng.choice(np.asarray(cell_group)))
        else:
            cell_group = cell_group[0]
        similarities = {}
        count = 0
        for key, cluster_count in cluster_counts.copy().items():
            if cluster_count != cell_count:
                continue
            count += 1
            sim = abs(cluster_count - cell_count)
            if sim == 0:
                cluster_to_ctype[key] = cell_group
                del ctype_deconv[cell_group]
                del cluster_counts[key]
                break
    return cluster_to_ctype


# Naive Model
def random_classifier(group, barcode_dict: dict, random_state=0):
    """Assign randomly cell_type to individual cells
    Arguments:
        group {pd.DataFrame} -- instance of groupby object
        barcode_dict {dict} -- dictionary with barcodes and associated cell types
    Returns:
        group {pd.DataFrame} -- modified group dataframe with clustered results
    """

    # Get celltype labels for this spot
    # based on the cell counts for each celltype
    barcode = get_barcode(group)
    deconvolved_spots = list(barcode_dict.keys())
    ctype_deconv = barcode_dict.get(barcode)

    if barcode not in deconvolved_spots:
        cell_group = "not deconvolved"
        group["spot_clusters"] = 0
        group["cell_type"] = cell_group

        return group

    if ctype_deconv is None:
        cell_group = "not deconvolved"
        group["spot_clusters"] = 0
        group["cell_type"] = cell_group
        return group

    celltype_ls = []
    for celltype, count in barcode_dict[barcode].items():
        if celltype not in ["cell_n", "distinct_cell_types"]:
            labels = [celltype] * int(count)
            celltype_ls.extend(labels)

    # Assign celltype labels to random cell ids
    np.random.seed(random_state)

    cell_ids = np.array(list(group.cell_ids))
    shuffled_ids = cell_ids.copy()
    np.random.shuffle(shuffled_ids)
    shuffled_ids = list(shuffled_ids)
    celltype_dict = dict(zip(shuffled_ids, celltype_ls))
    group["cell_type"] = group["cell_ids"].map(celltype_dict)
    cl_dict = dict(
        zip(
            list(group["cell_type"].unique()),
            range(len(list(group["cell_type"].unique()))),
        )
    )
    group.insert(
        group.shape[1] - 1,
        "spot_clusters",
        group.replace({"cell_type": cl_dict})["cell_type"].to_list(),
    )
    return group


def assign_random_celltype(
    cell_df_path, reduced_df, barcode_dict: dict, random_state=0
) -> pd.DataFrame:
    """Assign cell_type to individual cells
    Arguments:
        cell_df_path -- path to visium spot information table
        reduced_df {np.ndarray} -- transformed extracted feature table
        barcode_dict {dict} -- dictionary with barcodes and associated cell types
    Returns:
        barcodes {pd.DataFrame} -- modified dataframe with clustered results
    """
    spot_table = pd.read_csv(cell_df_path)
    feature_table = pd.DataFrame(reduced_df)
    feature_table.reset_index(inplace=True)
    feature_table.rename(columns={"index": "cell_ids"}, inplace=True)
    spot_table = spot_table.drop(["cell_x", "cell_y", "cell_polygons"], axis=1)

    if np.min(spot_table.cell_ids) > np.min(
        feature_table.cell_ids
    ):  # TODO: There seems to be a bug here. Need to better see how to handle it
        # TODO: Need to revise if it's just an indexing start issue (0 vs 1) or if it has further implciations
        logger.warn("spot_table and spot cell ids are disaligned.")
        dif = np.min(spot_table.cell_ids) - np.min(feature_table.cell_ids)
        feature_table.cell_ids = feature_table.cell_ids + dif
        feature_table.index = feature_table.index + dif
        logger.info("spot_table, spot cell and feature table have been aligned")

    spot = pd.merge(
        feature_table, spot_table[["cell_ids", "barcode"]], on="cell_ids", how="left"
    )

    logger.info(
        "spot_tablecells_id values - Max:"
        + str(np.max(spot_table.cell_ids))
        + " Min: "
        + str(np.min(spot_table.cell_ids))
    )
    logger.info(
        "spot cells_id values (result from merge) - Max:"
        + str(np.max(spot.cell_ids))
        + " Min: "
        + str(np.min(spot.cell_ids))
    )
    logger.info(
        "Feature table cells_id values- Max:"
        + str(np.max(feature_table.cell_ids))
        + " Min: "
        + str(np.min(feature_table.cell_ids))
    )

    barcodes = spot.groupby(by=["barcode"], as_index=False, group_keys=False).apply(
        random_classifier, barcode_dict=barcode_dict, random_state=random_state
    )

    return barcodes


# Global-based model
def probability_classifier(
    group, spot_probability_pd, cluster_probability_pd, random_state=0
):
    """Assign cell_type to individual cells using simulated annealing
    Arguments:
        group {pd.DataFrame} -- instance of groupby object
        random_state
    Returns:
        group {pd.DataFrame} -- modified group dataframe with clustered results
    """

    # Compute the combined cell type probability for the spot
    barcode = get_barcode(group)
    cell_probabilities = get_combined_probabilities(
        barcode, group, spot_probability_pd, cluster_probability_pd
    )

    # Count the number of cells in the spot
    n_cells = group.deconvolution_results.count()

    # Get the best allocation for the spot using simulated annealing
    # This part is computationally expensive
    if n_cells >= 3:
        new_allocation = get_best_allocation_simanel(
            group,
            cell_probabilities,
            repeats=min(25, n_cells * 2),
            T_start=n_cells * 10,
        )
        group["cell_type"] = new_allocation

    cl_dict = dict(
        zip(
            list(group["cell_type"].unique()),
            range(len(list(group["cell_type"].unique()))),
        )
    )
    group.insert(
        group.shape[1] - 1,
        "spot_clusters",
        group.replace({"cell_type": cl_dict})["cell_type"].to_list(),
    )
    return group


def global_clustering(
    reduced_df: np.ndarray, n_clusters, batch_size, random_state
) -> np.ndarray:
    """Clusterng on the dimensionality reduced feature table
    Arguments:
        reduced_df {np.ndarray} -- Dimensonality reduced feature table
        n_clusters {int} : number of clusters
        batch_size {int} : size to batch process clustering
    Returns:
        np.ndarray -- result of clustering
        np.ndarray -- unique clusters
    """
    model = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=random_state,
        batch_size=batch_size,
        n_init="auto",
    )
    result = model.fit_predict(reduced_df)
    clusters = unique(result)

    return result, clusters


def assign_probabilistic_celltype(
    cell_df_path,
    reduced_df,
    barcode_dict: dict,
    n_clusters=20,
    batch_size=1024,
    random_state=0,
) -> pd.DataFrame:
    """Assign cell_type to individual cells
    Arguments:
        cell_df_path -- path to visium spot information table
        reduced_df {np.ndarray} -- transformed extracted feature table
        barcode_dict {dict} -- dictionary with barcodes and associated cell types
    Returns:
        barcodes {pd.DataFrame} -- modified dataframe with clustered results
    """
    spot_table = pd.read_csv(cell_df_path)
    feature_table = pd.DataFrame(reduced_df)

    # Add morphological clusters
    clusters = global_clustering(
        reduced_df,
        n_clusters=n_clusters,
        batch_size=batch_size,
        random_state=random_state,
    )
    feature_table["cluster_name"] = clusters[0]

    feature_table.reset_index(inplace=True)
    feature_table.rename(columns={"index": "cell_ids"}, inplace=True)
    spot_table = spot_table.drop(["cell_x", "cell_y", "cell_polygons"], axis=1)

    if np.min(spot_table.cell_ids) > np.min(
        feature_table.cell_ids
    ):  # TODO: There seems to be a bug here. Need to better see how to handle it
        # TODO: Need to revise if it's just an indexing start issue (0 vs 1) or if it has further implciations
        logger.warn("spot_table and spot cell ids are disaligned.")
        dif = np.min(spot_table.cell_ids) - np.min(feature_table.cell_ids)
        feature_table.cell_ids = feature_table.cell_ids + dif
        feature_table.index = feature_table.index + dif
        logger.info("spot_table, spot cell and feature table have been aligned")

    spot = pd.merge(
        feature_table, spot_table[["cell_ids", "barcode"]], on="cell_ids", how="left"
    )

    logger.info(
        "spot_tablecells_id values - Max:"
        + str(np.max(spot_table.cell_ids))
        + " Min: "
        + str(np.min(spot_table.cell_ids))
    )
    logger.info(
        "spot cells_id values (result from merge) - Max:"
        + str(np.max(spot.cell_ids))
        + " Min: "
        + str(np.min(spot.cell_ids))
    )
    logger.info(
        "Feature table cells_id values- Max:"
        + str(np.max(feature_table.cell_ids))
        + " Min: "
        + str(np.min(feature_table.cell_ids))
    )

    # Generate random allocation
    random_alloc = (
        spot.drop("cluster_name", axis=1)
        .groupby(by=["barcode"], as_index=False, group_keys=False)
        .apply(random_classifier, barcode_dict=barcode_dict, random_state=random_state)
    )
    random_alloc = random_alloc.rename(columns={"cell_type": "deconvolution_results"})
    clust_alloc = spot.merge(
        random_alloc[["cell_ids", "barcode", "deconvolution_results"]],
        left_on=["cell_ids", "barcode"],
        right_on=["cell_ids", "barcode"],
    )

    # Create column to store new basic reassignment
    clust_alloc["cell_type"] = clust_alloc["deconvolution_results"]

    # Estimate probability distributions
    spot_probability_pd = estimate_spot_probability(
        load_cells_in_spots(clust_alloc)
    ).reset_index()
    cluster_probability_pd = estimate_cluster_probability(
        load_cells_in_spots(clust_alloc)
    ).reset_index()

    logger.info("Starting simulated annealing")
    barcodes = clust_alloc.groupby(
        by=["barcode"], as_index=False, group_keys=False
    ).apply(
        probability_classifier,
        spot_probability_pd,
        cluster_probability_pd,
        random_state=random_state,
    )
    barcodes.drop(["cluster_name", "deconvolution_results"], axis=1, inplace=True)
    return barcodes
