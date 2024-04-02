"""
Global assignment
"""

import random
from itertools import permutations

import frigidum
import pandas as pd

# from tqdm.auto import tqdmjob

### Basic functions to load data


def load_spot_cells(spots_pd):
    """Returns a dataframe with the cells contained in each spot"""
    return spots_pd[["barcode", "num_contained_cells", "deconvolution_results"]]


def load_cells_in_spots(cells_pd):
    """Returns cells in spots dataframe with no null barcodes"""
    return cells_pd[cells_pd.barcode.notnull()]


### Basic functions to compute probabilities


def estimate_spot_probability(cells_in_spot_pd):
    """Returns a dataframe with the probability of each cell type in each spot"""

    cells_per_spot_pd = (
        cells_in_spot_pd.groupby(["barcode"], as_index=False)
        .agg({"deconvolution_results": ["count"]})
        .droplevel(1, axis=1)
    )
    cells_per_spot_pd.columns = ["barcode", "ncells"]

    grouped_cells_in_spot_pd = (
        cells_in_spot_pd.groupby(["barcode", "deconvolution_results"], as_index=True)
        .agg({"deconvolution_results": ["count"]})
        .droplevel(1, axis=1)
    )
    grouped_cells_in_spot_pd.columns = ["cell_type"]
    grouped_cells_in_spot_pd = grouped_cells_in_spot_pd.reset_index()

    prob_by_spot_pd = grouped_cells_in_spot_pd.merge(cells_per_spot_pd, how="left")
    prob_by_spot_pd["cell_proportion"] = (
        prob_by_spot_pd.cell_type / prob_by_spot_pd.ncells
    )
    prob_by_spot_pd.columns = [
        "barcode",
        "cell_type",
        "cell_count",
        "total_cells_in_spot",
        "cell_proportion",
    ]

    prob_by_spot_pd = prob_by_spot_pd[["barcode", "cell_type", "cell_proportion"]]
    prob_by_spot_pd = pd.pivot_table(
        prob_by_spot_pd,
        values="cell_proportion",
        index=["barcode"],
        columns=["cell_type"],
        fill_value=0,
    )

    return prob_by_spot_pd  # .reset_index())


def estimate_cluster_probability(cells_pd):
    """Returns a dataframe with the probability of each cell type in each cluster"""

    # Get probabilities at spot level
    prob_cell = estimate_spot_probability(load_cells_in_spots(cells_pd)).reset_index()

    # propagate to clustering gropup
    cluster_probability = cells_pd.merge(prob_cell, on="barcode")
    celltype_index = list(cluster_probability.columns).index("deconvolution_results")
    celltypes = list(cluster_probability.columns)[celltype_index + 1 :]
    prob_cluster = cluster_probability.groupby("cluster_name")[celltypes].sum()
    prob_cluster = prob_cluster.div(prob_cluster.sum(axis=1), axis=0)

    return prob_cluster


def get_combined_probabilities(
    barcode, list_of_cells, spot_probability_pd, cluster_probability_pd
):
    """Returns a dataframe with the combined probabilities of each cell type in each cluster
    based on the cluster probability and the deconvolution probability of the spot"""

    spot_cells_pd = list_of_cells[list_of_cells.barcode == barcode]

    spot_prob = spot_probability_pd[spot_probability_pd.barcode == barcode]
    spot_cluster_prob = cluster_probability_pd[
        cluster_probability_pd.cluster_name.isin(spot_cells_pd.cluster_name.unique())
    ]

    celltypes_of_interest = [
        spot_prob.columns[i]
        for i in range(1, len(spot_prob.columns))
        if (spot_prob.iloc[0][spot_prob.columns[i]] > 0)
    ]

    spot_cluster_prob = spot_cluster_prob[["cluster_name"] + celltypes_of_interest]
    spot_cluster_prob_aux = spot_cluster_prob[celltypes_of_interest]
    spot_cluster_prob_aux = spot_cluster_prob_aux.div(
        spot_cluster_prob_aux.sum(axis=1), axis=0
    )
    spot_cluster_prob[celltypes_of_interest] = spot_cluster_prob_aux

    # Compute combined probabilities

    rows = []
    for i in range(0, len(spot_cluster_prob)):
        rows.append(
            dict(
                spot_cluster_prob.iloc[i][celltypes_of_interest]
                * spot_prob.iloc[0][celltypes_of_interest]
            )
        )
    combined_probabilities = pd.DataFrame(rows)
    combined_probabilities = combined_probabilities.div(
        combined_probabilities.sum(axis=1), axis=0
    )
    combined_probabilities.insert(
        0, "cluster_name", list(spot_cluster_prob["cluster_name"])
    )

    return combined_probabilities


### Basic functions for computing the best allocation of cells to clusters using brute force


def compute_allocation_score(allocation, spotcells, cell_probabilities):
    """Computes the score of an allocation of cells to clusters based on the probability of each cell type in each cluster"""
    score = 0
    for i in range(0, len(allocation)):
        cluster_type = spotcells.iloc[i].cluster_name
        score = (
            score
            + cell_probabilities[cell_probabilities.cluster_name == cluster_type][
                allocation[i]
            ].iloc[0]
        )
    return score


def get_best_allocation(spotcells, cell_probabilities):
    """Returns the best allocation of cells to clusters based on the probability of each cell type in each cluster using brute force"""
    allocation_permutations = list(set(permutations(spotcells.deconvolution_results)))
    best_allocation = []
    best_score = -1
    for allocation in allocation_permutations:
        score = compute_allocation_score(allocation, spotcells, cell_probabilities)
        if score > best_score:
            best_score = score
            best_allocation = allocation

    return best_allocation


### Functions needed to run simulated annealing


class alloc_obj:
    """Object to store the allocation, the cell probabilities and the spotcells. This will be a state in simulated annealing"""

    allocation = []
    cell_probabilities = None
    spotcells = None

    def __init__(self, allocation=[], cell_probabilities=None, spotcells=None):
        """Constructor for the allocation object"""
        self.allocation = allocation
        self.cell_probabilities = cell_probabilities
        self.spotcells = spotcells

    def copy(alloc_obj_a):
        """Returns a copy of the allocation object"""
        return alloc_obj(
            alloc_obj_a.allocation.copy(),
            alloc_obj_a.cell_probabilities,
            alloc_obj_a.spotcells,
        )


### Neighbour functions for simulated annealing


def random_small_step(alloc):
    """Returns a new allocation object with a random swap of two cells"""
    swaps = random.sample(range(len(alloc.allocation)), 2)
    alloc.allocation[swaps[0]], alloc.allocation[swaps[1]] = (
        alloc.allocation[swaps[1]],
        alloc.allocation[swaps[0]],
    )
    return alloc


def random_big_step(alloc):
    """Returns a new allocation object with two random swaps of two cells"""
    swaps = random.sample(range(len(alloc.allocation)), 2)
    alloc.allocation[swaps[0]], alloc.allocation[swaps[1]] = (
        alloc.allocation[swaps[1]],
        alloc.allocation[swaps[0]],
    )
    swaps = random.sample(range(len(alloc.allocation)), 2)
    alloc.allocation[swaps[0]], alloc.allocation[swaps[1]] = (
        alloc.allocation[swaps[1]],
        alloc.allocation[swaps[0]],
    )
    return alloc


def invert_step(alloc):
    """Returns a new allocation object with the allocation reversed"""
    alloc.allocation.reverse()
    return alloc


def objective_min_score(alloc):
    """Returns the score of the allocation object.
    The score is the number of cells minus the score of the allocation.
    This is done to maximize the score of the allocation"""
    return len(alloc.allocation) - compute_allocation_score(
        alloc.allocation, alloc.spotcells, alloc.cell_probabilities
    )


def get_best_allocation_simanel(
    spotcells, cell_probabilities, T_start=1000, T_stop=0.05, repeats=50
):
    """Returns the best allocation of cells types based on the probability of each cell type in each cluster using simulated annealing"""

    spot_allocation = alloc_obj()
    spot_allocation.spotcells = spotcells
    spot_allocation.allocation = list(spotcells.deconvolution_results)
    spot_allocation.cell_probabilities = cell_probabilities

    # Temporal funciton used to initialize the state
    def init_state():
        return spot_allocation

    local_opt = frigidum.sa(
        random_start=init_state,
        neighbours=[random_small_step, random_big_step, invert_step],
        objective_function=objective_min_score,
        T_start=T_start,
        T_stop=T_stop,
        repeats=repeats,
        copy_state=frigidum.annealing.copy,
    )

    return local_opt[0].allocation
