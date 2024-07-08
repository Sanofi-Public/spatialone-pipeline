"""
Cell type assignment confidence metrics
"""

import numpy as np
import pandas as pd
from scipy.stats import entropy
from sklearn.metrics import silhouette_samples


def calculate_sil_score(features, labels):
    """Compute Silhouette score for a particular spot
    based on cell assignment and morph features
    Arguments:
        features {np.array} -- cell morphological features
        labels {np.array} -- sequence of cell type labels
    Returns:
        sil_score -- Silhouette coefficient score for spot
    """
    n_labels = len(np.unique(labels))
    n_cells = len(labels)
    if 1 < n_labels < n_cells:
        sil_score = silhouette_samples(features, labels)
    else:
        sil_score = 1
    return sil_score


def calculate_entropy(labels):
    """Compute Shannon entropy of spot
    Arguments:
        labels {np.array} -- sequence of cell type labels
    Returns:
        shannon_entropy -- Shannon entropy
        normalized_entropy -- normalized Shannon entropy
    """
    # Calculate probability
    values, counts = np.unique(list(labels), return_counts=True)
    probabilities = counts / len(labels)
    # Calculate Shannon entropy
    shannon_entropy = entropy(probabilities)
    # Calculate normalized entropy
    if len(values) == 1:
        normalized_entropy = shannon_entropy
    else:
        # Calculate maximum entropy (ln of the number of unique elements)
        max_entropy = np.log(len(values))
        # normalize
        normalized_entropy = shannon_entropy / max_entropy
    return shannon_entropy, normalized_entropy


def compute_confidence_metrics(group):
    """Compute confidence metrics
    Arguments:
        group {pd.DataFrame} -- instance of groupby object
    Returns:
        group {pd.DataFrame} -- modified group dataframe with confidence metrics
    """
    X = group.iloc[:, 1:-4].to_numpy()  # morphological features
    labels = group["spot_clusters"].to_numpy()  # cell_type labels
    sil_score = calculate_sil_score(X, labels)
    shannon_entropy, normalized_entropy = calculate_entropy(labels)
    group["sil_score"] = sil_score
    group["normalized_entropy"] = normalized_entropy
    group["confidence_score"] = group.apply(
        lambda row: row["sil_score"] * (1 - row["normalized_entropy"]), axis=1
    )
    return group
