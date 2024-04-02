"""
transform functions for clustering
"""

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def scaling(master_table: pd.DataFrame) -> pd.DataFrame:
    """Dataframe with scaled features
    Arguments:
        master_table {pd.DataFrame} -- Dataframe of features
    Returns:
        pd.DataFrame -- Dataframe with scaled features
    """
    training_data = master_table.to_numpy()
    scaler = StandardScaler()
    scaled_training_data = scaler.fit_transform(training_data)
    return scaled_training_data


def dimensionality_reduction(scaled_training_data: pd.DataFrame):
    """Dataframe with dimensionality reduced features
    Arguments:
        scaled_training_data {pd.DataFrame} -- Dataframe with scaled features
    Returns:
        sklearn.decomposition._pca.PCA -- PCA model
        np.ndarray -- Dimensonality reduced feature table
    """
    pca = PCA(10, random_state=0)
    reduced_df = pca.fit_transform(scaled_training_data)
    # add path to save reduced_df
    return pca, reduced_df
