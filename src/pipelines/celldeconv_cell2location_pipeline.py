"""
inference pipeline for cell deconvolution
model: cell2location
"""
# pylint:disable=no-member

import os
from datetime import datetime

import anndata
import cell2location
import hdf5plugin
import hydra
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes
from dotenv import load_dotenv
from omegaconf import OmegaConf

from src.pipelines.dataio_pipeline import DataIO
from src.utils.logger import Logger

load_dotenv()
logger = Logger()


class Cell2Location:
    """
    Python class for cell deconvolution using cell2location algorithm
    """

    def __init__(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
        hydra.initialize(config_path="../../conf")
        self.pipeline_name = "celldeconv"

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

        if "params" in self.configs and "seed" in self.configs["params"]:
            scivi_seed = self.configs["params"]["seed"]
        else:
            scivi_seed = 2023
            logger.warning(
                f"Seed not defined in parameters file. {scivi_seed} will be used"
            )
        scvi.settings.seed = scivi_seed
        logger.info(f"scivi.seed set to {scivi_seed}")

    def load_data(self, prep_dir, results_dir, reference_dir, mode="prod"):
        ### read cell segment file
        logger.info("load tissue position table")
        tissue_position = pd.read_csv(f"{results_dir + 'spots_df.csv'}")
        if "barcode" in tissue_position.columns:
            tissue_position.set_index("barcode", inplace=True, drop=False)
        # ### filter in tissue and number of cell over 0
        self.tissue_position = tissue_position.query(
            "in_tissue == 1 and num_contained_cells > 0"
        )
        ### read st data
        logger.info("load st data")
        self.adata_st = sc.read_10x_h5(f"{prep_dir + 'filtered_feature_bc_matrix.h5'}")
        self.adata_st.var_names_make_unique()
        self.adata_st.var["SYMBOL"] = self.adata_st.var_names
        ## change parameter for unit test use only
        if mode == "unit test":
            self.configs["params"]["sc_use_gpu"] = False
            self.configs["params"]["sc_max_epoches"] = 10
            self.configs["params"]["st_max_epoches"] = 100
        logger.info(f'mode is {mode}, use gpu {self.configs["params"]["sc_use_gpu"]}')

    def train_sig(self, reference_dir, mode="prod"):
        adata_st = sc.read_h5ad(
            reference_dir + self.configs["params"]["atlas_type"] + "_cell_atlas.h5ad"
        )
        adata_st.var_names_make_unique()
        adata_ref = adata_st.copy()
        logger.info(f"reference sc data shape: {str(adata_st.shape)}")
        ### use untransformed and unnormalised count matrix

        selected = filter_genes(
            adata_ref,
            cell_count_cutoff=self.configs["params"][
                "sc_cell_count_cutoff"
            ],  ### count > 0 in over 20 cells
            cell_percentage_cutoff2=self.configs["params"][
                "sc_cell_percentage_cutoff2"
            ],  ### count > 0 in many cells >5%
            nonz_mean_cutoff=self.configs["params"]["sc_nonz_mean_cutoff"],
        )  ## mean expression across non-zero cells slightly larger than 1 (> 1.12)
        logger.info(f"filter genes for sc data, now gene length is {len(selected)}")
        # filter the object
        adata_ref = adata_ref[:, selected].copy()
        # prepare anndata for the regression model
        logger.info("Start training cell signature")
        cell2location.models.RegressionModel.setup_anndata(
            adata=adata_ref,
            # 10X reaction / sample / batch
            batch_key=self.configs["params"]["sc_batch_key"],
            # cell type, covariate used for constructing signatures
            labels_key=self.configs["params"]["sc_label_key"],
            # multiplicative technical effects (platform, 3' vs 5', donor effect)
            categorical_covariate_keys=self.configs["params"][
                "sc_categorical_covariate_keys"
            ],
        )
        # create the regression model
        mod = RegressionModel(adata_ref)
        # view anndata_setup as a sanity check
        mod.view_anndata_setup()
        if adata_ref.shape[0] > 10000:
            batch_size = 2500
        else:
            batch_size = int(np.floor(sc_exp.shape[0] / 2))
        mod.train(
            max_epochs=self.configs["params"]["sc_max_epoches"],
            use_gpu=self.configs["params"]["sc_use_gpu"],
            batch_size=batch_size,
            lr=self.configs["params"]["sc_lr"],
        )
        return mod, adata_ref

    def save_train_sig(self, reference_dir, mode="prod"):
        # In this section, we export the estimated cell abundance (summary of the posterior distribution).
        logger.info(f"mode is: {mode}")
        mod, adata_ref = self.train_sig(reference_dir, mode)
        adata_ref = mod.export_posterior(
            adata_ref,
            sample_kwargs={
                "num_samples": 1000,
                "batch_size": 2500,
                "use_gpu": self.configs["params"]["sc_use_gpu"],
            },
        )
        adata_file = (
            f'{reference_dir + self.configs["params"]["atlas_type"] + "_cell_sig.h5ad"}'
        )
        adata_ref.write_h5ad(adata_file, compression=hdf5plugin.FILTERS["zstd"])
        logger.info(f"reference signatures model save to: {adata_file}")
        return adata_ref

    def train_location_model(self, cell_sig):
        # extract estimated cell abundance
        adata_ref = cell_sig.copy()
        adata_ref.var_names_make_unique()
        adata_vis = self.adata_st.copy()
        # export estimated expression in each cluster
        if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
            inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
                [
                    f"means_per_cluster_mu_fg_{i}"
                    for i in adata_ref.uns["mod"]["factor_names"]
                ]
            ].copy()
        else:
            inf_aver = adata_ref.var[
                [
                    f"means_per_cluster_mu_fg_{i}"
                    for i in adata_ref.uns["mod"]["factor_names"]
                ]
            ].copy()
        inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
        # find shared genes and subset both anndata and reference signatures
        if np.intersect1d(adata_vis.var_names, inf_aver.index).shape[0] == 0:
            ## adjust index for spatial when reference index is not match up with spatial ones
            adata_vis.var.set_index("gene_ids", drop=True, inplace=True)
            logger.info(
                f"adjust index for spatial, intersect gene number is: {np.intersect1d(adata_vis.var_names, inf_aver.index).shape[0]}"
            )
            logger.info(
                f"spatial data gene format: {adata_vis.var_names[0]}, sc data gene num: {inf_aver.index[0]}"
            )
        intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
        adata_vis = adata_vis[:, intersect]
        adata_vis.var_names = intersect
        inf_aver = inf_aver.loc[intersect, :].copy()
        logger.info(
            f"spatial data gene num: {adata_vis.var_names.nunique()}, sc data gene num: {inf_aver.index.nunique()}"
        )
        logger.info(
            f"spatial data gene format: {adata_vis.var_names[0]}, sc data gene num: {inf_aver.index[0]}"
        )

        # training for stRNA-seq
        logger.info(f"Start mapping spatial data to get cell aboundance")
        adata_st = self.adata_st.copy()
        cell2location.models.Cell2location.setup_anndata(adata_vis)
        mod_st = cell2location.models.Cell2location(
            adata_vis,
            cell_state_df=inf_aver,
            N_cells_per_location=int(self.configs["params"]["st_N_cells_per_location"]),
            detection_alpha=self.configs["params"]["st_detection_alpha"],
        )

        cell2location.models.Cell2location.train(
            mod_st,
            max_epochs=int(self.configs["params"]["st_max_epoches"]),
            train_size=1.0,
            use_gpu=self.configs["params"]["sc_use_gpu"],
        )
        parameters_setting_st = {
            "batch_size": adata_vis.shape[1],
            "use_gpu": self.configs["params"]["sc_use_gpu"],
        }
        self.adata_st = cell2location.models.Cell2location.export_posterior(
            mod_st, adata_vis, sample_kwargs=parameters_setting_st
        )

    def get_cell_type_proportion(self):
        # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
        # to adata.obs with nice names for plotting
        adata_vis = self.adata_st.copy()
        logger.info(f"representing confident cell abundance")
        adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = adata_vis.obsm[
            "q05_cell_abundance_w_sf"
        ]
        ### healthy cell type aboundance
        # cell_type = cell_type.iloc[:, 7:-2]
        cell_type = adata_vis.obs[adata_vis.uns["mod"]["factor_names"]].copy()
        cell_type_stad = cell_type.copy()
        logger.info(
            f"representing determined cell abundance by set cell aboundance threshold "
        )
        ### cell aboundance threshold set at 0.1, cell abundance >0.1 was used as a gold standard label
        threshold = self.configs["params"]["cell_aboundance_threshold"]
        cell_type_stad[cell_type_stad < threshold] = 0
        cell_type_stad.columns = cell_type.columns
        ### measure cell type aboundace proportion per spot
        logger.info(f"{cell_type_stad.shape}, representing refined cell proportion")
        cell_type_matrix = np.apply_along_axis(
            lambda x: x / np.sum(x) if np.sum(x) > 0 else np.zeros(x.shape),
            1,
            cell_type_stad,  # fix to return an array of 0s
        )
        cell_type_stad = pd.DataFrame(
            cell_type_matrix, index=cell_type.index, columns=cell_type.columns
        )
        #### healthy cell proportion = original proporition * (1- tumor rate), which means healthy cell + tumor proportion equal to 1
        adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = cell_type_stad[
            cell_type.columns
        ]
        return adata_vis

    def props_to_counts(self, props, cell_num):
        bins = len(props)  # number of bins
        dist = np.zeros(bins)  # initial distribution filled with zeros
        proportion = props.copy()  # create a copy
        i = 0
        while i < int(cell_num):
            total = max(np.sum(dist), 1)  # how many items total
            prop = dist / total  # current distribution
            error = proportion - prop  # errors with the desired distribution
            idx = np.argmax(error)  # which bin has a max error
            dist[idx] = dist[idx] + 1  # increment it
            i = i + 1
        return dist

    def save_data(self, adata_st, results_dir):
        tissue_clone_label_pos = self.tissue_position
        assert (
            len(adata_st.obs_names & tissue_clone_label_pos.index) > 0
        ), "spatial coor index is not match with tissue position index"
        ### map cell count per spot
        logger.info(
            f"representing determined cell count by take total cell count per spot into consideration"
        )
        adata_st.obs.loc[
            adata_st.obs_names & tissue_clone_label_pos.index, "cell_n"
        ] = tissue_clone_label_pos.loc[
            adata_st.obs_names & tissue_clone_label_pos.index, "num_contained_cells"
        ].values
        adata_st.obs["cell_n"].fillna(0, inplace=True)
        ## fill cell aboundance of un-deconvplved cell type as 0
        adata_st.obs[adata_st.uns["mod"]["factor_names"]].fillna(0, inplace=True)
        ### get healthy cell count
        cell_list = list(adata_st.uns["mod"]["factor_names"])
        ### get healthy cell count by increment 1 cell count to max props repeartly
        cell_props = adata_st.obs[cell_list].copy()
        cell_num = adata_st.obs["cell_n"].copy()
        cell_type_count = pd.DataFrame(
            columns=dict.fromkeys(adata_st.uns["mod"]["factor_names"])
        )
        for row in range(adata_st.shape[0]):
            round_cell = self.props_to_counts(
                cell_props.loc[cell_num.index[row]].values, cell_num[row]
            )
            row_dict = pd.DataFrame(
                round_cell.reshape(1, -1),
                index=[cell_num.index[row]],
                columns=cell_list,
            )
            cell_type_count = pd.concat([cell_type_count, row_dict])
        self.cell_type_count = cell_type_count
        self.cell_num = cell_num
        ### save result to cache
        cell_type_count.reset_index(names="barcode", inplace=True)
        cell_type_count.to_csv(
            f"{results_dir + 'cell2location_results.csv'}", index=False
        )
        logger.info(
            f"save cell2location result to {results_dir + 'cell2location_results.csv'}"
        )
        self.cell2location_cell_count_result = adata_st


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
    cell_deconv = Cell2Location()
    cell_deconv.load_model_configs_from_flow(
        pipeline_configs_from_flow=config_flow[cell_deconv.pipeline_name]
    )

    # logger.info(f"[{cell_deconv.pipeline_name}] Fetch input files")
    # data_io.fetch_input_files(pipeline_configs=cell_deconv.configs)

    # Check if trained <atlas_type>_cell_signature exists in S3
    # search_result_cell_sig = data_io.check_if_file_exists(
    #     "prep",
    #     source_base_path=data_io.reference_dir,
    #     parent_folder="cell_sig",
    #     file_name_ext=cell_deconv.configs["params"]["atlas_type"] + "_cell_sig.h5ad",
    # )
    # If it cell_sig exists and there is no retraining of cell_sig requried, then load cell_sig from S3

    if cell_deconv.configs["params"]["retrain_cell_sig"] == False:
        logger.info(
            f"[{cell_deconv.pipeline_name}] Fetch trained cell signature data from reference directory"
        )
        # data_io.fetch_data(
        #     cell_deconv.configs,
        #     data_io.reference_path,
        #     fetch_parent_folder="cell_sig",
        #     target_base_path=data_io.cache_dir,
        #     file_name_prefix=cell_deconv.configs["params"]["atlas_type"],
        # )  ## fetch sc data
        cell_sig = sc.read_h5ad(
            data_io.reference_dir
            + cell_deconv.configs["params"]["atlas_type"]
            + "_cell_sig.h5ad",
        )

    else:
        logger.info(
            f"[{cell_deconv.pipeline_name}] Fetch single cell atlas data from s3 for training cell signature"
        )
        cell_sig = cell_deconv.save_train_sig(
            results_dir=data_io.reference_dir
        )  ### train sc data and save to reference directory

    logger.info(f"[{cell_deconv.pipeline_name}] Start pipeline")
    logger.info(f"[{data_io.results_dir}] cache direction")
    cell_deconv.load_data(
        prep_dir=data_io.prep_dir,
        results_dir=data_io.results_dir,
        reference_dir=data_io.reference_dir,
    )  ### load st data

    cell_deconv.train_location_model(cell_sig)
    adata_st = cell_deconv.get_cell_type_proportion()
    cell_deconv.save_data(adata_st, results_dir=data_io.results_dir)
