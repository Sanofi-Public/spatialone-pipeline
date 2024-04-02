"""
Python class containing methods to perform comparative analysis.
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

from src.utils.logger import Logger

logger = Logger()


class ComparativeAnalysis:
    """Python class containing methods to perform comparative analysis."""

    def __init__(self, configs=None):
        """Initializing ComparativeAnalysis instance variables.

        Args:
            configs (dict, optional): Dictionary containing SpatialAnalysis configs. Defaults to None.
        """
        self.configs = configs

    def compare_regions_fc(self, spots1_df, spots2_df, genes_list, level="tissue"):
        """Compare the gene expression in two regions of the tissue.
        Computes the changes in the gene expression of two regions of tissue
        and perform statistical tests to check if the differences are significant.
        Statistical tests performed: two sample t-test, wilcoxon test and mann-whitney test.
        Correction for multiple testing: Benjamini-Hochberg correction.

        Args:
            spots1_df (pd.DataFrame): Dataframe with the spots of the first region
            spots2_df (pd.DataFrame): List of barcodes of the cells in the second region
            genes_list (list): List of genes to compare

        Returns:
            results_df (pd.DataFrame): Dataframe with the results of the comparative analyses.
        """
        # define dataframe to store results
        results_df = []
        for gene in genes_list:
            # calculate mean expression levels for each region
            mean1 = spots1_df[gene].mean()
            mean2 = spots2_df[gene].mean()

            # calculate fold change
            if mean2 == 0:
                fc = np.inf
            elif mean1 > 0:
                fc = mean1 / mean2
            elif mean2 > 0:
                fc = np.inf
            else:
                fc = np.nan
            # perform two-sample t-test
            ttest = stats.ttest_ind(spots1_df[gene], spots2_df[gene])

            # perform Mann-Whitney U test
            mannwhitney = stats.mannwhitneyu(spots1_df[gene], spots2_df[gene])

            # add results to dataframe
            results_df.append(
                {
                    "Gene": gene,
                    "Mean1": mean1,
                    "Mean2": mean2,
                    "FC": fc,
                    "T-test stat": ttest.statistic,
                    "T-test pvalue": ttest.pvalue,
                    "Mann-Whitney stat": mannwhitney.statistic,
                    "Mann-Whitney pvalue": mannwhitney.pvalue,
                }
            )
        # Converting to a Pandas dataframe
        results_df = pd.DataFrame(results_df)
        # apply Benjamini-Hochberg correction to p-values
        corr_t_label = "Corrected T-test pvalue"
        corr_mw_label = "Corrected MW pvalue"
        significance_label = "Is Significant"
        results_df[corr_t_label] = multipletests(
            results_df["T-test pvalue"].fillna(1), method="fdr_bh"
        )[1]
        results_df[corr_mw_label] = multipletests(
            results_df["Mann-Whitney pvalue"].fillna(1), method="fdr_bh"
        )[1]

        # check if gene is significant based on corrected p-values and fold change
        pvalue = self.configs["params"][level]["diff_pval"]
        fc_threshold = self.configs["params"][level]["diff_fc_thresh"]
        results_df[significance_label] = (results_df[corr_t_label] < pvalue) & (
            results_df[corr_mw_label] < pvalue
        )
        results_df["Direction"] = "not significant"
        results_df.loc[
            (results_df.FC >= fc_threshold) & results_df[significance_label],
            "Direction",
        ] = "upregulated"
        results_df.loc[
            (results_df.FC <= 1 / fc_threshold) & results_df[significance_label],
            "Direction",
        ] = "downregulated"
        #'upregulated' if results_df["FC"] >= fc_threshold else 'downregulated' if results_df["FC"] <= 1/fc_threshold else 'not significant'
        return results_df

    def compare_gene_exp(self, regions_dict, genes_list, level="tissue"):
        """Compares gene expression levels for the regions defined in regions_dict.

        Args:
            regions_dict (dict): Dictionary of regions. See SpatialAnalysis.load_region_level_data()
             for details.
            genes_list (list): gene types to compare.

        Returns:
            pandas.DataFrame: Dataframe with the gene expression differences and fold changes for
        """
        logger.info("<compare_gene_exp> Comparing gene expression in tissue regions")
        if not self.configs["params"][level]["spot_diff_exp"] or (
            len(regions_dict) == 0
        ):
            return None
        region_data = {}
        for region_type, sub_regions in regions_dict.items():
            sub_spots_df = pd.concat(
                [val["spots_df"] for _, val in sub_regions.items()]
            )
            region_data[region_type] = sub_spots_df
        gene_diff_df = pd.DataFrame()
        for region1, region1_df in region_data.items():
            for region2, region2_df in region_data.items():
                if len(region1_df) < 2 or len(region2_df) < 2:
                    continue
                results_df = self.compare_regions_fc(
                    region1_df, region2_df, genes_list, level=level
                )
                results_df["Region A"] = region1
                results_df["Region B"] = region2
                # append results_df to gene_diff_df
                gene_diff_df = pd.concat([gene_diff_df, results_df])
        # perform log2 of column FC
        if not gene_diff_df.empty:
            gene_diff_df["log2FC"] = np.log2(gene_diff_df.FC)
            gene_diff_df["-log10pvalue"] = -np.log10(
                gene_diff_df["Corrected T-test pvalue"]
            )
        return gene_diff_df

    def test_infiltration(
        self, counts_df, cell_type, pval_thresh=0.05, alt_hyp="two-sided", buffer=False
    ):
        """Runs a Z-test to compare cell proportions inside vs outside a region

        Args:
            counts_df (pandas.DataFrame): Pandas dataframe with the cell types and their abundance at each level
            cell_type (str): cell type to study
            pval_thresh (float): p-value threshold
            alt_hyp (str, optional): Alternative hypothesis [‘two-sided’, ‘smaller’, ‘larger’]. Defaults to "larger".

        Returns:
            float: (stat) Z-score
            float: (pval) p-value
            bool: (sig_difference) Significant difference or not
        """
        df = counts_df[counts_df.cell_type == cell_type].reset_index(drop=True)
        if buffer:
            if len(df) < 2:
                return 0, 1, False
            buffer_level = df.level[0] - df.level[1]
        else:
            buffer_level = 0
        num_target_cells_within = df[df.level <= -buffer_level].abundance.sum()
        num_target_cells_outside = df[df.level > buffer_level].abundance.sum()
        num_cells_within = counts_df[counts_df.level <= -buffer_level].abundance.sum()
        num_cells_outside = counts_df[counts_df.level > buffer_level].abundance.sum()
        count = np.array([num_target_cells_within, num_target_cells_outside])
        nobs = np.array([num_cells_within, num_cells_outside])
        if np.any(nobs == 0):
            return 0, 1, False
        # null hypothesis: p(within) == p(outside)
        # alternative hypothesis: p(within) > p(outside)
        stat, pval = proportions_ztest(count, nobs, value=0, alternative=alt_hyp)
        sig_difference = pval < pval_thresh
        return stat, pval, sig_difference
