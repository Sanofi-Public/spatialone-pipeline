"""
Python class containing methods to generate plots for the spatial structure analysis report.
"""

import itertools
import os

import dash_bio as dashbio
import holoviews as hv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
import squidpy as sq
from distinctipy import distinctipy
from holoviews import dim, opts
from plotly.subplots import make_subplots
from sklearn.preprocessing import LabelEncoder

from src.spatial_report_generator.comparative_analysis import ComparativeAnalysis
from src.spatial_report_generator.report_generator import ReportGenerator
from src.utils.logger import Logger

logger = Logger()


class PlotGenerator:
    """
    Python class containing methods to generate plots for the spatial structure analysis report.
    """

    def __init__(
        self,
        cells_adata,
        spots_adata,
        configs=None,
        level="tissue",
        cache_path="",
        report_template_path="",
        scale_factors=None,
    ):
        """Initializing PlotGenerator instance variables.

        Args:
            cells_adata (anndata.AnnData): Cells AnnData object.
            spots_adata (anndata.AnnData): Spots AnnData object.
            configs (dict, optional): Dictionary containing SpatialAnalysis configs. Defaults to None.
        """
        self.pipeline_name = "PlotGenerator"
        self.cells_adata = cells_adata
        self.spots_adata = spots_adata
        self.configs = configs
        self.level = level
        self.cache_path = cache_path
        self.scale_factors = scale_factors
        self.params = self.configs["params"][level]
        self.width = self.params["plot_width"]
        self.height = self.params["plot_height"]
        self.fig_configs = self.configs["params"]["fig_configs"]
        self.inline_figs = []
        self.report_gen = ReportGenerator(report_template_path)

    def convert_micron(self, pix_distance, spot_size=65.0):
        """
        Converts the pixel distance to microns using the fact that spot_diameter_fullres in
        scale_factors refers to the number of pixels that span the diameter of a theoretical
        65µm spot in the original, full-resolution image.

        Returns:
            distance (float): distance in microns
        """
        distance = pix_distance * (
            spot_size / self.scale_factors["spot_diameter_fullres"]
        )
        return distance

    @classmethod
    def generate_distinct_colors(cls, n):
        """
        Generates distinct colors (hex format) for the different cell types.

        Args:
            n (int): An integer describing the number of distinct colours to generate.

        Returns:
            colors_hex (list): List of n hex colors
        """
        colors = distinctipy.get_colors(n)
        colors_hex = [
            "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
            for r, g, b in colors
        ]
        return colors_hex

    def generate_gene_distrib_collage(
        self, save_location, adata, markers, n, n_rows, n_cols
    ):
        """Generates a collage of the distribution plots for the markers specified

        Args:
            save_location (str): plot temporary save location
            adata (AnnData): spots anndata object containing the markers
            markers (list): list of gene markers to plot
            n_rows (int): number of genes to show
            n_rows (int): number of rows
            n_cols (int): number of columns

        Returns:
            matplotlib figure: figure containing the gene distruibution plots
        """
        adata.obsm["spatial"] = np.array(adata.obsm["spatial"])
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(self.width // 100, self.height // 100)
        )
        axes = axes.reshape(n_rows, n_cols)
        i = 0
        for row in range(n_rows):
            for col in range(n_cols):
                if i == n:
                    break
                ax = axes[row, col]
                marker = markers[i]
                save_path = save_location.format(marker=marker, i=i)
                save_path = os.path.join(os.getcwd(), save_path)
                sq.pl.spatial_scatter(adata, color=[marker], save=save_path)
                gene_exp_plot = plt.imread(save_path)
                ax.imshow(gene_exp_plot)
                ax.set_axis_off()
                ax.set_aspect("equal")
                os.remove(save_path)
                i += 1
        fig.tight_layout(pad=0.1, h_pad=0.1, w_pad=0.1)
        return fig

    def add_qc_report_button(self, exp_id, flow_id, task_name="qc_report"):
        """Adds button to link to QC report from space ranger

        Args:
            exp_id (str): Experiment name
        """
        if self.params[task_name]:
            configs = dict(self.fig_configs[task_name])
            button_html = f"""
                <button onclick=\"window.open('{exp_id}_web_summary.html?path={exp_id}|{flow_id}', '_blank')\">
                SpaceRanger Report
                </button>"""
            configs["figure"] = self.report_gen.format_to_html(button_html)
            self.inline_figs.append(configs)

    def cell_summary_table(self, df, task_name="cell_summary"):
        """Creates a cell summary table from the provided statistics dataframe.

        Args:
            df (pandas.DataFrame): spot statistics dataframe (SpatialAnalysis.spot_stats)

        Returns:
            plotly.graph_objects.Figure: Table of cell summary statistics
        """
        logger.info("<cell_summary_table> Generating cell summary table")
        if self.params[task_name]:
            cell_summary_table = go.Figure(
                data=[
                    go.Table(
                        header=dict(
                            values=list(df.columns),
                            fill_color="paleturquoise",
                            align="left",
                        ),
                        cells=dict(
                            values=[df[col] for col in df.columns],
                            fill_color="lavender",
                            align="left",
                        ),
                    )
                ]
            )
            cell_summary_table.update_layout(width=self.width)
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(cell_summary_table)
            configs["data"] = df
            configs["name"] = task_name
            self.inline_figs.append(configs)
            return cell_summary_table

    def cell_count_plot(
        self,
        cells_df,
        title="Total number of deconvolved cells",
        ignore_not_deconvolved=True,
        logscale=False,
        task_name="cell_counts",
    ):
        """
        Creates a barplot of the number of cells per cell type using plotly.

        Args:
            cells_df (pandas.DataFrame): A DataFrame containing the information of each cell.
            title (str): Title of the plot.
            ignore_not_deconvolved (bool): If True, the cells that were not deconvolved will be ignored.
            logscale (bool): If True, the y-axis will be in log scale.

        Returns:
            plotly.graph_objects.Figure: A barplot of the number of cells per cell type.
        """
        if self.params[task_name]:
            logger.info("<cell_count_plot> Generating cell count figure")
            if ignore_not_deconvolved:
                cells_df = cells_df[cells_df.cell_type != "undefined"]
                cells_df = cells_df[cells_df.cell_type != "Undefined"]
                cells_df = cells_df[cells_df.cell_type != "not deconvolved"]

            cells_df_grouped = (
                cells_df.groupby("cell_type")
                .count()["barcode"]
                .sort_values(ascending=False)
            )

            # Create a plotly barplot of the number of cells per cell type
            fig = px.bar(
                cells_df_grouped,
                x=cells_df_grouped.index,
                y=cells_df_grouped.values,
                color=cells_df_grouped.index,
                log_y=logscale,
                height=self.height,
                width=self.width,
            )
            # Add title
            fig.update_layout(
                title=f"{title}: {len(cells_df)}",
                xaxis=dict(title="Cell type"),
                yaxis=dict(title=f"Cell count"),
            )
            cell_count_fig = fig
            cell_count_fig.update_xaxes(tickangle=45)
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(cell_count_fig)
            self.inline_figs.append(configs)
            return cell_count_fig

    def cells_within_spot_plot(self, spot_summary, task_name="cell_per_spot"):
        """
        Creates a violin plot of the number of cells per cell type in each spot using plotly.

        Args:
            spot_stats (pandas.DataFrame): A DataFrame containing stats on the cells contained within spots.
                Created by SpatialAnalysis.summarize_cells_within_spots
            spot_cells (pandas.DataFrame): A DataFrame containing the cell types contained within each spot.
                Created by SpatialAnalysis.summarize_cells_within_spots

        Returns:
            plotly.graph_objects.Figure: A box plot of the number of cells within spots.
        """
        if self.params[task_name]:
            logger.info(
                "<cells_within_spot_plot> Generating cell counts within spot figure"
            )
            spot_stats, spot_cells = spot_summary
            cell_order = spot_stats.sort_values(
                "cell_count", ascending=False
            ).cell_type.tolist()
            spot_cells_flat = spot_cells.melt(value_vars=cell_order)
            fig = px.box(
                spot_cells_flat,
                x="cell_type",
                y="value",
                height=self.height,
                width=self.width,
                points=False,
            )
            fig.update_xaxes(
                tickangle=45, categoryorder="array", categoryarray=cell_order
            )
            fig.update_layout(
                xaxis_title="Cell type", yaxis_title="No. of cells within a spot"
            )
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            self.inline_figs.append(configs)
            return fig

    def avg_gene_exp_plot(
        self, spots_df, logscale=False, task_name="spot_avg_gene_counts"
    ):
        """Generates average gene expression plots

        Args:
            spots_df (_type_): _description_
            logscale (bool, optional): _description_. Defaults to False.
            task_name (str, optional): _description_. Defaults to "spot_avg_gene_counts".

        Returns:
            _type_: _description_
        """
        adata = self.spots_adata
        if self.params[task_name]:
            logger.info("<avg_gene_exp_plot> Generating average gene expression plots")
            genes = list(adata.uns["gene_columns"])
            gene_data = spots_df[genes]
            n = np.minimum(self.params["n_genes_bar"], len(genes))
            avg_gene_exp = gene_data.mean(axis=0).sort_values(ascending=False).head(n)
            # Create a plotly barplot of the number of cells per cell type
            fig = px.bar(
                avg_gene_exp,
                x=avg_gene_exp.index,
                y=avg_gene_exp.values,
                log_y=logscale,
                height=self.height,
                width=self.width,
            )
            y_label = "Log Avg. Gene Count" if logscale else "Avg. Gene Count"
            fig.update_xaxes(tickangle=45)
            fig.update_layout(xaxis_title="Gene", yaxis_title=y_label)
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            configs["data"] = pd.DataFrame(avg_gene_exp)
            configs["name"] = task_name
            self.inline_figs.append(configs)
            return fig

    def neighborhood_enrichment_plots(self, task_name="cell_net_matrix"):
        """
        Plots neighborhood enrichment analysis.

        Args:
            adata (anndata.AnnData): cells anndata object.

        Returns:
            plotly.graph_objs._figure.Figure: Clustergram with the neighborhood enrichment results.
        """
        adata = self.cells_adata
        net_clustermap_fig = None
        if self.params[task_name]:
            cluster_key = self.configs["params"]["cell_cluster_key"]
            if f"{cluster_key}_nhood_enrichment" in adata.uns.keys():
                logger.info(
                    "<neighborhood_enrichment_plots> Generating neighborhood enrichment plots"
                )
                net_results = adata.uns["cell type_nhood_enrichment"]["zscore"]
                cell_types = adata.obs[cluster_key].cat.categories
                df_dash = pd.DataFrame(
                    net_results, columns=cell_types, index=cell_types
                )
                min_val = net_results.min()
                max_val = net_results.max()
                clustermap = dashbio.Clustergram(
                    data=df_dash,
                    column_labels=list(df_dash.columns.values),
                    row_labels=list(df_dash.index),
                    height=int(0.7 * self.width),
                    width=self.width,
                    color_map=[
                        [0, "rgb(255, 0, 0)"],
                        [(0 - min_val) / (max_val - min_val), "rgb(255,254,224)"],
                        [1, "rgb(0, 0, 255)"],
                    ],
                )
                clustermap.data[-1].colorbar.x = 1.5
                net_clustermap_fig = clustermap
                net_clustermap_fig.update_xaxes(tickangle=45)
                configs = dict(self.fig_configs[task_name])
                configs["figure"] = self.report_gen.format_to_html(net_clustermap_fig)
                configs["data"] = df_dash
                configs["name"] = task_name
                self.inline_figs.append(configs)
        return net_clustermap_fig

    def chord_plot(self, spot_stats, task_name="cell_net_chord"):
        """
        Generates a Chord plot for the neighborhood enrichment results.

        Args:
            adata (anndata.AnnData): cells anndata object.
            cluster_key (str): feature name in obs/var to be plotted.

        Returns:
            str: html for the chord plot
        """
        adata = self.cells_adata
        net_chord_fig = None
        n_cell_types = len(spot_stats.cell_type)
        cluster_key = self.configs["params"]["cell_cluster_key"]
        key_name = f"{cluster_key}_nhood_enrichment"
        if self.params[task_name] and n_cell_types > self.params["chord_n_cells"]:
            if key_name in adata.uns.keys():
                logger.info(
                    "<chord_plot> Generating Chord plot for neighborhood enrichment"
                )
                net_results = adata.uns[key_name]["zscore"]
                cell_types = adata.obs[cluster_key].cat.categories
                df = pd.DataFrame(net_results, columns=cell_types, index=cell_types)

                # Defining the edges of the Chord plot
                chord_df = df.unstack()
                chord_df = chord_df.reset_index().sort_values(
                    ["level_0", 0], ascending=[True, False]
                )
                chord_df.columns = ["source_name", "target_name", "value"]
                # create a LabelEncoder instance
                le = LabelEncoder()
                chord_df["cell_pair"] = chord_df.apply(
                    lambda row: "".join(
                        sorted([row["source_name"], row["target_name"]])
                    ),
                    axis=1,
                )
                chord_df = chord_df.drop_duplicates(subset=["cell_pair"]).drop(
                    columns=["cell_pair"]
                )
                chord_df["source"] = le.fit_transform(chord_df["source_name"])
                chord_df["target"] = le.transform(chord_df["target_name"])
                chord_df = chord_df[chord_df.source != chord_df.target]
                chord_df = chord_df[["source", "target", "value"]]

                # Defining the node labels
                nodes_df = pd.DataFrame(columns=["cell", "group"])
                nodes_df["group"] = range(0, len(df.index))
                nodes_df["cell"] = df.index
                nodes_cell = hv.Dataset(nodes_df, "index")

                # Creating Chord plot
                hv.extension("bokeh")
                hv.output(size=300)
                # Generating distinct colors for the differenc cell types
                colors = self.generate_distinct_colors(len(df.index))
                try:
                    chord_plot = hv.Chord((chord_df, nodes_cell)).select(
                        value=(self.params["net_cutoff"], None)
                    )
                except ValueError:
                    logger.info("<chord_plot> Failed to create Chord plot")
                    return None
                chord_plot.opts(
                    opts.Chord(
                        cmap=colors,
                        edge_cmap=colors,
                        edge_color=dim("source").str(),
                        labels="cell",
                        node_color=dim("index").str(),
                        tools=["hover", "box_select", "tap", "save", "reset"],
                        active_tools=["box_select"],
                    )
                )
                net_chord_fig = hv.renderer("bokeh").html(chord_plot)
                configs = dict(self.fig_configs[task_name])
                configs["figure"] = self.report_gen.format_to_html(net_chord_fig)
                self.inline_figs.append(configs)
        return net_chord_fig

    def cooccurrence_plots(self, task_name="cell_cooccur"):
        """
        Plots co-occurence analysis on a Plotly figure. Prerequisite: run SpatialStructAnalysis.cooccurence_analysis().

        Args:
            adata (anndata.AnnData): cells anndata object.

        Returns:
            list: list of cooccurence plots (plotly.graph_objs._figure.Figure).
            list: cell types.
        """
        adata = self.cells_adata
        cluster_key = self.configs["params"]["cell_cluster_key"]
        key_name = f"{cluster_key}_co_occurrence"
        cell_types = []
        # Creating a Plotly figure
        fig = go.Figure()
        categories = []
        if self.params[task_name] and key_name in adata.uns.keys():
            logger.info("<cooccurence_plots> Generating co-occurrence plots")
            # Getting intervals and occurence data
            occurrence_data = adata.uns[f"{cluster_key}_co_occurrence"]
            out = occurrence_data["occ"]
            interval = occurrence_data["interval"][:-1]
            distances = [self.convert_micron(i) for i in interval]
            # Getting cell types
            cell_types = adata.obs[cluster_key].cat.categories

            cooccur_list = []
            for target_cell in cell_types:
                # Getting index of the target cell type
                target_cell_idx = np.where(cell_types == target_cell)[0][0]
                # Getting subset of relevant co-occurence data for the target cell
                relevant_occ_data = out[target_cell_idx, :, :]
                # Looping through all other cell types and adding traces to the figure
                for other_cell in cell_types:
                    cell_type_idx = np.where(cell_types == other_cell)[0][0]
                    cell_occ_data = relevant_occ_data[cell_type_idx, :]
                    trace = go.Scatter(
                        x=distances,
                        y=cell_occ_data,
                        name=other_cell,
                    )
                    fig.add_trace(trace)
                    categories.append((target_cell, other_cell))
                    cooccur_list.append(
                        cell_occ_data.tolist() + [target_cell, other_cell]
                    )
            drop_down_entries = self.create_drop_down(
                categories, mode="cooccurence", showlegend=True
            )
            # Adding figure title
            fig.update_layout(
                title=rf"$\frac{{p(exp|targetcell)}}{{p(exp)}}$",
                xaxis_title="Distance (microns)",
                yaxis_title="Value",
                width=self.width,
                height=self.height,
                xaxis_range=[0, 1000],
            )
            if drop_down_entries != []:
                default_plot_ix = 0
                fig.update_traces(visible=False)  # Hide all traces initially
                fig.data[default_plot_ix].visible = True  # Show the default trace
                fig.update_layout(
                    showlegend=False,
                    updatemenus=[
                        go.layout.Updatemenu(
                            active=default_plot_ix, buttons=drop_down_entries
                        )
                    ],
                )
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            configs["data"] = pd.DataFrame(
                cooccur_list, columns=distances + ["target_cell", "other_cell"]
            )
            configs["name"] = task_name
            self.inline_figs.append(configs)
        return fig

    def morans_i_exp_plots(self, task_name="morans_gene_heatmap", mode="gene"):
        """
        Plots the expression of the top n genes that show the highest Moran's I.

        Returns:
            adata (anndata.AnnData): Spots ann data object.
        """
        adata = self.spots_adata
        fig = None
        pval = self.params["moran_pval"]
        n_cols = self.params["n_cols"]
        data_label = f"{mode}_moranI"
        morans_df = pd.DataFrame()
        if self.params[task_name] and data_label in adata.uns.keys():
            logger.info("<morans_i_plots> Generating Moran's I plots")
            morans_df = adata.uns[data_label]
            morans_df = morans_df[morans_df["pval_norm"] < pval]
        if not morans_df.empty:
            n = np.minimum(self.params["n_genes_exp"], len(morans_df))
            save_location = self.cache_path + self.configs["params"]["morans_save_path"]
            fig = self.generate_gene_distrib_collage(
                save_location=save_location,
                adata=adata,
                markers=morans_df.index.tolist(),
                n=n,
                n_rows=int(np.ceil(n / n_cols)),
                n_cols=n_cols,
            )
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            configs["data"] = morans_df
            configs["name"] = task_name
            self.inline_figs.append(configs)
        return fig

    def morans_i_bar_plot(self, task_name="morans_gene_bar", mode="gene"):
        """
        Generates Moran's I bar plot

        Returns:
            adata (anndata.AnnData): Spots ann data object.
        """
        adata = self.spots_adata
        bar_fig = None
        pval = self.params["moran_pval"]
        data_label = f"{mode}_moranI"
        morans_df = pd.DataFrame()
        if self.params[task_name] and data_label in adata.uns.keys():
            logger.info("<morans_i_bar_plot> Generating Moran's I bar plot")
            morans_df = adata.uns[data_label]
            morans_df = morans_df[morans_df["pval_norm"] < pval]
        if not morans_df.empty:
            n = np.minimum(self.params["n_genes_bar"], len(morans_df))
            # Create a plotly barplot of the number of cells per cell type
            df = morans_df.head(n)
            bar_fig = px.bar(
                df, x=df.index, y=df.I, height=self.height, width=self.width
            )
            # Add title
            bar_fig.update_layout(
                title="Moran's I",
                xaxis=dict(title="Gene"),
                yaxis=dict(title=f"Morans I"),
            )
            bar_fig.update_xaxes(tickangle=45)
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(bar_fig)
            self.inline_figs.append(configs)
        return bar_fig

    def spatialde_heatmap_plot(
        self, task_name="spatialde_heatmap", min_count_thresh=1000
    ):
        """Generates heatmaps to visualize top spatially variable genes
        Returns:
            adata (anndata.AnnData): Spots ann data object.
        """
        adata = self.spots_adata
        fig = None
        pval = self.params["spatialde_pval"]
        n_cols = self.params["n_cols"]
        spatialde_df = pd.DataFrame()
        data_label = "spatialde"
        if self.params[task_name] and data_label in adata.uns.keys():
            logger.info(
                "<spatialde_heatmap_plot> Generating average gene expression plots"
            )
            spatialde_df = adata.uns[data_label]
            spatialde_raw_df = spatialde_df.copy()
            # Only showing the siginficantly spatially variable genes
            significance_flag = spatialde_df["qval"] <= pval
            spatialde_df = spatialde_df[significance_flag]

            # Only showing genes with atleast min_count_thresh counts
            valid_counts_flag = spatialde_df["total_counts"] >= min_count_thresh
            spatialde_df = spatialde_df[valid_counts_flag]
        if not spatialde_df.empty:
            n = np.minimum(self.params["n_genes_exp"], len(spatialde_df))
            # Sort by significance and FSV (Fraction of variance explained by spatial variation)
            top_genes_by_sig = spatialde_df.sort_values(
                by=["qval", "FSV"], ascending=[True, False], na_position="first"
            )
            # Showing the top n genes by spatial variance
            genes_to_render = top_genes_by_sig["g"].tolist()[:n]

            save_location = (
                self.cache_path + self.configs["params"]["spatialde_save_path"]
            )
            fig = self.generate_gene_distrib_collage(
                save_location=save_location,
                adata=adata,
                markers=genes_to_render,
                n=n,
                n_rows=int(np.ceil(n / n_cols)),
                n_cols=n_cols,
            )
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            configs["data"] = spatialde_raw_df
            configs["name"] = task_name
            self.inline_figs.append(configs)
        return fig

    def spatialde_bar_plot(self, task_name="spatialde_bar", min_count_thresh=1000):
        """
        Generates SpatialDE bar plot
        Returns:
            adata (anndata.AnnData): Spots ann data object.
        """
        adata = self.spots_adata
        bar_fig = None
        pval = self.params["spatialde_pval"]
        spatialde_df = pd.DataFrame()
        data_label = "spatialde"
        if self.params[task_name] and data_label in adata.uns.keys():
            logger.info("<spatialde_bar_plot> Generating SpatialDE bar plot")
            spatialde_df = adata.uns[data_label]
            spatialde_raw_df = spatialde_df.copy()
            # Only showing the siginficantly spatially variable genes
            significance_flag = spatialde_df["qval"] <= pval
            spatialde_df = spatialde_df[significance_flag]

            # Only showing genes with atleast min_count_thresh counts
            valid_counts_flag = spatialde_df["total_counts"] >= min_count_thresh
            spatialde_df = spatialde_df[valid_counts_flag]
        if not spatialde_df.empty:
            n = np.minimum(self.params["n_genes_bar"], len(spatialde_df))
            # Sort by significance and FSV (Fraction of variance explained by spatial variation)
            top_genes_by_sig = spatialde_df.sort_values(
                by=["qval", "FSV"], ascending=[True, False], na_position="first"
            )

            # Create a plotly barplot of the number of cells per cell type
            df = top_genes_by_sig.head(n)
            bar_fig = px.bar(df, x=df.g, y=df.FSV, height=self.height, width=self.width)
            # Add title
            bar_fig.update_layout(
                title="SpatialDE Results",
                xaxis=dict(title="Gene"),
                yaxis=dict(title=f"FSV"),
            )
            bar_fig.update_xaxes(tickangle=45)
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(bar_fig)
            self.inline_figs.append(configs)
        return bar_fig

    def create_drop_down(self, entries, mode="volcano", sort=True, showlegend=False):
        """Generates plotly drop down menu fields

        Args:
            entries (_type_): _description_

        Returns:
            _type_: _description_
        """
        options = []
        if entries == []:
            return options
        elif type(entries[0]) is not tuple:
            entries = np.array(list(zip(entries, entries)))
        else:
            entries = np.array(entries)
        visible = [False] * len(entries)
        # Applying unique function on array
        res, ind = np.unique(entries[:, 0], return_index=True)
        # Sorting indices
        drop_down_entries = res[np.argsort(ind)]
        if sort:
            drop_down_entries = sorted(drop_down_entries)
        for entry in drop_down_entries:
            visibility_field = visible.copy()
            trace_ix = np.where(entries[:, 0] == entry)[0]
            default_ix = np.where((entries == np.array([entry, entry])).all(axis=1))[0]
            visibility_field = [
                "legendonly" if ix in trace_ix else False for ix in range(len(entries))
            ]
            visibility_field[default_ix[0]] = True
            if mode == "cooccurence":
                title = (rf"$\frac{{p(exp|{entry})}}{{p(exp)}}$",)
            else:
                title = entry
            options.append(
                dict(
                    label=entry,
                    method="update",
                    args=[
                        {"visible": visibility_field},
                        {"title": title, "showlegend": showlegend},
                    ],
                )
            )
        return options

    def volcano_plots(self, task_name="diff_exp_annotations"):
        """Generates a grid of volcano plots for comparing gene expression.

        Args:
            adata (anndata.AnnData): Dataframe of gene statistics.

        Returns:
            list: (figures) list of figures.
            list: (categories) list of categories being compared.
        """
        adata = self.spots_adata
        fc_threshold = self.params["diff_fc_thresh"]
        color_discrete_map = self.params["color_discrete_map"]
        pval = self.params["diff_pval"]
        stats_df = pd.DataFrame()
        categories = []
        if self.params[task_name] and adata.uns[task_name] is not None:
            logger.info("<volcano_plots> Generating volcano plots.")
            stats_df = adata.uns[task_name]
        fig = go.Figure()

        if not stats_df.empty:
            combinations = (
                stats_df[["Region A", "Region B"]].drop_duplicates().to_numpy()
            )
            for cat_a, cat_b in combinations:
                if cat_a == cat_b:
                    continue
                df = stats_df[
                    (stats_df["Region A"] == cat_a) & (stats_df["Region B"] == cat_b)
                ]
                trace_name = f"{cat_a}-{cat_b}"
                color_list = [
                    color_discrete_map[direction] for direction in df["Direction"]
                ]
                scatter_trace = go.Scatter(
                    x=df["log2FC"],
                    y=df["-log10pvalue"],
                    mode="markers",
                    marker=dict(color=color_list),
                    customdata=df[
                        [
                            "Gene",
                            "Mean1",
                            "Mean2",
                            "FC",
                            "Corrected T-test pvalue",
                            "Corrected MW pvalue",
                        ]
                    ],
                    hovertemplate="<br>".join(
                        [
                            "Gene: %{customdata[0]}",
                            "Mean1: %{customdata[1]}",
                            "Mean2: %{customdata[2]}",
                            "FC: %{customdata[3]}",
                            "Corrected T-test pvalue: %{customdata[4]}",
                            "Corrected MW pvalue: %{customdata[5]}",
                        ]
                    ),
                    name=trace_name,
                )
                fig.add_trace(scatter_trace)
                categories += [trace_name]

            title = f"Volcano Plot of Gene Expression Differences Between {cat_a} and {cat_b} - Fold Change Threshold = {fc_threshold} - corrected p-value = {pval}"
            formatted_title = f'<span style="font-size: 14px;">{title}</span>'
            fig.update_layout(
                title_text=formatted_title, height=self.height, width=self.width
            )
            fig.add_hline(y=-np.log10(pval), line_dash="dash", line_color="grey")
            fig.add_vline(x=np.log2(fc_threshold), line_dash="dash", line_color="grey")
            fig.add_vline(
                x=np.log2(1 / fc_threshold), line_dash="dash", line_color="grey"
            )
            drop_down_entries = self.create_drop_down(categories)
            if drop_down_entries != []:
                default_plot_ix = 0
                fig.update_traces(visible=False)  # Hide all traces initially
                fig.data[default_plot_ix].visible = True  # Show the default trace
                fig.update_layout(
                    showlegend=False,
                    updatemenus=[
                        go.layout.Updatemenu(
                            active=default_plot_ix, buttons=drop_down_entries
                        )
                    ],
                )
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            configs["data"] = stats_df
            configs["name"] = task_name
            self.inline_figs.append(configs)
        return fig

    def abundance_vs_cells_plot(
        self, task_name="infilt_comparing_cell_types", infiltrating_cells=[]
    ):
        """Creates plots to study cell proportions across the different infiltration levels plots

        Args:
            infilt_results (pandas.DataFrame): Results of infiltration analysis. Generated
                by SpatialAnalysis.infiltration_analysis()
            scale_factors (dict): Scale factors dictionary

        Returns:
            list: (figures) list of figures
            list: (categories) list of intervals of study.
        """
        adata = self.cells_adata
        categories = []
        color_discrete_map = self.params["color_discrete_map"]
        fig = go.Figure()
        if self.params[task_name] and "infilt_analysis" in adata.uns.keys():
            logger.info(
                "<abundance_vs_cells_plot> Generating cell abundance plots at infiltration levels"
            )
            infilt_results = self.cells_adata.uns["infilt_analysis"]
            num_cells_per_level = infilt_results.groupby(["infilt_level", "cell_type"])[
                "cell_x"
            ].count()
            cell_types = sorted(infilt_results.cell_type.unique())
            levels = sorted(infilt_results.infilt_level.unique())[::-1]
            step_size = int(self.scale_factors["spot_diameter_fullres"])
            for i, level in enumerate(levels):
                counts = pd.DataFrame()
                counts["cell_type"] = cell_types
                counts = counts.merge(
                    num_cells_per_level[level], on="cell_type", how="left"
                )
                counts = counts.rename(columns={"cell_x": "abundance"})
                counts["proportion"] = counts["abundance"] / counts["abundance"].sum()
                counts = counts.fillna(0)
                if i == 0:
                    counts = counts.sort_values(
                        "abundance", ascending=False
                    ).reset_index(drop=True)
                    cell_types = counts.cell_type.tolist()
                upper_bound = int(
                    self.convert_micron(
                        level,
                    )
                )
                lower_bound = int(self.convert_micron(level - step_size))
                counts["infiltration"] = "not infiltrating"
                counts.loc[
                    counts["cell_type"].isin(infiltrating_cells), "infiltration"
                ] = "infiltrating"
                color_list = [
                    color_discrete_map[direction]
                    for direction in counts["infiltration"]
                ]
                trace_name = f"{lower_bound} to {upper_bound} um"
                trace = go.Bar(
                    x=counts["cell_type"],
                    y=counts["abundance"],
                    marker=dict(color=color_list),
                    name=trace_name,
                )
                fig.add_trace(trace)
                categories += [trace_name]
            fig.update_xaxes(tickangle=45)
            fig.update_layout(
                xaxis=dict(title="Cell types"),
                yaxis=dict(title=f"Cell abundance"),
            )
            drop_down_entries = self.create_drop_down(categories, sort=False)
            if drop_down_entries != []:
                default_plot_ix = 0
                fig.update_traces(visible=False)  # Hide all traces initially
                fig.data[default_plot_ix].visible = True  # Show the default trace
                fig.update_layout(
                    showlegend=False,
                    updatemenus=[
                        go.layout.Updatemenu(
                            active=default_plot_ix, buttons=drop_down_entries
                        )
                    ],
                )
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            self.inline_figs.append(configs)
        return fig

    def abundance_vs_level(self, task_name="infilt_comparing_levels", buffer=True):
        """Creates plots to study cell proportions across the different infiltration levels plots

        Args:
            infilt_results (pandas.DataFrame): Results of infiltration analysis. Generated
                by SpatialAnalysis.infiltration_analysis()
            scale_factors (dict): Scale factors dictionary

        Returns:
            list: (figures) list of figures
            list: (cell_types) list of cell types
        """
        adata = self.cells_adata
        cell_types = []
        labels = []
        infiltrating_cells = []
        comp_analysis = ComparativeAnalysis(configs=self.configs)
        fig = go.Figure()
        if self.params[task_name] and "infilt_analysis" in adata.uns.keys():
            logger.info(
                "<abundance_vs_level> Generating cell proportions across infiltration levels plots"
            )
            infilt_results = self.cells_adata.uns["infilt_analysis"]
            num_cells_per_level = infilt_results.groupby(["infilt_level", "cell_type"])[
                "cell_x"
            ].count()
            levels = sorted(infilt_results.infilt_level.unique())[::-1]
            # Getting proportions at each level
            counts_df = pd.DataFrame()
            step_size = int(self.scale_factors["spot_diameter_fullres"])
            step_size_micron = int(self.convert_micron(step_size))
            for level in levels:
                counts = pd.DataFrame(num_cells_per_level[level]).reset_index()
                counts = counts.rename(columns={"cell_x": "abundance"})
                counts = counts.sort_values("abundance", ascending=False)
                counts["proportion"] = counts["abundance"] / counts["abundance"].sum()
                counts["level"] = level
                counts_df = pd.concat([counts_df, counts])
            counts_df["upper_bound"] = [
                int(self.convert_micron(l)) for l in counts_df.level
            ]
            counts_df["lower_bound"] = [
                int(self.convert_micron(l - step_size)) for l in counts_df.level
            ]
            counts_df["mid_point"] = (
                counts_df["upper_bound"] + counts_df["lower_bound"]
            ) / 2
            cell_types = sorted(infilt_results.cell_type.unique())
            for cell_type in cell_types:
                df = counts_df[counts_df.cell_type == cell_type]
                df = df.sort_values("level")
                _, pval, sig_difference = comp_analysis.test_infiltration(
                    counts_df,
                    cell_type,
                    pval_thresh=self.params["diff_pval"],
                    buffer=buffer,
                )
                label = cell_type
                if sig_difference:
                    label += "*"
                    infiltrating_cells += [cell_type]
                labels += [label]
                trace = go.Bar(
                    x=df["mid_point"],
                    y=df["proportion"],
                    name=label,
                    text=f"p value: {pval:.3E}",
                )
                fig.add_trace(trace)

            # Add title
            fig.update_layout(
                xaxis=dict(title="Distance from region boundary"),
                yaxis=dict(title=f"Cell proportion"),
            )
            fig.update_xaxes(
                tickangle=45,
                tickmode="array",
                tickvals=counts_df["mid_point"],
                ticktext=counts_df["lower_bound"].astype(str)
                + " to "
                + counts_df["upper_bound"].astype(str)
                + " um",
            )
            fig.update_traces(width=50)
            fig.add_vline(x=0, line_dash="dash", line_color="grey")
            if buffer:
                # Highlight buffer region
                fig.add_vrect(
                    x0=-step_size_micron,
                    x1=step_size_micron,
                    line_width=0,
                    fillcolor="grey",
                    opacity=0.2,
                )
            fig.update_xaxes(
                range=[
                    counts_df["lower_bound"].min(),
                    counts_df["upper_bound"].max(),
                ]
            )
            drop_down_entries = self.create_drop_down(labels, sort=True)
            if drop_down_entries != []:
                default_plot_ix = 0
                fig.update_traces(visible=False)  # Hide all traces initially
                fig.data[default_plot_ix].visible = True  # Show the default trace
                fig.update_layout(
                    showlegend=False,
                    updatemenus=[
                        go.layout.Updatemenu(
                            active=default_plot_ix, buttons=drop_down_entries
                        )
                    ],
                )
            configs = dict(self.fig_configs[task_name])
            configs["figure"] = self.report_gen.format_to_html(fig)
            self.inline_figs.append(configs)
        return fig, infiltrating_cells

    def get_report(self, exp_id, report_name):
        """Populates the report

        Args:
            exp_id (_type_): Experiment ID
            report_name (_type_): Report name

        Returns:
            str: html report
        """
        return self.report_gen.populate_report(
            image_id=exp_id,
            report_name=report_name,
            inline_figs=self.inline_figs,
        )
