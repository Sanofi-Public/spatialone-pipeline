"""
Python class for generation a .tmap file for viewing results on TissUUmaps
"""
import json
import os

import pandas as pd
import yaml
from distinctipy import distinctipy


class tmapGenerator:
    """
    Class to to create the .tmap file
    """

    def __init__(self):
        """
        Creating instance variables
        """
        self.templates_path = "templates/"
        self.tmap_template_path = self.templates_path + "tmap_template.tmap"
        self.marker_configs_path = self.templates_path + "marker_templates.json"
        self.tmap_template = self.load_tmap()
        self.marker_configs = self.load_marker_configs()
        self.mode = "source-over"
        self.markers_to_show = {
            "num_cells_per_spot": False,
            "morphological_clustering": False,
            "cell_assignment": False,
            "cell_proportions": False,
            "graphical_clustering": False,
            "kmeans_clustering": False,
            "spot_qc": False,
            "primary_gene_dropdown": True,
            "secondary_gene_dropdown": True,
        }
        self.tmap_params = {
            "piechart_col_name": "",
            "gene_names": [],
            "qc_metrics": [],
            "graphical_clustering": [],
            "kmeans_clustering": [],
        }

    def load_tmap(self, tmap_path=None):
        """
        Load the tmap template
        """
        if tmap_path is None:
            tmap_path = self.tmap_template_path
        with open(tmap_path, "r") as stream:
            return yaml.safe_load(stream)

    def load_marker_configs(self):
        """
        Load the marker configs template
        """
        with open(self.marker_configs_path, "r") as f_contents:
            return json.load(f_contents)

    def save_tmap(self, yaml_file, tmap_path):
        """
        Save the tmap template file
        """
        with open(tmap_path, "w") as outfile:
            outfile.write(json.dumps(yaml_file, indent=4))

    def create_piechart_df(self, adata, save_dir):
        """Creates a dataframe for rendering piecharts on TissUUmaps

        Args:
            adata (AnnData): spots AnnData object
            save_dir (str): directory to save the piechart dataframe

        Returns:
            str: string describing the piechart column name
        """
        piechart_col_name = ""
        if "deconv_results" in adata.obsm:
            adata.obsm["deconv_results"] = adata.obsm["deconv_results"].fillna(0)
            deconv_results = adata.obsm["deconv_results"]
            piechart_col_name = ";".join(deconv_results.columns)
            adata.obsm["deconv_piechart"] = pd.DataFrame(
                deconv_results.round(2).astype(str).agg(";".join, axis=1),
                columns=[piechart_col_name],
            )
            # Creating piechart table
            pd.concat(
                [
                    adata.obsm.get("spatial", pd.DataFrame()),
                    adata.obsm.get("deconv_piechart", pd.DataFrame()),
                    adata.obsm.get("cell2spot", pd.DataFrame()),
                    adata.obs.get("collection_id", pd.DataFrame()),
                ],
                axis=1,
            ).to_csv(os.path.join(save_dir, "piechart_df.csv"))
        return piechart_col_name

    def get_tmap_params(self, spots_adata, cells_adata, piechart_col_name):
        """Function to get the tmap file parameters and identify which
        sections of the .tmap file to enable.

        Args:
            spots_adata (AnnData): spots AnnData object
            spots_adata (AnnData): cells AnnData object
            piechart_col_name (str): string describing the piechart column name

        Returns:
            dict: dictionary dictating which sections of the .tmap to enable
            dict: dictionaru of tmap parameters
        """

        self.tmap_params.update({"gene_names": spots_adata.var_names.tolist()})
        self.tmap_params.update({"piechart_col_name": piechart_col_name})

        if piechart_col_name != "":
            self.markers_to_show["cell_proportions"] = True

        if "qc" in spots_adata.obsm:
            self.tmap_params.update(
                {"qc_metrics": spots_adata.obsm["qc"].columns.tolist()}
            )
            self.markers_to_show["spot_qc"] = True

        if "graphical_clustering" in spots_adata.obsm:
            self.tmap_params.update(
                {
                    "graphical_clustering": sorted(
                        spots_adata.obsm["graphical_clustering"].columns.tolist()
                    )
                }
            )
            self.markers_to_show["graphical_clustering"] = True

        if "kmeans_clustering" in spots_adata.obsm:
            self.tmap_params.update(
                {
                    "kmeans_clustering": sorted(
                        spots_adata.obsm["kmeans_clustering"].columns.tolist()
                    )
                }
            )
            self.markers_to_show["kmeans_clustering"] = True

        if "cell_assignment" in cells_adata.obsm:
            self.markers_to_show["cell_assignment"] = True
        if "morphological_cluster" in cells_adata.obsm:
            self.markers_to_show["morphological_clustering"] = True
        return self.markers_to_show, self.tmap_params

    def update_tmap_layers(self, tmap_yaml, run_ids):
        """Updating the tmap file layers section. Controls each layer's contrast, opacity, etc.

        Args:
            tmap_yaml (dict): tmap file as a dictionary
            run_ids (list): list of runs to visualize

        Returns:
            dict : tmap file as a dictionary
        """
        layer_filters = {}
        layer_opacities = {}
        layer_visibilties = {}
        layers = []
        file_name = ""
        filter_defaults = [
            {"name": "Saturation", "value": "0"},
            {"name": "Brightness", "value": "0"},
            {"name": "Contrast", "value": "1"},
        ]
        for i, run in enumerate(run_ids):
            layer_filters[str(i)] = filter_defaults
            layer_opacities[str(i)] = "1"
            layer_visibilties[str(i)] = True
            layers.append(
                {
                    "name": f"../../prep/{run}/wsi.tif",
                    "tileSource": f"../../prep/{run}/wsi.tif.dzi",
                }
            )
            if i == 0:
                file_name = run
            else:
                file_name += f", {run}"
        if len(run_ids) == 1:
            additional_layers = ["cells_layer.png", "spots_layer.png"]
            for i, layer_name in enumerate(additional_layers):
                ix = i + 1
                layer_filters[str(ix)] = filter_defaults
                layer_opacities[str(ix)] = "0.3"
                layer_visibilties[str(ix)] = False
                layers.append({"name": layer_name, "tileSource": f"{layer_name}.dzi"})
        tmap_yaml.update(
            {
                "filename": file_name,
                "layerFilters": layer_filters,
                "layerOpacities": layer_opacities,
                "layerVisibilities": layer_visibilties,
                "layers": layers,
            }
        )
        return tmap_yaml

    def update_tmap_markers(self, tmap_yaml, markers, tmap_params, mode="collection"):
        """Updates the tmap file markers section

        Args:
            tmap_yaml (dict): tmap file as a dictionary
            markers (dict): dictionary dictating which sections of the .tmap to enable
            tmap_params (dict): dictionary of tmap parameters

        Returns:
            dict : tmap file as a dictionary
        """
        marker_configs = []
        if markers["cell_proportions"]:
            piechart_col_name = tmap_params["piechart_col_name"]
            celltype_colors = self.get_cell_color_dict(piechart_col_name)
        else:
            celltype_colors = {}

        for i, marker_name in enumerate(markers):
            if markers[marker_name]:
                configs = self.marker_configs[marker_name]
                configs["fromButton"] = i
                if mode == "collection":
                    configs["expectedRadios"]["collectionItem_col"] = True
                else:
                    configs["expectedRadios"]["collectionItem_col"] = False
                if marker_name == "cell_proportions":
                    configs["expectedHeader"]["pie_col"] = piechart_col_name
                    configs["expectedHeader"]["pie_dict"] = celltype_colors
                if marker_name == "cell_assignment":
                    configs["expectedHeader"]["cb_gr_dict"] = celltype_colors
                if marker_name == "graphical_clustering":
                    graphical_clustering = tmap_params["graphical_clustering"]
                    drop_down_options = []
                    for cluster_type in graphical_clustering:
                        drop_down_entry = dict(configs["dropdownOptions"][0])
                        drop_down_entry[
                            "expectedHeader.cb_col"
                        ] = f"/obsm/graphical_clustering/{cluster_type}"
                        drop_down_entry[
                            "expectedHeader.gb_col"
                        ] = f"/obsm/graphical_clustering/{cluster_type}"
                        drop_down_entry[
                            "expectedHeader.sortby_col"
                        ] = f"/obsm/graphical_clustering/{cluster_type}"
                        drop_down_entry["name"] = f"Gene Clustering: {cluster_type}"
                        drop_down_entry["optionName"] = cluster_type
                        drop_down_options.append(drop_down_entry)
                    configs["dropdownOptions"] = drop_down_options
                if marker_name == "kmeans_clustering":
                    kmeans_clustering = tmap_params["kmeans_clustering"]
                    drop_down_options = []
                    for cluster_type in kmeans_clustering:
                        drop_down_entry = dict(configs["dropdownOptions"][0])
                        drop_down_entry[
                            "expectedHeader.cb_col"
                        ] = f"/obsm/kmeans_clustering/{cluster_type}"
                        drop_down_entry[
                            "expectedHeader.gb_col"
                        ] = f"/obsm/kmeans_clustering/{cluster_type}"
                        drop_down_entry[
                            "expectedHeader.sortby_col"
                        ] = f"/obsm/kmeans_clustering/{cluster_type}"
                        drop_down_entry["name"] = f"Gene Clustering: {cluster_type}"
                        drop_down_entry["optionName"] = cluster_type
                        drop_down_options.append(drop_down_entry)
                    configs["dropdownOptions"] = drop_down_options
                if marker_name == "spot_qc":
                    qc_metrics = tmap_params["qc_metrics"]
                    drop_down_options = []
                    for qc_metric in qc_metrics:
                        drop_down_entry = dict(configs["dropdownOptions"][0])
                        drop_down_entry[
                            "expectedHeader.cb_col"
                        ] = f"/obsm/qc/{qc_metric}"
                        drop_down_entry[
                            "expectedHeader.sortby_col"
                        ] = f"/obsm/qc/{qc_metric}"
                        drop_down_entry["name"] = f"QC Metric: {qc_metric}"
                        drop_down_entry["optionName"] = qc_metric
                        drop_down_options.append(drop_down_entry)
                    configs["dropdownOptions"] = drop_down_options
                if marker_name in ["primary_gene_dropdown", "secondary_gene_dropdown"]:
                    gene_names = tmap_params["gene_names"]
                    drop_down_options = []
                    for i, gene in enumerate(gene_names):
                        drop_down_entry = dict(configs["dropdownOptions"][0])
                        drop_down_entry["expectedHeader.cb_col"] = f"/X;{i}"
                        drop_down_entry["expectedHeader.sortby_col"] = f"/X;{i}"
                        drop_down_entry["name"] = f"Gene expression: {gene}"
                        drop_down_entry["optionName"] = gene
                        drop_down_options.append(drop_down_entry)
                    configs["dropdownOptions"] = drop_down_options
                marker_configs.append(self.marker_configs[marker_name])
        tmap_yaml["markerFiles"] = marker_configs
        return tmap_yaml

    def get_cell_color_dict(self, piechart_col_name):
        """
        Builds a colors dictionary to ensure color consistency between
        the piechart tab and cell assignment tab
        """
        cell_classes = piechart_col_name.split(";")
        colors = distinctipy.get_colors(len(cell_classes))
        colors_hex = [
            "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
            for r, g, b in colors
        ]
        class_colors = dict(zip(cell_classes, colors_hex))
        return class_colors

    def create_tmap_file(
        self, run_ids, spots_adata, cells_adata, results_dir, annotation_file
    ):
        """
        Creates a tmap file for the given experiment
        """
        tmap_yaml = self.tmap_template.copy()
        tmap_yaml["compositeMode"] = self.mode
        # Adding regions file to .tmap
        tmap_yaml["regionFiles"] = [
            {
                "path": f"../../prep/{run_ids[0]}/{annotation_file}",
                "title": "Load regions",
            }
        ]
        # Getting .tmap file parameters
        piechart_col_name = self.create_piechart_df(spots_adata, results_dir)
        markers_to_show, tmap_params = self.get_tmap_params(
            spots_adata, cells_adata, piechart_col_name
        )
        # Adding images to .tmap file
        tmap_yaml = self.update_tmap_layers(tmap_yaml=tmap_yaml, run_ids=run_ids)
        # Adding buttons and dropdowns to .tmap file
        tmap_yaml = self.update_tmap_markers(
            tmap_yaml=tmap_yaml,
            markers=markers_to_show,
            tmap_params=tmap_params,
            mode=self.mode,
        )
        # Saving the .tmap file
        save_path = os.path.join(results_dir, "experiment_tmap.tmap")
        self.save_tmap(yaml_file=tmap_yaml, tmap_path=save_path)
        print(
            f"[SPATIALONE] Successfully created .tmap file for {run_ids} under folder named {save_path}"
        )
        return tmap_yaml
