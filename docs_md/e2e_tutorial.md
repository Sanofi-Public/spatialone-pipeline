# Analysis of 10x's Lung Cancer Squamous Cell Carcinoma sample
This tutorial provides a step-by-step guide on how to analyze a [10x sample dataset](https://www.10xgenomics.com/datasets/human-lung-cancer-11-mm-capture-area-ffpe-2-standard) consisting on a FFPE human lung cancer tissue diagnosed with _Neuroendocrine Carcinom_a.

#### 1. Clone the SpatialOne git repository from Github and move into the spatialone-pipeline folder


 ```sh
    git clone https://github.com/Sanofi-Public/spatialone-pipeline.git
    cd spatialone-pipeline
 ```
This will download all the necessary source code to run the analysis

#### 2. Download the sample data.

The following script will download all the necessarily data for running the analysis from zenodo. The script also checks the file checksum to check data integrity, and then unzips it.
```sh
./download_experiment_data.sh
```

Alternatively, you can run the following commands to directly download and unzip the data:
```sh
curl -L -o SpatialOne_Data.tar.gz "https://zenodo.org/records/12605154/files/SpatialOne_Data.tar.gz?download=1"
tar -xzvf SpatialOne_Data.tar.gz
```

Make sure that the sample data has been stored under the **prep** folder.

#### 3. Set up the system environment.

Edit your **.env** file to reflect your system setup.

- The _**HOST_DATA_PATH**_ variable should reflect the path where the data downloaded in the previous step has been stored.
  If you are following this tutorial, this variable corresponds to the on where you have cloned the github project. After running step 2, the The HOST_DATA_PATH folder should contain the _prep_, _conf_ and _reference_ folders.

- The **_GPU_DEVICE_ID_** variable should reflect the ID of your GPU (0 if only 1 GPU is available in the system). This parameter is not required if the analysis is run using only cpu.

For instance:
```sh
# Update your configurations for data and GPU
HOST_DATA_PATH=/Users/demo_user/Documents/spatialone-pipeline/ ## This path should point to the 'data' path where your experiment data is stored
GPU_DEVICE_ID=0 # select GPU ID access for docker image, 0 if you only have 1 GPU.

# Project variables, optional
REPOSITORY_NAME="spatialone-pipeline"
METAFLOW_USER="user@email.com"
METAFLOW_DEFAULT_DATASTORE="local"
PROJECT_NAME="spatialone"
```

#### 4. Build or retrieve the SpatialOne docker container.

You can generate the SpatialOne docker container by executing the following command:
```sh
make build
```
This process may take up to 20 minutes.

Alternatively, if you are using an amd architecture, you can retrieve the docker container from dockerhub:
```sh
docker pull albertpla/spatialone_amd:latest
docker tag albertpla/spatialone_amd:latest spatialone-pipeline:latest # This will rename the docker image to be coherent with the spatialone setup
```

#### 5. Setting up the reference data for cell deconvolution

The reference dataset that will be used during the cell deconvolution step needs to be stored under the __reference__ folder using the following nomenclature **reference_name_cell_atlas.h5ad**. In this tutorial we'll use the _luca_cell_atlas.h5ad_ dataset that was derived from the (Single-cell Lung Cancer Atlas)[https://luca.icbi.at/]. The reference data need to be anndata packaged following the guidelines provided by the scverse and [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html#Loading-Visium-and-scRNA-seq-reference-data)

The reference anndata file used in this example follows this structure:
- Raw gene expression counts from the atlas are stored in /X
- Data is gzip compressed
- Single cell barcode is "_experiment_", Expeirment/batch tag is "_batch_", and cell annotation is defined in "_cell_type_"

#### 6. Set up the configuration of the analysis
Spatial one uses yaml files to set up its configuraiton. Configuration files should be stored in the _conf_ folder; _visium_config_flow.yaml_ is the file name used by default. We reccomend keeping this config filename so the analysis can be run by simply running _make run_

In this analysis we will use _cellpose_ to segment cells in the H&E image, _cell2location_ to estimate the cell type proportion at each step, the local assignment method to estimate cell types of each segemented cell, and we'll run a standard downstream spatial analyisis.

The first configuration block corresponds to the experiment metadata

```yaml
user.id: "user@email.com" # for internal metadata handling
run.name: "SpatialOne Tutorial Example"
run.description: "Squamous Cell Carcinoma analyzed with Cellpose & Cell2location (LuCA)"
run.summary: "SCC analyiss"
experiment.ids: [
    # List of experiments to analyze
    # They will be analyzed concurrently using the same configuration params.
    "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma",    # This name should match the folder name under 'prep' that contains the sample data
    ]
```

The second block corresponds to the pipelines that will be executed. As we plan to run and end-to-end analysis all will be set to *True**
```yaml
pipelines.enabled:
    imgseg: True
    cell2spot: True
    celldeconv: True
    cluster: True
    assign: True
    qc: True
    datamerge: True
    spatialanalysis: True
```

For the image segmentation block it is key to specify that we want to use **cellpose** and the **nuclei** model. The _flow_threshold_ parameter determines how aggressive the segmentation will be, a value of **0.8** will result in most of the cells being segmented. For a deatailed description fo the parameters please check [here](docs_md/parameters.md).
```yaml
imgseg:
    image.resolution: "" #TBD
    image.magnification: "" #TBD
    model:
        name: "cellpose"
        version: "2.1.1" # Shouldn't this be gone?
        params:
            patch_size: 512
            overlap: 496
            downsample_factor: 1
            n_channels: 3
            channels: [1,0]
            model_type: "nuclei"
            batch_size: 100
            diameter: 16
            flow_threshold: 0.8
```
For the cell deconvolution module, we need to specify **cell2location** as the deconvolution method. It is essential that the _atlas_type_ parmeter is aligned with the reference dataset name defined in step 5 (in this case, *luca*). Some of the key parameters to be customized are _st_max_epoches_, which will determine the maxium number of training epochs, and _sc_use_gpu_ which determines if SpatialOne will use a GPU or a cpu for the training.
For a deatailed description fo the parameters please check [here](docs_md/parameters.md).
```yaml
celldeconv:
    model:
        name: "cell2location"
        version: "0.1.3" # Shouldn't this be gone?
        params:
            seed : 2023
            #params for clean sc datasets
            sc_cell_count_cutoff : 20 #a int variable for gene filter. A gene should be detected over minimum of cells counts e.g. should be detected in over 20 cells
            sc_cell_percentage_cutoff2 : 0.05 #(0,1) float variable for gene filter. A gene should be detected over minimum percent of cells.
            sc_nonz_mean_cutoff : 1.12 # (1, ) float variable for gene filter. A gene should have mean expression across non-zero cells slightly larger than 1
            #params for training sc to get st signatures
            sc_batch_key: "batch" #  single cell data batch category, e.g. 10X reaction / sample / batch
            sc_label_key: "cell_type" # cell type, covariate used for constructing signatures
            sc_categorical_covariate_keys: [] # multiplicative technical effects (platform, 3' vs 5', donor effect)
            sc_max_epoches : 50
            sc_lr : 0.02
            sc_use_gpu : True
            #parames for training st to get cell aboundance
            st_N_cells_per_location : 20
            st_detection_alpha : 200
            st_max_epoches : 25000
            cell_aboundance_threshold : 0.1 #cell aboundance threshold set, reduce aboundance below threshold to 0. Default,  cell abundance >0.1 was used as a gold standard label
            atlas_type : 'luca'
            mode : ### model
            retrain_cell_sig: True ## if need to re-train a signature or a treained signature is exist for using
```
For the morphological cell clustering block it is key to define the number of expected cell types (e.g. *20*). In order to allow the *local* cell type estimation method to infer the cell types for the segmented cells, it is important to enable the _spot_clustering_ step, as it will cluster both the morphological in the whole slide as in the spot.
```yaml
cluster:
    model:
        name: "gaussian"
        version: "X.X.X"
        params:
            n_clusters: 20
            spot_clustering: True
assign:
    model:
        name: "local"
        version: "X.X.X" # Shouldn't this be gone?
```

The cell2spot, QC, and spatialanalysis modules  do not require any specific configuration:
```yaml
cell2spot:
    model:
        name: "default"
        version: "X.X.X"
qc:
    model:
        name: "default"
        version: "X.X.X"
spatialanalysis:
    model:
        name: "default"
        version: "X.X.X"
```
Finally, in the datamerge step the user can define if geojson annotations are available, as well as the number of top genes that will be considered in the downstream analyisis.
```yaml
datamerge:
    run_data_merge_only: True
    model:
        name: "default"
        version: "X.X.X"
        params:
            annotation_file: "annotations.geojson"
            spot_index: "barcode"
            cell_index: "cell_ids"
            target_genes: []
            n_top_expressed_genes: 750      # Forces to include the top 500 most expressed genes in the reporting
            n_top_variability_genes: 750
```

#### 7. Running the analysis
If the user has mantained the default config file name "_visium_config_flow.yaml_" (reccomended) the analysis can simply triggered by running the following command which will run the SpatialOne light version:
```sh
make run
```
or
```sh
make run-cpu
```

If the file name has been changed or the user wants to run the analysis using a customized set up the analyis can be triggered by running:
```sh
set -a
source .env
set +a
docker run --gpus device=${GPU_DEVICE_ID} -it -v ${HOST_DATA_PATH}:/app/data spatialone-pipeline
```

Alternatively, you can directly run the analysis without loading the .env environment file by replacing `GPU_DEVICE_ID` and `HOST_DATA_PATH` by their valures:
```sh
docker run --gpus device=0 -it -v /Users/demo_user/Documents/spatialone-pipeline/app/data spatialone-pipeline
```



[Back to README](../README.md)