# Running a SpatialOne Analysis using CARD & Hovernet

For this tutorial, we assume the user is using the following path and has Docker properly set up: `/home/user/card_demo`. It is also assumed that the user has access to a machine with a GPU and the right memory set up.

## 1. Retrieve the Code
Donwload the code from the github repository
```bash
git clone https://github.com/Sanofi-Public/spatialone-pipeline.git
cd spatialone-pipeline
```

## 2. Download Required Data

### 2.1. Sample Data
Download the two cancer samples, the reference data for deconvolution, and the default configuration files.
```bash
./download_experiment_data.sh
```

### 2.2. Downlaod Hovernet Weights
Hovernet can work with different weights, as described in their [official repository]([https://drive.google.com/drive/folders/17IBOqdImvZ7Phe0ZdC5U1vwPFJFkttWp](https://github.com/vqdang/hover_net?tab=readme-ov-file#model-weights)). To use the CONSEP weights and labels described in the SpatialOne paper you should download them on the reference folder.

You can manually download them from the [author's site](https://drive.google.com/drive/folders/17IBOqdImvZ7Phe0ZdC5U1vwPFJFkttWp), from Zenodo ([weights](https://zenodo.org/records/12801948/files/hovernet_original_consep_type_tf2pytorch.tar?download=1), [labels](https://zenodo.org/records/12801948/files/type_info.json?download=1)) or you can execute the following two instructions:
```bash
curl -L -o reference/hovernet_original_consep_type_tf2pytorch.tar "https://zenodo.org/records/12801948/files/hovernet_original_consep_type_tf2pytorch.tar?download=1"
curl -L -o reference/type_info.json "https://zenodo.org/records/12801948/files/type_info.json?download=1"
```


## 3. Set Up the Environment
Prepare your environment to build the container images if necessary. To that end, edit the `.env` file:
```env
HOST_DATA_PATH = /home/user/card_demo/spatialone-pipeline
GPU_DEVICE_ID = 0
```

Export proxies only if necessary:
```bash
export HTTP_PROXY=http://your.proxy.url:port
export HTTPS_PROXY=http://your.proxy.url:port
```

## 4. Build the Docker Image

Build the docker containers:
```bash
make docker-build
```
Note that this operation may take up to 20 minutes

## 5. Edit the Configuration Files

Edit the config file located at `$HOST_DATA_PATH/conf/visium_config_flow.yaml` (in our case `./conf/visium_config_flow.yaml`) to reflect an appropriate configuration. You'll need to define the pipelines to run, and set up Hovernet & CARD.

### 5.1 Enable All Pipelines

Ensure that all the pipelines are enabled:
```yaml
pipelines.enabled:
    imgseg: True     # cell segmentation
    cell2spot: True  # matching cells to visium spots
    celldeconv: True # cell deconvoluiton
    cluster: True    # morphological clustering
    assign: True     # cell assignment integrates celldeconvolution with cell segmentation
    qc: True         # QC metrics generation
    datamerge: True  # To visualize in Tissuumaps enable "datamerge: true"
    spatialanalysis: True # Spatial analysis reporting
```

### 5.2 Enable Hovernet Configuration

Remove the `cellpose` configuration section and uncomment the `hovernet` section:
```yaml
imgseg:
    image.resolution: "400dpi"
    image.magnification: "20x"
    model:
        name: "hovernet"
        version: "0.0.1"
        params:
            gpu: '0'
            nr_types: 5
            batch_size: 16
            model_mode: "original"
            nr_inference_workers: 1
            nr_post_proc_workers: 1
            wsi: True
            tile: False
            help: False
            model_name: "hovernet_original_consep_type_tf2pytorch.tar"
            type_filename: "type_info.json"
```
Please note that if you plan to use a different set of weights than _hovernet_original_consep_type_tf2pytorch_ you should update the `model_name` and the `type_filename` parameters.

### 5.3 Set Up CARD Configuration

Remove the `cell2location` section and uncomment the `CARD` block:
```yaml
celldeconv:
    model:
        name: "card"
        version: "X.X.X"
        params:
            min_count_gene: 100
            min_count_spot: 5
            ct_varname: 'cellType'
            ct_select: 'NULL'
            atlas_type: 'luca'
            sc_label_key: 'cell_type'
            sc_sample_key: 'batch'
```
Ensure to update the `atlas_type` parameter to `luca` to match the reference dataset with the type of tissue analyzed in the example.

## 6. Run the Analysis

Run the analysis:
```bash
make docker-start
```

By following these steps, you will set up and run a SpatialOne analysis using CARD and Hovernet.
