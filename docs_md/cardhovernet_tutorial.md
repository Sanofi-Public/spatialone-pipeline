# Running a SpatialOne Analysis using CARD & Hovernet

For this tutorial, we assume the user is using the following path and has Docker properly set up: `/home/user/card_demo`. It is also assumed that the user has access to a machine with a GPU and the right memory set up.

## 1. Retrieve the Code

```bash
git clone https://github.com/Sanofi-Public/spatialone-pipeline.git
cd spatialone-pipeline
```

## 2. Download Sample Data

```bash
./download_experiment_data.sh
```

## 3. Set Up the Environment

Edit the `.env` file if necessary:
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

Build Docker containers:
```bash
make docker-build
```

## 5. Edit the Configuration Files

Edit the config file located at `./conf/visium_config_flow.yaml` to reflect an appropriate configuration.

### 5.1 Enable All Pipelines

Ensure that all the pipelines are enabled:
```yaml
pipelines.enabled:
    imgseg: False     # cell segmentation
    cell2spot: False  # matching cells to visium spots
    celldeconv: False # cell deconvoluiton
    cluster: False    # morphological clustering
    assign: False     # cell assignment integrates celldeconvolution with cell segmentation
    qc: False         # QC metrics generation
    datamerge: False  # To visualize in Tissuumaps enable "datamerge: true"
    spatialanalysis: False # Spatial analysis reporting
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