# Spatial One Flow

This repository contains the code and dockefiles to replicate results present in the paper: Spatial One. It can be downloaded and used to run end-end spatial transcriptomics analysis on a local or remote machine (cloud agnostic).

# Hardware requirements

The requirements on the hardware are based on the type and size of input data available. For the data presented in the publication, we present the following recommended hardware specifications:

16 CPU, 120GB RAM, 1 GPU

# Software requirements

The requirements for the packages and environments are present in the respective requirements and docker files for each container. The containers are set to operate in R and Python and GPU access is enabled via CUDA 11.

# Data preparation
Choose a single folder, where all your data and configurations will stay, as well as future outputs. We recommend the following structure and exact file naming format to ensure compliance with the main pipeline:

```
data/
├── conf/
│   └── visium_config_flow.yaml
├── prep/
│   ├── sample1/
│   │   └── #input files
│   └── sample2/
│       └── #input files
├── reference/
│   ├── # input reference data
│   ├── luca.h5ad
│   └── ...
└── results/
    └── #populated by the pipeline
```

**visium_config_flow.yaml** - contains all configuration parameters such as sample_id, pipelines to run or specific run parameters. Example can be found under the conf/ folder of this repository.

**prep/** - this is where all sample specific input data lies. Please separate each sample with own folder.

**reference/** - this is where non-sample specific input data lies, for example reference datasets or pretrained models. Pipeline expects the following data to be present:

 - reference atlas
 - cell signature files - (can be added automatically when retrain_cell_sig: True is set)
 - pretrained model (if using hovernet) - default model can be found [here for download](https://zenodo.org/records/10854151/files/hovernet_original_consep_type_tf2pytorch.tar?download=1)
 - label info (if using hovernet) - default model can be found [here for download](https://zenodo.org/records/10854151/files/type_info.json?download=1)

# Docker setup (GPU only)
Run following commands to setup the docker environment:

(Optional) In case you are using proxies export them as environmental variables:

    export HTTP_PROXY
    export HTTPS_PROXY

To build the needed docker containers:

    export HOST_DATA_PATH = {path to your folder with data}

    make docker-build

or (for cpu only machines)

**Disclaimer**: hovernet image segmentation not supported in this mode

    make docker-build-cpu

This step can take up to 20 minutes.

## Run an experiment
Make sure the config in HOST_DATA_PATH/conf/visium_config_flow is correctly setup.

To trigger the execution, run:

    make docker-start

or (for cpu oncly machines)
    make docker-start-cpu

You will then see the progress logs in the terminal and outputs will be produced to HOST_DATA_PATH/results/{sample_id}.

### Run single container spatialone only (without CARD or Hovernet) - supports CPU

edit variables set in makefile if needed, then run:

      make build
      make run #gpu version
      make run-cpu #cpu only version
