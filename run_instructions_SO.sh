#!/bin/bash
export HOST_DATA_PATH=$(pwd)
export GPU_DEVICE_ID=0

# Clone and set up the repository
#git clone https://github.com/Sanofi-Public/spatialone-pipeline

#cd spatialone-pipeline
#git checkout staging
make build
#cd ..

export HOST_DATA_PATH=$(pwd)

# Initialize flag variables
download_data=1
use_gpu=1  # Default to GPU

# Parse script arguments for flags
for arg in "$@"
do
    case $arg in
        --no_data_download)
            download_data=0
            ;;
        --cpu)
            use_gpu=0
            ;;
        --gpu)
            use_gpu=1
            ;;
    esac
done

# Download the data unless --no_data_download is specified
if [ "$download_data" -eq 1 ]; then
    ./download_experiment_data.sh
    # Check the exit status of the previous command
    if [ $? -ne 0 ]; then
        echo "The download and check script failed. Stopping the execution."
        exit 27
    fi

fi

# Set the Docker run command based on the processing unit choice
if [ "$use_gpu" -eq 0 ]; then
    DOCKER_RUN_CMD="docker run -it -v ${HOST_DATA_PATH}:/app/data"
else
    DOCKER_RUN_CMD="docker run --gpus device=${GPU_DEVICE_ID} -it -v ${HOST_DATA_PATH}:/app/data"
fi

# Running the pipeline with the selected configuration
#configs=("lung_scc" "lung_nec" "mouse_brain" "mouse_kidney_v1" "mouse_kidney_FFPE")
configs=("mouse_brain" "mouse_kidney_v1" "mouse_kidney_FFPE" "lung_scc" "lung_nec")
for config in "${configs[@]}"; do
    echo "Starting execution of ${configs[@]}"
    if [ "$use_gpu" -eq 0 ]; then
        cp "./conf/cpu/visium_config_flow_${config}.yaml" ./conf/visium_config_flow.yaml
    else
        cp "./conf/gpu/visium_config_flow_${config}.yaml" ./conf/visium_config_flow.yaml
    fi
    $DOCKER_RUN_CMD spatialone-pipeline 2>&1 | tee -a "log_${config}.txt"
    echo "Finished execution of ${configs[@]}"
done

# Save results to a folder
mkdir -p ./results/ref_atlases
cp ./reference/*sig* ./results/ref_atlases/
mv "results" "results_$(date +%s)"
