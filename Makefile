include .env
export

setup-env:
	conda install -y black
	conda install -y pre-commit
	pre-commit install
	pip install scalene

install-local:
	sudo apt-get update && sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev -y
	sudo apt-get install -y libgl1-mesa-glx
	sudo apt-get install -y libglib2.0-0
	python3 -m pip install torch==2.3.0 torchvision==0.18.0 torchaudio==2.3.0 --index-url https://download.pytorch.org/whl/cu118
	python3 -m pip install -r requirements.txt --no-cache-dir --user

install:
	apt-get update && apt-get install libatlas-base-dev liblapack-dev libblas-dev -y
	apt-get install -y libgl1-mesa-glx
	apt-get install -y libglib2.0-0
	python3 -m pip install torch==2.3.0 torchvision==0.18.0 torchaudio==2.3.0 --index-url https://download.pytorch.org/whl/cu118
	python3 -m pip install -r requirements.txt --no-cache-dir --user

install-arm:
	apt-get update && apt-get install libatlas-base-dev liblapack-dev libblas-dev -y
	apt-get install -y libgl1-mesa-glx
	apt-get install -y libglib2.0-0
	python3 -m pip install -r requirements-arm.txt --no-cache-dir --user


install-conda:
	sudo apt-get update -y
	sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev -y
	python3 -m pip install -r requirements.txt --no-cache-dir --user

pre-commit:
	poetry run pre-commit run --all-files

lint:
	black src/
	isort src/
	pylint --fail-under=7 src/

checks:
	black --check src/
	isort src/* --check-only
	pylint --fail-under=7 src/

build:
	DOCKER_BUILDKIT=1 docker build \
		--build-arg http_proxy=${HTTP_PROXY} \
		--build-arg https_proxy=${HTTPS_PROXY} \
		--build-arg no_proxy=${NO_PROXY} \
		--build-arg HTTP_PROXY=${HTTP_PROXY} \
		--build-arg HTTPS_PROXY=${HTTPS_PROXY} \
		--build-arg NO_PROXY=${NO_PROXY} \
		--platform linux/amd64 \
		-t spatialone-pipeline .
build-amd:
	DOCKER_BUILDKIT=1 docker build \
		--file Dockerfile_amd \
		--build-arg http_proxy=${HTTP_PROXY} \
		--build-arg https_proxy=${HTTPS_PROXY} \
		--build-arg no_proxy=${NO_PROXY} \
		--build-arg HTTP_PROXY=${HTTP_PROXY} \
		--build-arg HTTPS_PROXY=${HTTPS_PROXY} \
		--build-arg NO_PROXY=${NO_PROXY} \
		--platform linux/amd64 \
		-t spatialone-pipeline .
build-arm:
	DOCKER_BUILDKIT=1 docker build \
		--file Dockerfile_arm \
		--build-arg http_proxy=${HTTP_PROXY} \
		--build-arg https_proxy=${HTTPS_PROXY} \
		--build-arg no_proxy=${NO_PROXY} \
		--build-arg HTTP_PROXY=${HTTP_PROXY} \
		--build-arg HTTPS_PROXY=${HTTPS_PROXY} \
		--build-arg NO_PROXY=${NO_PROXY} \
		--platform linux/arm64 \
		-t spatialone-pipeline .

run:
	docker run --gpus device=${GPU_DEVICE_ID} -it -v ${HOST_DATA_PATH}:/app/data -it spatialone-pipeline

run-cpu:
	docker run -it -v ${HOST_DATA_PATH}:/app/data -it spatialone-pipeline

test:
	PYTHONPATH=${PYTHONPATH}:.:src:tests \
	pytest --cov=src/ --cov-fail-under=70 tests/unit_tests

clean:
	find ./ -name "*~" | xargs rm -v || :
	find ./ -name "*.pyc" | xargs rm -v || :
	find ./ -name "__pycache__" | xargs rm -rf || :

run-spatialone:
	python3 -m src.pipelines.visium_flow --no-pylint run

docker-build:
# needs HOST_DATA_PATH variable
	docker compose build

docker-start:
	docker compose up

docker-build-cpu:
	docker compose -f docker-compose-cpu.yml build

docker-start-cpu:
	docker compose -f docker-compose-cpu.yml up
