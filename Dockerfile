FROM nvidia/cuda:11.8.0-base-ubuntu20.04

# Install Python 3.9
RUN apt-get clean && rm -rf /var/lib/apt/lists/* && rm -rf /var/cache/apt/archives/lock && rm -rf /var/lib/dpkg/lock \
    && apt-get update \
    # && apt-get install libatlas-base-dev liblapack-dev libblas-dev -y \
    && apt-get install -y python3.9 python3.9-venv python3.9-dev \
    && apt-get install -y python3-pip

# Upgrade pip to the latest version
RUN python3.9 -m pip install --upgrade pip

# Optionally set Python 3.9 as the default Python version
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.9 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.9 1 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

#Set the working directory to "/app" in the Docker container
WORKDIR /app
# Set directory for /data to access conf/,prep/,results/,reference/ data from external volume
RUN mkdir /data
COPY requirements.txt /app/
COPY Makefile /app/
COPY .env /app/
RUN make install
COPY src/ /app/src/
COPY conf/ /app/conf/
COPY logging.yml /app/
COPY templates/ /app/templates/

# Specify what port are exposed in your application
# For further documentation, refer https://docs.docker.com/engine/reference/builder/#expose
EXPOSE 80

# Execute the uvicorn command to run the application on port 80
ENV PYTHONPATH "${PYTHONPATH}:/app"
CMD ["python3", "-m", "src.pipelines.visium_flow", "run"]
