""" api.py
"""
from typing import Dict

from fastapi import FastAPI
from pydantic import BaseModel, Field

from src.hovernet import HoverNet
from src.utils.dataio_pipeline import DataIO
from src.utils.logger import Logger

logger = Logger()

app = FastAPI(debug=True)

# to run in terminal: python -m uvicorn src.api:app --reload  --proxy-headers --host "0.0.0.0"  --port 80 --root-path .


class PredictionRequest(BaseModel):
    """PredictionRequest Class with parameters for prediction

    Args:
        BaseModel (_type_): _description_
    """

    configs: Dict
    exp_id: str
    input_path: str
    reference_path: str
    output_path: str
    model_filename: str = "hovernet_original_consep_type_tf2pytorch.tar"
    type_filename: str = "type_info.json"


@app.get("/hello")
async def root() -> Dict[str, str]:
    """Root route

    Returns:
        Dict[str, str]: Return simple hello world message
    """
    return {"message": "Hello World"}


@app.post("/predict")
async def predict(request: PredictionRequest) -> Dict[str, str]:
    """Route to Predict

    Returns:
        str: status message
    """
    # logger.info(f"Request: {request}")
    logger.info(f"configs: {request.configs}")
    logger.info(f"exp_id: {request.exp_id}")
    logger.info(request)
    exp_id = request.exp_id
    image_seg = HoverNet(request.configs)
    data_io = DataIO(exp_id)
    logger.info(f"[{image_seg.pipeline_name}] Start pipeline")

    logger.info(f"Running HoverNet Inference...")

    image_seg.predict(
        exp_id,
        request.input_path,
        request.output_path,
        data_io.cache_dir,
        request.reference_path,
        request.model_filename,
        request.type_filename,
    )

    # Return the prediction as a JSON response
    return {"prediction status": "completed"}
