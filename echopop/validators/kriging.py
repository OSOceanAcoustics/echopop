from typing import Any, Callable, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from pydantic import ConfigDict, Field, field_validator, model_validator
from pydantic_core import PydanticCustomError

from ..core.validators import BaseDictionary
from . import spatial
from .variogram import VariogramModelParameters


class KrigingParameters(BaseDictionary):
    aspect_ratio: float = Field(default=1e-3, gt=0.0, le=1.0, allow_inf_nan=False)
    k_min: int = Field(default=5, ge=2)
    k_max: int = Field(default=20)
    search_radius: float = Field(gt=0.0, allow_inf_nan=False)
    model_config = ConfigDict(title="kriging model parameters")

    @model_validator(mode="after")
    @classmethod
    def validate_k_interval(cls, values):
        # Get `k_min` and `k_max`
        k_min = getattr(values, "k_min", 5)
        k_max = getattr(values, "k_max", 20)

        # Ensure that the interval is sensible
        if k_min > k_max:
            # ---- Raise Error
            raise ValueError(
                f"Define 'k_max' ({k_max}) must be greater than or equal to 'k_min' ({k_min})."
            )

        return values


class VariogramKrigeModelParameters(VariogramModelParameters):
    model: Union[str, List[str]] = Field(union_mode="left_to_right")
    model_config = ConfigDict(title="theoretical variogram model parameters")

class ValidateKrigingClass(BaseDictionary):
    mesh: pd.DataFrame
    kriging_params: KrigingParameters
    variogram_params: VariogramKrigeModelParameters
    coordinate_names: Tuple[str, str]
    model_config = ConfigDict(title="kriging analysis parameters")

    @field_validator("mesh", mode="after")
    def validate_mesh(cls, v):
        # Validate with pandera
        return spatial.MeshDF.validate(v)


class ValidateMeshCropArgs(BaseDictionary):
    crop_function: Callable
    coordinate_names: Tuple[str, str]
    kwargs: Dict[str, Any]
    model_config = ConfigDict(title="mesh cropping parameters for kriging")
