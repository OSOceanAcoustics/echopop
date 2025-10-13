from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from pydantic import Field, field_validator, model_validator
from pydantic_core import PydanticCustomError

from ..core.validators import BaseDictionary
from . import spatial


class KrigingParameters(
    BaseDictionary, arbitrary_types_allowed=True, title="kriging model parameters"
):
    aspect_ratio: float = Field(default=1e-3, gt=0.0, le=1.0, allow_inf_nan=False)
    k_min: int = Field(default=5, ge=2)
    k_max: int = Field(default=20)
    search_radius: float = Field(gt=0.0, allow_inf_nan=False)

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


class VariogramKrigeModelParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="theoretical variogram model parameters",
):
    model: Union[str, List[str]] = Field(union_mode="left_to_right")
    lags: np.ndarray = Field(default=None)
    sill: Optional[float] = Field(default=None, gt=0.0, allow_inf_nan=False)
    nugget: Optional[float] = Field(default=None, ge=0.0, allow_inf_nan=False)
    hole_effect_range: Optional[float] = Field(default=None, ge=0.0, allow_inf_nan=False)
    correlation_range: Optional[float] = Field(default=None, gt=0.0, allow_inf_nan=False)
    enhance_semivariance: Optional[bool] = Field(default=None)
    decay_power: Optional[float] = Field(default=None, gt=0.0, allow_inf_nan=False)

    @field_validator("lags", mode="before")
    def validate_lags(cls, v):

        # If a list, coerce to an array
        if isinstance(v, list):
            v = np.array(v)

        # If an array, check if coercible to an array of floats
        if isinstance(v, np.ndarray):
            try:
                return v.astype(float)
            except:
                raise PydanticCustomError(
                    "invalid_array_dtype",
                    "Values within the array must all be coercible to floats",
                )
        return v


class ValidateKrigingClass(
    BaseDictionary, arbitrary_types_allowed=True, title="kriging analysis parameters"
):
    mesh: pd.DataFrame
    kriging_params: KrigingParameters
    variogram_params: VariogramKrigeModelParameters
    coordinate_names: Tuple[str, str]

    @field_validator("mesh", mode="after")
    def validate_mesh(cls, v):
        # Validate with pandera
        return spatial.MeshDF.validate(v)


class ValidateMeshCropArgs(
    BaseDictionary, arbitrary_types_allowed=True, title="mesh cropping parameters for kriging"
):
    crop_function: Callable
    coordinate_names: Tuple[str, str]
    kwargs: Dict[str, Any]
