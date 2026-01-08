from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd
from lmfit import Parameters
from pydantic import ConfigDict, Field, field_validator, model_validator

from ..core.validators import BaseDictionary
from ..geostatistics.variogram_models import get_variogram_arguments
from . import spatial


class VariogramModelParameters(BaseDictionary):
    correlation_range: Optional[float] = Field(default=None, gt=0.0, allow_inf_nan=False)
    decay_power: Optional[float] = Field(default=None, gt=0.0, le=2.0, allow_inf_nan=False)
    enhance_semivariance: Optional[bool] = Field(default=None)
    hole_effect_range: Optional[float] = Field(default=None, ge=0.0, allow_inf_nan=False)
    sill: Optional[float] = Field(default=None, gt=0.0, allow_inf_nan=False)
    nugget: Optional[float] = Field(default=None, ge=0.0, allow_inf_nan=False)
    smoothness_parameter: Optional[float] = Field(
        default=None, gt=0.0, le=10.0, allow_inf_nan=False
    )
    shape_parameter: Optional[float] = Field(default=None, gt=0.0, le=100.0, allow_inf_nan=False)
    power_exponent: Optional[float] = Field(default=None, gt=0.0, lt=2.0, allow_inf_nan=False)
    model_config = ConfigDict(title="theoretical variogram model parameters")

    @field_validator("decay_power", mode="before")
    def validate_decay_power(cls, v):
        if v is None:
            return v

        if not 0 < v <= 2:
            raise ValueError(
                f"decay_power must be in interval (0, 2]. Got {v}. Values > 2 create "
                f"non-stationary processes."
            )
        return v

    @field_validator("power_exponent", mode="before")
    def validate_power_exponent(cls, v):
        if v is None:
            return v

        if not 0.0 < v < 2.0:
            raise ValueError(
                f"power_exponent must be in interval (0, 2). Got {v}. Must be strictly less than 2 "
                f"for valid spatial covariance structure."
            )
        return v

    @field_validator("smoothness_parameter", mode="before")
    def validate_smoothness_parameter(cls, v):
        if v is None:
            return v

        if not 0.0 < v <= 10.0:
            raise ValueError(
                f"smoothness_parameter must be in interval (0, 10]. Got {v}. Common values: 0.5 "
                f"(exponential), 1.5 (once differentiable), 2.5 (twice differentiable)."
            )
        return v

    @field_validator("shape_parameter", mode="before")
    def validate_shape_parameter(cls, v):
        if v is None:
            return v

        if not 0.0 < v <= 100.0:
            raise ValueError(f"shape_parameter must be in interval (0, 100]. Got {v}.")
        return v

    @model_validator(mode="after")
    def validate_sill_nugget_relationship(self):
        if self.sill is not None and self.nugget is not None:
            if self.sill <= self.nugget:
                raise ValueError(
                    f"sill ({self.sill}) must be greater than nugget ({self.nugget}). "
                    f"The partial sill (sill - nugget) must be positive."
                )
        return self


class ValidateVariogramClass(BaseDictionary):
    coordinate_names: Tuple[str, str]
    lag_resolution: float = Field(gt=0.0, allow_inf_nan=False)
    n_lags: int = Field(gt=0)
    model_config = ConfigDict(title="variogram analysis parameters")


class ValidateEmpiricalVariogramArgs(BaseDictionary):
    azimuth_angle_threshold: float = Field(ge=0.0, le=180.0, allow_inf_nan=None)
    azimuth_filter: bool
    coordinate_names: Tuple[str, str]
    data: pd.DataFrame
    force_lag_zero: bool
    variable: str
    model_config = ConfigDict(title="empirical variogram parameters")

    @field_validator("data", mode="after")
    def validate_transects(cls, v):
        # Validate with pandera
        return spatial.TransectsDF.validate(v)

    @model_validator(mode="after")
    def validate_df_columns(self):
        # Get the mesh and transects DataFrames
        coords = self.coordinate_names
        data = self.data
        variable = self.variable

        # Initialize error message
        error_msg = "The input DataFrame is missing the defined column(s) for "
        coord_flag = False
        variable_flag = False

        # Check for coordinate columns
        coord_check = set(coords) <= set(data.columns)
        # ---- Add to error message
        if not coord_check:
            coord_flag = True
            error_msg += f"coordinates {coords}"

        # Check for the variable column
        variable_check = variable in data.columns
        # ---- Add to error message
        if not variable_check:
            variable_flag = True
            if coord_flag:
                error_msg += " and "
            error_msg += f"variable ('{variable}')"

        # Raise error, if needed
        if coord_flag or variable_flag:
            raise KeyError(error_msg + ".")

        return self


class ValidateFitVariogramArgs(BaseDictionary):
    model: Union[str, List[str]] = Field(union_mode="left_to_right")
    model_parameters: Parameters
    optimizer_kwargs: Dict[str, Any]
    model_config = ConfigDict(
        title="theoretical variogram fitting parameters", protected_namespaces=()
    )

    @field_validator("model", mode="before")
    def validate_model(cls, v):
        # Check applicable models
        try:
            _ = get_variogram_arguments(v)
        except Exception as e:
            e.args = (
                str(e) + " Check the documentation for `variogram` for allowed variogram models.",
            )
            raise
        return v
