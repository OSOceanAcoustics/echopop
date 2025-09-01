from typing import Literal

import numpy as np
from pydantic import Field

from ..core.validators import BaseDictionary


class ValidatePCDWBAParams(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="PCDWBA model parameters",
):
    """
    Validation schema for Phase-Corrected Distorted Wave Born Approximation (PCDWBA) model
    parameters.
    """

    g: float = Field(gt=0.0, le=np.inf, allow_inf_nan=True)
    h: float = Field(gt=0.0, le=np.inf, allow_inf_nan=True)
    length_mean: float = Field(gt=0.0, le=np.inf, allow_inf_nan=True)
    length_radius_ratio: float = Field(gt=0.0, le=np.inf, allow_inf_nan=True)
    length_sd_norm: float = Field(ge=0.0, le=np.inf, allow_inf_nan=True)
    number_density: float = Field(gt=0.0, le=np.inf, allow_inf_nan=True)
    radius_of_curvature_ratio: float = Field(gt=0.0, le=np.inf, allow_inf_nan=True)
    theta_mean: float = Field(ge=-180.0, le=180.0, allow_inf_nan=False)
    theta_sd: float = Field(ge=0.0, le=np.inf, allow_inf_nan=True)


class DistributionParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="distribution parameters",
):
    """
    Configuration parameters for probability distributions used in scattering models.
    """

    bins: int = Field(default=30)
    family: Literal["gaussian", "uniform"] = Field(default="gaussian")


class ValidatePCDWBASettings(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="PCDWBA model settings and additional arguments",
):
    """
    Validation schema for PCDWBA model computational settings and configuration.
    """

    frequency_interval: float = Field(gt=0.0, allow_inf_nan=False)
    length_distribution: DistributionParameters = Field(default_factory=DistributionParameters)
    n_integration: int = Field(default=50, gt=0)
    n_wavelength: int = Field(default=10, gt=0)
    orientation_distribution: DistributionParameters = Field(default_factory=DistributionParameters)
    taper_order: float = Field(default=10.0, gt=0.0, allow_inf_nan=False)
    type: str
