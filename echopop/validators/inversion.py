from pydantic import field_validator, Field
from typing import List, Optional
import numpy as np

from .base import BaseDictionary

class TSLRegressionParameters(BaseDictionary, title="TS-length regression parameters"):
    """
    Target strength - length regression parameters

    Parameters
    ----------
    slope : float
        TS-length regression slope.
    intercept : float
        TS-length regression intercept.
    """

    slope: float = Field(allow_inf_nan=False)
    intercept: float = Field(allow_inf_nan=False)

class ValidateLengthTS(
    BaseDictionary, 
    arbitrary_types_allowed=True, 
    title="TS-length inversion model parameters"
):
    """
    Validation model for TS-length inversion parameters used by InversionLengthTS.

    This Pydantic model validates and documents the configuration parameters required by the 
    InversionLengthTS class for acoustic inversion using length-TS regression.

    Parameters
    ----------
    ts_length_regression : utils.TSLRegressionParameters
        Regression parameters for converting fish length to target strength (TS).
    stratify_by : List[str]
        List of column names used for data stratification (e.g., 'stratum_ks').
    expected_strata : np.ndarray, optional
        Array of expected strata identifiers to process.
    impute_missing_strata : bool, default=True
        Whether to impute missing strata values during inversion.
    haul_replicates : bool, default=True
        Whether to use hauls as the statistical unit of replication (recommended to avoid
        pseudoreplication[1]_).

    Notes
    -----
    This model is intended for use with the InversionLengthTS class, which performs
    acoustic inversion by relating fish length to acoustic backscatter using empirical
    TS-length relationships.

    Using hauls as the unit of replication is recommended to avoid pseudoreplication when
    individuals within hauls are not independent[1]_.

    References
    ----------
    .. [1] Hurlbert, S.H. (1984). Pseudoreplication and the Design of Ecological Field Experiments.
       *Ecological Monographs*, 54(2), 187-211. https://doi.org/10.2307/1942661
    """

    ts_length_regression: TSLRegressionParameters
    stratify_by: List[str]
    expected_strata: Optional[np.ndarray[np.number]] = Field(default=None)
    impute_missing_strata: bool = Field(default=True)
    haul_replicates: bool = Field(default=True)

    @field_validator("stratify_by", mode="before")
    def validate_stratify_by(cls, v):
        if isinstance(v, str):
            v = [v]
        return v

    @field_validator("expected_strata", mode="before")
    def validate_expected_strata(cls, v):
        if v is None:
            return v

        if isinstance(v, list):
            return np.array(v)

        return v