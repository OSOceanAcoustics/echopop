from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional

import numpy as np
import pandas as pd
from pydantic import ConfigDict, Field, RootModel, ValidationError, field_validator, model_validator

from ..core.validators import BaseDataFrame, BaseDictionary

if TYPE_CHECKING:
    from ..inversion import InvParameters


class TSLRegressionParameters(BaseDictionary):
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
    model_config = ConfigDict(title="TS-length regression parameters")


class ValidateLengthTS(BaseDictionary):
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
    model_config = ConfigDict(title="TS-length inversion model parameters")

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


class ScatterDF(BaseDataFrame):
    """
    Validates MultiIndex DataFrame structure with required top-level columns and nested frequency
    index.

    Ensures the top-level column names include 'sv_mean', 'nasc', and 'thickness_mean', and that
    the second-level column names ('frequency') are numeric. All values under each frequency must
    be floats.
    """

    class Config(BaseDataFrame.Config):
        title = "acoustic backscatter export DataFrame"
        multiindex_strict = True
        multiindex_coerce = True
        unique_column_names = False

    @classmethod
    def pre_validate(cls, df: pd.DataFrame) -> pd.DataFrame:
        # Raise Error if not a MultiIndex DataFrame
        if not isinstance(df.columns, pd.MultiIndex):
            raise TypeError("Expected MultiIndex columns.")

        # Check for required top-level columns
        # ---- Get the top level indices
        top_level = set(df.columns.get_level_values(0))
        # ---- Identify missing columns
        missing_columns = {"sv_mean", "nasc", "thickness_mean"} - top_level
        # ---- Raise Error, if needed
        if missing_columns:
            # ---- Format
            missing_str = ", ".join(f"'{col}'" for col in missing_columns)
            raise KeyError(f"Missing required top-level columns: {missing_str}.")

        # Check for nested column 'frequency'
        if "frequency" not in df.columns.names:
            # ---- Get current column names
            current_names = list(df.columns.names)
            # ---- Format
            current_str = ", ".join(f"'{col}'" for col in current_names if col)
            # ---- Raise Error
            raise KeyError(
                f"Missing required nested column index 'frequency'. Current column indices are: "
                f"{current_str}."
            )

        return df


class EnvironmentParameters(BaseDictionary):
    """
    Environmental and medium parameters for acoustic scattering models.
    """

    sound_speed_sw: float = Field(default=1500.0, gt=0.0, allow_inf_nan=False)
    density_sw: float = Field(default=1025.0, gt=0.0, allow_inf_nan=False)
    model_config = ConfigDict(title="scattering model environment/medium parameters")


class SimulationParameters(BaseDictionary):
    """
    Configuration parameters for inversion simulation and optimization.
    """

    iter_cb: Optional[Callable] = Field(default=None)
    monte_carlo: Optional[bool] = Field(default=None)
    mc_realizations: Optional[int] = Field(default=None)
    mc_seed: Optional[int] = Field(default=None)
    scale_parameters: bool = Field(default=True)
    minimum_frequency_count: int = Field(default=2)
    model_config = ConfigDict(title="inversion simulation parameters")

    @model_validator(mode="after")
    def validate_monte_carlo_realizations(self):
        # Get `monte_carlo`
        monte_carlo = self.monte_carlo or False

        # Default n_realizations if monte_carlo
        if monte_carlo and self.mc_realizations is None:
            self.mc_realizations = 100

        # Always update monte_carlo
        self.monte_carlo = monte_carlo

        return self


class ValidateInversionMatrix(BaseDictionary):
    """
    Validation model for InversionMatrix class configuration.
    """

    data: pd.DataFrame
    simulation_settings: SimulationParameters = Field(default_factory=SimulationParameters)
    model_config = ConfigDict(title="scattering model inversion analysis")

    @field_validator("data", mode="after")
    def validate_data(cls, v):
        # Validate with pandera
        return ScatterDF.validate(v)


class ModelSettingsParameters(BaseDictionary):
    """
    Configuration parameters for acoustic scattering models.
    """

    type: str
    environment: EnvironmentParameters = Field(default_factory=EnvironmentParameters)
    model_config = ConfigDict(title="generate inversion model setting", extra="allow")


class ValidateBuildModelArgs(BaseDictionary):
    """
    Validation model for scattering model build arguments.
    """

    model_parameters: "InvParameters"
    model_settings: ModelSettingsParameters
    model_config = ConfigDict(
        title="scattering model preparation and assembly", protected_namespaces=()
    )

    @staticmethod
    def _check_variable(params: Dict[str, Any], variable: str, validator: Any):

        # Get the schema JSON
        properties = validator.model_json_schema()["properties"]

        # Get fill variable if needed
        if variable in ["min", "max"]:
            # ---- Format variable
            var = f"{variable}imum"
            var_ex = f"exclusive{variable.capitalize()}imum"
            # ---- Get variable defaults for boundaries
            var_defaults = {
                k: (v[var] if var in v else v[var_ex] + np.finfo(float).eps)
                for k, v in properties.items()
                if var in v or var_ex in v
            }
        else:
            var_defaults = {k: None for k in properties}

        # Set up the values
        param_slice = {
            k: v[variable] if variable in v else var_defaults[k] for k, v in params.items()
        }

        # Get the validation scheme
        return {k: {variable: v} for k, v in validator.create(**param_slice).items()}

    @model_validator(mode="after")
    def validate_model_parameterization(self):
        from ..inversion import InvParameters
        from .scattering_models import SCATTERING_MODEL_PARAMETERS

        # Check for model-type
        # ---- Dump the model
        model_settings = self.model_settings.model_dump()
        # ---- Check against reference
        if model_settings["type"] not in SCATTERING_MODEL_PARAMETERS:
            raise LookupError(f"Scattering model '{model_settings['type']}' could not be found.")

        # Get the model-specific validator
        model_validators = SCATTERING_MODEL_PARAMETERS[model_settings["type"]]

        # Validate `model_parameters`
        # ---- Dump the model [NOTE: this runs a check on the UNSCALED parameters]
        model_parameters = self.model_parameters._unscaled_parameters
        # ---- Check minimum values
        min_params = self._check_variable(model_parameters, "min", model_validators["parameters"])
        # ---- Set up values
        value_params = self._check_variable(
            model_parameters, "value", model_validators["parameters"]
        )
        # --- Set up maximum values
        max_params = self._check_variable(model_parameters, "max", model_validators["parameters"])
        # ---- Get the `vary` definitions
        vary_params = {k: {"vary": v["vary"]} for k, v in model_parameters.items()}
        # ---- Rebuild the model parameters
        validated_parameters = InvParameters(
            ModelInputParameters.create(
                **{
                    k: {
                        **min_params.get(k, {}),
                        **value_params.get(k, {}),
                        **max_params.get(k, {}),
                        **vary_params.get(k, {}),
                    }
                    for k in (
                        set(min_params) | set(value_params) | set(max_params) | set(vary_params)
                    )
                }
            )
        )
        # ---- Restore scaling, if required
        if self.model_parameters.is_scaled:
            validated_parameters.scale()
        # ---- Update the parameter set
        self.model_parameters = validated_parameters

        # Validate `model_settings`
        self.model_settings = ModelSettingsParameters.model_validate(
            model_validators["settings"].create(**self.model_settings.model_dump())
        )

        return self


class SingleParameter(BaseDictionary):
    min: Optional[float] = Field(default=-np.inf, allow_inf_nan=True)
    value: float = Field(allow_inf_nan=False)
    max: Optional[float] = Field(default=np.inf, allow_inf_nan=True)
    vary: bool = Field(default=False)
    model_config = ConfigDict(title="values for lmfit.Parameters class required for optimization")

    @model_validator(mode="after")
    def check_bounds(self):
        if self.min > self.max:
            raise ValueError(f"min ({self.min}) cannot be greater than max ({self.max}).")
        if not (self.min <= self.value <= self.max):
            raise ValueError(
                f"value ({self.value}) must be between min ({self.min}) and max ({self.max})."
            )
        return self


class ModelInputParameters(RootModel[Dict[str, "SingleParameter"]]):

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        title="scattering model parameters",
    )

    @model_validator(mode="before")
    @classmethod
    def prevalidator_trans(cls, data):
        # Use TYPE_CHECKING import
        # if TYPE_CHECKING:
        #     from ..inversion import InvParameters  # noqa: F401

        # Check type using string comparison to avoid import
        if type(data).__name__ == "InvParameters":
            from ..inversion import InvParameters  # noqa: F401

            return data.parameters
        return data

    # Validator method
    @classmethod
    def judge(cls, **kwargs):
        """
        Validator method
        """
        try:
            return cls(**kwargs)
        except ValidationError as e:
            raise e

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """
        return cls.judge(**kwargs).model_dump(exclude_none=True)
