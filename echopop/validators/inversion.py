from typing import Any, Callable, Dict, List, Optional

import numpy as np
import pandas as pd
from pydantic import ConfigDict, Field, field_validator, model_validator

from ..core.validators import BaseDataFrame, BaseDictionary
from ..inversion.api import SCATTERING_MODEL_PARAMETERS
from ..typing import InvParameters
from ..typing.inversion import ModelInputParameters


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
    BaseDictionary, arbitrary_types_allowed=True, title="TS-length inversion model parameters"
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


class ScatterDF(BaseDataFrame):
    """
    Validate hierarchical columns in a MultiIndex DataFrame.

    Ensures the top-level column names include 'sv_mean', 'nasc', and 'thickness_mean',
    and that the second-level column names ('frequency') are numeric. All values
    under each frequency must be floats.
    """

    class Config:
        multiindex_strict = True
        multiindex_coerce = True

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


class EnvironmentParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="scattering model environment/medium parameters",
):
    """
    Environmental and medium parameters for acoustic scattering models.

    This class validates and stores physical properties of the medium
    (typically seawater) that affect acoustic propagation and scattering.
    These parameters are essential for accurate scattering calculations.

    Attributes
    ----------
    sound_speed_sw : float, default=1500.0
        Sound speed in seawater in m/s. Must be positive.
        Typical range: 1450-1550 m/s depending on temperature,
        salinity, and pressure.
    density_sw : float, default=1025.0
        Seawater density in kg/m³. Must be positive.
        Typical range: 1020-1030 kg/m³ depending on temperature
        and salinity.

    Notes
    -----
    These parameters directly affect acoustic impedance and scattering
    cross-sections through the wavenumber and reflection coefficients.

    Sound speed variation with depth should be considered in deep-water
    applications, though single representative values are often adequate
    for survey-scale analyses.

    Examples
    --------
    >>> env = EnvironmentParameters(sound_speed_sw=1485.0, density_sw=1026.5)
    >>> env.sound_speed_sw
    1485.0
    """

    sound_speed_sw: float = Field(default=1500.0, gt=0.0, allow_inf_nan=False)
    density_sw: float = Field(default=1025.0, gt=0.0, allow_inf_nan=False)


class SimulationParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="inversion simulation parameters",
):
    """
    Configuration parameters for inversion simulation and optimization.

    This class manages settings that control the behavior of the inversion
    process, including Monte Carlo sampling, parameter scaling, and
    optimization callbacks.

    Attributes
    ----------
    iter_cb : Optional[Callable], default=None
        Callback function called during optimization iterations.
        Function signature: iter_cb(params, iter, resid, *args, **kwargs)
    monte_carlo : Optional[bool], default=None
        Whether to use Monte Carlo initialization for optimization.
        Enables warm-start strategy with multiple parameter realizations.
    mc_realizations : Optional[int], default=None
        Number of Monte Carlo realizations to generate.
        Automatically set to 100 if monte_carlo=True and not specified.
    mc_seed : Optional[int], default=None
        Random seed for reproducible Monte Carlo sampling.
        If None, uses system entropy for random initialization.
    scale_parameters : bool, default=True
        Whether to scale parameters to [0,1] range before optimization.
        Recommended for improved numerical conditioning.
    minimum_frequency_count : int, default=2
        Minimum number of frequencies required for inversion.
        Prevents under-determined optimization problems.

    Methods
    -------
    validate_monte_carlo_realizations()
        Model validator ensuring consistent Monte Carlo configuration

    Notes
    -----
    Monte Carlo initialization can significantly improve optimization
    convergence by providing multiple starting points. The best-performing
    realization is selected as the initialization for the main optimization.

    Parameter scaling transforms all parameters to the unit interval [0,1],
    which improves optimization performance when parameters have very
    different scales (e.g., length in mm vs density in kg/m³).

    Examples
    --------
    >>> sim_params = SimulationParameters(
    ...     monte_carlo=True,
    ...     mc_realizations=50,
    ...     scale_parameters=True
    ... )
    """

    iter_cb: Optional[Callable] = Field(default=None)
    monte_carlo: Optional[bool] = Field(default=None)
    mc_realizations: Optional[int] = Field(default=None)
    mc_seed: Optional[int] = Field(default=None)
    scale_parameters: bool = Field(default=True)
    minimum_frequency_count: int = Field(default=2)

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


class ValidateInversionMatrix(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="scattering model inversion analysis",
):
    """
    Validation model for InversionMatrix class configuration.

    This validator ensures that data and simulation settings are properly
    formatted and compatible for acoustic scattering inversion analysis.
    It enforces data structure requirements and validates simulation parameters.

    Attributes
    ----------
    data : pd.DataFrame
        MultiIndex DataFrame containing acoustic measurements with required
        columns: 'sv_mean', 'nasc', 'thickness_mean' at the top level,
        and 'frequency' as a nested index level.
    simulation_settings : SimulationParameters
        Configuration parameters for the inversion simulation including
        Monte Carlo settings, parameter scaling, and optimization callbacks.

    Methods
    -------
    validate_data(v)
        Field validator ensuring data conforms to ScatterDF schema

    Notes
    -----
    This validator is typically used internally by InversionMatrix during
    initialization to ensure all inputs are properly formatted and contain
    the required information for successful inversion.

    The data validation ensures the presence of essential acoustic quantities:
    - sv_mean: Volume backscattering strength (dB re 1 m^-1)
    - nasc: Nautical Area Scattering Coefficient (m²/nmi²)
    - thickness_mean: Mean layer thickness (m)

    Examples
    --------
    >>> validator = ValidateInversionMatrix(
    ...     data=acoustic_dataframe,
    ...     simulation_settings=SimulationParameters(monte_carlo=True)
    ... )
    """

    data: pd.DataFrame
    simulation_settings: SimulationParameters = Field(default_factory=SimulationParameters)

    @field_validator("data", mode="after")
    def validate_data(cls, v):
        # Validate with pandera
        return ScatterDF.validate(v)


class ModelSettingsParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="general model settings",
):
    """
    Configuration parameters for acoustic scattering models.

    This class stores model-specific settings including the scattering model
    type and environmental parameters. It uses flexible configuration to
    accommodate different scattering model requirements.

    Attributes
    ----------
    type : str
        Scattering model type identifier (e.g., 'pcdwba', 'spheres').
        Must match an available model in SCATTERING_MODEL_PARAMETERS.
    environment : EnvironmentParameters
        Environmental parameters including sound speed and seawater density.
        Defaults to standard oceanographic values if not specified.

    Notes
    -----
    This class uses Pydantic's extra="allow" configuration to permit
    model-specific parameters beyond the core type and environment settings.
    This flexibility allows each scattering model to specify its own
    additional configuration parameters.

    The model type determines which validation schema is applied during
    parameter checking and which forward model is used for scattering
    calculations.

    Examples
    --------
    >>> model_settings = ModelSettingsParameters(
    ...     type="pcdwba",
    ...     environment=EnvironmentParameters(sound_speed_sw=1485.0),
    ...     taper_order=10.0,  # PCDWBA-specific parameter
    ...     frequency_interval=2000.0
    ... )
    """

    type: str
    environment: EnvironmentParameters = Field(default_factory=EnvironmentParameters)
    model_config = ConfigDict(extra="allow")


class ValidateBuildModelArgs(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="scattering model preparation and assembly",
):
    """
    Validation model for scattering model build arguments.

    This validator ensures that model parameters and settings are compatible
    with the specified scattering model type. It performs cross-validation
    between parameter specifications and model requirements.

    Attributes
    ----------
    model_parameters : InvParameters
        Container with biological and physical parameters for the scattering model.
        Parameters must be compatible with the specified model type.
    model_settings : ModelSettingsParameters
        Model configuration including type identifier and environmental parameters.
        The model type determines validation requirements for parameters.

    Methods
    -------
    _check_variable(params, variable, validator)
        Static method to validate parameter bounds against model schema
    validate_model_parameterization()
        Model validator ensuring parameter compatibility with model type

    Notes
    -----
    This validator performs sophisticated cross-validation between the
    specified scattering model type and the provided parameters. It ensures
    that all required parameters are present and that parameter bounds
    are within physically reasonable ranges for the model.

    The validation process:
    1. Checks that the model type exists in SCATTERING_MODEL_PARAMETERS
    2. Validates parameter bounds against model-specific schemas
    3. Ensures all required parameters are specified with appropriate ranges

    This prevents runtime errors during model execution by catching
    configuration problems at initialization time.

    Examples
    --------
    >>> validator = ValidateBuildModelArgs(
    ...     model_parameters=InvParameters(param_dict),
    ...     model_settings=ModelSettingsParameters(type="pcdwba")
    ... )
    """

    model_parameters: InvParameters
    model_settings: ModelSettingsParameters
    model_config = ConfigDict(protected_namespaces=())

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

        # Check for model-type
        # ---- Dump the model
        model_settings = self.model_settings.model_dump()
        # ---- Check against reference
        if model_settings["type"] not in SCATTERING_MODEL_PARAMETERS:
            raise LookupError(f"Scattering model '{model_settings["type"]}' could not be found.")

        # Get the model-specific validator
        model_validators = SCATTERING_MODEL_PARAMETERS[model_settings["type"]]

        # Validate `model_parameters`
        # ---- Dump the model
        model_parameters = self.model_parameters.parameters
        # ---- Check minima
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
        self.model_parameters = InvParameters(
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

        # Validate `model_settings`
        self.model_settings = ModelSettingsParameters.model_validate(
            model_validators["settings"].create(**self.model_settings.model_dump())
        )

        return self
