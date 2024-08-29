from typing import (
    Any,
    Callable,
    Dict,
    List,
    Literal,
    Mapping,
    Optional,
    Tuple,
    Type,
    TypedDict,
    Union,
)

import numpy as np

from echopop.spatial.variogram import create_optimization_options
from echopop.survey import Survey

file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
survey = Survey(init_config, file_config)
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis(exclude_age1=True)
self = survey
survey.fit_variogram()
self.analysis['settings']['variogram']['optimization']

##################
model = ["bessel", "exponential"]
n_lags = 30

# Get the default variogram parameters
default_variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
# ---- Update model, n_lags
default_variogram_parameters.update({
    "model": model,
    "n_lags": n_lags
    })

# Extract specific variogram parameters
# ---- Number of lags
n_lags = default_variogram_parameters["n_lags"]
# ---- Lag resolution
lag_resolution = default_variogram_parameters["lag_resolution"]

# Compute the lag distances and max range
# ---- Add to the `variogram_parameters` dictionary
default_variogram_parameters["distance_lags"] = np.arange(1, n_lags) * lag_resolution
# ---- Update the max range parameter, if necessary
default_variogram_parameters["max_range"] = lag_resolution * n_lags

# ---- Update `VariogramBase`
from echopop.utils.validate import VariogramBase

VariogramBase.DEFAULT_VALUES
VariogramBase.update_defaults(default_variogram_parameters)
VariogramBase.DEFAULT_VALUES
variogram_parameters = VariogramBase.create(**test_dict)
# ---- Create VariogramOptimize
optimization_parameters = VariogramOptimize.create(**{"max_fun_evaluations": 400, "trust_region_solver": "exact"})



# Define the parameter property TypedDict
class InitialValues(TypedDict):
    low: Optional[float]
    upper: Optional[float]
    initial: Optional[float]

    # Default values for parameters
    DEFAULT_STRUCTURE = {"low": 0.0, "upper": np.inf, "initial": 0.0}

    # Define expected types for parameter details
    EXPECTED_DTYPES = {
        "low": Optional[float],
        "upper": Optional[float],
        "initial": Optional[float],
    }

# Define the main TypedDict
class VariogramInitial(TypedDict, total=False):
    fit_parameters: Union[List[str], Dict[str, InitialValues]]

    def __init__(self,
                 fit_parameters: Union[List[str], Dict[str, InitialValues]]):
        # Initialize the instance
        self.fit_parameters = fit_parameters

    @classmethod
    def create(cls,
               fit_parameters: Union[List[str], Dict[str, InitialValues]]) -> 'VariogramInitial':
        """
        Create a ParameterConfig from a list of parameter names or a dictionary.
        """
        if isinstance(fit_parameters, list):
            fit_parameters = {param: InitialValues.DEFAULT_STRUCTURE.copy() for param in fit_parameters}
        else:
            for param, details in fit_parameters.items():
                if not isinstance(details, dict):
                    raise TypeError(f"Value for '{param}' must be a dictionary.")
                for key in InitialValues.DEFAULT_STRUCTURE:
                    if key not in details:
                        fit_parameters[param][key] = InitialValues.DEFAULT_STRUCTURE[key]
        # ---- Update
        cls.validate(fit_parameters)
        return cls(fit_parameters)

    @staticmethod
    def validate(params: Dict[str, InitialValues]) -> None:
        """
        Validate the input dictionary against the expected structure.
        """
        for param, details in params.items():
            for key, default_value in InitialValues.DEFAULT_STRUCTURE.items():
                if key not in details:
                    details[key] = default_value
                # Validate type
                if not isinstance(details[key], float):
                    expected_type = float
                    actual_type = type(details[key])
                    expected_description = describe_type(expected_type)
                    actual_description = describe_type(actual_type)
                    raise TypeError(
                        f"Value for '{key}' in parameter '{param}' (type: {actual_description}) "
                        f"does not match expected type {expected_description}."
                    )

test_params = {
    "nugget": {"low": 0.1, "initial": 1.0},
    "sill": {"low": 0.0, "initial": 2.0},
    "correlation_range": {"low": 0.5},
    "hole_effect_range": {"low": 1.0, "upper": 5},
    "decay_power": {"upper": np.inf, "initial": -0.5}
}

parameter_config = VariogramInitial(test_params)
VariogramInitial.create(test_params)

test_params_list = ["nugget", "sill", "correlation_range", "hole_effect_range", "decay_power"]
fit_parameters=VariogramInitial.create(test_params_list)
params = fit_parameters
fit_parameters = test_params.copy()
def test_fun(
        params1: VariogramBase,
        params2: VariogramOptimize,
        params3: VariogramInitial,
        verbose = True,
):
    print(VariogramBase.create(**params1))
    print(VariogramOptimize.create(**params2))
    print(VariogramInitial.create(params3))

test_fun(params1={"n_lags": 30}, params2={"cost_fun_tolerance": 1e-10}, params3=["nugget"])

config_params = survey.input["statistics"]["variogram"]["model_config"].copy()
input_variogram_params = VariogramBase.create(**{"n_lags": 50, "sill": 1.50})
VariogramOptimize.create(**{key: {"initial": value} for key, value in input_variogram_params.items()} )
{**config_params, **input_variogram_params}
{**input_variogram_params, **config_params}








def variogram_analysis(
    variogram_parameters: Dict[str, Union[float, bool]]
)

def fit_variogram(
    self,
    base_variogram_model: Dict[str, Union[str, List[str]]] = {
        "model": Union[str, List[str]] = ["bessel", "exponential"],
        "n_lags": int = 30,
        # ... other parameters
    },
    optimization_parameters: Dict[str, Union[List[str], List[Tuple[str, float]], Dict[str, Dict[str, float]], None]] = {
        "fit_parameters": [
            "nugget",
            "sill",
            "correlation_range",
            "hole_effect_range",
            "decay_power",
        ],
        # ... other optimization parameters
    },
    variable: Literal["biomass", "abundance"] = "biomass",
    verbose: bool = True,
    standardize_coordinates: bool = True,
    force_lag_zero: bool = True,
) -> None:
    # ... function body


def fit_variogram(
    self,
    base_variogram_model={
        "model": ["bessel", "exponential"],
        "n_lags": 30,
        "azimuth_range": 360.0,
        "lag_resolution": None,
        "max_range": None,
        "sill": None,
        "nugget": None,
        "hole_effect_range": None,
        "correlation_range": None,
        "enhance_semivariance": None,
        "decay_power": None,
    },
    optimization_parameters={
        "fit_parameters": [
            "nugget",
            "sill",
            "correlation_range",
            "hole_effect_range",
            "decay_power",
        ],
        "initial_values": None,
        "lower_bounds": [
            ("nugget", 0.0),
            ("sill", 0.0),
            ("correlation_range", 0.0),
            ("hole_effect_range", 0.0),
            ("decay_power", 0.0),
        ],
        "upper_bounds": None,
        "max_fun_evaluations": 500,
        "cost_fun_tolerance": 1e-6,
        "solution_tolerance": 1e-4,
        "gradient_tolerance": 1e-4,
        "finite_step_size": 1e-8,
        "trust_region_solver": "exact",
        "x_scale": "jacobian",
        "jacobian_approx": "forward",
        "force_lag_zero": True,
    },
    variable="biomass",
    verbose=True,
    standardize_coordinates=True,
    force_lag_zero=True,
):
    # Combine base_variogram_model and optimization_parameters into a single dict
    combined_parameters = {**base_variogram_model, **optimization_parameters}

    # Use combined_parameters in your function logic
    # ...

    # Access individua
