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
survey.load_acoustic_data(ingest_exports="echoview")
survey.load_survey_data()
survey.transect_analysis(exclude_age1=True)
self = survey
survey.fit_variogram()
self.analysis['settings']['variogram']['optimization']


def fit_variogram(
        self,
        variogram_parameters: Dict[str, Union[List[str], int, float, None]] = {
            "model": Union[str, List[str]] = ["bessel", "exponential"]
        }):
    pass


def describe_type(expected_type: Any) -> str:
    """
    Convert a type hint into a human-readable string.
    """

    if hasattr(expected_type, "fail_state"):
        return expected_type.fail_state

    if hasattr(expected_type, "__origin__"):
        origin = expected_type.__origin__

        if origin is Literal:
            return f"one of {expected_type.__args__}"

        if origin is Union:
            args = expected_type.__args__
            descriptions = [describe_type(arg) for arg in args if arg is not type(None)]
            if not descriptions:
                return "any value"
            return " or ".join(descriptions)

        if origin is np.ndarray:
            return "numpy.ndarray of floats"

        if origin is list:
            item_type = expected_type.__args__[0]
            return f"list of '{describe_type(item_type)}'"

    if isinstance(expected_type, type):
        return f"'{expected_type.__name__}'"

    # return str(expected_type)
    return f"'{expected_type}'"


def validate_literal(value: Any, literals: List[Any]) -> bool:
    """
    Validate if the value matches one of the literal values.
    """
    return value in literals


def validate_type(value: Any, expected_type: Any) -> bool:
    """
    Validate if the value matches the expected type.
    """

    # Handle base `type` case
    if isinstance(expected_type, type):
        # ---- Handle `posint`
        if expected_type == posint:
            try:
                posint(value)
                return True
            except (TypeError, ValueError):
                return False
        else:
            return isinstance(value, expected_type)

    # Handle Union
    if hasattr(expected_type, "__origin__") and expected_type.__origin__ is Union:
        return any(validate_type(value, arg) for arg in expected_type.__args__)

    # Handle Literal
    if hasattr(expected_type, "__origin__") and expected_type.__origin__ is Literal:
        return validate_literal(value, expected_type.__values__)

    # Handle List
    if hasattr(expected_type, "__origin__") and expected_type.__origin__ is list:
        if not isinstance(value, list):
            return False
        item_type = expected_type.__args__[0]
        return all(validate_type(item, item_type) for item in value)

    # Handle numpy.ndarray
    if hasattr(expected_type, "__origin__") and expected_type.__origin__ is np.ndarray:
        return isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.number) and value.dtype == np.float64

    return False

# def validate_type(value: Any, expected_type: Type) -> bool:
#     """
#     Validate if the value matches the expected type.
#     """
#     if isinstance(expected_type, type):
#         return isinstance(value, expected_type)
#     elif hasattr(expected_type, "__origin__"):
#         origin = expected_type.__origin__
#         if origin is list:
#             if not isinstance(value, list):
#                 return False
#             # Check if all items in the list match the expected type
#             item_type = expected_type.__args__[0]
#             return all(validate_type(item, item_type) for item in value)
#         # Handle Tuple types
#         elif origin is tuple:
#             if not isinstance(value, tuple):
#                 return False
#             if len(value) != len(expected_type.__args__):
#                 return False
#             return all(validate_type(v, t) for v, t in zip(value, expected_type.__args__))

#         elif origin is Union:
#             return any(validate_type(value, arg) for arg in expected_type.__args__)
#     return False
# def validate_type(value: Any, expected_type: Type) -> bool:
#     """
#     Validate if the value matches the expected type.
#     """
#     if isinstance(expected_type, type):
#         return isinstance(value, expected_type)

#     # Handle Literal types
#     elif hasattr(expected_type, "__origin__") and expected_type.__origin__ is Literal:
#         allowed_values = expected_type.__args__
#         return value in allowed_values

#     # Handle Union types
#     elif hasattr(expected_type, "__origin__") and expected_type.__origin__ is Union:
#         return any(validate_type(value, arg) for arg in expected_type.__args__)

#     # Handle List types
#     elif hasattr(expected_type, "__origin__") and expected_type.__origin__ is list:
#         if not isinstance(value, list):
#             return False
#         item_type = expected_type.__args__[0]
#         return all(validate_type(item, item_type) for item in value)

#     # Handle Tuple types
#     elif hasattr(expected_type, "__origin__") and expected_type.__origin__ is tuple:
#         if not isinstance(value, tuple):
#             return False
#         if len(value) != len(expected_type.__args__):
#             return False
#         return all(validate_type(v, t) for v, t in zip(value, expected_type.__args__))

#     return False
# def validate_dict(data: Dict[str, Any], expected_types: Dict[str, Type]) -> None:
#     """
#     Validate a dictionary against expected types.
#     """
#     for key, expected_type in expected_types.items():
#         if key in data:
#             if not validate_type(data[key], expected_type):
#                 raise TypeError(f"Value for '{key}' does not match expected type {expected_type}.")

# def validate_dict(data: Dict[str, Any], expected_types: Dict[str, Type]) -> None:
#     """
#     Validate a dictionary against expected types.
#     """
#     for key, expected_type in expected_types.items():
#         if key in data:
#             if not validate_type(data[key], expected_type):
#                 raise TypeError(f"Value for '{key}' does not match expected type {expected_type}.")

def validate_dict(data: Dict[str, Any], expected_types: Dict[str, Any]) -> None:
    """
    Validate a dictionary against expected types.
    """
    for key, expected_type in expected_types.items():
        if key in data:
            if not validate_type(data[key], expected_type):
                expected_description = describe_type(expected_type)
                actual_description = type(data[key]).__name__
                if hasattr(expected_type, "fail_state"):
                    raise TypeError(f"Value for '{key}' (type: '{actual_description}') {expected_description}.")
                else:
                    raise TypeError(f"Value for '{key}' (type: '{actual_description}') does not match expected type {expected_description}.")

from typing import Annotated


class posint(int):

    fail_state = "must be a non-negative integer"

    def __new__(cls, value):
        if not isinstance(value, int) or value < 0:
            raise ValueError("Value must be a non-negative integer.")
        return super().__new__(cls, value)

expected_type = PositiveInt
isinstance(expected_type, PositiveInt)
type(expected_type)

posint = Annotated[PositiveInt, "must be a non-negative integer"]

def validate_positive_int(value: int) -> bool:
    """
    Validate if the value is a non-negative integer.
    """
    return isinstance(value, int) and value >= 0

class VariogramBase(TypedDict, total=False):
    # Provide the type-hints
    model: Union[str, List[str]]
    # n_lags: int
    n_lags: posint
    lag_resolution: Optional[float]
    max_range: Optional[float]
    sill: Optional[float]
    nugget: Optional[float]
    hole_effect_range: Optional[float]
    correlation_range: Optional[float]
    enhance_semivariance: Optional[bool]
    decay_power: Optional[float]

    # Define default values
    DEFAULT_VALUES = {
        "model": ["bessel", "exponential"],
        "n_lags": 30,
        "lag_resolution": None,
        "max_range": None,
        "sill": None,
        "nugget": None,
        "hole_effect_range": None,
        "correlation_range": None,
        "enhance_semivariance": None,
        "decay_power": None
    }

    # Define the expected datatypes
    EXPECTED_DTYPES = {
        "model": Union[str, List[str]],
        # "n_lags": int,
        "n_lags": posint,
        "lag_resolution": Optional[float],
        "max_range": Optional[float],
        "sill": Optional[float],
        "nugget": Optional[float],
        "hole_effect_range": Optional[float],
        "correlation_range": Optional[float],
        "enhance_semivariance": Optional[bool],
        "decay_power": Optional[float],
    }

    @classmethod
    def create(cls,
               model: Union[str, List[str]] = ["bessel", "exponential"],
            #    n_lags: int = 30,
               n_lags: posint = 30,
               lag_resolution: Optional[float] = None,
               max_range: Optional[float] = None,
               sill: Optional[float] = None,
               nugget: Optional[float] = None,
               hole_effect_range: Optional[float] = None,
               correlation_range: Optional[float] = None,
               enhance_semivariance: Optional[bool] = None,
               decay_power: Optional[float] = None,
               **kwargs
            ) -> "VariogramBase":
        """
        Base variogram model parameters

        Parameters
        ----------
        model: Union[str, List[str]]
            A string or list of model names. A single name represents a single family model. Two
            inputs represent the desired composite model (e.g. the composite J-Bessel and
            exponential model). Defaults to: ['bessel', 'exponential']. Available models and their
            required arguments can be reviewed in the :fun:`echopop.spatial.variogram.variogram`
            function.
        n_lags: int
            See the `variogram_parameters` argument in
            :fun:`echopop.spatial.variogram.empirical_variogram` for more details on
            `n_lags`.
        sill: Optional[float]
            See the description of `sill` in
            :fun:`echopop.spatial.variogram.variogram`.
        nugget: Optional[float]
            See the description of `nugget` in
            :fun:`echopop.spatial.variogram.variogram`.
        correlation_range: Optional[float]
            See the description of `correlation_range` in
            :fun:`echopop.spatial.variogram.variogram`.
        hole_effect_range: Optional[float]
            See the description of `hole_effect_range` in
            :fun:`echopop.spatial.variogram.variogram`.
        decay_power: Optional[float]
            See the description of `decay_power` in
            :fun:`echopop.spatial.variogram.variogram`.
        enhanced_semivariance: Optional[bool]
            See the description of `enhanced_semivariance` in
            :fun:`echopop.spatial.variogram.variogram`.
        max_range: Optional[float]
            An optional input defining the maximum lag distance range that will be computed for
            fitting the theoretical variogram parameters.

        Returns
        ----------
        VariogramBase: A validated dictionary with the user-defined variogram parameter values and
        default values for any missing parameters/keys.
        """

        # User-defined parameters
        inputs = {
            "model": model,
            "n_lags": n_lags,
            "lag_resolution": lag_resolution,
            "max_range": max_range,
            "sill": sill,
            "nugget": nugget,
            "hole_effect_range": hole_effect_range,
            "correlation_range": correlation_range,
            "enhance_semivariance": enhance_semivariance,
            "decay_power": decay_power
        }

        # Merge the default values with those provided by the user
        params = {**cls.DEFAULT_VALUES, **inputs}

        # Filter the parameter keys
        filtered_params = {key: params[key] for key in cls.EXPECTED_DTYPES if key in params}

        # Validate the parameter datatypes
        cls.validate(filtered_params)

        return filtered_params

    # Create validation method
    @staticmethod
    def validate(data: Dict[str, Any]):
        """
        Validate the input dictionary against the `VariogramBase` class definition/schema.

        Parameters
        ----------
        data: Dict[str, Any]
            A dictionary containing the parameters that will be validated.

        Raises
        ----------
        TypedError:
            If any value does not match the expected datatype.
        """

        # Define expected datatypes
        validate_dict(data, VariogramBase.EXPECTED_DTYPES)
        # FOR DEBUGGING
        # --------
        # print("Validate passed.")

test_dict =  {"enhance_semivariance": True, "model": ["exponential", "bessel"], "correlation_range": 0.008, "n_lags": 50}
result = VariogramBase.create(**test_dict)
result
VariogramBase.create(**{"enhance_semivariance": True, "bababooey": float})
invalid_dict = {"enhance_semivariance": 100}
result = VariogramBase.create(**invalid_dict)
result

class VariogramOptimize(TypedDict):
    max_fun_evaluations: int
    cost_fun_tolerance: float
    gradient_tolerance: float
    finite_step_size: float
    trust_region_solver: Literal["base", "exact"]
    x_scale: Union[Literal["jacobian"], np.ndarray[float]]
    jacobian_approx: Literal["forward", "central"] = "forward"

    # Define default values
    DEFAULT_VALUES = {
        "max_fun_evaluations": 500,
        "cost_fun_tolerance": 1e-6,
        "solution_tolerance": 1e-4,
        "gradient_tolerance": 1e-4,
        "finite_step_size": 1e-8,
        "trust_region_solver": "exact",
        "x_scale": "jacobian",
        "jacobian_approx": "forward",
        }

    # Define the expected datatypes
    EXPECTED_DTYPES = {
        "max_fun_evaluations": int,
        "cost_fun_tolerance": float,
        "solution_tolerance": float,
        "gradient_tolerance": float,
        "finite_step_size": float,
        "trust_region_solver": Literal["base", "exact"],
        "x_scale": Union[Literal["jacobian"], np.ndarray[float]],
        "jacobian_approx":Literal["forward", "central"],
    }

    @classmethod
    def create(cls,
               max_fun_evaluations: int = 500,
               cost_fun_tolerance: float = 1e-6,
               solution_tolerance: float = 1e-4,
               gradient_tolerance: float = 1e-4,
               finite_step_size: float = 1e-8,
               trust_region_solver: Literal["exact", "base"] = "exact",
               x_scale: Union[Literal["jacobian"], np.ndarray[float]] = "jacobian",
               jacobian_approx: Literal["forward", "central"] = "forward",
               **kwargs
            ) -> "VariogramOptimize":
        """
        Base variogram model parameters

        Parameters
        ----------
        max_fun_evaluations: int
            The maximum number of evaluations. Defaults to 500.
        cost_fun_tolerance: float
            Threshold used for determining convergence via incremental changes of the cost function.
            Defaults to 1e-6.
        solution_tolerance: float
            Threshold used for determining convergence via change of the independent variables.
            Defaults to 1e-8.
        gradient_tolerance: float
            Threshold used for determining convergence via the gradient norma. Defaults to 1e-8.
        finite_step_size: float
            The relative step sizes used for approximating the Jacobian via finite differences.
        trust_region_solver: Literal["exact", "float"]
            The method used for solving the trust-region problem by either using the Jacobian
            computed from the first iteration (`"base"`) or via singular value decomposition
            (`"exact"`). Defaults to "exact".
        x_scale: Union[Literal["jacobian"], np.ndarray[float]]
            When `x_scale="jacobian"`, the characteristic scale is updated across numerical
            iterations via the inverse norms of the Jacobian matrix. Otherwise, a `np.ndarray`
            of the same length as `fit_parameters` can provide a constant scaling factor.
        jacobian_approx: Literal["forward", "central"]
            Indicates whether forward differencing (`"forward"`) or central differencing
            (`"central"`) should be used to approximate the Jacobian matrix.

        Returns
        ----------
        VariogramOptimize: A validated dictionary with the user-defined variogram optimizatization
        parameter values and default values for any missing parameters/keys.
        """

        # User-defined parameters
        inputs = {
            "max_fun_evaluations": max_fun_evaluations,
            "cost_fun_tolerance": cost_fun_tolerance,
            "solution_tolerance": solution_tolerance,
            "gradient_tolerance": gradient_tolerance,
            "finite_step_size": finite_step_size,
            "trust_region_solver": trust_region_solver,
            "x_scale": x_scale,
            "jacobian_approx": jacobian_approx
        }

        # Merge the default values with those provided by the user
        params = {**cls.DEFAULT_VALUES, **inputs}

        # Filter the parameter keys
        filtered_params = {key: params[key] for key in cls.EXPECTED_DTYPES if key in params}

        # Validate the parameter datatypes
        cls.validate(filtered_params)

        return filtered_params

    # Create validation method
    @staticmethod
    def validate(data: Dict[str, Any]):
        """
        Validate the input dictionary against the `VariogramBase` class definition/schema.

        Parameters
        ----------
        data: Dict[str, Any]
            A dictionary containing the parameters that will be validated.

        Raises
        ----------
        TypedError:
            If any value does not match the expected datatype.
        """

        # Define expected datatypes
        validate_dict(data, VariogramOptimize.EXPECTED_DTYPES)
        # FOR DEBUGGING
        # --------
        # print("Validate passed.")

test_dict_invalid = {
    "jacobian_approx": 1  # Invalid Literal value
}
result = VariogramOptimize.create(**test_dict_invalid)

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
    "hole_effect_range": {"low": 1.0, "upper": 5.0},
    "decay_power": {"upper": np.inf, "initial": 0.5}
}
parameter_config = VariogramInitial(fit_parameters=test_params["fit_parameters"])
VariogramInitial.create(test_params)


test_params_list = ["nugget", "sill", "correlation_range", "hole_effect_range", "decay_power"]
VariogramInitial.create(fit_parameters=test_params_list)

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
