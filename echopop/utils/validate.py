"""
Validation functions.
"""

# TODO: Compile all package validators here since they may not belong elsewhere
from typing import Any, Dict, List, Literal, Optional, TypedDict, Union, get_args

import numpy as np


# CLASS-SPECIFIC CORE API
class posint(int):
    """Positive-only integer (includes 0)"""

    __failstate__ = "must be a non-negative integer"

    def __new__(cls, value):
        if not isinstance(value, int) or value < 0:
            raise ValueError("Value must be a non-negative integer.")
        return super().__new__(cls, value)


class posfloat(float):
    """Positive-only float (includes 0.0)"""

    __failstate__ = "must be a non-negative float"

    def __new__(cls, value):
        if not isinstance(value, (float, int)) or value < 0:
            raise ValueError("Value must be a non-negative float.")
        return super().__new__(cls, value)


class realposfloat(posfloat):
    """Real number positive-only float (includes 0.0)"""

    __failstate__ = "must be a non-negative real number"

    def __new__(cls, value):
        if not isinstance(value, (float, int)) or np.isinf(value):  # Check if value is infinity
            raise ValueError(f"Value {cls.__failstate__}.")
        return super().__new__(cls, value)


class realcircle(realposfloat):
    """Real number in a unit circle"""

    __failstate__ = "must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees"

    def __new__(cls, value):
        if not isinstance(value, (float, int)) or (value < 0.0 or value > 360.0):
            raise ValueError(f"Value {cls.__failstate__}.")
        return super().__new__(cls, value)


# Validation functions
def describe_type(expected_type: Any) -> str:
    """
    Convert a type hint into a human-readable string.
    """

    if hasattr(expected_type, "__failstate__"):
        return expected_type.__failstate__

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


def validate_type(value: Any, expected_type: Any) -> bool:
    """
    Validate if the value matches the expected type.
    """

    # Handle numpy.ndarray
    if hasattr(expected_type, "__origin__") and expected_type.__origin__ is np.ndarray:
        return (
            isinstance(value, np.ndarray)
            and np.issubdtype(value.dtype, np.number)
            # and value.dtype == np.float64
        )

    # Handle base `type` case
    if isinstance(expected_type, type):
        # ---- Handle `posint`
        if expected_type in [posint, posfloat, realposfloat, realcircle]:
            try:
                expected_type(value)
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
        # ---- Get allowed values
        allowed_values = get_args(expected_type)
        if isinstance(value, np.ndarray):
            allowed_values_array = np.array(list(allowed_values), dtype=object)
            value_array = np.array(value, dtype=object)
            return np.array_equal(np.sort(value_array), np.sort(allowed_values_array))
        return value in allowed_values

    # Handle List
    if hasattr(expected_type, "__origin__") and expected_type.__origin__ is list:
        if not isinstance(value, list):
            return False
        item_type = expected_type.__args__[0]
        return all(validate_type(item, item_type) for item in value)

    # # Handle numpy.ndarray
    # if hasattr(expected_type, "__origin__") and expected_type.__origin__ is np.ndarray:
    #     return (
    #         isinstance(value, np.ndarray)
    #         and np.issubdtype(value.dtype, np.number)
    #         # and value.dtype == np.float64
    #     )

    return False


def validate_typed_dict(data: Dict[str, Any], expected_types: Dict[str, Any]) -> None:
    """
    Validate a dictionary against expected types.
    """
    for key, expected_type in expected_types.items():
        if key in data:
            if not validate_type(data[key], expected_type):
                expected_description = describe_type(expected_type)
                actual_description = type(data[key]).__name__
                if hasattr(expected_type, "__failstate__"):
                    raise TypeError(
                        f"Value for '{key}' ({data[key]}, type: '{actual_description}') "
                        f"{expected_description}."
                    )
                else:
                    raise TypeError(
                        f"Value for '{key}' ({data[key]}, type: '{actual_description}') does not "
                        f"match expected type {expected_description}."
                    )


# Validation classes
class VariogramBase(TypedDict, total=False):
    model: Union[str, List[str]]
    n_lags: posint
    lag_resolution: Optional[realposfloat]
    max_range: Optional[realposfloat]
    sill: Optional[realposfloat]
    nugget: Optional[realposfloat]
    hole_effect_range: Optional[realposfloat]
    correlation_range: Optional[realposfloat]
    enhance_semivariance: Optional[bool]
    decay_power: Optional[realposfloat]

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
        "decay_power": None,
    }

    # Root default values -- used for reinitialization
    _ROOT_DEFAULT = {
        "model": ["bessel", "exponential"],
        "n_lags": 30,
        "lag_resolution": None,
        "max_range": None,
        "sill": None,
        "nugget": None,
        "hole_effect_range": None,
        "correlation_range": None,
        "enhance_semivariance": None,
        "decay_power": None,
    }

    # Define the expected datatypes
    EXPECTED_DTYPES = {
        "model": Union[str, List[str]],
        "n_lags": posint,
        "lag_resolution": Optional[realposfloat],
        "max_range": Optional[realposfloat],
        "sill": Optional[realposfloat],
        "nugget": Optional[realposfloat],
        "hole_effect_range": Optional[realposfloat],
        "correlation_range": Optional[realposfloat],
        "enhance_semivariance": Optional[bool],
        "decay_power": Optional[realposfloat],
    }

    @classmethod
    def create(
        cls,
        model: Union[str, List[str]] = None,
        n_lags: posint = None,
        lag_resolution: Optional[realposfloat] = None,
        max_range: Optional[realposfloat] = None,
        sill: Optional[realposfloat] = None,
        nugget: Optional[realposfloat] = None,
        hole_effect_range: Optional[realposfloat] = None,
        correlation_range: Optional[realposfloat] = None,
        enhance_semivariance: Optional[posfloat] = None,
        decay_power: Optional[realposfloat] = None,
        **kwargs,
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
        n_lags: posint
            See the `variogram_parameters` argument in
            :fun:`echopop.spatial.variogram.empirical_variogram` for more details on
            `n_lags`.
        sill: Optional[realposfloat]
            See the description of `sill` in
            :fun:`echopop.spatial.variogram.variogram`.
        nugget: Optional[realposfloat]
            See the description of `nugget` in
            :fun:`echopop.spatial.variogram.variogram`.
        correlation_range: Optional[realposfloat]
            See the description of `correlation_range` in
            :fun:`echopop.spatial.variogram.variogram`.
        hole_effect_range: Optional[realposfloat]
            See the description of `hole_effect_range` in
            :fun:`echopop.spatial.variogram.variogram`.
        decay_power: Optional[realposfloat]
            See the description of `decay_power` in
            :fun:`echopop.spatial.variogram.variogram`.
        enhanced_semivariance: Optional[bool]
            See the description of `enhanced_semivariance` in
            :fun:`echopop.spatial.variogram.variogram`.
        max_range: Optional[realposfloat]
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
            "decay_power": decay_power,
        }

        # Drop missing `inputs`
        input_filtered = {key: value for key, value in inputs.items() if value is not None}

        # Merge the default values with those provided by the user
        # ---- Copy defaults
        params = cls.DEFAULT_VALUES.copy()
        # ---- Update
        params.update(input_filtered)

        # Filter the parameter keys
        filtered_params = {key: params[key] for key in cls.EXPECTED_DTYPES if key in params}

        # Validate the parameter datatypes
        cls.validate(filtered_params)

        return filtered_params

    # Create updating method
    @classmethod
    def update_defaults(cls, new_defaults: Dict[str, Any]):
        """
        Update the DEFAULT_VALUES attribute with new default values.

        Parameters
        ----------
        new_defaults: Dict[str, Any]
            A dictionary containing the new default values.

        Raises
        ----------
        ValueError:
            If any new default value does not match the expected datatype.
        """

        # Filter the parameter keys
        filtered_new_defaults = {
            key: new_defaults[key] for key in cls.EXPECTED_DTYPES if key in new_defaults
        }

        # Validate the parameter datatypes
        cls.validate(filtered_new_defaults)

        # Update `DEFAULT_VALUES` attribute
        cls.DEFAULT_VALUES.update(filtered_new_defaults)

    # Create default value restoration method
    @classmethod
    def restore_defaults(cls):
        """
        Restore DEFAULT_VALUES attribute to original pre-updated state.
        """

        # Restore original state
        cls.DEFAULT_VALUES.update(cls._ROOT_DEFAULT)

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
        validate_typed_dict(data, VariogramBase.EXPECTED_DTYPES)
        # FOR DEBUGGING
        # --------
        # print("Validate passed.")
        # --------


class VariogramEmpirical(TypedDict, total=False):
    # Provide the type-hints
    azimuth_range: realcircle
    force_lag_zero: bool
    standardize_coordinates: bool

    # Define default values
    DEFAULT_VALUES = {
        "azimuth_range": 360.0,
        "force_lag_zero": True,
        "standardize_coordinates": True,
    }

    # Define the expected datatypes
    EXPECTED_DTYPES = {
        "azimuth_range": realcircle,
        "force_lag_zero": bool,
        "standardize_coordinates": bool,
    }

    @classmethod
    def create(
        cls,
        azimuth_range: realcircle = 360.0,
        force_lag_zero: bool = True,
        standardize_coordinates: bool = True,
        **kwargs,
    ) -> "VariogramEmpirical":
        """
        Empirical variogram parameters

        Parameters
        ----------
        azimuth_range: float
            The total azimuth angle range that is allowed for constraining
            the relative angles between spatial points, particularly for cases where a high degree
            of directionality is assumed.
        force_lag_zero: bool
            See the `variogram_parameters` argument in
            :fun:`echopop.spatial.variogram.empirical_variogram` for more details on
            `force_lag_zero`.
        standardize_coordinates: bool
            When set to `True`, transect coordinates are standardized using reference coordinates.

        Returns
        ----------
        VariogramEmpirical: A validated dictionary with the user-defined empirical variogram
        parameter values and default values for any missing parameters/keys.
        """

        # User-defined parameters
        inputs = {
            "azimuth_range": azimuth_range,
            "force_lag_zero": force_lag_zero,
            "standardize_coordinates": standardize_coordinates,
        }

        # Drop missing `inputs`
        input_filtered = {key: value for key, value in inputs.items() if value is not None}

        # Merge the default values with those provided by the user
        # ---- Copy defaults
        params = cls.DEFAULT_VALUES.copy()
        # ---- Update
        params.update(input_filtered)

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
        validate_typed_dict(data, VariogramEmpirical.EXPECTED_DTYPES)
        # FOR DEBUGGING
        # --------
        # print("Validate passed.")
        # --------


class VariogramOptimize(TypedDict):
    max_fun_evaluations: posint
    cost_fun_tolerance: realposfloat
    gradient_tolerance: realposfloat
    finite_step_size: realposfloat
    trust_region_solver: Literal["base", "exact"]
    x_scale: Union[Literal["jacobian"], np.ndarray[realposfloat]]
    jacobian_approx: Literal["forward", "central"]

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
        "max_fun_evaluations": posint,
        "cost_fun_tolerance": realposfloat,
        "solution_tolerance": realposfloat,
        "gradient_tolerance": realposfloat,
        "finite_step_size": realposfloat,
        "trust_region_solver": Literal["base", "exact"],
        "x_scale": Union[Literal["jacobian"], np.ndarray[realposfloat]],
        "jacobian_approx": Literal["forward", "central"],
    }

    @classmethod
    def create(
        cls,
        max_fun_evaluations: posint = 500,
        cost_fun_tolerance: realposfloat = 1e-6,
        solution_tolerance: realposfloat = 1e-4,
        gradient_tolerance: realposfloat = 1e-4,
        finite_step_size: realposfloat = 1e-8,
        trust_region_solver: Literal["exact", "base"] = "exact",
        x_scale: Union[Literal["jacobian"], np.ndarray[realposfloat]] = "jacobian",
        jacobian_approx: Literal["forward", "central"] = "forward",
        **kwargs,
    ) -> "VariogramOptimize":
        """
        Base variogram model parameters

        Parameters
        ----------
        max_fun_evaluations: posint
            The maximum number of evaluations. Defaults to 500.
        cost_fun_tolerance: realposfloat
            Threshold used for determining convergence via incremental changes of the cost function.
            Defaults to 1e-6.
        solution_tolerance: realposfloat
            Threshold used for determining convergence via change of the independent variables.
            Defaults to 1e-8.
        gradient_tolerance: realposfloat
            Threshold used for determining convergence via the gradient norma. Defaults to 1e-8.
        finite_step_size: realposfloat
            The relative step sizes used for approximating the Jacobian via finite differences.
        trust_region_solver: Literal["exact", "base"]
            The method used for solving the trust-region problem by either using the Jacobian
            computed from the first iteration (`"base"`) or via singular value decomposition
            (`"exact"`). Defaults to "exact".
        x_scale: Union[Literal["jacobian"], np.ndarray[realposfloat]]
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
            "jacobian_approx": jacobian_approx,
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
        validate_typed_dict(data, VariogramOptimize.EXPECTED_DTYPES)
        # FOR DEBUGGING
        # --------
        # print("Validate passed.")


class InitialValues(TypedDict):
    """Initial values type-class"""

    min: Optional[posfloat]
    value: Optional[realposfloat]
    max: Optional[posfloat]

    # Default values for parameters
    DEFAULT_STRUCTURE = {"min": 0.0, "value": 0.0, "max": np.inf}

    # Define expected types for parameter details
    EXPECTED_DTYPES = {
        "min": Optional[posfloat],
        "value": Optional[realposfloat],
        "max": Optional[posfloat],
    }


class VariogramInitial(TypedDict, total=False):
    fit_parameters: Union[List[str], Dict[str, InitialValues]]

    # Define valid parameters
    VALID_PARAMETERS = ["correlation_range", "decay_power", "hole_effect_range", "nugget", "sill"]

    def __init__(self, fit_parameters: Union[List[str], Dict[str, InitialValues]]):
        # Initialize the instance
        self.fit_parameters = fit_parameters

    @classmethod
    def create(
        cls, fit_parameters: Union[List[str], Dict[str, InitialValues]]
    ) -> "VariogramInitial":
        """
        Create a ParameterConfig from a list of parameter names or a dictionary.
        """

        # Fill missing keys
        if isinstance(fit_parameters, list):
            fit_parameters = {
                param: InitialValues.DEFAULT_STRUCTURE.copy() for param in fit_parameters
            }
        else:
            for param, details in fit_parameters.items():
                if not isinstance(details, dict):
                    raise TypeError(f"Value for '{param}' must be a dictionary.")
                for key in InitialValues.DEFAULT_STRUCTURE:
                    if key not in details:
                        fit_parameters[param][key] = InitialValues.DEFAULT_STRUCTURE[key]

        # Validate
        cls.validate(fit_parameters)

        # Update order
        fit_parameters = {
            param: {
                key: fit_parameters[param][key] for key in InitialValues.DEFAULT_STRUCTURE.keys()
            }
            for param in fit_parameters.keys()
        }

        return cls(fit_parameters)

    @staticmethod
    def validate(params: Dict[str, InitialValues]) -> None:
        """
        Validate the input dictionary against the expected structure.
        """

        # Subset `fit_parameters` to only include the appropriate parameters
        unexpected_param = set(params.keys()).difference(VariogramInitial.VALID_PARAMETERS)
        # ---- Raise Error
        if unexpected_param:
            # ---- Create list
            unexpected_list = ", ".join(f"'{par}'" for par in unexpected_param)
            # ---- Format expected list
            expected_list = ", ".join(f"'{par}'" for par in VariogramInitial.VALID_PARAMETERS)
            # ---- Print Error
            raise ValueError(
                f"Unexpected parameter(s): {unexpected_list}. Initial values ('min', 'value' "
                f"and 'max') for variogram optimization are only accepted/valid for the "
                f"following parameters: {expected_list}."
            )

        # Validate type
        for param, details in params.items():
            try:
                validate_typed_dict(details, InitialValues.EXPECTED_DTYPES)
            except Exception as e:
                raise TypeError(f"Parameter '{param}': initial {str(e)[0].lower() + str(e)[1:]}")

        # Ensure that values are ordered appropriately
        invalid_values = [
            key
            for key, value in params.items()
            if not (value["min"] <= value["value"] <= value["max"])
        ]
        # ---- Raise Error
        if invalid_values:
            # ---- Add outer apostrophes
            invalid_list = ", ".join(f"'{par}'" for par in invalid_values)
            # ---- Print
            raise TypeError(
                f"Invalid initial values for: {invalid_list}. Values must satisfy the logic: "
                f"`min` <= value` <= `max`."
            )

        # FOR DEBUGGING
        # --------
        # print("Validate passed.")
        # --------
