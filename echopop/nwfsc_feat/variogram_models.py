import inspect
import warnings
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from lmfit import Minimizer, Parameters
from scipy import special

from ..utils.validate_dict import VariogramBase, VariogramInitial, VariogramOptimize
from .mesh import griddify_lag_distances

# Set warnings filter
warnings.simplefilter("always")


# Single family models
# ---- Exponential
def exponential(distance_lags: np.ndarray, sill: float, nugget: float, correlation_range: float):
    """
    Calculates the exponential semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-(distance_lags / correlation_range))

    # Compute the exponential semivariogram
    return partial_sill * decay + nugget


# ---- Gaussian
def gaussian(distance_lags: np.ndarray, sill: float, nugget: float, correlation_range: float):
    """
    Calculates the Gaussian semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-(distance_lags**2 / correlation_range**2.0))

    # Compute the Gaussian semivariogram
    return partial_sill * decay + nugget


# ---- J-Bessel
def jbessel(distance_lags: np.ndarray, sill: float, nugget: float, hole_effect_range: float):
    """
    Calculates the J-Bessel semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - special.j0(hole_effect_range * distance_lags)

    # Compute the J-Bessel semivariogram
    return partial_sill * decay + nugget


# ---- K-Bessel
def kbessel(distance_lags: np.ndarray, sill: float, nugget: float, hole_effect_range: float):
    """
    Calculates the K-Bessel semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Avoid the case where `hole_effect_range` = 0.0
    if hole_effect_range == 0.0:
        return np.where(distance_lags == 0.0, 0.0, nugget)

    # Compute the cyclical term
    cycle = np.where(
        distance_lags / hole_effect_range < 1e-4,
        0.0,
        special.kv(1.0, (distance_lags / hole_effect_range)),
    )

    # Calculate the spatial decay term
    decay = np.where(
        distance_lags / hole_effect_range < 1e-4,
        0.0,
        1.0 - (distance_lags / hole_effect_range) * cycle,
    )

    # Compute the J-Bessel semivariogram
    return partial_sill * decay + nugget


# ---- Linear
def linear(distance_lags: np.ndarray, sill: float, nugget: float):
    """
    Calculates the linear semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Compute the linear semivariogram
    return partial_sill * distance_lags + nugget


# ---- Nugget
def nugget(distance_lags: np.ndarray, sill: float, nugget: float):
    """
    Calculates the semivariogram model comprising only a nugget effect

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    """

    # Sum together except at lag == 0.0
    return np.where(distance_lags == 0.0, 0.0, sill + nugget)


# ---- Sinc
def sinc(distance_lags: np.ndarray, sill: float, nugget: float, hole_effect_range: float):
    """
    Calculates the sine cardinal ('sinc') semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    """

    # Get machine epsilon
    eps = np.finfo(float).eps

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.sin(hole_effect_range * (distance_lags + eps)) / (distance_lags + eps)

    # Compute the sinc semivariogram
    return partial_sill * decay + nugget


# ---- Spherical
def spherical(distance_lags: np.ndarray, sill: float, nugget: float, correlation_range: float):
    """
    Calculates the spherical semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = (3.0 * distance_lags) / (2.0 * correlation_range) - distance_lags**3.0 / (
        2.0 * correlation_range**3.0
    )

    # Compute the spherical semivariogram
    return np.where(
        distance_lags < correlation_range,
        partial_sill * decay + nugget,
        sill + nugget,
    )


# Composite family models (i.e. hole-effects)
# ---- J-Bessel and Gaussian
def bessel_gaussian(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
):
    """
    Calculates the composite J-Bessel and Gaussian semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** 2))

    # Calculate the hole effect
    hole_effect = special.j0(hole_effect_range * distance_lags)

    # Compute the composite J-Bessel and Gaussian semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- J-Bessel and exponential
def bessel_exponential(
    distance_lags: np.ndarray,
    nugget: float,
    sill: float,
    correlation_range: float,
    decay_power: float,
    hole_effect_range: float,
):
    """
    Calculates the composite J-Bessel and exponential semivariogram model at defined lagged
    distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    decay_power: float
        An exponential term that is used in certain generalized exponential (or related)
        semivariogram models that modulates the ascending rate for a semivariogram.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** decay_power))

    # Calculate the hole effect
    hole_effect = special.j0(hole_effect_range * distance_lags)

    # Compute the composite J-Bessel and exponential semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- cosine and exponential
def cosine_exponential(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
    enhance_semivariance: bool,
):
    """
    Calculates the composite cosine and exponential semivariogram model at defined lagged distances
    with and without semivariance enhancement

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    enhance_semivariance: bool
        A boolean term that determines whether the correlation decay in certain cosine-related
        variogram models are enhanced (or not) are further lag distances.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay_modifier = -1.0 if enhance_semivariance is True else 1.0
    decay = decay_modifier * np.exp(-(distance_lags / correlation_range))

    # Calculate the hole effect
    hole_effect = np.cos(hole_effect_range * distance_lags)

    # Compute the composite cosine and exponential semivariogram
    return partial_sill * (1.0 - decay * hole_effect) + nugget


# ---- cosine and Gaussian
def cosine_gaussian(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
):
    """
    Calculates the composite cosine and Gaussian semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = np.exp(-((distance_lags / correlation_range) ** 2))

    # Calculate the hole effect
    hole_effect = np.cos(hole_effect_range * distance_lags)

    # Compute the composite cosine and Gaussian semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- exponential and linear
def exponential_linear(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
    decay_power: float,
):
    """
    Calculates the composite exponential and linear semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    decay_power: float
        An exponential term that is used in certain generalized exponential (or related)
        semivariogram models that modulates the ascending rate for a semivariogram.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** decay_power))

    # Calculate the hole effect
    hole_effect = 1.0 - hole_effect_range * distance_lags**decay_power

    # Compute the composite exponential and linear semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- Gaussian and linear
def gaussian_linear(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
):
    """
    Calculates the composite Gaussian and linear semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    sill: float
        The asymptotic value as lags approach infinity.
    nugget: float
        The semivariogram y-intercept that corresponds to variability at lag distances shorter than
        the lag resolution.
    correlation_range: float
        The ascending rate for the semivariogram.
    hole_effect_range: float
        The (normalized) length scale/range that 'holes' are observed, which represent 'null' (or
        very small) points compared to their neighboring lags.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** 2))

    # Calculate the hole effect
    hole_effect = 1.0 - hole_effect_range * distance_lags**2

    # Compute the composite Gaussian and linear semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# Variogram function API
VARIOGRAM_MODELS = {
    "single": {
        "exponential": exponential,
        "gaussian": gaussian,
        "jbessel": jbessel,
        "kbessel": kbessel,
        "linear": linear,
        "nugget": nugget,
        "sinc": sinc,
        "spherical": spherical,
    },
    "composite": {
        ("bessel", "exponential"): bessel_exponential,
        ("bessel", "gaussian"): bessel_gaussian,
        ("cosine", "exponential"): cosine_exponential,
        ("cosine", "gaussian"): cosine_gaussian,
        ("exponential", "linear"): exponential_linear,
        ("gaussian", "linear"): gaussian_linear,
    },
}


# Variogram wrapper function
def variogram(
    distance_lags: np.ndarray,
    variogram_parameters: Optional[Dict[str, float]] = None,
    model: Optional[Union[str, List[str]]] = None,
    **kwargs,
):
    """
    Compute the theoretical semivariogram

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: Optional[Dict[str, float]]
        An optional dictionary that can contain values for variogram model parameters (see the
        below table associated with the argument `model`). Alternatively, these parameters can be
        entered directly and are contained within `kwargs`. Possible parameters include:
            - `sill` (Sill): The asymptotic value as lags approach infinity.
            - `nugget` (Nugget): The semivariogram y-intercept that corresponds to variability
            at lag distances shorter than the lag resolution.
            - `correlation_range` (Correlation length scale/range): The ascending rate for the
            semivariogram.
            - `hole_effect_range` (Hole effect range): The (normalized) length scale/range that
            'holes' are observed, which represent 'null' (or very small) points compared to their
            neighboring lags.
            - `decay_power` (Decay term exponent): An exponential term that is used in certain
            generalized exponential (or related) semivariogram models that modulates the ascending
            rate for a semivariogram.
            - `enhance_semivariance` (Semivariance enhancement): A boolean term that determines
            whether the correlation decay in certain cosine-related variogram models are enhanced
            (or not) are further lag distances.
    model: Optional[Union[ str , list ]]
        A string or list of model names. A single name represents a single family model. Two inputs
        represent the desired composite model (e.g. the composite J-Bessel and exponential model).
        Available variogram models and their respective arguments include (alongside
        `distance_lags`):

        +----------------------------+-----------------+--------------------------+
        | :fun:`variogram`           | Input           | Parameters               |
        | model                      |                 |                          |
        +============================+=================+==========================+
        |  :fun:`exponential`        | 'exponential    | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`gaussian`           | 'gaussian'      | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`jbessel`            | 'jbessel'       | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`kbessel`            | 'kbessel'       | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`linear`             | 'linear'        | - `nugget`               |
        |                            |                 | - `sill`                 |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`nugget`             | 'nugget'        | - `nugget`               |
        |                            |                 | - `sill`                 |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`sinc`               | 'sinc'          | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`spherical`          | 'spherical'     | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`bessel_exponential` | ['bessel',      | - `sill`                 |
        |                            |  'exponential'] | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `decay_power`          |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`bessel_gaussian`    | ['bessel',      | - `sill`                 |
        |                            |  'gaussian']    | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `decay_power`          |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`cosine_exponential` | ['cosine',      | - `sill`                 |
        |                            |  'exponential'] | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        |                            |                 | - `enhance_semivariance` |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`cosine_gaussian`    | ['cosine',      | - `sill`                 |
        |                            |  'gaussian']    | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`exponential_linear` | ['exponential', | - `sill`                 |
        |                            |  'linear']      | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        |                            |                 | - `decay_power`          |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`gaussian_linear`    | ['gaussian',    | - `sill`                 |
        |                            |  'linear']      | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+

    Returns
    ----------
    variogram: np.ndarray
        An array containing the (normalized) semivariance for each lag bin.
    """

    # Determine model source
    if model is not None:
        # ---- Get the variogram arguments and function from `model`
        model_source = model
    elif variogram_parameters is not None:
        # ---- Get the variogram arguments and function from `variogram_parameters`
        model_source = variogram_parameters["model"]
    else:
        raise ValueError("Argument `model` is missing.")

    # Extract function signatures
    variogram_args, variogram_function = get_variogram_arguments(model_source)

    # Evaluate whether required function parameters are present
    if variogram_parameters is not None:
        # ---- Get input arguments
        input_args = variogram_parameters
        # ---- Use `variogram_parameters` as source
        arg_diff = set(list(variogram_args)) - set(input_args) - set(["distance_lags"])
    else:
        # ---- Get input arguments
        input_args = kwargs
        # ---- Use `kwargs` as source
        arg_diff = set(list(variogram_args)) - set(input_args) - set(["distance_lags"])

    # Raise error if any are missing
    if len(arg_diff) > 0:
        raise ValueError(
            f"The following variogram parameters are missing: {', '.join(list(arg_diff))}"
        )

    # Filter out only the variogram parameters required for the model
    required_args = dict((k, input_args[k]) for k in input_args if k in list(variogram_args))

    # Pipe the parameters into the appropriate variogram function
    return variogram_function["model_function"](distance_lags, **required_args)

def get_variogram_arguments(model_name: Union[str, List[str]]):
    """
    Get the variogram function arguments and model function for a given model.

    Parameters
    ----------
    model_name : Union[str, List[str]]
        A string or list of model names. A single name represents a single family model.
        Two inputs represent the desired composite model (e.g. the composite J-Bessel and
        exponential model).

    Returns
    -------
    Tuple[inspect.Signature.parameters, Dict[str, Any]]
        A tuple containing the function signature parameters and a dictionary with the
        model function.
    """

    # Convert to lowercase to match reference model dictionary
    if isinstance(model_name, str):
        model_input = model_name.lower()
    elif isinstance(model_name, list) & len(model_name) == 1:
        model_input = "".join(model_name).lower()
    else:
        model_input = [name.lower() for name in model_name]
        # ---- Alphabetic sort
        model_input.sort()

    # Parse user input from reference model dictionary
    # ---- Check against VARIOGRAM_MODELS API to ensure model exists
    if isinstance(model_input, list) and (tuple(model_input) in VARIOGRAM_MODELS["composite"]):
        # if (len(model_input) > 1) & (tuple(model_input) in VARIOGRAM_MODELS["composite"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["composite"][tuple(model_input)]
    # elif (len([model_input]) == 1) & (model_input in VARIOGRAM_MODELS["single"]):
    elif not isinstance(model_input, list) and (model_input in VARIOGRAM_MODELS["single"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["single"][model_input]
    else:
        raise LookupError(
            f"The model input ({model_name}) could not be matched to an"
            f" existing variogram method."
        )

    # Check input parameters against the required function arguments
    # ---- Get the function signature
    function_signature = inspect.signature(model_function)
    # ---- Create ordered dictionary of required arguments
    return function_signature.parameters, {"model_function": model_function}
