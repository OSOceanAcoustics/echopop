from typing import Union

import numpy as np
from scipy import special


# Single family models
# ---- J-Bessel
def bessel(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the J-Bessel semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - special.j0(variogram_parameters["hole_effect_range"] * distance_lags)

    # Compute the J-Bessel semivariogram
    return partial_sill * decay + variogram_parameters["nugget"]


# ---- Exponential
def exponential(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the exponential semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-(distance_lags / variogram_parameters["correlation_range"]))

    # Compute the exponential semivariogram
    return partial_sill * decay + variogram_parameters["nugget"]


# ---- Gaussian
def gaussian(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the Gaussian semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-(distance_lags**2 / variogram_parameters["correlation_range"] ** 2.0))

    # Compute the Gaussian semivariogram
    return partial_sill * decay + variogram_parameters["nugget"]


# ---- Linear
def linear(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the linear semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Compute the linear semivariogram
    return partial_sill * distance_lags + variogram_parameters["nugget"]


# ---- Sinc
def sinc(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the sine cardinal ('sinc') semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Get machine epsilon
    eps = np.finfo(float).eps

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.sin(
        (distance_lags * eps) * variogram_parameters["hole_effect_range"] / (distance_lags * eps)
    )

    # Compute the sinc semivariogram
    return partial_sill * decay + variogram_parameters["nugget"]


# ---- Spherical
def spherical(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the spherical semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = (3.0 * distance_lags) / (
        2.0 * variogram_parameters["correlation_range"]
    ) - distance_lags**3.0 / (2.0 * variogram_parameters["correlation_range"] ** 3.0)

    # Compute the spherical semivariogram
    return np.where(
        distance_lags <= variogram_parameters["correlation_range"],
        partial_sill * decay + variogram_parameters["nugget"],
        partial_sill,
    )


# Composite family models (i.e. hole-effects)
# ---- J-Bessel and Gaussian
def bessel_gaussian(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the composite J-Bessel and Gaussian semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / variogram_parameters["correlation_range"]) ** 2))

    # Calculate the hole effect
    hole_effect = special.j0(variogram_parameters["hole_effect_range"] * distance_lags)

    # Compute the composite J-Bessel and Gaussian semivariogram
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


# ---- J-Bessel and exponential
def bessel_exponential(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the composite J-Bessel and exponential semivariogram model at defined lagged
    distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(
        -(
            (distance_lags / variogram_parameters["correlation_range"])
            ** variogram_parameters["decay_power"]
        )
    )

    # Calculate the hole effect
    hole_effect = special.j0(variogram_parameters["hole_effect_range"] * distance_lags)

    # Compute the composite J-Bessel and exponential semivariogram
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


# ---- cosine and exponential
def cosine_exponential(distance_lags: np.ndarray, variogram_parameters: dict, **kwargs):
    """
    Calculates the composite cosine and exponential semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay_modifier = -1.0 if kwargs["enhance_semivariance"] is True else 1.0
    decay = decay_modifier * np.exp(-(distance_lags / variogram_parameters["correlation_range"]))

    # Calculate the hole effect
    hole_effect = np.cos(variogram_parameters["hole_effect_range"] * distance_lags)

    # Compute the composite cosine and exponential semivariogram
    return partial_sill * (1.0 - decay * hole_effect) + variogram_parameters["nugget"]


# ---- cosine and Gaussian
def cosine_gaussian(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the composite cosine and Gaussian semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = np.exp(-((distance_lags / variogram_parameters["correlation_range"]) ** 2))

    # Calculate the hole effect
    hole_effect = np.cos(variogram_parameters["hole_effect_range"] * distance_lags)

    # Compute the composite cosine and Gaussian semivariogram
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


# ---- exponential and linear
def exponential_linear(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the composite exponential and linear semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(
        -(
            (distance_lags / variogram_parameters["correlation_range"])
            ** variogram_parameters["decay_power"]
        )
    )

    # Calculate the hole effect
    hole_effect = (
        1.0
        - variogram_parameters["hole_effect_range"]
        * distance_lags ** variogram_parameters["decay_power"]
    )

    # Compute the composite exponential and linear semivariogram
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


# ---- Gaussian and linear
def gaussian_linear(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the composite exponential and linear semivariogram model at defined lagged distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / variogram_parameters["correlation_range"]) ** 2))

    # Calculate the hole effect
    hole_effect = 1.0 - variogram_parameters["hole_effect_range"] * distance_lags**2

    # Compute the composite Gaussian and linear semivariogram
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


# Variogram function API
VARIOGRAM_MODELS = {
    "single": {
        "bessel": bessel,
        "exponential": exponential,
        "gaussian": gaussian,
        "linear": linear,
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
    variogram_parameters: dict,
    model: Union[str, list] = ["bessel", "exponential"],
    **kwargs
):
    """
    Compute the theoretical semivariogram

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    model: Union[ str , list ]
        A string or list of model names. A single name represents a single family model. Two inputs
        represent the desired composite model (e.g. the composite J-Bessel and exponential model)
    """

    # Convert to lowercase to match reference model dictionary
    if isinstance(model, str):
        model_input = model.lower()
    elif isinstance(model, list) & len(model) == 1:
        model_input = "".join(model).lower()
    else:
        model_input = [name.lower() for name in model]
        # ---- Alphabetic sort
        model_input.sort()

    # Parse user input from reference model dictionary
    # ---- Check against VARIOGRAM_MODELS API to ensure model exists
    if (len(model_input) > 1) & (tuple(model_input) in VARIOGRAM_MODELS["composite"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["composite"][tuple(model_input)]
    elif (len([model_input]) == 1) & (model_input in VARIOGRAM_MODELS["single"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["single"][model_input]
    # ---- Pass the additional user arguments (kwargs) to the child function
    return model_function(distance_lags, variogram_parameters=variogram_parameters, **kwargs)
