from typing import Union

import numpy as np
from scipy import special


# Single-family models
def bessel(distance_lags: np.ndarray, variogram_parameters, **kwargs):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = 1.0 - special.j0(variogram_parameters["hole_effect_range"] * distance_lags)

    #
    return partial_sill * decay + variogram_parameters["nugget"]


def exponential(distance_lags: np.ndarray, variogram_parameters, **kwargs):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = 1.0 - np.exp(-(distance_lags / variogram_parameters["correlation_range"]))
    #
    return partial_sill * decay + variogram_parameters["nugget"]


def gaussian(distance_lags: np.ndarray, variogram_parameters, **kwargs):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = 1.0 - np.exp(-(distance_lags**2 / variogram_parameters["correlation_range"] ** 2.0))

    #
    return partial_sill * decay + variogram_parameters["nugget"]


def linear(distance_lags: np.ndarray, variogram_parameters, **kwargs):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    return partial_sill * distance_lags + variogram_parameters["nugget"]


def spherical(distance_lags: np.ndarray, variogram_parameters, **kwargs):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = (3.0 * distance_lags) / (
        2.0 * variogram_parameters["correlation_range"]
    ) - distance_lags**3.0 / (2.0 * variogram_parameters["correlation_range"] ** 3.0)

    #
    return np.where(
        distance_lags <= variogram_parameters["correlation_range"],
        partial_sill * decay + variogram_parameters["nugget"],
        partial_sill,
    )


# Composite-family models
def bessel_gaussian(distance_lags: np.ndarray, variogram_parameters, **kwargs):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = 1.0 - np.exp(-((distance_lags / variogram_parameters["correlation_range"]) ** 2))

    #
    hole_effect = special.j0(variogram_parameters["hole_effect_range"] * distance_lags)

    #
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


def bessel_exponential(
    distance_lags: np.ndarray, variogram_parameters, decay_power: np.float64 = 1.0, **kwargs
):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = 1.0 - np.exp(
        -((distance_lags / variogram_parameters["correlation_range"]) ** decay_power)
    )

    #
    hole_effect = special.j0(variogram_parameters["hole_effect_range"] * distance_lags)

    #
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


def cosine_exponential(
    distance_lags: np.ndarray, variogram_parameters, enhance_semivariance=False, **kwargs
):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay_modifier = -1.0 if enhance_semivariance else 1.0
    decay = decay_modifier * np.exp(-(distance_lags / variogram_parameters["correlation_range"]))

    #
    hole_effect = np.cos(variogram_parameters["hole_effect_range"] * distance_lags)

    #
    return partial_sill * (1.0 - decay * hole_effect) + variogram_parameters["nugget"]


def cosine_gaussian(distance_lags: np.ndarray, variogram_parameters, **kwargs):

    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = np.exp(-((distance_lags / variogram_parameters["correlation_range"]) ** 2))

    #
    hole_effect = np.cos(variogram_parameters["hole_effect_range"] * distance_lags)

    #
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


def exponential_linear(
    distance_lags: np.ndarray, variogram_parameters, decay_power: np.float64 = 1.0, **kwargs
):
    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay_function = 1.0 - np.exp(
        -((distance_lags / variogram_parameters["correlation_range"]) ** decay_power)
    )

    #
    hole_effect = 1.0 - variogram_parameters["hole_effect_range"] * distance_lags**decay_power

    #
    return partial_sill * (decay_function * hole_effect) + variogram_parameters["nugget"]


def gaussian_linear(
    distance_lags: np.ndarray, variogram_parameters, decay_power: np.float64 = 1.0, **kwargs
):
    #
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    #
    decay = 1.0 - np.exp(-((distance_lags / variogram_parameters["correlation_range"]) ** 2))

    #
    hole_effect = 1.0 - variogram_parameters["hole_effect_range"] * distance_lags**decay_power

    #
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


# Dictionary containing available variogram models for user input
VARIOGRAM_MODELS = {
    "single": {
        "bessel": bessel,
        "exponential": exponential,
        "gaussian": gaussian,
        "linear": linear,
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


def variogram(
    distance_lags,
    variogram_parameters,
    model: Union[str, list] = ["bessel", "exponential"],
    **kwargs
):

    # Convert to lowercase to match reference model dictionary
    # And then convert to lowercase
    if isinstance(model, str):
        model_input = model.lower()
    elif isinstance(model, list) & len(model) == 1:
        model_input = "".join(model).lower()
    else:
        model_input = [name.lower() for name in model]

        # Alphabetic sort
        model_input.sort()

    # Parse user input from reference model dictionary
    # Check against VARIOGRAM_MODELS options to ensure model exists
    if (len(model_input) > 1) & (tuple(model_input) in VARIOGRAM_MODELS["composite"]):

        # Parse model function
        model_function = VARIOGRAM_MODELS["composite"][tuple(model_input)]

    elif (len([model_input]) == 1) & (model_input in VARIOGRAM_MODELS["single"]):

        # Parse model function
        model_function = VARIOGRAM_MODELS["single"][model_input]

    # Pass the additional user arguments (kwargs) to the child function
    return model_function(distance_lags, variogram_parameters=variogram_parameters, **kwargs)
