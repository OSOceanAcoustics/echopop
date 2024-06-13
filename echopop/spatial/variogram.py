from typing import Union, Optional
import inspect 

import pandas as pd
import numpy as np
from scipy import special

from .mesh import griddify_lag_distances

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
def bessel_exponential(distance_lags: np.ndarray, nugget: np.ndarray, sill: np.ndarray,
                       correlation_range: np.ndarray, decay_power:np.ndarray,
                       hole_effect_range: np.ndarray):
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
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(
        -(
            (distance_lags / correlation_range)
            ** decay_power
        )
    )

    # Calculate the hole effect
    hole_effect = special.j0(hole_effect_range * distance_lags)

    # Compute the composite J-Bessel and exponential semivariogram
    return partial_sill * (decay * hole_effect) + nugget


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

def empirical_variogram(transect_data: pd.DataFrame,
                        variogram_parameters: dict,
                        settings_dict: dict):
    """
    Compute the empirical variogram from transect data.
    """
    
    # Convert the estimate column to an array
    # ---- estimates
    estimates = transect_data[settings_dict['variable']].to_numpy()
    
    # Extract relevant variogram parameters
    # ---- Lag resolution ['lcsl']
    lag_resolution = variogram_parameters["lag_resolution"]
    # ---- Number of lags 
    n_lags = variogram_parameters["n_lags"]
    # ---- Compute the lags ['h']
    lags = np.arange(1, n_lags) * lag_resolution

    # Calculate the lag distance matrix among transect data
    transect_distance_matrix, transect_azimuth_matrix = (
        griddify_lag_distances(transect_data, transect_data, angles = True)
    )
    # ---- Convert to lags
    lag_matrix = (np.round(transect_distance_matrix / lag_resolution).astype(int) + 1)

    # Pre-allocate vectors/arrays that will be iteratively filled 
    # ---- Counts for each lag
    lag_counts = np.zeros(n_lags - 1)
    # ---- Summed deviation for each lag 
    lag_deviations = np.zeros_like(lag_counts)
    # ---- Summmed estimates for each lag
    lag_estimates = np.zeros_like(lag_counts)
    # ---- Summed squared estimates for each lag
    lag_estimates_squared = np.zeros_like(lag_counts)
    # ---- Head (or first) observation indices
    head_index = np.zeros((len(estimates), n_lags - 1))

    # Create a triangle mask with the diaganol offset to the left by 1
    # ---- Initial mask
    triangle_mask = np.tri(len(estimates), k = -1, dtype = bool)
    # ---- Vertically and then horizontally flip to force the 'True' and 'False' positions
    triangle_mask_flp = np.flip(np.flip(triangle_mask), axis = 1)

    # Tally the counts of each lag present
    # ---- Filter the lag matrix and convert to a 1D array
    equivalent_lags = variogram_matrix_filter(lag_matrix, 
                                              triangle_mask_flp, 
                                              transect_azimuth_matrix, 
                                              variogram_parameters["azimuth_range"])
    # ---- Compute the binned sum
    lag_counts = np.bincount(equivalent_lags)[1:n_lags]
    
    # Compute the summed estimates per lag
    # ---- Filter the field estimates
    estimates_filtered = variogram_matrix_filter(estimates, 
                                                 triangle_mask_flp, 
                                                 transect_azimuth_matrix, 
                                                 variogram_parameters["azimuth_range"])
    # ---- Compute the binned sum
    lag_estimates = np.bincount(equivalent_lags, weights=estimates_filtered)[1:n_lags]
    # ---- Compute the binned squared-sum
    lag_estimates_squared = np.bincount(equivalent_lags, weights=estimates_filtered ** 2)[1:n_lags]
    
    # Calculate the deviations within each lag
    # ---- Subset the lag array via a boolean bitmap
    lag_bitmap = (equivalent_lags < n_lags)
    # ---- Create a dummy array that produces the row indices for the estimate matrix/array and then
    # ---- apply the triangle mask and azimuth filter
    estimate_rows = variogram_matrix_filter(np.arange(len(estimates))[:, np.newaxis], 
                                            triangle_mask_flp, 
                                            transect_azimuth_matrix, 
                                            variogram_parameters["azimuth_range"])
    # ---- Calculate the deviations between indexed estimates and lag-specific ones
    deviations = (estimates[estimate_rows][lag_bitmap] - estimates_filtered[lag_bitmap]) ** 2
    # ---- Sum for each lag bin
    lag_deviations = np.bincount(equivalent_lags[lag_bitmap], weights=deviations)[1:n_lags]

    # Compute the mean and standard deviation of the head estimates for each lag bin
    # ---- Apply a mask using the triangle bitmap
    head_mask = np.where(triangle_mask_flp, lag_matrix, -1)
    # ---- Helper function for computing the binned summations for each row
    def bincount_row(row, n_lags):
        return np.bincount(row[row != -1], minlength=n_lags)[1:n_lags]
    # ---- Find the head indices of each lag for each row
    head_index = np.apply_along_axis(bincount_row, axis=1, arr=head_mask, n_lags=n_lags)

    # Compute the standardized semivariance [gamma(h)]
    gamma_h, lag_covariance = semivariance(estimates,
                                           lag_estimates,
                                           lag_estimates_squared,
                                           lag_counts,
                                           lag_deviations,
                                           head_index)
    
    # Prepend a 0.0 and force the nugget effect to be 0.0, if necessary
    # ---- Return the computed lags, empirical variogram estimate [gamma(h)], and lag counts
    if variogram_parameters["force_lag_zero"]:
        return (
            np.concatenate([[0], lags]),
            np.concatenate([[0.0], gamma_h]),
            np.concatenate([[len(estimates) - 1], lag_counts]),
            lag_covariance
        )
    else:
        return lags, gamma_h, lag_counts, lag_covariance


def variogram_matrix_filter(data_matrix: np.ndarray,
                            mask_matrix: np.ndarray,
                            azimuth_matrix: np.ndarray,
                            azimuth_range: float):
    """
    Apply a triangle and azimuth filter to a data matrix required for computing the empirical 
    variogram
    """

    # Convert array to matrix, if needed
    if data_matrix.ndim == 1:
        data_matrix = np.tile(data_matrix, (len(data_matrix), 1))
    else:
        if data_matrix.shape != azimuth_matrix.shape:
            # ---- Determine which dimension is mismatched
            dimension_diff = np.where(np.array(data_matrix.shape) != np.array(mask_matrix.shape))[0]
            if dimension_diff == 0:
                data_matrix = np.tile(data_matrix, (len(data_matrix), 1))
            else:
                data_matrix = np.tile(data_matrix, (1, len(data_matrix)))

    # Define azimuth angle threshold
    azimuth_threshold = 0.5 * azimuth_range

    # Replace any azimuth NaN values with 0's, if necessary
    azimuth_matrix[np.isnan(azimuth_matrix)] = 0.0

    # Create the azimuth angle bitmap
    azimuth_bitmap = (azimuth_matrix >= - azimuth_threshold) & (azimuth_matrix < azimuth_threshold)

    # Mask the data matrix and broadcast out into a 1D array
    return data_matrix[mask_matrix & azimuth_bitmap]

def semivariance(estimates: np.ndarray,
                 lag_estimates: np.ndarray,
                 lag_estimates_squared: np.ndarray,
                 lag_counts: np.ndarray,
                 lag_deviations: np.ndarray,
                 head_index: np.ndarray):
    """
    Compute the standardized semivariance.
    """
    
    # Calculate the mean head estimate per lag bin
    mean_head = (estimates[:, np.newaxis] * (head_index / lag_counts)).sum(axis = 0)

    # Calculate the standard deviation of head values per lag
    sigma_head = np.sqrt(
        ((estimates[:, np.newaxis] - mean_head) ** 2 * (head_index / lag_counts)).sum(axis = 0)
    )

    # Calculate the global mean and variance for each lag bin
    # ---- Mean
    lag_means = lag_estimates / lag_counts
    # ---- Variance
    lag_variance = lag_estimates_squared / lag_counts - lag_means ** 2

    # Estimate the standard deviation of tail estimates
    sigma_tail = np.sqrt(np.abs(lag_variance))

    # Calculate the semivariance
    # ---- Compute the partial sill that is applied as a weighted calculation
    partial_sill = sigma_tail * sigma_head
    # ---- Semivariance [gamma(h)] and cross-sill estimate
    gamma_h = 0.5 * lag_deviations / (lag_counts * partial_sill)

    # Calculate the mean lag distance covariance
    # ---- Find non-zero head and tail variances
    non_zero_variance = np.where((sigma_head > 0.0) & (sigma_tail > 0.0))[0]
    # ---- Mean lag covariance
    mean_lag_covariance = (sigma_head[non_zero_variance] * sigma_tail[non_zero_variance]).mean()
    return gamma_h, mean_lag_covariance

def initialize_variogram_parameters(gamma_h: np.ndarray,
                                    lags: np.ndarray):
    
    """
    Approximate appropriate bounds for each parameter that are used for optimization.
    """
    
    # Compute empirical boundary estimates for the sill
    # ---- Lower
    sill_min = 0.7 * gamma_h.max()
    # ---- Upper
    sill_max = 1.2 * gamma_h.max() 

    # Compute empirical boundary for the nugget
    # ---- Lower
    nugget_min = 0.0
    # ---- Upper
    nugget_max = sill_max

    # Estimate the length scale
    # ---- Find the maximum value within the first half of the variogram
    semivariogram_max = np.argmax(gamma_h[:int(np.round(0.5 * len(gamma_h)))])
    # ---- Find the lag distance corresponding to the 6 dB attenuation from the maximum
    lag_6db = np.argmin(np.abs(gamma_h[:semivariogram_max] - 0.3))
    # ---- Assign the length scale
    if lags[lag_6db] <= 0.01 * lags.max():
        # ---- Minimum value
        length_scale = 0.1 * lags.max()
    else:
        length_scale = lags[lag_6db]

    # Return bounds
    return (sill_min, sill_max), (nugget_min, nugget_max), length_scale

def map_args_to_variogram(func, variogram_parameters, init_parameters: Optional[dict] = None, **kwargs):
    """
    Map user input dictionary to function arguments based on the variogram function's signature.

    Parameters:
    - func: The target function.
    - variogram_parameters: A dictionary containing user-provided arguments.

    Returns:
    - A dictionary of arguments mapped to the function's parameters.
    """

    # Parse the variogram parameters dictionary for the model name
    model_name = variogram_parameters["model"]
    # ---- Search for correct model
    if len(model_name) > 1:
        args = inspect.signature(VARIOGRAM_MODELS["composite"][tuple(model_name)])
    else:
        args = inspect.signature(VARIOGRAM_MODELS["single"][model_name])
    
    # Get argument values/names
    arg_names = args.parameters.values()

    # Validate that all required parameters are present
    missing_args = []
    # ---- Iterate through
    for param in arg_names:
        # ---- Ignore "distance_lags"
        if param.name != "distance_lags":
            if param.default is param.empty and param.name not in variogram_parameters:
                if param.name not in kwargs:
                    if init_parameters is not None and param.name not in init_parameters:
                        missing_args.append(param.name)
                    elif init_parameters is None:
                        missing_args.append(param.name)
    # ---- Generate error if any are missing
    if missing_args:
        raise ValueError(f"Missing variogram parameters: {', '.join(missing_args)}.")
    
    # Gather parameter values
    output = {}
    
    for param in arg_names:
        if init_parameters is not None and param.name in init_parameters:
            output[param.name] = init_parameters[param.name]
        elif param.name in kwargs:
            output[param.name] = kwargs[param.name]
        else:
            output[param.name] = variogram_parameters[param.name]
    
    # Update variogram parameters


    # Get the function's signature
    sig = inspect.signature(func)
    
    # Create a dictionary of arguments to pass to the function
    kwargs = {}
    
    for param in sig.parameters.values():
        if param.name in user_input:
            kwargs[param.name] = user_input[param.name]
        elif param.default is param.empty:
            raise ValueError(f"Missing required argument: {param.name}")

    return kwargs
