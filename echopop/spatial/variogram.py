from typing import Union, Optional, List
import inspect 
import warnings

import pandas as pd
import numpy as np
from scipy import special
from lmfit import Parameters, Minimizer, minimize

from typing import Dict, List, Tuple, Union, Optional, Literal

from .mesh import griddify_lag_distances

# Set warnings filter
warnings.simplefilter("always")

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

    # Get the variogram arguments and function
    variogram_args, variogram_function = get_variogram_arguments(model)

    # Extract the appropriate variogram arguments from `variogram_parameters`
    parameter_values
    {k: v for k, v in 

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
                        settings_dict: dict) -> tuple:
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
    lags = variogram_parameters["distance_lags"]

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
                            azimuth_range: float) -> np.ndarray:
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
                 head_index: np.ndarray) -> tuple:
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
                                    lags: np.ndarray) -> tuple:
    
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

def get_variogram_arguments(model_name: Union[str, List[str]]):
    """
    Get the variogram function arguments
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
    if (len(model_input) > 1) & (tuple(model_input) in VARIOGRAM_MODELS["composite"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["composite"][tuple(model_input)]
    elif (len([model_input]) == 1) & (model_input in VARIOGRAM_MODELS["single"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["single"][model_input]
    else:
        raise LookupError(f"The model input ({model_name}) could not be matched to an" 
                          f" existing variogram method.")
    
    # Check input parameters against the required function arguments
    # ---- Get the function signature
    function_signature = inspect.signature(model_function)
    # ---- Create ordered dictionary of required arguments
    return function_signature.parameters, {f"model_function": model_function}


def create_optimization_options(fit_parameters: Union[str, List[str], Dict[str, Dict[str, float]]],
                                default_variogram_settings: Dict[str, float],
                                model: Union[str, List[str]],
                                initial_values: Optional[Union[List[Tuple[str, float]], Dict[str, Dict[str, float]]]] = None,
                                lower_bounds: Optional[Union[List[Tuple[str, float]], Dict[str, Dict[str, float]]]] = None,
                                upper_bounds: Optional[Union[List[Tuple[str, float]], Dict[str, Dict[str, float]]]] = None,
                                max_fun_evaluations: Optional[int] = None,
                                cost_fun_tolerance: Optional[float] = None,
                                solution_tolerance: Optional[float] = None,
                                gradient_tolerance: Optional[float] = None,
                                finite_step_size: Optional[float] = None,
                                trust_region_solver: Optional[Literal["exact", "base"]] = None,
                                x_scale: Optional[Union[Literal["jacobian"], np.ndarray[float]]] = None,
                                jacobian_approx: Optional[Literal["forward", "central"]] = None) -> dict:
    """
    Construct the variogram optimization parameter dictionary
    """ 

    # Convert to a list if `parameter_values` is a single-parameter string
    if isinstance(fit_parameters , str):
        fit_parameters  = [fit_parameters]

    # Construct the parameter values
    # ---- If `parameter_values` is already a correctly formatted dict
    if isinstance(fit_parameters, dict):
        _parameters = fit_parameters 
    elif isinstance(fit_parameters, list):
        _parameters = {}
        # ---- Iteratively populate `_parameters` with the varied parameters
        _parameters = {param: {} for param in fit_parameters}
    else: 
        raise TypeError("Argument `fit_parameters` must either be a single string, a list of strings, or a pre-formatted"
                        " dictionary. Please see the documentation for `create_optimization_options` for futher details.")

    # Helper function for populating `_parameters`
    def populate_parameters(argument, key, name):
        # ---- Check argument existence
        if argument: 
            if isinstance(argument, dict):
                _argument = {param: {key: value} for param, value in (argument.items())}
            elif isinstance(argument, list):
                if (all(isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], str) 
                        and isinstance(item[1], float) for item in argument)):
                    _argument = {param: value for param, value in argument}
                else:
                    raise TypeError(f"The list of values for `{name}` must be a list of tuples (str, float).")
            else: 
                raise TypeError(f"Argument `{name} must either be either a dictionary or a list.")
        
            # Update `_parameters` with `_argument`
            for param, value in zip(fit_parameters, _argument.items()):
                if isinstance(value, tuple):
                    _parameters[param][key] = value[1]
                else:
                    _parameters[param][key] = value

    # Populate lower bounds
    populate_parameters(lower_bounds, "min", "lower_bounds")

    # Populate upper bounds
    populate_parameters(upper_bounds, "max", "upper_bounds")

    # Populate starting values
    if initial_values is None:
        # ---- Initialize empty list
        init_values = []
        # ---- Create empty fixed list
        fixed_parameters = []
        # ---- If no initial values are provided, use default values
        for param in fit_parameters:
            init_values.append((param, default_variogram_settings[param]))
    else:
        # ---- Create copy of `initial_values`
        init_values = initial_values.copy()
        # ---- Get the required parameters
        variogram_args, variogram_function = get_variogram_arguments(model)
        # ---- Get names
        arg_names = list(variogram_args.keys())
        # ---- Drop `distance_lags`
        arg_names.remove("distance_lags")
        # ---- Assign a list of parameters that will be kept fixed
        fixed_parameters = []
        # ---- Get the initial values present
        init_parameters = {param for param, value in initial_values}
        # ---- Find missing values that are contained within the `fit_parameters` list
        missing_parameters = set(arg_names) - set(init_parameters)
        # ---- Append default values to `inital_values`
        if missing_parameters:
            for param in missing_parameters:
                fixed_parameters.append(param)
                init_values.append((param, default_variogram_settings[param]))
    # ---- And populate
    populate_parameters(init_values, "value", "initial_values")

    # Populate fixed constants, if any
    if fixed_parameters:
        for param in fixed_parameters:
            _parameters[param]["vary"] = False

    # Initial a Parameters class from the `lmfit` package
    # ---- Initialize `paramaeters` object
    parameters = Parameters()
    # ---- Populate `parameters`
    for param in _parameters.keys():
        parameters.add(param, **_parameters[param])

    # Save parameter names that are present
    # ---- Parameter names
    optimization_parameters = ["max_fun_evaluations", "cost_fun_tolerance", "solution_tolerance",
                               "gradient_tolerance", "finite_step_size", "trust_region_solver", 
                               "x_scale", "jacobian_approx"]
    # ---- Parameter values
    optimization_values = [max_fun_evaluations, cost_fun_tolerance, solution_tolerance, gradient_tolerance,
                           finite_step_size, trust_region_solver, x_scale, jacobian_approx]
    # ---- Drop unused 
    unused_parameters_bitmap = [k is not None for k in optimization_values]
    used_parameters = [value for bit, value in zip(unused_parameters_bitmap, optimization_parameters) if bit]
    used_values = [value for bit, value in zip(unused_parameters_bitmap, optimization_values) if bit]

    # Internal API (using argument names for the actual underlying functions)
    _options = {
        "max_nfev": max_fun_evaluations,
        "ftol": cost_fun_tolerance,
        "xtol": solution_tolerance,
        "gtol": gradient_tolerance,
        "diff_step": finite_step_size
    }

    # Filter
    _options = {k: v for k, v in _options.items() if v is not None}

    # Add trust region solver definition
    if trust_region_solver is not None:
        _options.update({"tr_solver": None if trust_region_solver == "base" else "exact"})

    # Add x-scale
    if x_scale is not None:
        _options.update({"x_scale": "jac" if x_scale == "jacobian" else x_scale})
    
    # Add Jacobian approximation method
    if jacobian_approx is not None:
        _options.update({"jac": "2-point" if jacobian_approx == "forward" else "3-point"})

    # Create and return a dictionary representing the full optimization parameterization
    return ( {"parameters": parameters, "fixed_parameters": fixed_parameters, "config": _options, 
            "used_parameters": pd.DataFrame(zip(used_parameters, used_values), columns=["names", "values"])} )


def validate_variogram_parameters(variogram_parameters: dict,
                                  fit_parameters: Optional[list[str]] = None,
                                  optimization_settings: Optional[dict] = None) -> None:
    """
    Validate variogram parameters and optimization settings
    """

    # Parse the variogram parameters dictionary for the model name
    if "model" not in variogram_parameters:
        raise KeyError("No variogram model was defined in the `variogram_parameters` dictionary."
                       " Please see the documentation for `variogram()` for valid model inputs.")
    else:
        model_name = variogram_parameters["model"]    

    # Validate that all required parameters are present
    # ---- Extract
    function_arguments, _ = get_variogram_arguments(model_name)
    # ----Checks against the user-input (`variogram_parameters`) and the initialized 
    # ---- dictionary (`init_parameters`) when present
    init_parameters = optimization_settings["fixed_parameters"]
    missing_args = []
    # ---- Iterate through
    for param in function_arguments.values():
        # ---- Ignore "distance_lags"
        if param.name != "distance_lags":
            # ---- For cases where there is no default argument
            if param.default is param.empty:
                # ---- For cases where there is no user-input
                if param.name not in variogram_parameters:
                    # ---- For cases where a value is not initialized 
                    if init_parameters is not None and param.name not in init_parameters:
                        missing_args.append(param.name)
                    else:  
                        missing_args.append(param.name)
    # ---- Generate error if any are missing
    if missing_args:
        raise ValueError(f"Missing variogram parameters: {', '.join(missing_args)}.")

    # Evaluate whether arguments defined in `fit_parameters` exist
    # ---- Find arguments that do not match any arguments 
    if fit_parameters is not None:
        erroneous_args = set(fit_parameters).difference(list(function_arguments.keys()))
        # ---- Generate warning 
        if list(erroneous_args): 
            warnings.warn( 
                f"Unnecessary variogram parameters included in the argument `fit_parameters`:"
                f" {', '.join(list(erroneous_args))}. These will be ignored.",
                stacklevel = 1
            )

    # Validate the optimization settings entries, if necessary
    if optimization_settings is not None:
        # ---- Pivot parameters
        present_parameters = optimization_settings["used_parameters"].set_index("names").T
        # ---- Validate that `optimization_settings` is a dictionary
        if not isinstance(optimization_settings, dict):
            raise TypeError("The argument `optimization_settings` must be a dictionary. Please see"
                            " documentation for `fit_variogram()` for valid optimization parameters/"
                            "arguments.")
        # ---- Check for any erroneous entries
        erroneous_args = (
            set(present_parameters.columns) 
            - set(["max_fun_evaluations", "cost_fun_tolerance", "solution_tolerance", 
                "gradient_tolerance", "finite_step_size", "trust_region_solver", 
                "x_scale", "jacobian_approx"])
        )
        # ---- Generate warning if needed
        if list(erroneous_args):
            raise KeyError(f"User optimization arguments ({', '.join(list(erroneous_args))})"
                        f" are undefined. Please see documentation for `fit_variogram()`"
                        f" for valid optimization parameters/arguments.")
        # ---- Check `max_fun_evaluations`
        if "max_fun_evaluations" in present_parameters.columns:
            if not isinstance(present_parameters["max_fun_evaluations"]["values"], int):
                raise TypeError("The optimization parameter `max_fun_evaluations` must be an integer.")
            elif present_parameters["max_fun_evaluations"]["values"] <= 0:
                raise ValueError("The optimization parameter `max_fun_evaluations` must be a positive," 
                                " non-zero integer.")
        # ---- Check `cost_fun_tolerance`
        if "cost_fun_tolerance" in present_parameters.columns:
            if not isinstance(present_parameters["cost_fun_tolerance"]["values"], float):
                raise TypeError("The optimization parameter `cost_fun_tolerance` must be a float.")
        # ---- Check `solution_tolerance`
        if "solution_tolerance" in present_parameters.columns:
            if not isinstance(present_parameters["solution_tolerance"]["values"], float):
                raise TypeError("The optimization parameter `solution_tolerance` must be a float.")
        # ---- Check `gradient_tolerance`
        if "gradient_tolerance" in present_parameters.columns:
            if not isinstance(present_parameters["gradient_tolerance"]["values"], float):
                raise TypeError("The optimization parameter `gradient_tolerance` must be a float.")
        # ---- Check `finite_step_size`
        if "finite_step_size" in present_parameters.columns:
            if not isinstance(present_parameters["finite_step_size"]["values"], float):
                raise TypeError("The optimization parameter `finite_step_size` must be a float.")
        # ---- Check `trust_region_solver`
        if "trust_region_solver" in present_parameters.columns:
            if not isinstance(present_parameters["trust_region_solver"]["values"], str):
                raise TypeError("The optimization parameter `trust_region_solver` must be a string.")
            elif present_parameters["trust_region_solver"]["values"] not in ["base", "exact"]:
                raise ValueError(f"The optimization parameter `trust_region_solver` "
                                f"`{present_parameters['trust_region_solver']['values']}` is invalid. This "
                                f"value must be one of the following: [None, 'exact'].")
        # ---- Check `x_scale`
        if "x_scale" in present_parameters.columns:
            if not isinstance(present_parameters["x_scale"]["values"], (float, str)):
                raise TypeError(f"The optimization parameter `x_scale` must either be a single float"
                                f" or be a string.")
            elif isinstance(present_parameters["x_scale"]["values"], str):
                if present_parameters["x_scale"]["values"] != "jacobian":
                    raise ValueError(f"The optimization parameter `x_scale` ({present_parameters['x_scale']['values']})"
                                    f" is invalid. The only string value accepted is: 'jacobian'.")
        # ---- Check `jacobian_approx`
        if "jacobian_approx" in present_parameters.columns:
            if not isinstance(present_parameters["jacobian_approx"]["values"], str):
                raise TypeError("The optimization parameter `jacobian_approx` must be a string.")
            elif present_parameters["jacobian_approx"]["values"] not in ["forward", "central"]:
                raise ValueError(f"The optimization parameter `jacobian_approx` ({present_parameters['jacobian_approx']['values']})"
                                f" is invalid. The only accepted values are: ['forward', 'central'].")

def optimize_variogram(lag_counts: np.ndarray,
                       lags: np.ndarray,
                       gamma_h: np.ndarray,
                       variogram_parameters: dict,
                       optimization_settings: dict,
                       settings_dict: dict):

    # Compute the lag weights
    lag_weights = lag_counts / lag_counts.sum()

    # Vertically stack the lags, semivariance, and weights
    data_stack = np.vstack((lags, gamma_h, lag_weights))

    # Index lag distances that are within the parameterized range
    within_range = np.where(lags <= variogram_parameters["range"])[0]

    # Truncate the data stack
    truncated_stack = data_stack[:, within_range]

    # Get model name
    _, variogram_fun = get_variogram_arguments(variogram_parameters["model"])

    # Create helper cost function
    def cost_function(parameters, x, y, w):
        yr = variogram_fun["model_function"](x, **parameters)
        return (yr - y) * w
    # get initial fit
    init_fit = cost_function(optimization_settings["parameters"], x=truncated_stack[0], y = truncated_stack[1], w = truncated_stack[2])
    rmse_init = np.mean(np.abs(init_fit))

    minimizer = Minimizer(cost_function, optimization_settings["parameters"], fcn_args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]))
    result = minimizer.minimize(method="least_squares", **optimization_settings["config"])
    new_fit = np.mean(np.abs(result.residual))
    best_fit_params = result.params.valuesdict()
    
    if settings_dict["verbose"]:

        params = [k for k in optimization_settings["parameters"].keys()]

        # ---- Create updated values for printing
        df = pd.DataFrame({"parameter": params})
        df["initial"] = list(optimization_settings["parameters"].valuesdict().values())
        df["new"] = list(best_fit_params.values())
        lst = [f"{name.replace("_", " ").capitalize()}: {'{:.3}'.format(init)} -> {'{:.3}'.format(new)}" for name, init, new in zip(df["parameter"], df["initial"], df["new"])]
        app = "\n".join(item for item in lst)
        print(
            f"-----------------------------\n"
            f"VARIOGRAM OPTIMIZATION\n"
            f"-----------------------------\n"
            f"| See `self.analysis['settings']['variogram']['optimization']"
            f" for parameter settings.\n"
            f"-----------------------------\n"
            f"| Initial -> Updated\n"
            f"-----------------------------\n"
            f"Overall fit [MAD]: {'{:.3}'.format(rmse_init)} -> {'{:.3}'.format(new_fit)}\n"
            f"{app}"
        )
    
    # Return result
    return best_fit_params
