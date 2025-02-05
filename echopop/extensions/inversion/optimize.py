from typing import Any, Callable, Dict, Literal, Tuple, Union

import numpy as np
from lmfit import Minimizer, Parameters
from numpy.typing import ArrayLike

from .math import generate_frequency_interval, wavenumber
from .scatterer import compute_Sv, compute_ts


# def mae(
#     prediction: ArrayLike[float],
#     measurement: ArrayLike[float],
# ):
#     """
#     Mean absolute deviation (MAD) in logarithmic space (dB)
#     """
#     # == functions/cost_functionALL.m
#     pass


# def rmse(
#     prediction: ArrayLike[float],
#     measurement: ArrayLike[float],
# ):
#     """
#     Root mean square deviation (RMSE) in logarithmic space (dB)
#     """
#     # == functions/cost_functionALL.m
#     pass


def normalize_parameters(
    scattering_parameters: Dict[str, Any], inverse: bool = False
) -> Dict[str, Any]:
    """
    Normalize the optimization parameters

    Parameters
    ----------
    scattering_parameters: Dict[str, Any]
        Dictionary comprising acoustic scattering model parameters that are to be optimized.
    inverse: bool, default: False
        Boolean flag to indicate whether the parameters should be normalized via min-max
        normalization (`inverse=True`) or denormalized (i.e. inverse normalization,
        `inverse=False`).
    """

    # Min-max normalization
    if inverse is False:
        scattering_parameters_norm = {
            key: (
                {
                    **value,
                    "initial": (
                        (value["initial"] - value["low"]) / (value["high"] - value["low"])
                        if value["high"] != value["low"]
                        else 0.0
                    ),
                }
                if isinstance(value, dict)
                else value
            )  # Skip scalar entries
            for key, value in scattering_parameters.items()
        }
    # Min-max inverse normalization
    else:
        scattering_parameters_norm = {
            key: (
                {
                    **value,
                    "initial": (
                        value["initial"] * (value["high"] - value["low"]) + value["low"]
                        if value["high"] != value["low"]
                        else value["low"]
                    ),  # Handle edge cases
                }
                if isinstance(value, dict)
                else value
            )  # Skip scalar entries
            for key, value in scattering_parameters.items()
        }

    # Return dictionary
    return scattering_parameters_norm


def simulate_Sv(
    scattering_parameters: Parameters,
    processing_parameters: Dict[str, Any],
) -> np.ndarray[float]:

    # Extract parameter values from dictionary for parsing
    parameters_dict = scattering_parameters.valuesdict()

    # Compute acoustic property metrics
    # ---------------------------------
    # Normalize the length standard deviation
    if "length_sd_norm" not in parameters_dict and "length_deviation" in parameters_dict:
        length_sd_norm = parameters_dict["length_deviation"] / parameters_dict["length_mean"]
    elif "length_sd_norm" not in parameters_dict and "length_sd_norm" in processing_parameters:
        length_sd_norm = processing_parameters["length_sd_norm"]
    else:
        length_sd_norm = parameters_dict["length_sd_norm"]

    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval(
        processing_parameters["center_frequencies"],
        length_sd_norm,
        processing_parameters["frequency_interval"],
    )

    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_center = wavenumber(
        processing_parameters["center_frequencies"], processing_parameters["water_sound_speed"]
    )
    # ---- Compute ka (center frequencies)
    ka_center = k_center * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k = wavenumber(frequencies, processing_parameters["water_sound_speed"])
    # -------- Now `ka`
    ka = k * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]

    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        parameters_dict["theta_mean"] - 3.1 * parameters_dict["theta_sd"],
        parameters_dict["theta_mean"] + 3.1 * parameters_dict["theta_sd"],
        processing_parameters["n_theta"],
    )
    theta_radians = theta_values * np.pi / 180.0

    # Compute over vector lengths
    length_values = np.linspace(
        parameters_dict["length_mean"] - 3 * (length_sd_norm * parameters_dict["length_mean"]),
        parameters_dict["length_mean"] + 3 * (length_sd_norm * parameters_dict["length_mean"]),
        processing_parameters["length_bin_count"],
    )

    # PCDWBA (TS modeling step)
    # ------
    fbs = compute_ts(
        processing_parameters["taper_order"],
        length_sd_norm,
        parameters_dict["length_mean"],
        parameters_dict["length_radius_ratio"],
        parameters_dict["radius_of_curvature_ratio"],
        theta_radians,
        k,
        ka,
        parameters_dict["g"],
        parameters_dict["h"],
        processing_parameters["n_integration"],
        processing_parameters["ni_wavelen"],
        model=processing_parameters["ts_model"],
    )

    # Compute S_V
    Sv_prediction = compute_Sv(
        parameters_dict["number_density"],
        theta_values,
        parameters_dict["theta_mean"],
        parameters_dict["theta_sd"],
        length_values,
        parameters_dict["length_mean"],
        parameters_dict["length_mean"] * length_sd_norm,
        fbs,
        ka,
        ka_center,
    )

    # Return array
    return Sv_prediction


def simulate_Sv_fit(
    scattering_parameters: Parameters,
    Sv_measured: np.ndarray[float],
    processing_parameters: Dict[str, Any],
) -> Dict[str, Any]:

    # Pre-allocate weight array
    wd = np.ones(processing_parameters["center_frequencies"].shape)

    # Compute S_V with the updated parameters
    Sv_prediction = simulate_Sv(scattering_parameters, processing_parameters)

    # Compute deviation
    error = Sv_prediction - Sv_measured

    # Compute the cost-function
    Q = np.sum(np.abs(error) * wd)

    # Return the summed absolute error
    return Q


def invert_population(
    best_fit_params: Dict[str, float],
    measurement_metadata: np.ndarray[float],
    processing_parameters: Dict[str, Any],
) -> np.ndarray[float]:  # or just a full DataFrame given the multiple estimates being calculated
    """
    Generate population estimates based on inverted TS model parameters
    """

    # Get the best-fit number density value
    nz = best_fit_params["number_density"]

    # Areal number density
    pa = nz * measurement_metadata["layer_thickness"]

    # Abundance
    abundance = nz * measurement_metadata["layer_volume"]

    # Density
    rho_z = processing_parameters["water_density"] * best_fit_params["g"]

    # Average radii
    a_ave = best_fit_params["length_mean"] / best_fit_params["length_radius_ratio"]

    # Get average volume [single animal]
    v_ave = np.pi * a_ave**2 * best_fit_params["length_mean"]

    # Compute average weight [single animal]
    wgt = v_ave * rho_z

    # Compute tonnage
    biomass_total = abundance * wgt * 1e-6  # mt

    # Compute the biomass density
    biomass_density = pa * wgt * 1e-3  # kg/m^2

    # Return array of relevant population estimates
    return np.array([nz, pa, abundance, biomass_total, biomass_density])


def prepare_optimization(
    scattering_parameters: Dict[str, Any],
    processing_parameters: Dict[str, Any],
    Sv_measured: np.ndarray[float],
) -> Tuple[Parameters, Minimizer]:
    """
    Prepare optimization settings
    """

    # Initialize the Parameters class from the `lmfit` package
    # ---- Initialize `parameters` object
    parameters = Parameters()
    # ---- Define the expected keymapping
    key_mapping = {"initial": "value", "low": "min", "high": "max", "vary": "vary"}

    # Iterate through to add to `parameters`
    _ = {
        parameters.add(param, **{key_mapping[k]: v for k, v in attrs.items() if k in key_mapping})
        for param, attrs in scattering_parameters.items()
        if isinstance(attrs, dict)
    }

    # Generate `Minimizer` function class required for bounded optimization
    minimizer = Minimizer(
        simulate_Sv_fit,
        parameters,
        fcn_args=(Sv_measured, processing_parameters),
    )

    # Return the objects
    return parameters, minimizer


def optimize_scattering_model(
    minimizer: Minimizer,
    Sv_measured: np.ndarray[float],
    optimization_parameters: Dict[str, Any],
    processing_parameters: Dict[str, Any],
) -> Tuple[Parameters, np.ndarray]:
    """
    Optimize scattering model parameters
    """

    # Minimize the cost-function to compute the best-fit/optimized variogram parameters
    parameters_optimized = minimizer.minimize(
        method="least_squares",
        **optimization_parameters,
    )

    # Predict the new S_V values
    best_fit_Sv = simulate_Sv(parameters_optimized.params, processing_parameters)

    # Compute updated fit
    Q = simulate_Sv_fit(parameters_optimized.params, Sv_measured, processing_parameters)

    # Return the best-fit parameters S_V value array, and fit
    return parameters_optimized, best_fit_Sv, Q
