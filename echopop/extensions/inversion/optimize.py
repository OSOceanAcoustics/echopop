from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
from lmfit import Minimizer, Parameters
import pandas as pd
import time
import warnings

from .math import generate_frequency_interval, inverse_normalize_series, wavenumber
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
    parameter_sets: Dict[str, Any], 
    inverse: bool = False,
    inverse_reference: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Normalize the optimization parameters

    Parameters
    ----------
    parameter_sets: Dict[str, Any]
        Dictionary comprising acoustic scattering model parameters that are to be optimized.
    inverse: bool, default: False
        Boolean flag to indicate whether the parameters should be normalized via min-max
        normalization (`inverse=True`) or denormalized (i.e. inverse normalization,
        `inverse=False`).
    """

    # Min-max normalization
    if inverse is False:
        return (
                {
                    realization: {
                        key: (
                            {
                                **value,
                                "initial": (
                                    (value["initial"] - value["low"]) 
                                    / (value["high"] - value["low"])
                                    if value["high"] != value["low"]
                                    else 0.0
                                ),
                                "low": 0.0, "high": 1.0,
                            }
                            if isinstance(value, dict)
                            else value
                        )  # Skip scalar entries
                        for key, value in sets.items()
                    }
                    for realization, sets in parameter_sets.items()
                }
        )
    # Min-max inverse normalization
    else:
        return (
            {
                key: (
                    value * (inverse_reference[key]["high"] - inverse_reference[key]["low"]) + inverse_reference[key]["low"]
                )    
                for key, value in parameter_sets.items()
            }
        )


def simulate_Sv(
    scattering_parameters: Parameters,
    config: Dict[str, Any],
) -> np.ndarray[float]:
    """
    Simulate the volumetric backscattering strength (Sv, dB re. m^-1)
    """

    # Extract parameter values from dictionary for parsing
    parameters_dict = scattering_parameters.valuesdict()

    # Rescale to the original scale
    if config["scale_parameters"]:
        parameters_dict = normalize_parameters(parameters_dict, 
                                               inverse=True, 
                                               inverse_reference=config["parameter_bounds"])

    # Compute acoustic property metrics
    # ---------------------------------
    # Normalize the length standard deviation
    if "length_sd_norm" not in parameters_dict and "length_deviation" in parameters_dict:
        length_sd_norm = parameters_dict["length_deviation"] / parameters_dict["length_mean"]
    elif "length_sd_norm" not in parameters_dict and "length_sd_norm" in config:
        length_sd_norm = config["length_sd_norm"]
    else:
        length_sd_norm = parameters_dict["length_sd_norm"]

    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval(
        config["center_frequencies"],
        length_sd_norm,
        config["frequency_interval"] * 1e3,
    )

    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_center = wavenumber(
        config["center_frequencies"], config["sound_speed_sw"]
    )
    # ---- Compute ka (center frequencies)
    ka_center = k_center * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k = wavenumber(frequencies, config["sound_speed_sw"])
    # -------- Now `ka`
    ka = k * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]

    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        parameters_dict["theta_mean"] - 3.1 * parameters_dict["theta_sd"],
        parameters_dict["theta_mean"] + 3.1 * parameters_dict["theta_sd"],
        config["orientation_bin_count"],
    )
    theta_radians = theta_values * np.pi / 180.0

    # Compute over vector lengths
    length_values = np.linspace(
        parameters_dict["length_mean"] - 3 * (length_sd_norm * parameters_dict["length_mean"]),
        parameters_dict["length_mean"] + 3 * (length_sd_norm * parameters_dict["length_mean"]),
        config["length_bin_count"],
    )

    # PCDWBA (TS modeling step)
    # ------
    fbs = compute_ts(
        config["taper_order"],
        length_sd_norm,
        parameters_dict["length_mean"],
        parameters_dict["length_radius_ratio"],
        parameters_dict["radius_of_curvature_ratio"],
        theta_radians,
        k,
        ka,
        parameters_dict["g"],
        parameters_dict["h"],
        config["n_integration"],
        config["n_wavelength"],
        model=config["ts_model"],
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
    config: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Cost-function representing the scattering model fit (Sv, dB re. 1 m^-1) 
    """

    # Pre-allocate weight array
    wd = np.ones(config["center_frequencies"].shape)

    # Compute S_V with the updated parameters
    Sv_prediction = simulate_Sv(scattering_parameters, config)

    # Compute deviation
    error = Sv_prediction - Sv_measured

    # Compute the cost-function
    Q = np.sum(np.abs(error) * wd)

    # Return the summed absolute error
    return Q

def parameterize(Sv_measured: pd.Series,
                 scattering_parameters: Dict[str, Any], 
                 simulation_parameters: Dict[str, Any],
                 center_frequencies: np.ndarray[float], 
                 sv_threshold: float,
                 processing_parameters: Dict[str, Any],
                 **kwargs) -> List[Minimizer]:
    """
    Generate `lmfit.Minimizer` object used for optimizing theroetical Sv (dB re. m^-1) based on
    measured values
    """
    
    # Prepare parameters and minimization procedure
    parameters = prepare_optimization(
        scattering_parameters, 
        simulation_parameters=simulation_parameters,
    )

    # Find which values are below the defined threshold
    valid_idx = np.argwhere(Sv_measured > sv_threshold).flatten()

    # Convert to a `numpy.ndarray`
    if not isinstance(Sv_measured, np.ndarray):
        Sv_measured = Sv_measured.to_numpy()[valid_idx]
    else:
        Sv_measured = Sv_measured[valid_idx]

    # Generate `Minimizer` function class required for bounded optimization
    # ---- This will generate per realization within a List
    return (
        [Minimizer(simulate_Sv_fit, 
                   parameters[idx], 
                   fcn_args=(Sv_measured, 
                             {**processing_parameters, 
                              **{"center_frequencies": center_frequencies[valid_idx]}}), 
                   nan_policy="omit") 
         for idx in parameters]
    )

def prepare_minimizer(
    data_df: pd.DataFrame,
    center_frequencies: np.ndarray[float],
    aggregate: str,
    sv_threshold: float,
    scattering_parameters: Dict[str, Any], 
    simulation_parameters: Dict[str, Any],    
    processing_parameters: Dict[str, Any],    
    **kwargs,
) -> pd.DataFrame:
    """
    Wrapper function that creates a list of `lmfit.Minimizer` objects for configured realizations 
    that will be used to predict Sv (dB re. m^-1)
    """

    # Initialize the minimizer column
    data_df["minimizer"] = np.nan

    # Convert to an `object` for flexible datatyping
    data_df["minimizer"] = data_df["minimizer"].astype(object)

    # Create list of `Minimizer` objects depending on number of defined realizations
    data_df["minimizer"] = data_df.apply(parameterize, 
                                         axis=1, 
                                         args=(scattering_parameters,
                                               simulation_parameters,
                                               center_frequencies,
                                               sv_threshold,
                                               processing_parameters))    

    # Define label (for verbosity)
    if aggregate == "interval":
        data_df["label"] = data_df.index.map(lambda tup: f"Transect: {tup[0]}; Interval: {tup[1]};")
    else:
        data_df["label"] = data_df.index.map(lambda tup: f"Transect: {tup};")

    # Return the dataset 
    return data_df

def estimate_population(
    inverted_data: pd.DataFrame,
    coordinates_df: np.ndarray[float],
    nasc_df: pd.DataFrame,
    layer_height_df: pd.DataFrame,
    density_sw: float,
    aggregate: Literal["interval", "transect"],
    reference_frequency: Optional[float] = None,
    **kwargs,    
) -> pd.DataFrame: 
    """
    Generate population estimates based on inverted TS model parameters
    """

    # Compute the average body radius [single animal]
    radius_mean = inverted_data["length_mean"] / inverted_data["length_radius_ratio"]
    
    # Estimate the average volume assuming an uniformly bent cylinder [single animal]
    body_volume = np.pi * radius_mean**2 * inverted_data["length_mean"]

    # Compute the animal body density [single animal]
    body_density = density_sw * inverted_data["g"]

    # Compute average body weight [single animal]
    weight_mean = body_volume * body_density

    # Subset the layer thicknesses for the correct reference frequency
    layer_thickness = layer_height_df[reference_frequency].copy()
    # ---- Drop the empty cells
    layer_thickness_filled = layer_thickness[layer_thickness > 0.].reset_index(name="thickness")

    # Apportion values if inversion done on a transect-by-transect basis
    if aggregate == "transect":
        # ---- Subset for the correct reference frequency
        nasc_df = nasc_df[reference_frequency].copy()
        # ---- Compute the summed NASC per transect
        transect_nasc = nasc_df.reset_index(name="nasc").groupby(["transect_num"])["nasc"].sum()
        # ---- Reset the index
        nasc_intervals = nasc_df.reset_index(name="nasc").set_index(["transect_num"])
        # ---- Reindex transect sums
        transect_nasc_sum = transect_nasc.reindex(nasc_intervals.index)
        # ---- Compute NASC weights
        nasc_intervals["weight"] = nasc_intervals["nasc"] / transect_nasc_sum
        # ---- Set the index
        nasc_intervals.set_index(["interval"], append=True, inplace=True)
        # ---- Calculate the summed layer height per interval
        interval_thickness = layer_thickness_filled.groupby(
            ["transect_num", "interval"]
        )["thickness"].sum().reset_index()
        # ---- Calculate the average thickness per transect
        average_thickness = interval_thickness.groupby(
            "transect_num"
        )["thickness"].mean()
        # ---- Get number of unique intervals per transect
        n_intervals = layer_thickness_filled.groupby(["transect_num"])["interval"].nunique()
        # ---- Compute the adjusted areal number density (animals nmi^-2)
        areal_number_density = (
            inverted_data["number_density"] 
            * nasc_intervals["weight"] * average_thickness * n_intervals
        ).fillna(0.0) * 1852 ** 2
    else:
        # ---- Index the filled layer thicknesses by 'interval'
        layer_thickness_filled.set_index(["transect_num", "interval", "layer"], inplace=True)
        # ---- Compute the areal number density for each layer
        layer_density = inverted_data["number_density"] * layer_thickness_filled["thickness"]
        # ---- Reset index
        layer_density = layer_density.copy().fillna(0.).reset_index(name="number_density")
        # ---- Vertically integrate across intervals
        areal_number_density = layer_density.groupby(
            ["transect_num", "interval"]
        )["number_density"].sum()
        # ---- Compute the adjusted areal number density (animals nmi^-2)
        areal_number_density = areal_number_density * 1852 ** 2

    # Convert to areal biomass density (kg nmi^-2)
    areal_biomass_density = (areal_number_density * weight_mean).fillna(0.)

    # Prepare the output
    output_df = coordinates_df.copy()
    # ---- Areal number density
    output_df["number_areal_density"] = areal_number_density
    # ---- Biomass density
    output_df["biomass_areal_density"] = areal_biomass_density
    # ---- Fill NaN values
    output_df.fillna(0., inplace=True)
    
    # Return the output DataFrame with the index reset
    return output_df.sort_index().reset_index()

def prepare_optimization(
    scattering_parameters: Dict[str, Any],
    simulation_parameters: Dict[str, Any],
) -> Parameters:
    """
    Prepare optimization settings
    """

    # Generate parameter sets
    parameter_sets = monte_carlo_initialization(scattering_parameters, simulation_parameters)    

    # Scale parameters, if defined
    if simulation_parameters["scale_parameters"]:
        # ---- Normalize
        parameter_sets = normalize_parameters(
            parameter_sets
        )
    
    # Format the dictionary into the required `lmfit.Parameters` object
    parameters = {
        realization: format_parameters(values)
        for realization, values in parameter_sets.items()
    }

    # Return the objects
    return parameters

def format_parameters(parameters_dict: Dict[str, Any]):
    
    # Initialize the Parameters class from the `lmfit` package
    # ---- Initialize `parameters` object
    parameters = Parameters()
    # ---- Define the expected keymapping
    key_mapping = {"initial": "value", "low": "min", "high": "max", "optimize": "vary"}

    # Iterate through to add to `parameters`
    _ = {
        parameters.add(param, **{key_mapping[k]: v for k, v in attrs.items() if k in key_mapping})
        for param, attrs in parameters_dict.items()
        if isinstance(attrs, dict)
    }

    # Return the `lmfit.Parameters` object
    return parameters

def monte_carlo_initialization(scattering_parameters: Dict[str, Any],
                               simulation_parameters: Dict[str, Any]) -> Dict[str, Any]:
    """
    Monte Carlo simulation of initial values for scattering parameters
    """

    # Create parameter sets for the defined number of realizations
    if simulation_parameters["monte_carlo"]:
        parameter_sets = dict(
            map(
                lambda i: (
                    i,
                    {
                        key: (
                            {
                                "initial": np.random.uniform(value["low"], value["high"]),
                                "low": value["low"],
                                "high": value["high"]
                            }
                            if value["optimize"]
                            else value
                        )
                        for key, value in scattering_parameters.items()
                    },                    
                ),
                range(simulation_parameters["n_realizations"]),
            )
        )
    # Otherwise, produce just the single-set
    else:
        parameter_sets = {0: scattering_parameters}
    # Return the parameter sets
    return parameter_sets

def optimize_scattering_model(
    minimizer: Minimizer,
    Sv_measured: np.ndarray[float],
    optimization_parameters: Dict[str, Any],
    config: Dict[str, Any],
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
    best_fit_Sv = simulate_Sv(parameters_optimized.params, config)

    # Compute updated fit
    Q = simulate_Sv_fit(parameters_optimized.params, Sv_measured, config)

    # Return the best-fit parameters S_V value array, and fit
    return parameters_optimized, best_fit_Sv, Q

def group_optimizer(
    Sv_measured: pd.Series,
    processing_parameters: Dict[str, Any],
    optimization_parameters: Dict[str, Any],
    simulation_parameters: Dict[str, Any],
    center_frequencies: np.ndarray[float], 
    sv_threshold: float,
    verbose: bool = True,
    **kwargs
):
    """
    Optimize the scattering model parameters across grouped by a defined index within a DataFrame
    """

    # Catch start time in case `verbose=True`
    start_time = time.time()
    
    # Find which values are below the defined threshold
    valid_idx = np.argwhere(Sv_measured[:-2] > sv_threshold).flatten()

    # Only run if the correct number of frequencies are valid
    if len(valid_idx) >= processing_parameters["minimum_frequency_count"]:

        # Assign message string
        frequency_msg = ""

        # Convert to a numpy array
        Sv = Sv_measured[:-2].to_numpy()[valid_idx]

        # Optimize over all realizations [initialize]
        # ---- Fit errors
        fit_errors = np.ones(simulation_parameters["n_realizations"]) * np.nan
        # ---- Parameter sets (column names/structure)
        parameter_sets = pd.DataFrame(Sv_measured.minimizer[0].params.valuesdict(), 
                                      index=range(simulation_parameters["n_realizations"])) * np.nan        

        # Iterate through the realizations
        for realization in range(simulation_parameters["n_realizations"]):
            with warnings.catch_warnings():
                warnings.filterwarnings(action="ignore", category=RuntimeWarning)
                # ---- Optimize
                best_fit_params, _, fit_errors[realization] = optimize_scattering_model(
                    Sv_measured["minimizer"][realization], 
                    Sv, 
                    optimization_parameters, 
                    {**processing_parameters, 
                    **{"center_frequencies": center_frequencies[valid_idx]}}
                )
                # ---- Create `pandas.Series` and update `parameter_sets` DataFrame
                parameter_sets.loc[realization] = pd.Series(best_fit_params.params.valuesdict())

        # Find the parameter set with the lowest Q
        best_fit_set = parameter_sets.loc[np.argmin(fit_errors)]
        
        # Add error, Q, to series
        best_fit_set["Q"] = np.min(fit_errors)
    else:
        # Assign message string
        frequency_msg = (
            f"\nWARNING: The number of frequencies with Sv > {sv_threshold} dB [{len(valid_idx)}] "
            f"was fewer than the minimum frequency count "
            f"[{processing_parameters["minimum_frequency_count"]}]. Values were not optimized."
        )
        
        # Create `pandas.Series`
        best_fit_set = pd.Series(Sv_measured.minimizer[0].params.valuesdict()) * np.nan
       
        # Add error, Q, to series
        best_fit_set["Q"] = np.nan

    # Inverse transformation if values are scaled
    if processing_parameters["scale_parameters"]:
        best_fit_set = inverse_normalize_series(best_fit_set, 
                                                processing_parameters["parameter_bounds"])

    # Catch end time in case `verbose=True`
    end_time = time.time()

    # Print results
    if verbose:
        # ---- Row label
        row = f"{Sv_measured['label']}"
        # ---- Get error value
        error_value = f" Sv error (Q): {np.round(best_fit_set['Q'], 3)} dB "
        # ---- Parameter values
        parameter_values = f"\n{best_fit_set[:-1].to_frame().T}"
        # ---- Get elapsed time (s)
        elapsed_time = f" Elapsed time: {np.round(end_time - start_time, 2)} s;" 
        # ---- Number of frequencies        
        valid_freq = (
            "[" 
            + "/".join(f"{freq * 1e-3}" for freq in center_frequencies[valid_idx]) + " kHz" 
            + "]"
        )
        
        # Print out     
        print(row + elapsed_time + error_value + valid_freq + frequency_msg + parameter_values) 
    
    # Return
    return best_fit_set
    