import abc
import numpy as np
import pandas as pd
from typing import Union, Dict, List, Optional, Any
from functools import reduce
import pytest

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs
from echopop.nwfsc_feat import utils
from typing import Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate

#######
azimuth_range = 360.
n_lags = 30
lag_resolution = 0.002
sill = 0.91
hole_effect_range = 0.
correlation_range = 0.007
decay_power = 1.5
transect_df = df_nasc_all_ages.copy()
estimates = transect_df["biomass"].to_numpy()
#######

#######
spatial.filter_lag_matrix = filter_lag_matrix
#######
transect_df: pd.DataFrame = df_nasc_all_ages.copy()
variable: str = "biomass_density"
coordinates: List[str] = ["x", "y"]
azimuth_filter: bool = True
azimuth_angle_threshold: float = 180.
lag_resolution = 0.002
n_lags = 30
force_lag_zero: bool = True


##################################################
def empirical_variogram(
    transect_df: pd.DataFrame,
    n_lags: int,
    lag_resolution: float,
    azimuth_filter: bool,
    azimuth_angle_threshold: float,
    variable: str = "biomass_density",
    coordinates: Tuple[str, str] = ("x", "y"),
    force_lag_zero: bool = True    
) -> Tuple[np.ndarray[float], np.ndarray[float], np.ndarray[int], float]:
    """
    Compute the empirical variogram from transect data.

    Parameters
    ----------
    transect_df : pd.DataFrame
        A dataframe containing georeferenced coordinates associated with a particular variable (e.g.
        biomass). This DataFrame must have at least two valid columns comprising the overall 2D 
        coordinates (e.g. 'x' and 'y').
    n_lags : int
        The number of lags used for computing the (semi)variogram.
    lag_resolution : float
        The distance interval represented by each lag interval. 
    azimuth_filter : bool
        When True, a 2D array of azimuth angles are generated. This subsequent array represents the 
        relative azimuth angles between spatial points, and can serve as a filter for case where 
        a high degree of directionality is assumed. This accompanies the argument 
        'azimuth_angle_threshold' that defines the threshold azimuth angle.
    azimuth_angle_threshold : float
        This threshold is used for filtering the azimuth angles. 
    variable : str, default = 'biomass_density'
        The variable used for computing the empirical variogram (e.g. 'biomass_density'), which 
        must exist as a column in 'transect_df'.
    coordinates : Tuple[str, str], default = ('x', 'y')
        A tuple containing the 'transect_df' column names defining the coordinates. The order of 
        this input matters where they should be defined as the (horizontal axis, vertical axis).
    force_lag_zero : bool, default = True 
        When True, the nugget effect is assumed to be 0.0 for the empirical variogram. This adds 
        lag 0 to the subsequent array outputs where semivariance (or 'gamma_h') is also equal to 
        0.

    Returns
    ----------
    Tuple[np.ndarray[float], np.ndarray[float], np.ndarray[int], float]
        A tuple containing arrays with the lag intervals, semivariance, and lag counts. The mean 
        lag covariance between the head and tail points computed for all input data is also 
        provided.
        
    """
    # Initialize lag array
    lags = np.concatenate([np.arange(1, n_lags) * lag_resolution])

    # Calculate the distance (and azimuth) matrix
    distance_matrix, azimuth_matrix = lag_distance_matrix(
        coordinates_1=transect_df,
        coordinate_names = coordinates,
        self=True,
        azimuth_matrix=azimuth_filter
    )

    # Convert lag distances to lag numbers
    lag_matrix = np.round(distance_matrix / lag_resolution).astype(int) + 1

    # Extract estimates column
    estimates = transect_df[variable].to_numpy()

    # Create a triangle mask with the diaganol offset to the left by 1
    # ---- Initial mask
    triangle_mask = np.tri(len(estimates), k=-1, dtype=bool)
    # ---- Vertically and then horizontally flip to force the 'True' and 'False' positions
    triangle_mask_flp = np.flip(np.flip(triangle_mask), axis=1)

    # Quantize the lags
    lag_counts, lag_estimates, lag_estimates_squared, lag_deviations = quantize_lags(
        estimates, lag_matrix, triangle_mask_flp, azimuth_matrix, n_lags, azimuth_angle_threshold
    )

    # Compute the mean and standard deviation of the head estimates for each lag bin
    # ---- Apply a mask using the triangle bitmap
    head_mask = np.where(triangle_mask_flp, lag_matrix, -1)

    # Helper function for computing the binned summations for each row
    def bincount_row(row, n_lags):
        return np.bincount(row[row != -1], minlength=n_lags)[1:n_lags]

    # Pre-allocate vectors/arrays that will be iteratively filled
    head_index = np.zeros((len(estimates), n_lags - 1))
    # ---- Find the head indices of each lag for each row
    head_index = np.apply_along_axis(bincount_row, axis=1, arr=head_mask, n_lags=n_lags)

    # Compute the standardized semivariance [gamma(h)]
    gamma_h, lag_covariance = semivariance(
        estimates, lag_estimates, lag_estimates_squared, lag_counts, lag_deviations, head_index
    )

    # Prepend a 0.0 and force the nugget effect to be 0.0, if necessary
    # ---- Return the computed lags, empirical variogram estimate [gamma(h)], and lag counts
    if force_lag_zero:
        return (
            np.concatenate([[0], lags]),
            np.concatenate([[0.], gamma_h]),
            np.concatenate([[len(estimates) - 1], lag_counts]),
            lag_covariance
        )
    else:
        return lags, gamma_h, lag_counts, lag_covariance

gamma_h = np.concatenate([[0], gamma_h])
lag_counts = np.concatenate([[len(estimates) - 1], lag_counts])

####
from echopop.spatial.variogram import initialize_initial_optimization_values, get_variogram_arguments
from lmfit import Minimizer, Parameters
initialization_variogram = ['nugget', 'sill', 'correlation_range', 'hole_effect_range', 'decay_power']
dict_optimization = {'max_nfev': 500, 'ftol': 1e-06, 'gtol': 0.0001, 'xtol': 1e-06, 'diff_step': 1e-08, 'tr_solver': 'exact', 'x_scale': 'jac', 'jac': '3-point'}
parameters = initialize_initial_optimization_values(
    initialization_variogram, 
    {"model": ["exponential", "bessel"],
     "n_lags": 30, "nugget": 0., "sill": 0.91, "hole_effect_range": 0.0, 
     "correlation_range": 0.007, "decay_power": 1.5},
)
model = ["exponential", "bessel"]
range = 0.06
####

# Compute the lag weights
lag_weights = lag_counts / lag_counts.sum()

# Vertically stack the lags, semivariance, and weights
data_stack = np.vstack((lags, gamma_h, lag_weights))

# Index lag distances that are within the parameterized range
# within_range = np.where(lags <= variogram_parameters["range"])[0]
within_range = np.where(lags <= range)[0]

# Truncate the data stack
truncated_stack = data_stack[:, within_range]

# Get model name
# _, variogram_fun = get_variogram_arguments(variogram_parameters["model"])
_, variogram_fun = get_variogram_arguments(model)

# Create helper cost-function that is weighted using the kriging weights (`w`), lag
# distances (`x`), and empirical semivariance (`y`)
def cost_function(parameters, x, y, w):
    yr = variogram_fun["model_function"](x, **parameters)
    return (yr - y) * w

# Compute the initial fit based on the pre-optimized parameter values
initial_fit = cost_function(
    parameters,
    x=truncated_stack[0],
    y=truncated_stack[1],
    w=truncated_stack[2],
)

# Compute the initial mean absolute deviation (MAD)
mad_initial = np.mean(np.abs(initial_fit))

# Generate `Minimizer` function class required for bounded optimization
minimizer = Minimizer(
    cost_function,
    parameters,
    fcn_args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]),
)

# Minimize the cost-function to compute the best-fit/optimized variogram parameters
parameters_optimized = minimizer.minimize(
    method="least_squares", **dict_optimization
)

# Calculate the optimized MAD
mad_optimized = np.mean(np.abs(parameters_optimized.residual))

# Extract the best-fit parameter values
best_fit_params = parameters_optimized.params.valuesdict()

return (
    best_fit_params,
    (
        list(optimization_settings["parameters"].keys()),
        list(optimization_settings["parameters"].valuesdict().values()),
        mad_initial,
    ),
    (list(best_fit_params.keys()), list(best_fit_params.values()), mad_optimized),
)

####################################################################################################
from echopop.survey import Survey
from echopop.spatial.transect import correct_transect_intervals
from echopop.acoustics import aggregate_sigma_bs, nasc_to_biomass
# from echopop.biology import (
#     # age1_metric_proportions,
#     # distribute_length_age,
#     # filter_species,
#     # fit_length_weight_relationship,
#     # fit_length_weights,
#     # impute_kriged_values,
#     # # number_proportions,
#     # # partition_transect_age,
#     # quantize_number_counts,
#     # quantize_weights,
#     # reallocate_kriged_age1,
#     # weight_proportions,
# )
from echopop.spatial.krige import kriging
from echopop.spatial.mesh import crop_mesh, mesh_to_transects, stratify_mesh
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import (
    edit_transect_columns,
    save_transect_coordinates,
    summarize_transect_strata,
    transect_spatial_features,
)
from echopop.spatial.variogram import (
    empirical_variogram,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
)
from echopop.statistics import stratified_transect_statistic
from echopop.utils.validate_dict import (
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
)

survey = Survey(init_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis()
survey.fit_variogram()
survey.analysis["settings"]["variogram"][""]
self = survey
input_dict, analysis_dict, configuration_dict, settings_dict = self.input, self.analysis["transect"], self.config, self.analysis["settings"]
# Extract the necessary correct strata mean sigma_bs
sigma_bs_strata = analysis_dict["acoustics"]["sigma_bs"]["strata_mean_df"]

# Pull out the length-weight conversion for each stratum
length_weight_strata = analysis_dict["biology"]["weight"]["weight_stratum_df"]

# Get the name of the stratum column
stratum_col = settings_dict["transect"]["stratum_name"]

# Get group-specific columns
age_group_cols = settings_dict["transect"]["age_group_columns"]

# Extract the correct strata dataframe
# ---- Define `strata_df` if KS
if settings_dict["transect"]["stratum"] == "ks":
    strata_df = input_dict["spatial"]["strata_df"].copy()
# Define `inpfc_strata_df` if INPFC
elif settings_dict["transect"]["stratum"] == "inpfc":
    strata_df = input_dict["spatial"]["inpfc_strata_df"].copy()

# Get group-specific column names and create conversion key
name_conversion_key = {age_group_cols["haul_id"]: "haul_num", age_group_cols["nasc_id"]: "nasc"}
# ---- Update if the stratum is not equal to INPFC
if settings_dict["transect"]["stratum"] != "inpfc":
    name_conversion_key.update({age_group_cols["stratum_id"]: stratum_col})

# Rename columns
# ---- Extract NASC data
nasc_data = input_dict["acoustics"]["nasc_df"].copy()
# ---- Change names
nasc_data.rename(columns=name_conversion_key, inplace=True)

# Correct the acoustic survey transect intervals
nasc_interval_df = correct_transect_intervals(nasc_data)

distributions_dict, proportions_dict, TS_L_parameters, settings_dict = (
    input_dict["biology"]["distributions"],
    analysis_dict["biology"]["proportions"],
    configuration_dict["TS_length_regression_parameters"]["pacific_hake"],
    settings_dict
)