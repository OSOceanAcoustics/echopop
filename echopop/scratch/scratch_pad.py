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


def load_isobath_data(
    isobath_filepath: Union[str, Path], 
    sheet_name: str, 
    column_name_map: Dict[str, str] = {}
) -> pd.DataFrame:
    """
    Load isobath data from an Excel file.

    Parameters
    ----------
    isobath_filepath : str or Path
        Path to the Excel file containing mesh data
    sheet_name : str
        Name of the sheet containing the mesh data
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names

    Returns
    -------
    pd.DataFrame
        DataFrame containing the isobath data with longitude and latitude

    Examples
    --------
    >>> isobath_df = load_isobath_data("isobath_file.xlsx", "Sheet1")
    >>> isobath_df = load_isobath_data("isobath_file.xlsx", column_name_map={"latitude_200m":
    "latitude"})
    """
    mesh_filepath = Path(isobath_filepath)

    if not mesh_filepath.exists():
        raise FileNotFoundError(f"Mesh data file not found: {isobath_filepath}")

    # Read Excel file into memory
    df = pd.read_excel(isobath_filepath, sheet_name=sheet_name, index_col=None, header=0)

    # Force the column names to be lower case
    df.columns = df.columns.str.lower()

    # Rename columns if mapping is provided
    if column_name_map:
        df.rename(columns=column_name_map, inplace=True)

    return df

def standardize_coordinates(
    data_df: pd.DataFrame,
    longitude_offset: float = 0.,
    latitude_offset: float = 0.,
    reference_df: Optional[pd.DataFrame] = None,
    delta_longitude: Optional[float] = None,
    delta_latitude: Optional[float] = None,
) -> Tuple[pd.DataFrame, Union[float, None], Union[float, None]]:
    """
    Standardizes the longitude and latitude of a dataset
    
    Parameters
    ----------
    data_df: pd.DataFrame
        DataFrame with longitude and latitude coordinates
    longitude_offset: float, default=0.
        Offset to apply to the longitude coordinates
    latitude_offset: float, default=0. 
        Offset to apply to the latitude coordinates
    reference_df: pd.DataFrame, optional
        Reference DataFrame with longitude and latitude coordinates for interpolation that is 
        used as an additional offset to longitude
    delta_longitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates
    delta_latitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates
        
    Returns
    -------
    pd.DataFrame
        DataFrame with the new standardized coordinates 'x' and 'y'.
    float or None
        Total longitudinal distance (degrees) used for standardizing coordinates that can be used 
        for the transformation of other georeferenced datasets.
    float or None
        Total latitudinal distance (degrees) used for standardizing coordinates that can be used 
        for the transformation of other georeferenced datasets.
    """

    # Create interpolation function from reference grid coordinates (to interpolate longitude)
    if reference_df is not None:
        reference_interp = interpolate.interp1d(
            reference_df["latitude"], reference_df["longitude"], kind="linear", bounds_error=False
        )
    else:
        reference_interp = 0.

    # Transform longitude
    transformed_longitude = (
        data_df["longitude"] - reference_interp(data_df["latitude"]) + longitude_offset
    )

    # Calculate the geospatial distances along the longitudinal and latitudinal axes [if missing]
    # ---- Longitude
    if delta_longitude is None:
        delta_longitude = transformed_longitude.max() - transformed_longitude.min()
    # ---- Latitude 
    if delta_latitude is None:
        delta_latitude = data_df.latitude.max() - data_df.latitude.min()

    # Standardize the x- and y-coordinates
    # ---- longitude --> x
    data_df["x"] = (
        np.cos(np.pi / 180.0 * data_df["latitude"])
        * (
            transformed_longitude
            - longitude_offset
        )
        / delta_longitude
    )
    # ---- latitude --> y
    data_df["y"] = (
        data_df["latitude"] - latitude_offset
    ) / delta_latitude

    # Return the output tuple
    return (data_df, delta_longitude, delta_latitude)

def lag_distance_matrix(
    coordinates_1: Union[pd.DataFrame, np.ndarray],
    coordinates_2: Optional[Union[pd.DataFrame, np.ndarray]] = None,
    self: bool = False,
    azimuth_matrix: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Calculate the lag distance matrix between two sets of coordinates.
    Parameters
    ----------
    coordinates_1 : pd.DataFrame or np.ndarray
        First set of coordinates (x, y) as a DataFrame or numpy array.
    coordinates_2 : pd.DataFrame or np.ndarray, optional
        Second set of coordinates (x, y) as a DataFrame or numpy array. If None, uses coordinates_1.
    self : bool, default=False
        If True, calculates distances within the same set of coordinates.
    azimuth_matrix : bool, default=False
        If True, returns both distance and azimuth angles; if False, returns only distances.
    
    Returns
    -------
    Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]
        If angles is True, returns a tuple of distance matrix and azimuth angles.
        If angles is False, returns only the distance matrix.
    """
    
    # Set reference to self if 'coordinates_2' is not defined
    if coordinates_2 is None:
        coordinates_2 = coordinates_1

    # Get the distance array coordinates in the 'x' and 'y' directions
    # ---- Case: Internal distances (as arrays) and self is True
    if self and all(isinstance(x, np.ndarray) for x in [coordinates_1, coordinates_2]):
        # ---- x-coordinates
        x_coords = (coordinates_1, coordinates_1)
        # ---- y-coordinates
        y_coords = (coordinates_2, coordinates_2)  
    # ---- Case: Distances (as arrays) and self is False
    elif not self and all(isinstance(x, np.ndarray) for x in [coordinates_1, coordinates_2]):
        # ---- x-coordinates
        x_coords = (coordinates_1, coordinates_2)
        # ---- y-coordinates
        y_coords = (coordinates_1, coordinates_2)     
    # ---- Case: DataFrames
    elif all(isinstance(x, pd.DataFrame) for x in [coordinates_1, coordinates_2]):
        # ---- x-coordinates
        x_coords = (coordinates_1["x"].to_numpy(), coordinates_2["x"].to_numpy())
        # ---- y-coordinates
        y_coords = (coordinates_1["y"].to_numpy(), coordinates_2["y"].to_numpy())  
    #  Resolve the distances
    # ---- x-distance
    x_distance = np.subtract.outer(*x_coords)
    # ---- y-distance
    y_distance = np.subtract.outer(*y_coords)
        
    # Get the azimuth angle matrix, if required
    if azimuth_matrix:
        # ---- Create copies of 'x_distance' and 'y_distance'
        x_angles = x_distance.copy(); y_angles = y_distance.copy()
        # ---- Replace the self-points with 'NaN'
        np.fill_diagonal(x_angles, np.nan); np.fill_diagonal(y_angles, np.nan)
        # ---- Calculate the azimuth angle grid
        azimuth_grid = (
            np.arctan(
                np.divide(y_angles, x_angles, where=(x_angles != 0.) & (~np.isnan(x_angles)))
            ) * 180. / np.pi + 180. % 180.
        )
        # ---- Return the resulting tuple of matrices for Euclidean distances and azimuth angles
        return np.sqrt(x_distance * x_distance + y_distance * y_distance), azimuth_grid
    else:
        # ---- Return Euclidean distance matrix
        return np.sqrt(x_distance * x_distance + y_distance * y_distance)

azimuth_range = 360.
n_lags = 30
lag_resolution = 0.002
sill = 0.91
hole_effect_range = 0.
correlation_range = 0.007
decay_power = 1.5
transect_df = df_nasc_all_ages.copy()
estimates = transect_df["biomass"].to_numpy()
######
# estimates, lag_matrix, triangle_mask_flp, azimuth_matrix, **variogram_parameters
######
lags = np.concatenate([[0], np.arange(1, n_lags) * lag_resolution])

# QUANTIZE LAGS
# Convert the lag distance matrix into number lags (based on lag resolution)
lag_matrix = np.round(lag_distance_matrix / lag_resolution).astype(int) + 1



# Create a triangle mask with the diaganol offset to the left by 1
# ---- Initial mask
triangle_mask = np.tri(len(estimates), k=-1, dtype=bool)
# ---- Vertically and then horizontally flip to force the 'True' and 'False' positions
triangle_mask_flp = np.flip(np.flip(triangle_mask), axis=1)

# Validate that `estimates` is a 1D array
if estimates.ndim > 1:
    raise ValueError("Estimates array ('estimates') must be a 1D array.")

# Validate that dimensions are 2D
if lag_matrix.ndim < 2 or triangle_mask_flp.ndim < 2 or azimuth_matrix.ndim < 2:
    error = (
        "The function `quantize_lags` requires arrays to be 2D. The following 1D arrays have "
        "invalid shapes: "
    )
    invalid_arrays = []
    if lag_matrix.ndim < 2:
        invalid_arrays += ["'lag_matrix'"]
    if triangle_mask_flp.ndim < 2:
        invalid_arrays += ["'mask_matrix'"]
    if azimuth_matrix.ndim < 2:
        invalid_arrays += ["'azimuth_matrix'"]
    raise ValueError(error + ", ".join(invalid_arrays))

#####
# lag_matrix, mask_matrix, azimuth_matrix, azimuth_range
#####
# variogram_matrix_filter <<<<<
data_matrix = lag_matrix.copy()
mask_matrix = triangle_mask_flp.copy()
azimuth_angle_threshold = 180.

def filter_lag_matrix(
    data_matrix: np.ndarray[int],
    mask_matrix: np.ndarray[bool],
    azimuth_matrix: Optional[np.ndarray[float]] = None,
    azimuth_angle_threshold: Optional[float] = None,
) -> np.ndarray[int]:
    """
    Filter the lag matrix based on a boolean mask and optional azimuth angle matrix
    
    Parameters
    ----------
    data_matrix : np.ndarray[int]
        The lag distance matrix to be filtered.
    mask_matrix : np.ndarray[bool]
        A boolean mask indicating which elements to keep in the lag matrix. This is typically a 
        triangle boolean matrix.
    azimuth_matrix : np.ndarray[float], optional
        An optional azimuth angle matrix that can be used to filter the lag matrix based on 
        azimuth angles.
    azimuth_angle_threshold : float, optional
        If provided, this threshold is used to filter the azimuth angles in the azimuth matrix. 
        This defines the The total azimuth angle range that is allowed for constraining the 
        relative angles between spatial points, particularly for cases where a high degree of 
        directionality is assumed.
    
    Returns
    -------
    np.ndarray[int]
        The filtered lag matrix, where elements not meeting the mask or azimuth criteria are set to 
        zero.
    """
    
    # Convert array to matrix, if needed
    if data_matrix.ndim == 1:
        data_matrix = np.tile(data_matrix, (len(data_matrix), 1))
    else:
        if data_matrix.shape != mask_matrix.shape:
            # ---- Determine which dimension is mismatched
            dimension_diff = np.where(np.array(data_matrix.shape) != np.array(mask_matrix.shape))[0]
            if dimension_diff == 0:
                data_matrix = np.tile(data_matrix, (len(data_matrix), 1))
            else:
                data_matrix = np.tile(data_matrix, (1, len(data_matrix)))
                
    # If 'azimuth_matrix' is supplied, then apply threshold as additional bitmap
    if azimuth_matrix is not None and azimuth_angle_threshold is not None:
        # ---- Replace any 'NaN' values with 0's
        azimuth_matrix[np.isnan(azimuth_matrix)] = 0.0
        # ---- Create bitmap
        azimuth_mask = (
            (azimuth_matrix >= -azimuth_angle_threshold) & 
            (azimuth_matrix < azimuth_angle_threshold)
        )
    else:
        # ---- Create empty azimuth mask
        azimuth_mask = np.ones_like(data_matrix, dtype=bool)
        
    # Mask the data matrix and broadcast into a 1D array
    return data_matrix[mask_matrix & azimuth_mask]

# Filter the lag matrix based on the mask and azimuth angle threshold
equivalent_lags = filter_lag_matrix(lag_matrix, triangle_mask_flp, azimuth_matrix, azimuth_angle_threshold)
# ---- Compute the binned sums
lag_counts = np.bincount(equivalent_lags)[1:n_lags]

# Sum the estimates for each lag
estimates_lagged = filter_lag_matrix(estimates, triangle_mask_flp, azimuth_matrix, azimuth_angle_threshold)
# ---- Compute the binned sums
lag_estimates = np.bincount(equivalent_lags, weights=estimates_lagged)[1:n_lags]

# Compute the binned squared-sum
lag_estimates_squared = np.bincount(equivalent_lags, weights=estimates_lagged**2)[1:n_lags]

# Generate a lag bitmap
lag_bitmap = equivalent_lags < n_lags

# Filter the estimates matrix
estimates_matrix = filter_lag_matrix(np.arange(len(estimates))[:, np.newaxis],
                                     triangle_mask_flp,
                                     azimuth_matrix, 
                                     azimuth_angle_threshold)

# Calculate the deviations between the indexed and lag-specific estimates
deviations = (estimates[estimates_matrix][lag_bitmap] - estimates_lagged[lag_bitmap]) ** 2

# Sum the deviations per lag bin
lag_deviations = np.bincount(equivalent_lags[lag_bitmap], weights=deviations)[1:n_lags]

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

# Calculate the mean head estimate per lag bin
mean_head = (estimates[:, np.newaxis] * (head_index / lag_counts)).sum(axis=0)

# Calculate the standard deviation of head values per lag
sigma_head = np.sqrt(
    ((estimates[:, np.newaxis] - mean_head) ** 2 * (head_index / lag_counts)).sum(axis=0)
)

# Calculate the global mean and variance for each lag bin
# ---- Mean
lag_means = lag_estimates / lag_counts
# ---- Variance
lag_variance = lag_estimates_squared / lag_counts - lag_means**2

# Estimate the standard deviation of tail estimates
sigma_tail = np.sqrt(np.abs(lag_variance))

# Calculate the semivariance
# ---- Compute the partial sill that is applied as a weighted calculation
partial_sill = sigma_tail * sigma_head
# ---- Semivariance [gamma(h)] and cross-sill estimate
# gamma_h = 0.5 * lag_deviations / (lag_counts * partial_sill)
with np.errstate(divide="ignore"):
    gamma_h = 0.5 * lag_deviations / (lag_counts * partial_sill)
# Calculate the mean lag distance covariance
# ---- Find non-zero head and tail variances
non_zero_variance = np.where((sigma_head > 0.0) & (sigma_tail > 0.0))[0]
# ---- Mean lag covariance
if non_zero_variance.size > 0:
    mean_lag_covariance = (sigma_head[non_zero_variance] * sigma_tail[non_zero_variance]).mean()
else:
    mean_lag_covariance = np.nan

lag_covariance = mean_lag_covariance

# Compute the standardized semivariance [gamma(h)]
# gamma_h, lag_covariance = semivariance(
#     estimates, lag_estimates, lag_estimates_squared, lag_counts, lag_deviations, head_index
# )

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