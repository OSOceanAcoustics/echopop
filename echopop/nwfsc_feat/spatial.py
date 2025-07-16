import warnings
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import geopandas as gpd
import numpy as np
import pandas as pd
from lmfit import Minimizer, Parameters
from scipy import interpolate
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from ..spatial.variogram import variogram
from .projection import wgs84_to_utm

# Set warnings filter
warnings.simplefilter("always")


def standardize_coordinates(
    data_df: pd.DataFrame,
    x_offset: float = 0.0,
    y_offset: float = 0.0,
    coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
    reference_df: Optional[pd.DataFrame] = None,
    delta_x: Optional[float] = None,
    delta_y: Optional[float] = None,
) -> Tuple[pd.DataFrame, Union[float, None], Union[float, None]]:
    """
    Standardize the longitude and latitude coordinates of a dataset.

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame with coordinates
    x_offset : float, default=0.
        Offset to apply to the x-coordinates that corresponds to `coordinate_names[0]`
    y_offset : float, default=0.
        Offset to apply to the y-coordinates that corresponds to `coordinate_names[0]`
    coordinate_names : Tuple[str, str], default=("longitude", "latitude")
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    reference_df : pd.DataFrame, optional
        Reference DataFrame with x and y coordinates for interpolation that is
        used as an additional offset to the x-axis.
    delta_x : float, optional
        Total x-axis distance used for standardizing coordinates
    delta_y : float, optional
        Total y-axis distance used for standardizing coordinates

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

    # Get the coordinate names
    x_coord, y_coord = coordinate_names

    # Create interpolation function from reference grid coordinates (to interpolate longitude)
    if reference_df is not None:
        reference_interp = interpolate.interp1d(
            reference_df[y_coord], reference_df[x_coord], kind="linear", bounds_error=False
        )
        reference_offset = reference_interp(data_df[y_coord])
    else:
        reference_offset = 0.0

    # Transform longitude
    transformed_x = data_df[x_coord] - reference_offset + x_offset

    # Calculate the geospatial distances along the x- and y-axes [if missing]
    # ---- Longitude
    if delta_x is None:
        delta_x = transformed_x.max() - transformed_x.min()
    # ---- Latitude
    if delta_y is None:
        delta_y = data_df[y_coord].max() - data_df[y_coord].min()

    # Standardize the x- and y-coordinates
    # ---- x
    data_df["x"] = np.cos(np.pi / 180.0 * data_df[y_coord]) * (transformed_x - x_offset) / delta_x
    # ---- y
    data_df["y"] = (data_df[y_coord] - y_offset) / delta_y

    # Return the output tuple
    return (data_df, delta_x, delta_y)


def lag_distance_matrix(
    coordinates_1: Union[pd.DataFrame, np.ndarray],
    coordinate_names: Optional[Tuple[str, str]] = None,
    coordinates_2: Optional[Union[pd.DataFrame, np.ndarray]] = None,
    self: bool = False,
    azimuth_matrix: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the lag distance matrix between two sets of coordinates.

    Parameters
    ----------
    coordinates_1 : pd.DataFrame or np.ndarray
        First set of coordinates (x, y) as a DataFrame or numpy array.
    coordinate_names : Tuple[str, str], optional
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    coordinates_2 : pd.DataFrame or np.ndarray, optional
        Second set of coordinates (x, y) as a DataFrame or numpy array. If None, uses coordinates_1.
    self : bool, default=False
        If True, calculates distances within the same set of coordinates.
    azimuth_matrix : bool, default=False
        If True, returns both distance and azimuth angles; if False, returns only distances.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        If azimuth_matrix is True, returns a tuple of (distance matrix, azimuth angles).
        If azimuth_matrix is False, returns (distance matrix, empty array).
    """

    # Get coordinate names
    if coordinate_names is not None:
        x_name, y_name = coordinate_names

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
        x_coords = (coordinates_1[x_name].to_numpy(), coordinates_2[x_name].to_numpy())
        # ---- y-coordinates
        y_coords = (coordinates_1[y_name].to_numpy(), coordinates_2[y_name].to_numpy())
    #  Resolve the distances
    # ---- x-distance
    x_distance = np.subtract.outer(*x_coords)
    # ---- y-distance
    y_distance = np.subtract.outer(*y_coords)

    # Get the azimuth angle matrix, if required
    if azimuth_matrix:
        # ---- Create copies of 'x_distance' and 'y_distance'
        x_angles = x_distance.copy()
        y_angles = y_distance.copy()
        # ---- Replace the self-points with 'NaN'
        np.fill_diagonal(x_angles, np.nan)
        np.fill_diagonal(y_angles, np.nan)
        # ---- Calculate the azimuth angle grid
        azimuth_grid = (
            np.arctan(
                np.divide(y_angles, x_angles, where=(x_angles != 0.0) & (~np.isnan(x_angles)))
            )
            * 180.0
            / np.pi
            + 180.0 % 180.0
        )
        # ---- Return the resulting tuple of matrices for Euclidean distances and azimuth angles
        return np.sqrt(x_distance * x_distance + y_distance * y_distance), azimuth_grid
    else:
        # ---- Return Euclidean distance matrix
        return np.sqrt(x_distance * x_distance + y_distance * y_distance), np.array([])


def filter_lag_matrix(
    data_matrix: np.ndarray[int],
    mask_matrix: np.ndarray[bool],
    azimuth_matrix: Optional[np.ndarray[float]] = None,
    azimuth_angle_threshold: Optional[float] = None,
) -> np.ndarray[int]:
    """
    Filter the lag matrix based on a boolean mask and optional azimuth angle matrix.

    Parameters
    ----------
    data_matrix : np.ndarray
        The lag distance matrix to be filtered.
    mask_matrix : np.ndarray[bool]
        A boolean mask indicating which elements to keep in the lag matrix. This is typically a
        triangle boolean matrix.
    azimuth_matrix : np.ndarray[float], optional
        An optional azimuth angle matrix that can be used to filter the lag matrix based on
        azimuth angles.
    azimuth_angle_threshold : float, optional
        If provided, this threshold is used to filter the azimuth angles in the azimuth matrix.
        This defines the total azimuth angle range that is allowed for constraining the
        relative angles between spatial points, particularly for cases where a high degree of
        directionality is assumed.

    Returns
    -------
    np.ndarray
        The filtered lag matrix, where elements not meeting the mask or azimuth criteria are
        extracted as a 1D array.
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
    if (
        azimuth_matrix is not None
        and len(azimuth_matrix) > 0
        and azimuth_angle_threshold is not None
    ):
        # ---- Replace any 'NaN' values with 0's
        azimuth_matrix[np.isnan(azimuth_matrix)] = 0.0
        # ---- Create bitmap
        azimuth_mask = (azimuth_matrix >= -azimuth_angle_threshold) & (
            azimuth_matrix < azimuth_angle_threshold
        )
    else:
        # ---- Create empty azimuth mask
        azimuth_mask = np.ones_like(data_matrix, dtype=bool)

    # Mask the data matrix and broadcast into a 1D array
    return data_matrix[mask_matrix & azimuth_mask]


def quantize_lags(
    estimates: np.ndarray[float],
    lag_matrix: np.ndarray[int],
    mask_matrix: np.ndarray[bool],
    azimuth_matrix: np.ndarray[float],
    n_lags: int,
    azimuth_angle_threshold: Optional[float] = None,
) -> Tuple[np.ndarray[int], np.ndarray[float], np.ndarray[float], np.ndarray[float]]:
    """
    Quantize lag counts, sums, square-sums, and deviations required for computing the empirical
    variogram

    Parameters
    ----------
    estimates : np.ndarray[float]
        A 1D array of field estimates, e.g. 'biomass'.
    lag_matrix : np.ndarray[int]
        A 2D array of integer lag values.
    mask_matrix : np.ndarray[bool]
        A boolean mask indicating which elements to keep in the lag matrix. This is typically a
        triangle boolean matrix.
    azimuth_matrix : np.ndarray[float]
        A 2D array of azimuth angle values. Alternatively, this may also be a 1D array of length 0
        for the case where the azimuth angle matrix is not defined.
    n_lags : int
        Number of lags used for the variogram analysis.
    azimuth_angle_threshold : float, optional
        If provided, this threshold is used to filter the azimuth angles in the azimuth matrix.
        This defines the total azimuth angle range that is allowed for constraining the
        relative angles between spatial points, particularly for cases where a high degree of
        directionality is assumed.

    Returns
    -------
    Tuple[np.ndarray[int], np.ndarray[float], np.ndarray[float], np.ndarray[float]]
        A tuple comprising the binned lag counts, summed estimates, square-summed estimates, and
        deviations.
    """

    # Validate that `estimates` is a 1D array
    if estimates.ndim > 1:
        raise ValueError("Estimates array ('estimates') must be a 1D array.")

    # Validate that dimensions are 2D
    if (
        lag_matrix.ndim < 2
        or mask_matrix.ndim < 2
        or (azimuth_matrix.ndim < 2 and len(azimuth_matrix) > 0)
    ):
        error = (
            "The function `quantize_lags` requires arrays to be 2D. The following 1D arrays have "
            "invalid shapes: "
        )
        invalid_arrays = []
        if lag_matrix.ndim < 2:
            invalid_arrays += ["'lag_matrix'"]
        if mask_matrix.ndim < 2:
            invalid_arrays += ["'mask_matrix'"]
        if azimuth_matrix.ndim < 2 and len(azimuth_matrix) > 0:
            invalid_arrays += ["'azimuth_matrix'"]
        raise ValueError(error + ", ".join(invalid_arrays))

    # Filter the lag matrix based on the mask and azimuth angle threshold
    equivalent_lags = filter_lag_matrix(
        lag_matrix, mask_matrix, azimuth_matrix, azimuth_angle_threshold
    )
    # ---- Compute the binned sums
    lag_counts = np.bincount(equivalent_lags)[1:n_lags]

    # Sum the estimates for each lag
    estimates_lagged = filter_lag_matrix(
        estimates, mask_matrix, azimuth_matrix, azimuth_angle_threshold
    )
    # ---- Compute the binned sums
    lag_estimates = np.bincount(equivalent_lags, weights=estimates_lagged)[1:n_lags]

    # Compute the binned squared-sum
    lag_estimates_squared = np.bincount(equivalent_lags, weights=estimates_lagged**2)[1:n_lags]

    # Generate a lag bitmap
    lag_bitmap = equivalent_lags < n_lags

    # Filter the estimates matrix
    estimates_matrix = filter_lag_matrix(
        np.arange(len(estimates))[:, np.newaxis],
        mask_matrix,
        azimuth_matrix,
        azimuth_angle_threshold,
    )

    # Calculate the deviations between the indexed and lag-specific estimates
    deviations = (estimates[estimates_matrix][lag_bitmap] - estimates_lagged[lag_bitmap]) ** 2

    # Sum the deviations per lag bin
    lag_deviations = np.bincount(equivalent_lags[lag_bitmap], weights=deviations)[1:n_lags]

    # Return the calculated quantities
    return lag_counts, lag_estimates, lag_estimates_squared, lag_deviations


def semivariance(
    estimates: np.ndarray[float],
    lag_estimates: np.ndarray[float],
    lag_estimates_squared: np.ndarray[float],
    lag_counts: np.ndarray[int],
    lag_deviations: np.ndarray[float],
    head_index: np.ndarray[int],
) -> Tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Compute the (standardized) semivariance

    Parameters
    ----------
    estimates : np.ndarray[float]
        A 1D array of field estimates, e.g. 'biomass'.
    lag_estimates : np.ndarray[float]
        Summed lag estimates weighted by values for the values within 'estimates'.
    lag_estimates_squared : np.ndarray[float]
        Square-summed lag estimates weighted by values for the values within 'estimates'.
    lag_counts : np.ndarray[int]
        Binned counts of values within each lag.
    lag_deviations : np.ndarray[float]
        Statistical deviation within each lag bin.
    head_index : np.ndarray[int]
        A 2D array containing the head indices of each lag for each row.

    Returns
    -------
    Tuple[np.ndarray[float], np.ndarray[float]]
        A tuple comprising 1D arrays with the semivariance estimates and mean lag covariance at
        each lag.
    """

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

    # Return the outputs
    return gamma_h, mean_lag_covariance


def empirical_variogram(
    transect_df: pd.DataFrame,
    n_lags: int,
    lag_resolution: float,
    azimuth_filter: bool,
    azimuth_angle_threshold: float,
    variable: str = "biomass_density",
    coordinates: Tuple[str, str] = ("x", "y"),
    force_lag_zero: bool = True,
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
    -------
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
        coordinate_names=coordinates,
        self=True,
        azimuth_matrix=azimuth_filter,
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
            np.concatenate([[0.0], gamma_h]),
            np.concatenate([[len(estimates) - 1], lag_counts]),
            lag_covariance,
        )
    else:
        return lags, gamma_h, lag_counts, lag_covariance


def fit_variogram(
    lags: np.ndarray[float],
    lag_counts: np.ndarray[int],
    gamma: np.ndarray[float],
    variogram_parameters: Parameters,
    model: Union[str, List[str]] = ["bessel", "exponential"],
    optimizer_kwargs: Dict[str, Any] = {},
) -> Tuple[Dict[str, Any], float, float]:
    """
    Compute the best-fit variogram parameters for input data

    Parameters
    ----------
    lags : np.ndarray[float]
        A 1D array of the lag distances.
    lag_counts : np.ndarray[int]
        A 1D array of the lag counts.
    gamma : np.ndarray[float]
        A 1D array comprising the semivariance estimates for each lag.
    variogram_parameters : Parameters
        A :class:`lmfit.parameter.Parameters` object containing the parameters required for the
        defined variogram model. See :func:`echopop.spatial.variogram.variogram` for more details.
        Valid parameters include:

            - **sill : float**
            The asymptotic value as lags approach infinity

            - **nugget : float**
            The semivariogram *y*-intercept that corresponds to variability at lag distances
            shorter than the lag resolution

            - **correlation_range : float**
            The relative length scale, or range, at which the autocorrelation between lag distances
            no longer increases and becomes asymptotic

            - **hole_effect_range : float**
            The (normalized) length scale/range that holes' are observed, which represent 'null'
            (or very small) points compared to their neighboring lags

            - **decay_power : float**
            An exponential term that is used in certain generalized exponential (or related)
            semivariogram models that modulates the ascending rate for a semivariogram

            - **enhance_semivariance : bool**
            A boolean term that determines whether the correlation decay in certain  cosine-related
            variogram models are enhanced (or not) with increasing lag distances

    model : Union[str, List[str]], default=['bessel', 'exponential']
        A string or list of model names. A single name represents a single family model. Two inputs
        represent the desired composite model (e.g. the composite J-Bessel and exponential model).
        Defaults to: ``model=["bessel", "exponential"]``. Available models and their required
        arguments can be reviewed in the :func:`echopop.spatial.variogram.variogram` function.
    optimizer_kwargs : Dict[str, Any], default={}
        A dictionary comprising the various function arguments used by
        :class:`lmfit.minimizer.Minimizer` for Least-Squares minimization that incorporates the
        Trust Reflective method. See :class:`lmfit.minimizer.Minimizer` for more details.

    Returns
    -------
    Tuple[Dict[str, Any], float, float]
        A tuple containing a dictionary with the best-fit keyword variogram parameter values,
        the mean absolute deviation (MAD) of the initial parameter values, and the MAD of the
        best-fit parameter values.

    See Also
    --------
    :func:`echopop.spatial.variogram.variogram` :
        Variogram model parameters.
    :class:`lmfit.parameter.Parameters` :
        Variogram parameter optimization leverages the ``Parameters`` class from ``lmfit`` for
        model optimization.
    :class:`lmfit.minimizer.Minimizer` :
        Optimization keyword arguments used for the Least-Squares fitting.
    """
    # Normalize the lag counts to get the lag weights
    lag_weights = lag_counts / lag_counts.sum()

    # Vertically stack the lags, semivariance, and weights
    data_stack = np.vstack((lags, gamma, lag_weights))

    # Recover the lag resolution
    delta_lag = np.diff(lags).mean()

    # Compute the range
    range = lags.max() + delta_lag

    # Index lag distances that are within the parameterized range
    within_range = np.where(lags <= range)[0]

    # Truncate the data stack
    truncated_stack = data_stack[:, within_range]

    # Create helper cost-function that is weighted using the kriging weights (`w`), lag
    # distances (`x`), and empirical semivariance (`y`)
    def cost_function(parameters, x, y, w, model):
        yr = variogram(x, {**parameters, **{"model": model}})
        return (yr - y) * w

    # Compute the initial fit based on the pre-optimized parameter values
    initial_fit = cost_function(
        variogram_parameters,
        x=truncated_stack[0],
        y=truncated_stack[1],
        w=truncated_stack[2],
        model=model,
    )

    # Compute the initial mean absolute deviation (MAD)
    mad_initial = np.mean(np.abs(initial_fit))

    # Generate `Minimizer` function class required for bounded optimization
    minimizer = Minimizer(
        cost_function,
        variogram_parameters,
        fcn_args=(truncated_stack[0], truncated_stack[1], truncated_stack[2], model),
    )

    # Minimize the cost-function to compute the best-fit/optimized variogram parameters
    parameters_optimized = minimizer.minimize(method="least_squares", **optimizer_kwargs)

    # Calculate the optimized MAD
    mad_optimized = np.mean(np.abs(parameters_optimized.residual))

    # Extract the best-fit parameter values
    best_fit_params = parameters_optimized.params.valuesdict()

    # Return the final tuple
    return best_fit_params, mad_initial, mad_optimized


def get_survey_western_extents(
    transect_df: pd.DataFrame,
    coordinate_names: Tuple[str, str],
    latitude_threshold: float,
) -> pd.DataFrame:
    """
    Get the western extents of each survey transect that can be used to constrain the adaptive
    nearest neighbors search algorithm incorporated into the kriging interpolation algorithm

    Parameters
    ----------
    transect_df : pd.DataFrame
        A dataframe containing georeferenced coordinates associated with a particular variable (e.g.
        biomass). This DataFrame must have at least two valid columns comprising the overall 2D
        coordinates (e.g. 'x' and 'y'). Furthermore, this function requires that `transect_df`
        also contain a column called 'latitude'.
    coordinates_names : Tuple[str, str], default = ('x', 'y')
        A tuple containing the 'transect_df' column names defining the coordinates. The order of
        this input matters where they should be defined as the (horizontal axis, vertical axis).
    latitude_threshold : float
        A threshold that is applied to the georeferenced coordinates that further constrains any
        extrapolation that occurs during the kriging analysis.

    Returns
    -------
    pd.DataFrame
        A DataFrame comprising three columns: 'transect_num', and the column names supplied by
        the argument `coordinate_names`.
    """

    # Apply the latitude filter
    transect_thresholded = transect_df.loc[transect_df["latitude"] < latitude_threshold]

    # Parse the western-most coordinate indices of each transect
    western_extent_idx = transect_thresholded.groupby(["transect_num"])[
        coordinate_names[0]
    ].idxmin()

    # Subset the DataFrame
    transect_western_extent = transect_thresholded.loc[western_extent_idx].reset_index(drop=True)

    # Return the reduced DataFrame
    return transect_western_extent.filter(["transect_num", *coordinate_names])


def uniform_search_strategy(
    sparse_radii: np.ndarray[int],
    valid_distances: np.ndarray[int],
    local_points: np.ndarray[float],
    k_min: int,
    k_max: int,
    search_radius: float,
    wr_indices: np.ndarray[int],
    oos_indices: np.ndarray[np.number],
    oos_weights: np.ndarray[float],
    **kwargs,
) -> Tuple[np.ndarray[np.number], np.ndarray[np.number], np.ndarray[np.number]]:
    """
    Uniform extrapolation search strategy for finding (and weighting) k-th nearest points
    (relative to a reference coordinate) required for computing the lagged semivariogram in an
    adaptive approach

    Parameters
    ----------
    sparse_radii : np.ndarray[int]
        Indices where there are fewer than `k_min` nearest neighbors.
    valid_distances : np.ndarray[int]
        The number of masked distance matrix values where extrapolation is required.
    local_points : np.ndarray[float]
        An array with the sorted distances (from nearest to furthest) relative to each point.
    k_min : int
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.
    k_max : int
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.
    search_radius : float
        The adaptive search radius that identifies the *k*-nearest neighbors around each
        georeferenced value that are subsequently kriged.
    wr_indices : np.ndarray[int]
        Indices of within-radius (WR) (i.e. < `k_max`) points.
    oos_indices : np.ndarray[np.number]
        Template array based on the size of the data input and `k_min` that will contain indices
        where extrapolation is required where there are fewer than `k_min` nearest neighbors.
    oos_weights : np.ndarray[float]
        Weights applied to extraplolated values.

    Returns
    -------
    Tuple[np.ndarray[np.number], np.ndarray[np.number], np.ndarray[np.number]]
        A tuple with updated values for `wr_indices`, `oos_indices`, and `oos_weights` via a
        search strategy that applies unconstrained and uniform extrapolation to out-of-sample (OOS)
        points.
    """

    # Index for areas with some valid points but fewer than k_min
    partial_indices = sparse_radii[valid_distances[sparse_radii] > 0]

    # Create boolean mask for within-range/sample points
    wr_mask = local_points[partial_indices, :k_max] < search_radius

    # Create temporary matrix for within-range samples
    wr_tmp = wr_indices[partial_indices].copy()
    # ---- Update temporary matrix by applying `wr_mask` for wr points
    wr_tmp[~wr_mask] = np.nan

    # Create temporary matrix for oos samples
    oos_tmp = wr_indices[partial_indices].copy()
    # ---- Update temporary matrix by applying `wr_mask` for oos points
    oos_tmp[wr_mask] = np.nan

    # Assign the OOS values to `oos_indices`
    oos_indices[partial_indices] = np.sort(oos_tmp[:, :k_min])

    # Apply the mask to the remaining `wr_indices` values
    wr_indices[partial_indices] = np.sort(wr_tmp[:, :k_max])

    # Get areas with no valid points within the search radius
    full_extrap_indices = sparse_radii[valid_distances[sparse_radii] == 0]
    if len(full_extrap_indices) > 0:
        # ---- Use all `k_min`-nearest neighbors for extrapolation
        oos_indices[full_extrap_indices] = wr_indices[full_extrap_indices, :k_min]
        wr_indices[full_extrap_indices] = np.nan
        # ---- Compute the OOS kriging weights
        oos_mean = np.apply_along_axis(np.nanmean, 1, local_points[full_extrap_indices, :k_min])
        # ---- Exponentiate the OOS mean
        oos_exp = np.exp(-oos_mean / search_radius)
        # ---- Update the OOS kriging weights
        oos_weights[full_extrap_indices] = oos_exp

    # Return Tuple
    return wr_indices, oos_indices, oos_weights


def western_boundary_search_strategy(
    kriging_mesh: pd.DataFrame,
    western_extent: pd.DataFrame,
    coordinate_names: Tuple[str, str],
    sparse_radii: np.ndarray[int],
    valid_distances: np.ndarray[int],
    local_points: np.ndarray[float],
    distance_matrix_masked: np.ndarray[float],
    nearby_indices: np.ndarray[int],
    k_min: int,
    k_max: int,
    search_radius: float,
    wr_indices: np.ndarray[int],
    oos_indices: np.ndarray[np.number],
    oos_weights: np.ndarray[float],
    **kwargs,
) -> Tuple[np.ndarray[np.number], np.ndarray[np.number], np.ndarray[np.number]]:
    """
    Search strategy that applies western boundary constraints for transect-based surveys

    Parameters
    ----------
    kriging_mesh : pd.DataFrame
        Kriging mesh used for interpolated data values via geostatistics.
    western_extent : pd.DataFrame
        DataFrame with the western-most extent of each transect line used for re-weighting the
        out-of-sample/extrapolated kriged values.
    coordinate_names : Tuple[str, str]
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    sparse_radii : np.ndarray[int]
        Indices where there are fewer than `k_min` nearest neighbors.
    valid_distances : np.ndarray[int]
        The number of masked distance matrix values where extrapolation is required.
    local_points : np.ndarray[float]
        An array with the sorted distances (from nearest to furthest) relative to each point.
    distance_matrix_masked : np.ndarray[float]
        An array with the search-radius-masked nearest neighbor distances.
    nearby_indices : np.ndarray[int]
        Indices of points that require extrapolation.
    k_min : int
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.
    k_max : int
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.
    search_radius : float
        The adaptive search radius that identifies the *k*-nearest neighbors around each
        georeferenced value that are subsequently kriged.
    wr_indices : np.ndarray[int]
        Indices of within-radius (WR) (i.e. < `k_max`) points.
    oos_indices : np.ndarray[np.number]
        Template array based on the size of the data input and `k_min` that will contain indices
        where extrapolation is required where there are fewer than `k_min` nearest neighbors.
    oos_weights : np.ndarray[float]
        Weights applied to extraplolated values.

    Returns
    -------
    Tuple[np.ndarray[np.number], np.ndarray[np.number], np.ndarray[np.number]]
        A tuple with updated values for `wr_indices`, `oos_indices`, and `oos_weights` via a
        search strategy that uses an extrapolation re-weighting based on transect extents.
    """

    # Parse ordered coordinate names
    x_name, y_name = coordinate_names

    # Index the mesh grid coordinates for bounding the search radius expansion/extrapolation
    # ---- y-coordinates with array transformation to access matrix operations
    mesh_y = kriging_mesh[y_name].to_numpy()[sparse_radii].reshape(-1, 1)
    # ---- x-coordinates
    mesh_x = kriging_mesh[x_name].to_numpy()[sparse_radii]

    # Calculate the mesh distance from the western boundary of the survey transects
    # ---- Find closest point
    mesh_western_distance = np.abs(mesh_y - western_extent["y"].to_numpy()).argmin(axis=1)
    # ---- Calculate the western limits (x-axis)
    western_limit = western_extent.iloc[np.ravel(mesh_western_distance)][x_name]
    # ---- Compute bounding threshold (for tapered extrapolation function)
    western_threshold = western_limit - search_radius
    # ---- Create a thresholded mask for lazy operations
    western_limit_mask = mesh_x < western_threshold

    # Adjust values that don't fall outside the western extent
    if np.any(~western_limit_mask):
        # ---- Grab all values that don't fall outside the western extent
        soft_extrapolation_index = sparse_radii[~western_limit_mask]
        # ---- Find the local points where there are at least some valid points
        partial_indices = soft_extrapolation_index[valid_distances[soft_extrapolation_index] > 0]
        # ---- Update the current values in `wr_indices`
        if len(partial_indices) > 0:
            # -------- Create boolean mask for within-range/sample points
            wr_mask = local_points[partial_indices, :k_max] < search_radius
            # -------- Create temporary matrix for within-range samples
            wr_tmp = wr_indices[partial_indices].copy()
            # -------- Create temporary matrix for oos samples
            oos_tmp = wr_indices[partial_indices].copy()
            # -------- Update temporary matrix by applying `wr_mask` for wr points
            wr_tmp[~wr_mask] = np.nan
            # -------- Update temporary matrix by applying `wr_mask` for oos points
            oos_tmp[wr_mask] = np.nan
            # -------- Assign the OOS values to `oos_indices`
            oos_indices[partial_indices] = np.sort(oos_tmp[:, :k_min])
            # -------- Apply the mask to the remaining `wr_indices` values
            wr_indices[partial_indices] = np.sort(wr_tmp[:, :k_max])

        # ---- Find the local points where there are no valid points within the search radius
        full_extrap_indices = soft_extrapolation_index[
            valid_distances[soft_extrapolation_index] == 0
        ]
        if len(full_extrap_indices) > 0:
            # -------- Update `oos_indices`
            oos_indices[full_extrap_indices] = wr_indices[full_extrap_indices, :k_min]
            # -------- Update `wr_indices`
            wr_indices[full_extrap_indices] = np.nan

    # Taper function for extrapolating values outside the search radius
    if np.any(western_limit_mask):
        # ---- Index these values
        extrapolation_index = sparse_radii[western_limit_mask]
        # ---- Compute the OOS kriging weights
        oos_mean = np.apply_along_axis(np.nanmean, 1, local_points[extrapolation_index, :k_min])
        # ---- Exponentiate the OOS mean
        oos_exp = np.exp(-oos_mean / search_radius)
        # ---- Update the OOS weights
        oos_weights[extrapolation_index] = oos_exp
        # ---- Get the outside indices that correspond to this tapered extrapolation
        sparse_extrapolation_index = nearby_indices[western_limit_mask].astype(float)
        # ---- Apply indices as a mask to the NaN-masked distance matrix
        extrapolated_distance = np.take_along_axis(
            distance_matrix_masked[sparse_radii][western_limit_mask],
            sparse_extrapolation_index.astype(int),
            axis=1,
        )
        # ---- Create NaN mask
        extrapolated_nan_mask = ~np.isnan(extrapolated_distance)
        # -------- Apply mask to indices
        sparse_extrapolation_index_nan = sparse_extrapolation_index.copy()
        sparse_extrapolation_index_nan[extrapolated_nan_mask] = np.nan
        # -------- Update `out_of_sample_indices` matrix
        oos_indices[extrapolation_index] = np.sort(sparse_extrapolation_index_nan)
        # ---- Get inside indices that apply to these points
        # -------- Create NaN mask for within-sample values
        interpolated_nan_mask = np.isnan(extrapolated_distance)
        # -------- Apply mask to indices
        sparse_interpolation_index_nan = sparse_extrapolation_index.copy()
        sparse_interpolation_index_nan[interpolated_nan_mask] = np.nan
        # -------- Pad NaN to match `within_sample_indices` matrix
        sparse_interpolation_pad = np.pad(
            sparse_interpolation_index_nan,
            [(0, 0), (0, k_max - k_min)],
            mode="constant",
            constant_values=np.nan,
        )
        # -------- Updated `within_sample_indices` matrix
        wr_indices[extrapolation_index] = np.sort(sparse_interpolation_pad)

    # Return Tuple
    return wr_indices, oos_indices, oos_weights


def search_radius_mask(distance_matrix: np.ndarray, search_radius: float) -> np.ndarray[float]:
    """
    Generate a mask for the search radius to identify points within the search distance.

    Parameters
    ----------
    distance_matrix : np.ndarray
        Distance matrix between points.
    search_radius : float
        Maximum search distance for identifying neighbors.

    Returns
    -------
    np.ndarray
        Masked distance matrix where distances beyond search_radius are set to NaN.
    """
    # Create a copy to avoid modifying the original
    masked_matrix = distance_matrix.copy()
    # Set distances beyond search radius to NaN
    masked_matrix[masked_matrix > search_radius] = np.nan
    return masked_matrix


def count_within_radius(distance_matrix_masked: np.ndarray) -> np.ndarray[int]:
    """
    Count the number of valid (non-NaN) distances for each point in the masked distance matrix.

    Parameters
    ----------
    distance_matrix_masked : np.ndarray
        Masked distance matrix where invalid distances are NaN.

    Returns
    -------
    np.ndarray
        Array containing the count of valid neighbors for each point.
    """
    # Count non-NaN values along each row
    return np.sum(~np.isnan(distance_matrix_masked), axis=1)


def adaptive_search_radius(
    distance_matrix: np.ndarray[float],
    k_min: int,
    k_max: int,
    search_radius: float,
    search_strategy: Callable = uniform_search_strategy,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray[np.number], np.ndarray[np.number], np.ndarray[float]]:
    """
    Find the indices of the k-th nearest points (relative to a reference coordinate) required
    for computing the lagged semivariogram

    Parameters
    ----------
    distance_matrix : np.ndarray[float]
        An array/matrix that includes the distances of each mesh points from every
        georeferenced along-transect interval.
    k_min : int
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.
    k_max : int
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.
    search_radius : float
        The adaptive search radius that identifies the *k*-nearest neighbors around each
        georeferenced value that are subsequently kriged.
    search_strategy : Callable, default=:func:`echopop.nwfsc_feat.spatial.uniform_search_strategy`
        A `Callable` function that defaults to using a uniform search strategy where out-of-sample
        points are extrapolated using equal weights. User-defined search strategies can be defined
        by parsing any of the internally computed variables:

        - **sparse_radii : np.ndarray[int]** \n
        Indices where there are fewer than `k_min` nearest neighbors.

        - **valid_distances : np.ndarray[int]** \n
        The number of masked distance matrix values where extrapolation is required.

        - **local_points : np.ndarray[float]**\n
        An array with the sorted distances (from nearest to furthest) relative to each point.

        - **distance_matrix_masked : np.ndarray[float]** \n
        An array with the search-radius-masked nearest neighbor distances.

        - **nearby_indices : np.ndarray[int]** \n
        Indices of points that require extrapolation.

        - **k_min : int** \n
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.

        - **k_max : int** \n
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.

        - **search_radius : float** \n
        The adaptive search radius that identifies the *k*-nearest neighbors around each
        georeferenced value that are subsequently kriged.

        - **wr_indices : np.ndarray[int]** \n
        Indices of within-radius (WR) (i.e. < `k_max`) points.

        - **oos_indices : np.ndarray[np.number]** \n
        Template array based on the size of the data input and `k_min` that will contain indices
        where extrapolation is required where there are fewer than `k_min` nearest neighbors.

        - **oos_weights : np.ndarray[float]**
        Weights applied to extraplolated values.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray[np.number], np.ndarray[np.number], np.ndarray[float]]
        A tuple with arrays representing the k-nearest distances, within-radius (WR) indices,
        out-of-sample (OOS) indices that indicate which points will be extrapolated, and the
        OOS weights applied to values that are extrapolated.
    """

    # Generate the search radius mask
    distance_matrix_masked = search_radius_mask(distance_matrix, search_radius)

    # Identify mesh points that require extrapolation
    # ---- Count the number of values within the search radius
    valid_distances = count_within_radius(distance_matrix_masked)
    # ---- Identify rows where the number of valid points are less than `k_min`
    sparse_radii = np.hstack(np.where(valid_distances < k_min))

    # Calculate the closest grid points and their indices
    # ---- Closest indices
    local_indices = distance_matrix.argsort(axis=1)
    # ---- Map the distance matrix to the local indices
    local_points = np.take_along_axis(distance_matrix, local_indices, axis=1)

    # Initialize matrices
    # ---- Within-radius (WR) samples
    wr_indices = local_indices[:, :k_max].astype(float)
    # ---- Out-of-sample (OOS) indices
    oos_indices = np.full((len(valid_distances), k_min), np.nan)
    # ---- OOS weights
    oos_weights = np.ones(len(valid_distances))

    # For points where there are fewer than `k_min` points within the search radius, extrapolate
    if np.size(sparse_radii) > 0:
        # ---- Extract the indices requiring extrapolation
        nearby_indices = local_indices[sparse_radii][:, :k_min]
        # ---- Update local points
        # -------- Fill NaN values
        local_points[sparse_radii, k_min:] = np.nan
        wr_indices[sparse_radii, k_min:] = np.nan
        # ---- Collect the available parameters for defining the search strategy
        search_params = {
            "sparse_radii": sparse_radii,
            "valid_distances": valid_distances,
            "local_points": local_points,
            "distance_matrix_masked": distance_matrix_masked,
            "nearby_indices": nearby_indices,
            "k_min": k_min,
            "k_max": k_max,
            "search_radius": search_radius,
            "wr_indices": wr_indices,
            "oos_indices": oos_indices,
            "oos_weights": oos_weights,
        }

        # ---- Apply search strategy to update the WR/OOS indices, and OOS weights
        if isinstance(search_strategy, Callable):
            wr_indices, oos_indices, oos_weights = search_strategy(**search_params)

    # Return output
    return local_points[:, :k_max], wr_indices, oos_indices, oos_weights


def ordinary_kriging_matrix(
    local_distance_matrix: np.ndarray[float], variogram_parameters: dict[str, Any]
) -> np.ndarray[float]:
    """
    Calculate the kriging covariance matrix

    Parameters
    ----------
    local_distance_matrix : np.ndarray[float]
        A 2D array with the numeric distances between points.
    variogram_parameters : Dict[str, Any]
        Dictionary describing the variogram model and parameters, such as:
        - 'model' (str or list of str): variogram model names (e.g., 'bessel', 'exponential'),
        - 'nugget' (float): nugget effect,
        - 'sill' (float): sill parameter,
        - 'correlation_range' (float): range parameter,
        - 'hole_effect_range' (float): hole effect range,
        - 'decay_power' (float): power parameter for power variogram.

    Returns
    -------
    np.ndarray[float]
        A 2D array representing the (n+1, n+1) kriging covariance matrix including the unbiasedness
        constraint.

    Notes
    -----
    Ordinary Kriging assumes a constant but unknown mean, requiring the
    kriging weights to sum to 1:

    .. math::

        \\sum_{i=1}^n \\lambda_i = 1

    This leads to the augmented linear system:

    .. math::

        \\begin{bmatrix}
        K & \\mathbf{1} \\\\
        \\mathbf{1}^\\top & 0
        \\end{bmatrix}
        \\begin{bmatrix}
        \\boldsymbol{\\lambda} \\\\
        \\mu
        \\end{bmatrix}
        =
        \\begin{bmatrix}
        \\mathbf{k} \\\\
        1
        \\end{bmatrix}

    Where:

    - :math:`K` is the :math:`n \\times n` covariance or variogram matrix between known points
    - :math:`\\mathbf{k}` is the vector of covariances (or variograms) between the known points and
      the interpolation location
    - :math:`\\boldsymbol{\\lambda}` is the vector of kriging weights
    - :math:`\\mu` is the Lagrange multiplier enforcing the unbiasedness constraint
    - :math:`\\mathbf{1}` is a column vector of ones

    The solution yields the kriging weights :math:`\\lambda_i` and the multiplier :math:`\\mu`.
    """

    # Calculate the covariance/kriging matrix (without the constant term)
    kriging_matrix_initial = variogram(
        distance_lags=local_distance_matrix, variogram_parameters=variogram_parameters
    )

    # Expand the covariance/kriging matrix with a constant
    # ---- In Ordinary Kriging, this should be '1'
    kriging_matrix = np.concatenate(
        [kriging_matrix_initial, np.ones((local_distance_matrix.shape[0], 1))], axis=1
    )

    # Add column and row of ones for Ordinary Kriging
    kriging_matrix = np.concatenate(
        [kriging_matrix, np.ones((1, local_distance_matrix.shape[0] + 1))], axis=0
    )

    # Diagonal fill (0.0)
    np.fill_diagonal(kriging_matrix, 0.0)

    return kriging_matrix


def kriging_lambda(
    anisotropy: float,
    lagged_semivariogram: np.ndarray,
    kriging_covariance_matrix: np.ndarray,
):
    """
    Solve the Kriging system (set of linear equations) using a truncated singular value
    decomposition (SVD) of the kriging covariance matrix.

    Parameters
    ----------
    anisotropy : float
        Directional threshold for retaining singular values in the decomposition, given as a
        fraction of the largest singular value. Singular values below this threshold are discarded
        to regularize the solution.
    lagged_semivariogram: np.array[float]
        Right-hand side vector of the kriging system comprising the local lagged semivariogram.
        This corresponds to the semivariogram estimates between prediction location and observed
        sample points, with an additional `1` appended for the unbiasedness constrain when using
        ordinary kriging.
    kriging_covariance_matrix : np.array[float]
        A 2D array representing the full kriging covariance matrix, which is typically of size
        (n+1, n+1) when using ordinary kriging. This includes the variogram matrix between sample
        points. When ordinary kriging is being used, an additional column and row of `1` and a
        bottom-right `0` are also included to enforce the unbiasedness constraint:

        .. math::

            \\begin{bmatrix}
            \\gamma_{ij} & 1 \\\\
            1^T & 0
            \\end{bmatrix}

    Returns
    -------
    np.ndarray[float]
        A 1D array of length :math:`n + 1`, where the first :math:`n` values are the kriging
        weights :math:`\\lambda_i`, and the final element is the Lagrange multiplier :math:`\\mu`
        enforcing the unbiasedness constraint.

    Notes
    -----
    This function performs a numerically stable inversion of the kriging system using a truncated
    SVD:

    .. math::

        K^+ = V \\Sigma^{-1} U^T

    where only the singular values satisfying:

    .. math::

        \\frac{\\sigma_i}{\\sigma_1} > \\text{anisotropy}

    are retained.

    The final weights are given by:

    .. math::

        \\boldsymbol{\\lambda} = K^+ \\cdot \\mathbf{k}

    where :math:`\\mathbf{k}` is the lagged semivariogram vector to the target location.
    """

    # Singular value decomposition (SVD)
    # ---- U: left singular vectors (directions of maximum variance)
    # ---- Sigma: singular values (amount of variance captured by each singular vector, U)
    # ---- VH: conjugate transpose of the right singular vectors
    U, Sigma, VH = np.linalg.svd(kriging_covariance_matrix, full_matrices=True)

    # Create Sigma mask informed by the ratio-threshold
    # ---- The ratio between all singular values and their respective
    # ---- maximum is used to apply a mask that informs which values
    # ---- are used to tabulate the kriging weights (aka lambda)
    Sigma_mask = np.abs(Sigma / Sigma[0]) > anisotropy

    # Inverse masked semivariogram (K)
    K_inv = np.matmul(
        np.matmul(VH.T[:, Sigma_mask], np.diag(1.0 / Sigma[Sigma_mask])), U[:, Sigma_mask].T
    )

    # Calculate kriging weights (lambda)
    return np.dot(K_inv, lagged_semivariogram)


def parse_stacked_kriging_array(
    stacked_array: np.ndarray[float], k_min: int, k_max: int, **kwargs
) -> Tuple[np.ndarray[float], np.ndarray[int], np.ndarray[float], np.ndarray[float]]:
    """
    Helper function for parsing the horizontally stacked array when interpolating points via kriging

    Parameters
    ----------
    stacked_array : np.ndarray[np.number]
        A 1D array that is encoded with various parameters including the out-of-sample weights,
        extrapolation and interpolation indices, local lagged semivariogram, and range distance
        estimates.
    k_min : int
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.
    k_max : int
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.

    Returns
    -------
    Tuple[np.ndarray[float], np.ndarray[int], np.ndarray[float], np.ndarray[float]]
        A tuple comprising the out-of-sample weights, composite indices for values that will be
        either interpolated or extrapolated, the local lagged semivariogram, and the range distance
        values.
    """

    # Break up the stacked index
    # ---- Extrapolation/out-of-sample (OOS) weight
    outside_weight = stacked_array[0]
    # ---- Get the Within- and out-of-sample indices
    interp_indices, extrap_indices = (
        stacked_array[s:e][~np.isnan(stacked_array[s:e])].astype(int)
        for s, e in [(1, k_max + 1), (k_max + 1, k_max + 1 + k_min)]
    )
    # ---- Combine the interpolated and extrapolated indices into a composite index
    composite = np.concatenate([interp_indices, extrap_indices])
    # ---- Extract the local lagged semivariogram and range distance values
    M2_vario, range_vals = (
        stacked_array[s:e][~np.isnan(stacked_array[s:e])]
        for s, e in [(-(2 * k_max) - 1, -k_max), (-k_max, None)]
    )

    # Return the values
    return outside_weight, composite, M2_vario, range_vals


def kriging_point_estimator(
    kriging_array: np.ndarray[np.number],
    data_array: np.ndarray[float],
    kriging_parameters: Dict[str, Any],
    variogram_parameters: Dict[str, Any],
) -> np.ndarray[float]:
    """
    Interpolate value at a specified point via kriging

    Parameters
    ----------
    kriging_array: np.ndarray
        Horizontally stacked 1D array containing kriging indices and local semivariogram values.
        Expected to include:
        - Outside weight scalar,
        - Within- and out-of-sample indices,
        - Lagged semivariogram values (M2),
        - Local range estimates.
    data_array: np.ndarray
        A 2D array of shape (n, 3), where columns represent [x-coordinates, y-coordinate,
        variable].
    kriging_parameters: dict
        Dictionary of kriging parameters, must contain keys:
        - 'k_min' (int): minimum number of neighbors,
        - 'k_max' (int): maximum number of neighbors,
        - 'search_radius' (float): radius to consider neighbors,
        - 'anisotropy' (float): truncation threshold for singular values in SVD solver.
    variogram_parameters : Dict[str, Any]
        Dictionary describing the variogram model and parameters, such as:
        - 'model' (str or list of str): variogram model names (e.g., 'bessel', 'exponential'),
        - 'nugget' (float): nugget effect,
        - 'sill' (float): sill parameter,
        - 'correlation_range' (float): range parameter,
        - 'hole_effect_range' (float): hole effect range,
        - 'decay_power' (float): power parameter for power variogram.

    Returns
    -------
    np.ndarray[float]
        A 1D array of length 3 containing:
        - point estimate (float): the kriging predicted value at the target location,
        - kriged variance (float): variance estimate associated with the kriging prediction,
        - sample variance (float): coefficient of variation based variance estimate,
          set to NaN if the point estimate is effectively zero to avoid division by zero.

    Notes
    -----
    - The function assumes ordinary kriging with an unbiasedness constraint.
    - The sample variance is computed as:

      .. math::

        \text{sample_variance} = \frac{
            \\sqrt{\\text{kriged_variance} \\times \\mathrm{var}(Z)}
        }{
            |\\hat{Z}|
        }

      where :math:`\\hat{Z}` is the kriging estimate and :math:`Z` the observed variable values.
    - If the point estimate is zero or extremely close to zero, the sample variance is
      undefined and returned as NaN.
    """

    # Break up the stacked index
    outside_weight, composite_index, M2_vario, range_vals = parse_stacked_kriging_array(
        kriging_array, **kriging_parameters
    )

    # Index the data array
    data_array_idx = data_array[composite_index]

    # Set extrapolated variable values to 0.0
    # Only apply the mask if range_vals has the same length as data_array_idx
    if len(range_vals) == len(data_array_idx):
        data_array_idx[range_vals > kriging_parameters["search_radius"], 2] = 0.0

    # Calculate the local distance matrix
    local_distance_matrix, _ = lag_distance_matrix(
        coordinates_1=data_array_idx[:, 0], coordinates_2=data_array_idx[:, 1], self=True
    )

    # Compute the covariance matrix using ordinary kriging
    kriging_covariance_matrix = ordinary_kriging_matrix(local_distance_matrix, variogram_parameters)

    # Calculate the kriging weights (lambda)
    kriging_weights = kriging_lambda(
        kriging_parameters["anisotropy"], M2_vario, kriging_covariance_matrix
    )

    # Calculate the point estimate
    point_estimate = (
        kriging_weights[: len(composite_index)] * data_array_idx[:, 2]
    ).sum() * outside_weight

    # Calculate the kriged variance
    kriged_variance = kriging_weights.dot(M2_vario)

    # Calculate the sample variance and CV
    if abs(point_estimate) < np.finfo(float).eps:
        sample_variance = np.nan
    else:
        sample_variance = np.sqrt(kriged_variance * np.var(data_array_idx[:, 2], ddof=1)) / np.abs(
            point_estimate
        )

    # Return the point estimates, kriged/sample variances
    return np.array([point_estimate, kriged_variance, sample_variance])


def krige(
    transect_df: pd.DataFrame,
    kriging_mesh: pd.DataFrame,
    coordinate_names: Tuple[str, str],
    variable: str,
    kriging_parameters: Dict[str, Any],
    variogram_parameters: Dict[str, Any],
    adaptive_search_strategy: Callable = uniform_search_strategy,
) -> np.ndarray[float]:
    """
    Use kriging to interoplate data

    Parameters
    ----------
    transect_df : pd.DataFrame
        Dataframe including georeferenced data
    kriging_mesh : pd.DataFrame
        Grid data that has been transformed
    coordinate_names : Tuple[str, str]
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    variable : str
        The variable used for computing the empirical variogram (e.g. 'biomass_density'), which
        must exist as a column in 'transect_df'.
    kriging_parameters: dict
        Dictionary of kriging parameters, must contain keys:
        - 'k_min' (int): minimum number of neighbors,
        - 'k_max' (int): maximum number of neighbors,
        - 'search_radius' (float): radius to consider neighbors,
        - 'anisotropy' (float): truncation threshold for singular values in SVD solver.
    variogram_parameters : Dict[str, Any]
        Dictionary describing the variogram model and parameters, such as:
        - 'model' (str or list of str): variogram model names (e.g., 'bessel', 'exponential'),
        - 'nugget' (float): nugget effect,
        - 'sill' (float): sill parameter,
        - 'correlation_range' (float): range parameter,
        - 'hole_effect_range' (float): hole effect range,
        - 'decay_power' (float): power parameter for power variogram.
    adaptive_search_strategy : Callable, \
        default=:func:`echopop.nwfsc_feat.spatial.uniform_search_strategy`
        A `Callable` function that defaults to using a uniform search strategy where out-of-sample
        points are extrapolated using equal weights. User-defined search strategies can be defined
        by parsing any of the internally computed variables:

        - **sparse_radii : np.ndarray[int]** \n
        Indices where there are fewer than `k_min` nearest neighbors.

        - **valid_distances : np.ndarray[int]** \n
        The number of masked distance matrix values where extrapolation is required.

        - **local_points : np.ndarray[float]**\n
        An array with the sorted distances (from nearest to furthest) relative to each point.

        - **distance_matrix_masked : np.ndarray[float]** \n
        An array with the search-radius-masked nearest neighbor distances.

        - **nearby_indices : np.ndarray[int]** \n
        Indices of points that require extrapolation.

        - **k_min : int** \n
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.

        - **k_max : int** \n
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.

        - **search_radius : float** \n
        The adaptive search radius that identifies the *k*-nearest neighbors around each
        georeferenced value that are subsequently kriged.

        - **wr_indices : np.ndarray[int]** \n
        Indices of within-radius (WR) (i.e. < `k_max`) points.

        - **oos_indices : np.ndarray[np.number]** \n
        Template array based on the size of the data input and `k_min` that will contain indices
        where extrapolation is required where there are fewer than `k_min` nearest neighbors.

        - **oos_weights : np.ndarray[float]**
        Weights applied to extraplolated values.

    Returns
    -------
    A 1D array of the same length as `transect_df` containing:
        - point estimate (float): the kriging predicted value at the target location,
        - kriged variance (float): variance estimate associated with the kriging prediction,
        - sample variance (float): coefficient of variation based variance estimate,
          set to NaN if the point estimate is effectively zero to avoid division by zero.
    """

    # Generate the distance matrix for all mesh points relative to all transect coordinates
    distance_matrix, _ = lag_distance_matrix(
        coordinates_1=kriging_mesh, coordinates_2=transect_df, coordinate_names=coordinate_names
    )

    # Run the adaptive search window to identify which points have to be re-weighted to account
    # for extrapolation
    range_grid, inside_indices, outside_indices, outside_weights = adaptive_search_radius(
        distance_matrix,
        search_strategy=adaptive_search_strategy,
        **kriging_parameters,
    )

    # Calculate the lagged semivariogram (M20)
    local_variogram = np.apply_along_axis(variogram, 1, range_grid, variogram_parameters)

    # Append 1.0 for the ordinary kriging assumptions (M2)
    local_variogram_M2 = np.sort(
        np.hstack((local_variogram, np.ones((local_variogram.shape[0], 1))))
    )

    # Stack all of the kriging parameters
    full_stack = np.hstack(
        (
            outside_weights.reshape(-1, 1).tolist(),
            inside_indices.tolist(),
            outside_indices.tolist(),
            local_variogram_M2.tolist(),
            range_grid.tolist(),
        )
    )

    # Extract x- and y-coordinates, and variable_data
    data_array = transect_df[list(coordinate_names) + [variable]].to_numpy()

    # Kriging point estimator (currently assumes ordinary kriging)
    kriged_values = np.apply_along_axis(
        kriging_point_estimator, 1, full_stack, data_array, kriging_parameters, variogram_parameters
    )

    # Return the 2D array
    return kriged_values


def project_kriging_results(
    kriged_estimates: np.ndarray,
    kriging_mesh: pd.DataFrame,
    transect_df: pd.DataFrame,
    variable: str,
    default_mesh_cell_area: Optional[float] = None,
) -> Tuple[pd.DataFrame, float]:
    """
    Project the kriged results over the defined mesh.

    Parameters
    ----------
    kriged_estimates : np.ndarray, shape (n, 3)
        Kriged values with columns: [estimate, kriging variance, sample variance]
    kriging_mesh : pd.DataFrame
        Grid data that has been transformed
    transect_df : pd.DataFrame
        Dataframe including georeferenced data
    coordinate_names : Tuple[str, str]
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    variable : str
        The variable used for computing the empirical variogram (e.g. 'biomass_density'), which
        must exist as a column in 'transect_df'.
    default_mesh_cell_area: Optional[float], default=None
        Default cell area for kriging mesh nodes when no area is provided within the kriging
        mesh DataFrame.

    Returns
    -------
    Tuple[pd.DataFrame, float]
        A tuple containing the output mesh DataFrame with columns including the kriged estimates
        and variance, sample variance, and cell coefficient of variation (CV). The other value is
        the overall CV computed for the entire kriging mesh.
    """

    # Create copy
    mesh_results = kriging_mesh.copy()

    # Get the adjusted areas
    if "fraction" in kriging_mesh.columns and default_mesh_cell_area is not None:
        mesh_results["area"] = mesh_results["fraction"] * default_mesh_cell_area
    elif "area" not in mesh_results.columns:
        # If no area column exists, create a default area of 1.0 for each cell
        if default_mesh_cell_area is not None:
            mesh_results["area"] = default_mesh_cell_area
        else:
            mesh_results["area"] = 1.0

    # Calculate the (global) survey variance
    survey_variance = transect_df[variable].var()

    # Distribute estimates over an area
    survey_estimate = np.nansum(kriged_estimates[:, 0] * mesh_results["area"])

    # Add the kriged array values
    # ---- Apply kriged estimateadjustment if any values are negative
    if any(kriged_estimates[:, 0] < 0.0):
        # -------- Get count
        N0 = np.sum(kriged_estimates[:, 0] < 0.0)
        # -------- Truncate
        kriged_estimates[kriged_estimates[:, 0] < 0.0, 0] = 0.0
        # -------- Print warning
        warnings.warn(
            f"{N0} invalid kriged estimates (< 0.0) found. These have been replaced with '0.0'; "
            f"however, this truncation may distort error distributions. If this is problematic, "
            f"adjust the choice of variogram and kriging algorithm parameterization.",
            stacklevel=2,
        )
    # ---- Kriged estimates
    mesh_results["estimate"] = kriged_estimates[:, 0]
    # ---- Kriged variance
    mesh_results["kriged_variance"] = kriged_estimates[:, 1]
    # ---- Sample variance
    mesh_results["sample_variance"] = kriged_estimates[:, 2]
    # ---- Mesh node CV
    mesh_results["cell_cv"] = (
        mesh_results["area"].mean()
        * np.sqrt(kriged_estimates[:, 1].dot(survey_variance))
        / survey_estimate
        * np.sqrt(len(kriged_estimates[:, 1]))
    )

    # Calculate the global CV
    global_cv = (
        np.sqrt(np.nansum(kriged_estimates[:, 1] * mesh_results["area"] ** 2) * survey_variance)
        / survey_estimate
    )

    # Return the DataFrame and global CV estimate
    return mesh_results, global_cv


def transect_coordinate_centroid(spatial_grouped: gpd.GeoSeries):
    """
    Calculate the centroid of a given spatial group.

    This function computes the geometric centroid of a collection of spatial points,
    which is useful for determining the center point of transect lines or other
    spatial groupings.

    Parameters
    ----------
    spatial_grouped: gpd.GeoSeries
        A GeoSeries comprising coordinates (i.e. points). Each element should be
        a Point geometry representing spatial locations.

    Returns
    -------
    Point
        A shapely Point object representing the centroid of all input coordinates.

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point
    >>> points = gpd.GeoSeries([Point(0, 0), Point(1, 1), Point(2, 0)])
    >>> centroid = transect_coordinate_centroid(points)
    >>> print(f"Centroid: ({centroid.x:.1f}, {centroid.y:.1f})")
    Centroid: (1.0, 0.3)

    Notes
    -----
    The function uses the union_all() method to combine all geometries before
    calculating the centroid, which ensures proper handling of the spatial
    reference system.
    """

    # Compute the union of all coordinates within `spatial_grouped`
    centroid_point = spatial_grouped.union_all().centroid

    # Return output
    return Point(centroid_point)


def transect_extent(transect_df: pd.DataFrame, projection: str, num_nearest_transects: int):
    """
    Compute the spatial extent of survey transects using convex hull generation.

    This function creates a polygon representing the spatial extent of survey transects
    by generating convex hulls around each transect and its nearest neighbors, then
    unioning all hulls to create the overall survey boundary.

    Parameters
    ----------
    transect_df : pd.DataFrame
        Dataframe containing survey transect data with columns:
        - 'longitude': Longitude coordinates
        - 'latitude': Latitude coordinates
        - 'transect_num': Transect identifier numbers
    projection : str
        EPSG projection code string (e.g., 'epsg:4326' for WGS84)
    num_nearest_transects : int
        Number of nearest neighbor transects to include when generating
        the convex hull around each transect

    Returns
    -------
    shapely.geometry.base.BaseGeometry
        A shapely geometry object representing the union of all transect convex hulls,
        defining the overall spatial extent of the survey area.

    Examples
    --------
    >>> import pandas as pd
    >>> transect_data = pd.DataFrame({
    ...     'longitude': [-125.0, -125.1, -125.2],
    ...     'latitude': [48.0, 48.1, 48.2],
    ...     'transect_num': [1, 2, 3]
    ... })
    >>> extent = transect_extent(transect_data, 'epsg:4326', 2)
    >>> print(f"Extent type: {type(extent)}")
    Extent type: <class 'shapely.geometry.polygon.Polygon'>

    Notes
    -----
    The function performs the following steps:
    1. Converts the DataFrame to a GeoDataFrame with point geometries
    2. Transforms coordinates from WGS84 to UTM for accurate distance calculations
    3. Calculates centroids for each transect
    4. For each transect, finds the nearest neighbor transects
    5. Generates convex hulls around each transect and its neighbors
    6. Returns the union of all convex hulls

    The resulting polygon can be used for spatial filtering, mesh cropping,
    or defining survey boundaries for analysis.
    """

    # Convert to GeoDataFrame
    transect_gdf = gpd.GeoDataFrame(
        transect_df,
        geometry=gpd.points_from_xy(transect_df["longitude"], transect_df["latitude"]),
        crs=projection,
    )

    # Convert from WGS84 to UTM
    wgs84_to_utm(transect_gdf)

    # Calculate the centroid of each transect line
    transect_centroid = transect_gdf.groupby("transect_num")["geometry"].apply(
        transect_coordinate_centroid
    )

    # Generate grouped polygons around each transect line
    # ---- Initialize polygon list
    transect_polygons = []
    # ---- Iterate through each transect
    for transect in transect_centroid.index:
        # ---- Extract coordinates of the transect
        coord_centroid = transect_centroid[transect]
        # ---- Extract all remaining centroids
        other_centroids = transect_centroid[transect_centroid.index != transect].to_frame()
        # ---- Handle case where there's only one transect (no other centroids)
        if len(other_centroids) == 0:
            # -------- Just use the current transect to create a polygon
            unique_transects = np.array([transect])
            transect_coords = transect_gdf[transect_gdf.transect_num.isin(unique_transects)]
            polygon = Polygon(list(transect_coords.geometry))
            transect_polygons.append(polygon.convex_hull)
            continue
        # ---- Calculate the distance between centroids
        other_centroids["distance_centroid"] = other_centroids.geometry.apply(
            lambda g: coord_centroid.distance(g)
        )
        # ---- Find the 'n' nearest transect centroids
        nearest_centroids = other_centroids.distance_centroid.nsmallest(num_nearest_transects)
        # ---- Filter the transect centroids
        nearest_transects = other_centroids[
            other_centroids.distance_centroid.isin(nearest_centroids)
        ]
        # ---- Parse the coordinates of the relevant transect numbers
        unique_transects = np.append(nearest_transects.index, transect)
        # ---- Get the full coordinates of the relevant transects
        transect_coords = transect_gdf[transect_gdf.transect_num.isin(unique_transects)]
        # ---- Generate the local polygon
        polygon = Polygon(list(transect_coords.geometry))
        # ---- Append the convex hull of the transect polygon to `transect_polygons`
        transect_polygons.append(polygon.convex_hull)

    # Merge the polygons via the union of each set
    return unary_union(transect_polygons)
