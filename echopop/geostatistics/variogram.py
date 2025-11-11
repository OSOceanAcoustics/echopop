import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from lmfit import Parameters
from pydantic import ValidationError

from .. import validators as val
from .variogram_models import fit_variogram

# Set warnings filter
warnings.simplefilter("always")


def lag_distance_matrix(
    coordinates_1: Union[pd.DataFrame, np.ndarray],
    coordinate_names: Optional[Tuple[str, str]] = None,
    coordinates_2: Optional[Union[pd.DataFrame, np.ndarray]] = None,
    self: bool = False,
    azimuth_matrix: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate Euclidean distance and optional azimuth matrices between coordinate sets.

    Computes pairwise distances between all combinations of coordinates from two datasets,
    with optional azimuth angle calculation for directional geostatistical analysis.

    Parameters
    ----------
    coordinates_1 : pd.DataFrame or np.ndarray
        First set of coordinates (x, y) as a DataFrame or numpy array.
    coordinate_names : Tuple[str, str], optional
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    coordinates_2 : pd.DataFrame or np.ndarray, optional
        Second set of coordinates (x, y) as a DataFrame or numpy array. If None, uses coordinates_1.
    self : bool, default=False
        If True, calculates distances within the same set of coordinates (ignored when
        coordinates_2 is None).
    azimuth_matrix : bool, default=False
        If True, returns both distance and azimuth angles; if False, returns only distances.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        If azimuth_matrix is True, returns (distance_matrix, azimuth_angles).
        If azimuth_matrix is False, returns (distance_matrix, empty_array).

    Notes
    -----
    The function computes Euclidean distances using:

    .. math::
        d_{ij} = \\sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}

    For azimuth calculations, angles are computed as:

    .. math::
        θ_{ij} = \\arctan\\left(\\frac{y_j - y_i}{x_j - x_i}\\right) \\cdot \\frac{180}{π} + 180°

    The azimuth matrix has NaN values on the diagonal to avoid undefined angles for identical
    coordinate pairs. Azimuth angles range from 0° to 360°.

    This function is optimized for variogram computation where distance matrices are fundamental
    for lag-based spatial correlation analysis.

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Isaaks, E.H. & Srivastava, R.M. (1989). Applied Geostatistics. Oxford University Press.
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
    Apply spatial filtering to distance matrices using boolean masks and azimuth constraints.

    This function extracts elements from 2D spatial matrices (typically distance or lag matrices)
    based on boolean masks and optional directional filtering. It's primarily used in variogram
    computation to ensure each spatial pair is counted once and to implement directional analysis.

    Parameters
    ----------
    data_matrix : np.ndarray
        The 2D matrix to be filtered (e.g., lag distances, estimates).
    mask_matrix : np.ndarray[bool]
        Boolean mask indicating which matrix elements to retain. Typically a triangular mask
        to avoid double-counting symmetric pairs.
    azimuth_matrix : np.ndarray[float], optional
        Matrix of azimuth angles between coordinate pairs. Used for directional filtering.
    azimuth_angle_threshold : float, optional
        Maximum angular deviation (degrees) from the reference direction. Only pairs within
        ±threshold are retained for directional variogram analysis.

    Returns
    -------
    np.ndarray
        1D array of filtered values from the original matrix.

    Notes
    -----
    The filtering process applies both geometric and directional constraints:

    **Geometric Filtering:**
    Uses boolean masks (typically triangular) to extract unique spatial pairs:

    .. math::
        \\text{mask}_{ij} = \\begin{cases}
        \\text{True} & \\text{if } i > j \\\\
        \\text{False} & \\text{otherwise}
        \\end{cases}

    **Directional Filtering:**
    For azimuth-based filtering, retains pairs where:

    .. math::
        |θ_{ij} - θ_0| ≤ \\text{threshold}

    where θ₀ is the reference direction (typically 0° for N-S analysis).

    This dual filtering enables:
    - Efficient variogram computation (avoiding redundant calculations)
    - Directional variogram analysis for anisotropic processes
    - Flexible spatial neighborhood definition

    References
    ----------
    .. [1] Journel, A.G. & Huijbregts, C.J. (1978). Mining Geostatistics. Academic Press.
    .. [2] Goovaerts, P. (1997). Geostatistics for Natural Resources Evaluation. Oxford University
           Press.
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
    Aggregate spatial data into lag bins for empirical variogram computation.

    Performs the quantization step of variogram analysis by grouping spatial data pairs into
    discrete lag distance bins and computing the statistical moments needed for semivariance
    estimation.

    Parameters
    ----------
    estimates : np.ndarray[float]
        1D array of field estimates at each spatial location (e.g., biomass density values).
    lag_matrix : np.ndarray[int]
        2D array of integer lag bin assignments for each spatial pair.
    mask_matrix : np.ndarray[bool]
        Boolean mask for extracting unique spatial pairs (typically triangular).
    azimuth_matrix : np.ndarray[float]
        2D array of azimuth angles between spatial pairs, or empty array if not used.
    n_lags : int
        Number of lag bins for the variogram analysis.
    azimuth_angle_threshold : float, optional
        Angular constraint for directional variogram analysis (degrees).

    Returns
    -------
    Tuple[np.ndarray[int], np.ndarray[float], np.ndarray[float], np.ndarray[float]]
        - lag_counts: Number of spatial pairs in each lag bin
        - lag_estimates: Sum of estimates for pairs in each lag bin
        - lag_estimates_squared: Sum of squared estimates for pairs in each lag bin
        - lag_deviations: Sum of squared differences between paired estimates

    Notes
    -----
    The quantization process computes lag-specific statistics for variogram estimation:

    **Lag Counts:**
    .. math::
        N(h_k) = |\\{(x_i, x_j) : h_{ij} \\in \\text{bin}_k\\}|

    **Lag Statistics:**
    .. math::
        S_1(h_k) = \\sum_{h_{ij} \\in \\text{bin}_k} Z(x_i)

    .. math::
        S_2(h_k) = \\sum_{h_{ij} \\in \\text{bin}_k} Z(x_i)^2

    .. math::
        D(h_k) = \\sum_{h_{ij} \\in \\text{bin}_k} [Z(x_i) - Z(x_j)]^2

    These statistics are the building blocks for computing standardized semivariance:

    .. math::
        γ(h_k) = \\frac{D(h_k)}{2N(h_k) \\cdot σ_{head}(h_k) \\cdot σ_{tail}(h_k)}

    The function uses `numpy.bincount` for efficient aggregation, making it suitable for large
    spatial datasets common in acoustic surveys.

    References
    ----------
    .. [1] Rivoirard, J., et al. (2000). Geostatistics for Estimating Fish Abundance. Blackwell
           Science.
    .. [2] Petitgas, P. (1993). Geostatistics for fish stock assessments. Reviews in Fish Biology
           and Fisheries, 3(4), 307-334.
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
    Compute standardized semivariance for empirical variogram estimation.

    Calculates the final step in empirical variogram computation using the standardized
    approach suitable for heteroscedastic spatial data. This method normalizes semivariance
    by local variance estimates to improve robustness for irregular sampling patterns.

    Parameters
    ----------
    estimates : np.ndarray[float]
        1D array of field estimates at spatial locations.
    lag_estimates : np.ndarray[float]
        Sum of estimates for each lag bin (from quantize_lags).
    lag_estimates_squared : np.ndarray[float]
        Sum of squared estimates for each lag bin (from quantize_lags).
    lag_counts : np.ndarray[int]
        Number of spatial pairs contributing to each lag bin.
    lag_deviations : np.ndarray[float]
        Sum of squared differences between spatial pairs in each lag bin.
    head_index : np.ndarray[int]
        2D array mapping estimates to lag bins for head variance calculation.

    Returns
    -------
    Tuple[np.ndarray[float], np.ndarray[float]]
        - gamma_h: Standardized semivariance values for each lag
        - mean_lag_covariance: Mean covariance across all lags

    Notes
    -----
    The standardized semivariance is computed as:

    .. math::
        γ_{std}(h) = \\frac{1}{2N(h)} \\cdot \\frac{\\sum [Z(x_i) - Z(x_j)]^2}
                     {σ_{head}(h) \\cdot σ_{tail}(h)}

    where:

    **Head Standard Deviation:**
    .. math::
        σ_{head}(h) = \\sqrt{\\frac{1}{N(h)} \\sum [Z(x_i) - \\bar{Z}_{head}(h)]^2}

    **Tail Standard Deviation:**
    .. math::
        σ_{tail}(h) = \\sqrt{\\frac{S_2(h)}{N(h)} - \\left(\\frac{S_1(h)}{N(h)}\\right)^2}

    **Mean Lag Covariance:**
    .. math::
        C_{lag} = \\frac{1}{K} \\sum_{k=1}^{K} σ_{head}(h_k) \\cdot σ_{tail}(h_k)

    This standardization approach provides several advantages:
    - Handles heteroscedastic spatial data more robustly
    - Reduces sensitivity to extreme values
    - Improves model fitting for irregular sampling patterns
    - Maintains interpretability of classical variogram theory

    The method is particularly well-suited for fisheries acoustic data where biomass estimates can
    vary dramatically in magnitude and spatial distribution.

    References
    ----------
    .. [1] Rivoirard, J., et al. (2000). Geostatistics for Estimating Fish Abundance. Blackwell
           Science.
    .. [2] Petitgas, P. (1993). Use of a disjunctive kriging to model areas of high pelagic fish
           density in acoustic fisheries surveys. Aquatic Living Resources, 6(3), 201-209.
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
    with np.errstate(divide="ignore"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            gamma_h = 0.5 * lag_deviations / (lag_counts * partial_sill)
    # ---- Create mask
    valid_mask = (lag_counts > 0) & (partial_sill > 0) & np.isfinite(gamma_h)
    # ---- Apply mask to produce only valid `gamma_h` estimates WITHOUT interpolation
    gamma_h = np.where(valid_mask, gamma_h, 0.0)
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
    data: pd.DataFrame,
    n_lags: int,
    lag_resolution: float,
    azimuth_filter: bool,
    azimuth_angle_threshold: float,
    variable: str = "biomass_density",
    coordinate_names: Tuple[str, str] = ("x", "y"),
    force_lag_zero: bool = True,
) -> Tuple[np.ndarray[float], np.ndarray[float], np.ndarray[int], float]:
    """
    Compute empirical variogram from spatial survey data using standardized semivariance.

    This function implements the complete empirical variogram computation workflow, from distance
    matrix calculation through standardized semivariance estimation. The method is optimized for
    fisheries acoustic survey data with irregular sampling patterns and heterogeneous variance
    structures.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing georeferenced coordinates and target variable. Must include
        columns specified in `coordinate_names` and `variable`.
    n_lags : int
        Number of lag bins for variogram computation. Typical values: 10-30.
    lag_resolution : float
        Distance interval for each lag bin (same units as coordinates).
    azimuth_filter : bool
        If True, computes azimuth angles for directional analysis capabilities.
    azimuth_angle_threshold : float
        Maximum angular deviation (degrees) for directional filtering. Use 180° for
        omnidirectional analysis.
    variable : str, default='biomass_density'
        Column name of the variable for variogram computation.
    coordinate_names : Tuple[str, str], default=('x', 'y')
        Column names for spatial coordinates. Format: (x_column, y_column).
    force_lag_zero : bool, default=True
        If True, prepends lag 0 with γ(0) = 0 to enforce zero nugget assumption.

    Returns
    -------
    Tuple[np.ndarray[float], np.ndarray[float], np.ndarray[int], float]
        - lags: Array of lag distances
        - gamma: Standardized semivariance values
        - lag_counts: Number of pairs contributing to each lag
        - lag_covariance: Mean covariance between head and tail points

    Notes
    -----
    The empirical variogram computation follows these steps:

    **1. Distance Matrix Calculation:**
    .. math::
        D_{ij} = \\sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}

    **2. Lag Quantization:**
    .. math::
        L_{ij} = \\left\\lfloor \\frac{D_{ij}}{\\Delta h} \\right\\rfloor + 1

    **3. Standardized Semivariance:**
    .. math::
        γ(h_k) = \\frac{1}{2N(h_k)} \\cdot \\frac{\\sum [Z(x_i) - Z(x_j)]^2}
                 {σ_{head}(h_k) \\cdot σ_{tail}(h_k)}

    **Triangular Masking:**
    Uses lower triangular matrix extraction to ensure each spatial pair is counted
    exactly once, avoiding redundant calculations in symmetric distance matrices.

    **Quality Indicators:**
    - `lag_counts`: Higher counts indicate more reliable lag estimates
    - `lag_covariance`: Overall measure of spatial correlation structure
    - Lags with few pairs (< 30) should be interpreted cautiously

    The standardized approach is particularly effective for:
    - Acoustic survey data with patchy biomass distributions
    - Irregular sampling grids typical in marine surveys
    - Data with strong heteroscedasticity

    References
    ----------
    .. [1] Rivoirard, J., et al. (2000). Geostatistics for Estimating Fish Abundance. Blackwell
           Science.
    .. [2] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [3] Petitgas, P. (1993). Geostatistics for fish stock assessments. Reviews in Fish Biology
           and Fisheries, 3(4), 307-334.
    """
    # Initialize lag array
    lags = np.concatenate([np.arange(1, n_lags) * lag_resolution])

    # Calculate the distance (and azimuth) matrix
    distance_matrix, azimuth_matrix = lag_distance_matrix(
        coordinates_1=data,
        coordinate_names=coordinate_names,
        self=True,
        azimuth_matrix=azimuth_filter,
    )

    # Convert lag distances to lag numbers
    lag_matrix = np.round(distance_matrix / lag_resolution).astype(int) + 1

    # Extract estimates column
    estimates = data[variable].to_numpy()

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


class Variogram:
    """ "
    Class for calculating and fitting variograms to quantify the spatial variance in population
    variables and other metrics.

    A variogram (or semivariogram) is a fundamental tool in geostatistics that describes
    the spatial correlation structure of a variable as a function of distance. The empirical
    semivariogram γ(h) is defined as:

    .. math::
        γ(h) = \\frac{1}{2N(h)} \\sum_{i=1}^{N(h)} [Z(x_i) - Z(x_i + h)]^2

    where Z(x_i) is the value at location x_i, h is the lag distance, and N(h) is the
    number of pairs separated by distance h.

    The class implements a standardized semivariogram approach where the semivariance
    is normalized by the product of head and tail standard deviations:

    .. math::
        γ_{std}(h) = \\frac{γ(h)}{σ_{head}(h) \\cdot σ_{tail}(h)}

    This standardization helps account for local variance heterogeneity and provides
    more robust variogram estimates for irregularly distributed data [1]_.

    Parameters
    ----------
    lag_resolution : float
        The distance interval represented by each lag. Must be positive and determines
        the spatial resolution of the variogram analysis.
    n_lags : int
        The number of lags used for computing the semivariogram. Must be positive.
        Typical values range from 10-30 depending on data density and spatial extent.
    coordinate_names : Tuple[str, str], default=('x', 'y')
        Column names for the spatial coordinates in input DataFrames.
        Format: (horizontal_axis, vertical_axis).

    Attributes
    ----------
    gamma : np.ndarray or None
        The computed semivariance values for each lag distance.
    lags : np.ndarray or None
        The lag distances used in the variogram analysis.
    lag_counts : np.ndarray or None
        The number of data pairs contributing to each lag estimate.
    lag_covariance : float or None
        The mean lag covariance between head and tail points.
    variogram_params_optimized : dict or None
        Optimized parameters from theoretical variogram model fitting.
    variogram_params_initial : dict or None
        Initial parameters used for theoretical variogram model fitting.
    variogram_fit_initial : float or None
        Mean absolute deviation of initial model fit.
    variogram_fit_optimized : float or None
        Mean absolute deviation of optimized model fit.

    Methods
    -------
    calculate_empirical_variogram(data, variable, azimuth_filter=True,
    azimuth_angle_threshold=180.0, force_lag_zero=True)
        Compute the empirical variogram from transect data.
    fit_variogram_model(model, model_parameters, optimizer_kwargs={})
        Fit a theoretical variogram model to the empirical variogram using weighted least squares.

    Examples
    --------
    Basic usage for computing an empirical variogram:

    >>> import pandas as pd
    >>> import numpy as np
    >>> from echopop.nwfsc_feat.variogram import Variogram
    >>>
    >>> # Create sample spatial data
    >>> data = pd.DataFrame({
    ...     'x': np.random.uniform(0, 100, 50),
    ...     'y': np.random.uniform(0, 100, 50),
    ...     'biomass_density': np.random.exponential(2, 50)
    ... })
    >>>
    >>> # Initialize variogram
    >>> vario = Variogram(lag_resolution=5.0, n_lags=15)
    >>>
    >>> # Compute empirical variogram
    >>> vario.calculate_empirical_variogram(data, 'biomass_density')
    >>>
    >>> # Fit theoretical model
    >>> from lmfit import Parameters
    >>> params = Parameters()
    >>> params.add('sill', value=2.0, min=0.1)
    >>> params.add('nugget', value=0.1, min=0.0)
    >>> params.add('correlation_range', value=20.0, min=1.0)
    >>>
    >>> best_fit = vario.fit_variogram_model('exponential', params)

    Notes
    -----
    The class implements several key applications for fisheries acoustic survey data:

    1. **Standardized Semivariogram**: Uses local variance normalization to handle
       heteroscedastic data common in biological surveys.

    2. **Azimuth Filtering**: Supports directional variogram analysis for anisotropic
       spatial patterns often observed in marine ecosystems.

    3. **Robust Lag Estimation**: Uses weighted binning approaches that account for
       irregular sampling patterns typical in acoustic transect surveys.

    The theoretical variogram models available include single-family models
    (exponential, Gaussian, spherical, Bessel functions) and composite models
    that can capture complex spatial structures with both short-range correlation
    and periodic patterns [2]_, [3]_.

    References
    ----------
    .. [1] Cressie, N.A.C. (1993). *Statistics for Spatial Data*. John Wiley & Sons.
    .. [2] Chilès, J.P., and Delfiner, P. (2012). *Geostatistics: Modeling Spatial
       Uncertainty*. 2nd ed. John Wiley & Sons.
    .. [3] Rivoirard, J., Simmonds, J., Foote, K.G., Fernandes, P., and Bez, N. (2000).
       *Geostatistics for Estimating Fish Abundance*. Blackwell Science.

    """

    def __new__(
        cls,
        lag_resolution: float,
        n_lags: int,
        coordinate_names: Tuple[str, str] = ("x", "y"),
    ):
        # Validate
        try:
            # ---- Check
            valid_args = val.ValidateVariogramClass.create(
                **dict(
                    lag_resolution=lag_resolution, n_lags=n_lags, coordinate_names=coordinate_names
                )
            )
        # Break creation
        except (ValidationError, Exception) as e:
            raise e from None

        # Create instance
        self = super().__new__(cls)

        # Update attributes
        self.__dict__.update(valid_args)

        # Generate
        return self

    def __init__(self, lag_resolution: float, n_lags: int, coordinate_names: Tuple[str, str]):

        # Initialize variables
        self.gamma = None
        self.lags = None
        self.lag_counts = None
        self.lag_covariance = None
        self.variogram_params_optimized = None
        self.variogram_params_inital = None
        self.variogram_fit_initial = None
        self.variogram_fit_optimized = None

    def calculate_empirical_variogram(
        self,
        data: pd.DataFrame,
        variable: str,
        azimuth_filter: bool = True,
        azimuth_angle_threshold: float = 180.0,
        force_lag_zero: bool = True,
    ) -> None:
        """
        Compute the empirical variogram from transect data.

        Calculates the standardized semivariogram using the method described in
        Rivoirard et al. (2000) [1]_, which is particularly suitable for fisheries
        acoustic data with irregular sampling patterns and heterogeneous variance.

        The empirical semivariogram is computed as:

        .. math::
            γ(h) = \\frac{1}{2N(h)} \\sum_{i=1}^{N(h)}
            \\frac{[Z(x_i) - Z(x_i + h)]^2}{σ_{head}(h) \\cdot σ_{tail}(h)}

        where the standardization by head and tail standard deviations accounts for
        local variance heterogeneity.

        Parameters
        ----------
        data : pd.DataFrame
            DataFrame containing georeferenced coordinates and the target variable. Must include
            columns specified in `coordinate_names` and `variable`.
        variable : str
            Column name of the variable for variogram computation (e.g., 'biomass_density').
        azimuth_filter : bool, default=True
            If True, computes azimuth angles between point pairs to enable directional filtering.
            Useful for detecting spatial anisotropy.
        azimuth_angle_threshold : float, default=180.0
            Maximum azimuth angle deviation (in degrees) for including point pairs. Values less
            than 180° enable directional variogram analysis.
        force_lag_zero : bool, default=True
            If True, forces the nugget effect to zero by prepending lag 0 with γ(0) = 0. This
            assumes no measurement error or micro-scale variation.

        Notes
        -----
        The method uses a triangular mask to ensure each point pair is counted only once, avoiding
        redundant calculations in the distance matrix. For acoustic survey data, typical
        `azimuth_angle_threshold` values of 45-90° can help identify along-track vs. across-track
        spatial patterns.

        The lag counts provide important information about the reliability of each lag estimate -
        lags with very few pairs should be interpreted cautiously.

        References
        ----------
        .. [1] Rivoirard, J., Simmonds, J., Foote, K.G., Fernandes, P., and Bez, N. (2000).
           *Geostatistics for Estimating Fish Abundance*. Blackwell Science.
        """

        # Validate
        try:
            # ---- Check
            valid_args = val.ValidateEmpiricalVariogramArgs.create(
                **dict(
                    azimuth_angle_threshold=azimuth_angle_threshold,
                    azimuth_filter=azimuth_filter,
                    coordinate_names=self.coordinate_names,
                    data=data,
                    force_lag_zero=force_lag_zero,
                    variable=variable,
                )
            )
        # Break creation
        except (ValidationError, Exception) as e:
            raise e from None

        # Compute the empirical variogram
        self.lags, self.gamma, self.lag_counts, self.lag_covariance = empirical_variogram(
            n_lags=self.n_lags,
            lag_resolution=self.lag_resolution,
            **valid_args,
        )

    def fit_variogram_model(
        self,
        model: Union[str, List[str]],
        model_parameters: Parameters,
        optimizer_kwargs: Dict[str, Any] = {},
    ) -> Dict[str, Any]:
        """
        Fit a theoretical variogram model to the empirical variogram using weighted least squares
        optimization.

        Fits parametric models of the form:

        .. math::
            γ(h) = C_0 + C_1 \\cdot f(h; θ)

        where C₀ is the nugget, C₁ is the partial sill, f(h; θ) is the correlation function, and θ
        represents model-specific parameters (range, shape, etc.).

        The optimization minimizes the weighted residual sum of squares:

        .. math::
            \\min_{θ} \\sum_{i=1}^{n} w_i [γ_{emp}(h_i) - γ_{model}(h_i; θ)]^2

        where w_i are weights proportional to the number of pairs at each lag.

        Parameters
        ----------
        model : str or List[str]
            Theoretical variogram model name(s). Single string for basic models (e.g.,
            'exponential', 'gaussian', 'spherical'). List of two strings for composite models
            (e.g., ['bessel', 'exponential']).
        model_parameters : lmfit.Parameters
            Initial parameter values and constraints for optimization. Required parameters depend
            on the chosen model. See `variogram_models` module for model-specific requirements.
        optimizer_kwargs : dict, default={}
            Additional arguments passed to `lmfit.minimize()`. Common options include 'max_nfev'
            (maximum function evaluations) and solver-specific parameters.

        Returns
        -------
        dict
            Dictionary of optimized parameter values with parameter names as keys.

        Notes
        -----
        **Available Models:**

        *Single Models:*
        - 'exponential': Exponential decay, suitable for continuous processes
        - 'gaussian': Gaussian (squared exponential), very smooth processes
        - 'spherical': Spherical model with finite range
        - 'jbessel': J-Bessel function, exhibits hole-effect patterns
        - 'linear': Linear growth (unbounded)

        *Composite Models:*
        - ['bessel', 'exponential']: Combines periodic and exponential decay
        - ['bessel', 'gaussian']: Combines periodic and Gaussian smoothness
        - ['cosine', 'exponential']: Cosine modulation with exponential decay

        **Parameter Guidelines:**
        - `nugget`: Usually 0-20% of total variance
        - `sill`: Total variance of the process
        - `correlation_range`: Distance where correlation becomes negligible
        - `hole_effect_range`: Controls periodicity in Bessel/cosine models

        The method uses Trust Region Reflective algorithm for bounded optimization, which handles
        parameter constraints robustly [1]_.

        References
        ----------
        .. [1] Branch, M.A., Coleman, T.F., and Li, Y. (1999). A subspace, interior, and conjugate
           gradient method for large-scale bound-constrained minimization problems. *SIAM Journal
           on Scientific Computing*, 21(1), 1-23.
        """

        # Validate
        try:
            # ---- Check
            valid_args = val.ValidateFitVariogramArgs.create(
                **dict(
                    model=model,
                    model_parameters=model_parameters,
                    optimizer_kwargs=optimizer_kwargs,
                )
            )
        # Break creation
        except (ValidationError, Exception) as e:
            raise e from None

        # Reconvert Parameters
        valid_args["model_parameters"] = model_parameters

        # Get the initial parameters and store
        self.variogram_params_inital = Parameters.valuesdict(valid_args["model_parameters"])

        # Fit the optimized theoretical variogram model
        (
            self.variogram_params_optimized,
            self.variogram_params_inital,
            self.variogram_fit_optimized,
        ) = fit_variogram(
            lags=self.lags,
            lag_counts=self.lag_counts,
            gamma=self.gamma,
            **valid_args,
        )

        # Return the best-fit variogram model parameters
        return self.variogram_params_optimized
