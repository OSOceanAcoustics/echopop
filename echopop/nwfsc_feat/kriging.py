import warnings
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np
import pandas as pd

from ..spatial.variogram import variogram
from .variogram import lag_distance_matrix

# Set warnings filter
warnings.simplefilter("always")

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

