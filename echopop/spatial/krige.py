import numpy as np
import pandas as pd

from ..spatial.mesh import griddify_lag_distances
from ..spatial.transect import define_western_extent
from ..spatial.variogram import variogram


def kriging(transect_data: pd.DataFrame, mesh_data: pd.DataFrame, settings_dict: dict):
    """
    Use kriging to interoplate data

    Parameters
    ----------
    transect_data: pd.DataFrame
        Dataframe including georeferenced data
    mesh_data: pd.DataFrame
        Grid data that has been transformed
    settings_dict: dict
        Kriging and variogram model parameters
    """

    # Extract biological variable values
    # ---- Define the variable name
    if "density" not in settings_dict["variable"]:
        variable_name = settings_dict["variable"] + "_density"
    else:
        variable_name = settings_dict["variable"]
    # ---- Extract array
    variable_data = transect_data[variable_name].to_numpy()

    # Generate the distance matrix for each mesh point relative to all transect coordinates
    distance_matrix = griddify_lag_distances(mesh_data, transect_data)

    # Calculate the western extent of the transect data
    western_extent = define_western_extent(transect_data)

    # Run the adaptive search window to identify which points have to be re-weighted to account
    # for extrapolation
    range_grid, inside_indices, outside_indices, outside_weights = adaptive_search_radius(
        distance_matrix, mesh_data, western_extent, settings_dict
    )

    # Calculate the lagged semivariogram (M20)
    local_variogram = np.apply_along_axis(
        variogram, 1, range_grid, settings_dict["variogram_parameters"]
    )

    # Append 1.0 for the ordinary kriging assumptions (M2)
    local_variogram_M2 = np.sort(
        np.hstack((local_variogram, np.ones((local_variogram.shape[0], 1))))
    )

    # Extract x- and y-coordinates, and variable_data
    # ---- x
    x_coordinates = transect_data["x"].to_numpy()
    # ---- y
    y_coordinates = transect_data["y"].to_numpy()

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

    # Ordinary kriging
    kriged_values = np.apply_along_axis(
        kriging_interpolation,
        1,
        full_stack,
        settings_dict["kriging_parameters"],
        settings_dict["variogram_parameters"],
        x_coordinates,
        y_coordinates,
        variable_data,
    )

    # Compute the coefficients of variation (CV)
    # ---- Compute the global/survey variance
    survey_variance = np.var(variable_data)
    # ---- Calculate the integrated variable when distributed over area
    # --------- Compute area
    area = settings_dict["kriging_parameters"]["A0"] * mesh_data["fraction_cell_in_polygon"]
    # -------- Drop erroneous negative values along edge
    kriged_values[kriged_values[:, 0] < 0, 0] = 0.0
    # -------- Distribute biological variable over area
    survey_estimate = np.nansum(kriged_values[:, 0] * area)
    # ---- Compute the georeferenced CV at each mesh node
    mesh_CV = (
        area.mean()
        * np.sqrt(kriged_values[:, 1].dot(survey_variance))
        / survey_estimate
        * np.sqrt(len(kriged_values[:, 1]))
    )
    # ---- Compute the global/survey CV
    survey_CV = (
        np.sqrt(np.nansum(kriged_values[:, 1] * area**2) * survey_variance) / survey_estimate
    )

    # Return a Tuple with the results
    # ---- Create DataFrame with the mesh node results
    # -------- Initialize
    mesh_results = mesh_data.copy()
    # -------- Add area
    mesh_results["area"] = area
    # -------- Add the kriged variable
    mesh_results["kriged_mean"] = kriged_values[:, 0]
    # -------- Add the kriged variance
    mesh_results["kriged_variance"] = kriged_values[:, 1]
    # -------- Add the sample variance
    mesh_results["sample_variance"] = kriged_values[:, 2]
    # -------- Add the sample CV
    mesh_results["sample_cv"] = mesh_CV
    # -------- Add the absolute kriged value (i.e. not normalized by area)
    mesh_results["biomass"] = kriged_values[:, 0] * area
    # -------- Extract only the necessary dataframe columns
    mesh_results = mesh_results.filter(regex="^(?!(fraction|x|y))")
    # ---- Create dictionary with survey-wide kriged results
    survey_results = {
        "variable": variable_name,
        "survey_mean": kriged_values[:, 0].mean(),
        "survey_estimate": survey_estimate,
        "survey_cv": survey_CV,
        "mesh_results_df": mesh_results,
    }
    # ---- Return output
    return survey_results


def kriging_interpolation(
    stacked_array: np.ndarray,
    kriging_parameters: dict,
    variogram_parameters: dict,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    variable_data: np.ndarray,
):
    """
    Interpolate data at georeferenced coordinates using ordinary kriging.

    Parameters
    ----------
    stacked_array: np.ndarray
        Horizontally stacked array that includes inside and outside weight indices, extrapolation
        weights, the local semivariogram (M2), and local range estimates.
    kriging_parameters: dict
        Kriging parameters.
    variogram_parameters: dict
        Variogram parameters.
    x_coordinates: np.array
        The x-axis coordinates
    y_coordinates: np.array
        The y-axis coordinates
    variable_data: np.ndarray
        An array of data that will be interpolated.
    """

    # Extract kriging parameter values
    # ---- Anisotropy
    anisotropy = kriging_parameters["anisotropy"]
    # ---- k_max
    k_max = kriging_parameters["kmax"]
    # ---- k_min
    k_min = kriging_parameters["kmin"]
    # ---- search radius
    search_radius = kriging_parameters["search_radius"]

    # Break up the stacked index
    # ---- Extrapolation/out-of-sample (OOS) weight
    outside_weight = stacked_array[0]
    # ---- Within-sample (IS) indices
    interp_indices = stacked_array[1 : k_max + 1]
    # -------- Drop NaN
    interp_indices = interp_indices[~np.isnan(interp_indices)].astype(int)
    # ---- OOS indices
    extrap_indices = stacked_array[k_max + 1 : (k_max + 1 + k_min)]
    # -------- Drop NaN
    extrap_indices = extrap_indices[~np.isnan(extrap_indices)].astype(int)
    # ---- Combine the IS and OOS indices
    composite = np.concatenate([interp_indices, extrap_indices])
    # ---- M2 lagged variogram
    M2_vario = stacked_array[-(k_max * 2) - 1 : -k_max]
    # -------- Drop NaN
    M2_vario = M2_vario[~np.isnan(M2_vario)]
    # ---- Range grid
    range_vals = stacked_array[-k_max:]
    # -------- Drop NaN
    range_vals = range_vals[~np.isnan(range_vals)]

    # Index the coordinates and variable data
    # ---- x
    x_indexed = x_coordinates[composite]
    # ---- y
    y_indexed = y_coordinates[composite]
    # ---- variable
    variable_indexed = variable_data[composite]

    # Set extrapolated variable values to 0.0
    variable_indexed[range_vals > search_radius] = 0.0

    # Compute the kriging covariance matrix
    kriging_covariance = kriging_matrix(x_indexed, y_indexed, variogram_parameters)

    # Compute the kriging weights (lambda)
    kriging_weights = kriging_lambda(anisotropy, M2_vario, kriging_covariance)

    # Calculate the point estimate
    point_estimate = (kriging_weights[: len(composite)] * variable_indexed).sum() * outside_weight

    # Calculate the kriged variance
    kriged_variance = kriging_weights.dot(M2_vario)

    # Calculate the sample variance and CV
    if abs(point_estimate) < np.finfo(float).eps:
        sample_variance = np.nan
    else:
        sample_variance = np.sqrt(kriged_variance * np.var(variable_indexed, ddof=1)) / np.abs(
            point_estimate
        )

    # Return output
    return np.array([point_estimate, kriged_variance, sample_variance])


def kriging_matrix(x_coordinates, y_coordinates, variogram_parameters):
    """
    Calculate the kriging covariance matrix

    Parameters
    ----------
    x_coordinates: np.array
        The x-axis coordinates
    y_coordinates: np.array
        The y-axis coordinates
    variogram_parameters: dict
        Dictionary containing variogram model parameters
    """

    # Calculate local distance matrix of within-range samples
    local_distance_matrix = griddify_lag_distances(x_coordinates, y_coordinates)

    # Calculate the covariance/kriging matrix (without the constant term)
    kriging_matrix_initial = variogram(
        distance_lags=local_distance_matrix, variogram_parameters=variogram_parameters
    )

    # Expand the covariance/kriging matrix with a constant
    # ---- In Ordinary Kriging, this should be '1'
    # Columns
    kriging_matrix = np.concatenate(
        [kriging_matrix_initial, np.ones((x_coordinates.size, 1))], axis=1
    )

    # Add column and row of ones for Ordinary Kriging
    kriging_matrix = np.concatenate([kriging_matrix, np.ones((1, x_coordinates.size + 1))], axis=0)

    # Diagonal fill (0.0)
    # ---- TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?  # noqa
    np.fill_diagonal(kriging_matrix, 0.0)

    return kriging_matrix


def search_radius_mask(distance_matrix: np.ndarray, search_radius: float):
    """
    Creates a NaN mask of values that fall beyond the search radius

    Parameters
    ----------
    distance_matrix: np.ndarray
        An array of lag distances between mesh points and every transect coordinate.
    search_radius: float
        The maximum lag distance allowed from each mesh point.
    """
    # Create copy of matrix
    matrix_copy = distance_matrix.copy()

    # Find values beyond the maximum search radius and mask (assign NaN) to those beyond
    # ---- Generate mask
    mask = matrix_copy > search_radius
    # ---- Set values outside the radius to NaN
    matrix_copy[mask] = np.nan

    # Return output
    return matrix_copy


def count_within_radius(distance_matrirx_masked: np.ndarray):
    """
    Counts the number of NaN-masked distance matrix values to determine points where extrapolation
    is required.

    Parameters
    ----------
    distance_matrix_masked: np.ndarray
        A NaN-masked array of lagged distances.
    """
    # Create copy of matrix
    matrix_copy = distance_matrirx_masked.copy()

    # Create boolean matrix of values outside search radius
    nan_mask = np.isnan(matrix_copy)

    # Count the number of values within the search radius and return output
    return np.sum(~nan_mask, axis=1)


def adaptive_search_radius(
    distance_matrix: np.ndarray,
    mesh_data: pd.DataFrame,
    western_extent: pd.DataFrame,
    settings_dict: dict,
):
    """
    Find the indices of the k-th nearest points (relative to a reference coordinate) required
    for computing the lagged semivariogram

    Parameters
    ----------
    distance_matrix: np.ndarray
        An array/matrix that includes the distances of each mesh points from every
        georeferenced along-transect interval
    mesh_data: pd.DataFrame
        Kriging mesh.
    western_extent: pd.DataFrame
        Coordinates of the western extent of transect lines.
    settings_dict:
        Dictionary that contains all of the analysis settings that detail specific algorithm
        arguments and user-defined inputs.
    """

    # Extract key search radius parameters
    # ---- k_min
    k_min = settings_dict["kriging_parameters"]["kmin"]
    # ---- k_max
    k_max = settings_dict["kriging_parameters"]["kmax"]
    # ---- Search radius (distance)
    search_radius = settings_dict["kriging_parameters"]["search_radius"]

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
        # Extract the indices requiring extrapolation
        nearby_indices = local_indices[sparse_radii][:, :k_min]
        # Index the mesh grid coordinates for bounding the search radius expansion/extrapolation
        # ---- y-coordinates with array transformation to access matrix operations
        mesh_y = mesh_data["y"].to_numpy()[sparse_radii].reshape(-1, 1)
        # ---- x-coordinates
        mesh_x = mesh_data["x"].to_numpy()[sparse_radii]

        # Update local points
        # ---- Fill NaN values
        local_points[sparse_radii, k_min:] = np.nan
        wr_indices[sparse_radii, k_min:] = np.nan

        # Calculate the mesh distance from the western boundary of the survey transects
        # ---- Find closest point
        mesh_western_distance = np.abs(mesh_y - western_extent["y"].to_numpy()).argmin(axis=1)
        # ---- Calculate the western limits (x-axis)
        western_limit = western_extent.iloc[np.ravel(mesh_western_distance)]["x"]
        # ---- Compute bounding threshold (for tapered extrapolation function)
        western_threshold = western_limit - search_radius
        # ---- Create a thresholded mask for lazy operations
        western_limit_mask = mesh_x < western_threshold

        # Adjust values that don't fall outside the western extent
        if np.any(~western_limit_mask):
            # ---- Grab all values that don't fall outside the western extent
            soft_extrapolation_index = sparse_radii[~western_limit_mask]

            # Find the local points where there are at least some valid points
            partial_indices = soft_extrapolation_index[
                valid_distances[soft_extrapolation_index] > 0
            ]
            # ---- Update the current values in `wr_indices`
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

            # Find the local points where there are no valid points within the search radius
            full_extrap_indices = soft_extrapolation_index[
                valid_distances[soft_extrapolation_index] == 0
            ]
            # ---- Update `oos_indices`
            oos_indices[full_extrap_indices] = wr_indices[full_extrap_indices, :k_min]
            # ----- Update `wr_indices`
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

        # Alert message (if verbose = True)
        if settings_dict["verbose"]:
            print(
                f"""Extrapolation applied to kriging mesh points ({len(sparse_radii)} of """
                f"""{wr_indices.shape[0]}):
            * {len(valid_distances[valid_distances == 0])} points had 0 valid range estimates"""
                f""" without extrapolation
            * {len(valid_distances[(valid_distances != 0) & (valid_distances < k_min)])} """
                f"""points had at least 1 valid point but fewer than {k_min} valid neighbors"""
            )

    # Return output
    return local_points[:, :k_max], wr_indices, oos_indices, oos_weights


def kriging_lambda(
    anisotropy: float,
    lagged_semivariogram: np.ndarray,
    kriging_matrix_input: np.ndarray,
):
    """
    Apply singular value decomposition (SVD) to compute kriging (lambda) weights

    Parameters
    ----------
    anisotropy: np.float64
        Anisotropy ratio.
    lagged_semivariogram: np.array
        Lagged semivariogram
    kriging_matrix_input: np.array
        Kriging matrix.
    """
    # Singular value decomposition (SVD)
    # ---- U: left singular vectors (directions of maximum variance)
    # ---- Sigma: singular values (amount of variance captured by each singular vector, U)
    # ---- VH: conjugate transpose of the right singular vectors
    U, Sigma, VH = np.linalg.svd(kriging_matrix_input, full_matrices=True)

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
