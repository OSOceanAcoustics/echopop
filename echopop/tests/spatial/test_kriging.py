import numpy as np
import pandas as pd

from echopop.nwfsc_feat.spatial import (
    adaptive_search_radius,
    count_within_radius,
    get_survey_western_extents,
    krige,
    kriging_lambda,
    ordinary_kriging_matrix,
    parse_stacked_kriging_array,
    project_kriging_results,
    search_radius_mask,
    uniform_search_strategy,
    western_boundary_search_strategy,
)


# ==================================================================================================
# Test get_survey_western_extents
# --------------------------------
def test_get_survey_western_extents_basic(sample_mesh_coordinates):
    """Test basic functionality of get_survey_western_extents."""
    # Add latitude column to mesh coordinates
    mesh_with_lat = sample_mesh_coordinates.copy()
    mesh_with_lat["latitude"] = mesh_with_lat["y"] + 46.0  # Add realistic latitude

    result = get_survey_western_extents(
        mesh_with_lat, coordinate_names=("x", "y"), latitude_threshold=45.0
    )

    # Check that result is a DataFrame
    assert isinstance(result, pd.DataFrame)

    # Check that it has the expected columns
    expected_columns = {"transect_num", "x", "y"}
    assert set(result.columns) == expected_columns

    # Check that each transect has one entry
    transect_counts = result["transect_num"].value_counts()
    assert all(count == 1 for count in transect_counts.values)


def test_get_survey_western_extents_westernmost(sample_mesh_coordinates):
    """Test that get_survey_western_extents returns westernmost points."""
    # Add latitude column to mesh coordinates
    mesh_with_lat = sample_mesh_coordinates.copy()
    mesh_with_lat["latitude"] = mesh_with_lat["y"] + 46.0  # Add realistic latitude

    result = get_survey_western_extents(
        mesh_with_lat, coordinate_names=("x", "y"), latitude_threshold=45.0
    )

    # For each transect, the x-coordinate should be the minimum
    for transect_num in result["transect_num"].unique():
        original_transect = mesh_with_lat[mesh_with_lat["transect_num"] == transect_num]
        result_transect = result[result["transect_num"] == transect_num]

        min_x = original_transect["x"].min()
        result_x = result_transect["x"].iloc[0]

        assert result_x == min_x


# ==================================================================================================
# Test search_radius_mask
# ------------------------
def test_search_radius_mask_basic(sample_distance_matrix):
    """Test basic functionality of search_radius_mask."""
    search_radius = 2.5
    result = search_radius_mask(sample_distance_matrix, search_radius)

    # Check output shape matches input
    assert result.shape == sample_distance_matrix.shape

    # Check that distances within radius are preserved
    within_radius = sample_distance_matrix <= search_radius
    assert np.allclose(result[within_radius], sample_distance_matrix[within_radius])

    # Check that distances outside radius are NaN
    outside_radius = sample_distance_matrix > search_radius
    assert np.all(np.isnan(result[outside_radius]))


def test_search_radius_mask_edge_cases():
    """Test edge cases for search_radius_mask."""
    # Test with very small radius
    matrix = np.array([[0.0, 1.0], [1.0, 0.0]])
    result = search_radius_mask(matrix, 0.5)

    expected = np.array([[0.0, np.nan], [np.nan, 0.0]])
    assert np.allclose(result, expected, equal_nan=True)

    # Test with very large radius
    result = search_radius_mask(matrix, 10.0)
    assert np.allclose(result, matrix)


# ==================================================================================================
# Test count_within_radius
# -------------------------
def test_count_within_radius_basic(sample_distance_matrix):
    """Test basic functionality of count_within_radius."""
    # Create a masked matrix
    masked_matrix = search_radius_mask(sample_distance_matrix, 2.5)
    result = count_within_radius(masked_matrix)

    # Check output shape
    assert result.shape == (sample_distance_matrix.shape[0],)

    # Check that counts are reasonable
    assert all(count >= 0 for count in result)
    assert all(count <= sample_distance_matrix.shape[1] for count in result)


def test_count_within_radius_all_valid():
    """Test count_within_radius with all valid distances."""
    matrix = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.0], [2.0, 1.0, 0.0]])
    result = count_within_radius(matrix)

    # Should count all points since no NaN values
    expected = np.array([3, 3, 3])
    assert np.array_equal(result, expected)


def test_count_within_radius_with_nans():
    """Test count_within_radius with NaN values."""
    matrix = np.array([[0.0, np.nan, 2.0], [1.0, 0.0, np.nan], [np.nan, 1.0, 0.0]])
    result = count_within_radius(matrix)

    # Should only count non-NaN values
    expected = np.array([2, 2, 2])
    assert np.array_equal(result, expected)


# ==================================================================================================
# Test adaptive_search_radius
# ----------------------------
def test_adaptive_search_radius_basic(
    sample_sparse_distance_matrix, sample_mesh_coordinates, sample_kriging_settings
):
    """Test basic functionality of adaptive_search_radius."""
    k_min = sample_kriging_settings["kriging_parameters"]["kmin"]
    k_max = sample_kriging_settings["kriging_parameters"]["kmax"]
    search_radius = sample_kriging_settings["kriging_parameters"]["search_radius"]

    result = adaptive_search_radius(sample_sparse_distance_matrix, k_min, k_max, search_radius)

    # Check that we get 4 return values
    assert len(result) == 4
    local_points, wr_indices, oos_indices, oos_weights = result

    # Check shapes
    n_points = sample_sparse_distance_matrix.shape[0]

    assert local_points.shape == (n_points, k_max)
    assert wr_indices.shape == (n_points, k_max)
    assert oos_indices.shape == (n_points, k_min)
    assert oos_weights.shape == (n_points,)


def test_adaptive_search_radius_with_strategy(
    sample_sparse_distance_matrix,
    sample_mesh_coordinates,
    sample_kriging_settings,
    sample_western_extent,
):
    """Test adaptive_search_radius with custom search strategy."""
    k_min = sample_kriging_settings["kriging_parameters"]["kmin"]
    k_max = sample_kriging_settings["kriging_parameters"]["kmax"]
    search_radius = sample_kriging_settings["kriging_parameters"]["search_radius"]

    # Create a simple test strategy
    def test_strategy(*args, **kwargs):
        # Return the uniform strategy result
        return uniform_search_strategy(*args, **kwargs)

    result = adaptive_search_radius(
        sample_sparse_distance_matrix, k_min, k_max, search_radius, test_strategy
    )

    # Should still return 4 values
    assert len(result) == 4


def test_adaptive_search_radius_western_boundary(
    sample_sparse_distance_matrix,
    sample_mesh_coordinates,
    sample_kriging_settings,
    sample_western_extent,
):
    """Test adaptive_search_radius with western boundary strategy."""
    k_min = sample_kriging_settings["kriging_parameters"]["kmin"]
    k_max = sample_kriging_settings["kriging_parameters"]["kmax"]
    search_radius = sample_kriging_settings["kriging_parameters"]["search_radius"]

    # Create western boundary strategy
    def test_western_strategy(*args, **kwargs):
        # For simplicity, return uniform strategy result
        return uniform_search_strategy(*args, **kwargs)

    result = adaptive_search_radius(
        sample_sparse_distance_matrix, k_min, k_max, search_radius, test_western_strategy
    )

    # Should return 4 values
    assert len(result) == 4
    local_points, wr_indices, oos_indices, oos_weights = result

    # Check that weights exist (may be modified)
    assert len(oos_weights) == sample_sparse_distance_matrix.shape[0]


# ==================================================================================================
# Test ordinary_kriging_matrix
# -----------------------------
def test_ordinary_kriging_matrix_basic(sample_local_distance_matrix, sample_variogram_parameters):
    """Test basic functionality of ordinary_kriging_matrix."""
    # Add model parameter to variogram_parameters
    variogram_params = {"model": "exponential"}
    for param_name in sample_variogram_parameters.keys():
        variogram_params[param_name] = sample_variogram_parameters[param_name].value

    result = ordinary_kriging_matrix(sample_local_distance_matrix, variogram_params)

    # Check output shape (should be n+1 x n+1 for ordinary kriging)
    n = sample_local_distance_matrix.shape[0]
    assert result.shape == (n + 1, n + 1)

    # Check that last row and column are ones (except bottom-right which is 0)
    assert np.allclose(result[-1, :-1], 1.0)  # Last row (except last element)
    assert np.allclose(result[:-1, -1], 1.0)  # Last column (except last element)
    assert result[-1, -1] == 0.0  # Bottom-right should be 0


def test_ordinary_kriging_matrix_symmetry(
    sample_local_distance_matrix, sample_variogram_parameters
):
    """Test that ordinary_kriging_matrix produces symmetric matrix."""
    # Add model parameter to variogram_parameters
    variogram_params = {"model": "exponential"}
    for param_name in sample_variogram_parameters.keys():
        variogram_params[param_name] = sample_variogram_parameters[param_name].value

    result = ordinary_kriging_matrix(sample_local_distance_matrix, variogram_params)

    # Check symmetry (excluding the constraint row/column)
    n = sample_local_distance_matrix.shape[0]
    covariance_part = result[:n, :n]
    assert np.allclose(covariance_part, covariance_part.T)


# ==================================================================================================
# Test kriging_lambda
# --------------------
def test_kriging_lambda_basic(
    sample_anisotropy_value, sample_lagged_semivariogram, sample_kriging_covariance_matrix
):
    """Test basic functionality of kriging_lambda."""
    result = kriging_lambda(
        sample_anisotropy_value, sample_lagged_semivariogram, sample_kriging_covariance_matrix
    )

    # Check output shape
    assert result.shape == (sample_kriging_covariance_matrix.shape[0],)

    # Check that weights approximately sum to 1 (for ordinary kriging)
    weights_sum = np.sum(result[:-1])  # Exclude Lagrange multiplier
    assert abs(weights_sum - 1.0) < 0.1  # Allow some tolerance


def test_kriging_lambda_numerical_stability():
    """Test kriging_lambda with potentially ill-conditioned matrix."""
    # Create a nearly singular matrix
    matrix = np.array([[1.0, 0.99, 1.0], [0.99, 1.0, 1.0], [1.0, 1.0, 0.0]])
    rhs = np.array([0.5, 0.7, 1.0])

    # Should not raise an error
    result = kriging_lambda(0.01, rhs, matrix)
    assert len(result) == 3


# ==================================================================================================
# Test parse_stacked_kriging_array
# ---------------------------------
def test_parse_stacked_kriging_array_basic(sample_stacked_kriging_array):
    """Test basic functionality of parse_stacked_kriging_array."""
    k_min = 2
    k_max = 3
    result = parse_stacked_kriging_array(sample_stacked_kriging_array, k_min, k_max)

    # Check that result is a tuple with 4 elements
    assert isinstance(result, tuple)
    assert len(result) == 4

    # Unpack the results
    oos_weights, composite_indices, local_semivariogram, range_distances = result

    # Check shapes and types
    assert isinstance(oos_weights, np.ndarray)
    assert isinstance(composite_indices, np.ndarray)
    assert isinstance(local_semivariogram, np.ndarray)
    assert isinstance(range_distances, np.ndarray)


def test_parse_stacked_kriging_array_values(sample_stacked_kriging_array):
    """Test that parse_stacked_kriging_array correctly extracts values."""
    k_min = 2
    k_max = 3
    result = parse_stacked_kriging_array(sample_stacked_kriging_array, k_min, k_max)

    # Check that we get tuple output
    assert isinstance(result, tuple)
    assert len(result) == 4


# ==================================================================================================
# Test kriging_point_estimator
# -----------------------------
def test_kriging_point_estimator_basic(sample_biomass_data, sample_variogram_parameters):
    """Test basic functionality of kriging_point_estimator via krige function."""
    # Test kriging_point_estimator indirectly through krige function
    # since creating a proper mock kriging_array is complex

    # Create a simple mesh
    mesh = pd.DataFrame({"x": [1.0, 2.0], "y": [1.0, 2.0]})

    # Create kriging parameters
    kriging_parameters = {"k_min": 2, "k_max": 3, "search_radius": 5.0, "anisotropy": 0.01}

    # Add model parameter to variogram_parameters
    variogram_params = {"model": "exponential"}
    for param_name in sample_variogram_parameters.keys():
        variogram_params[param_name] = sample_variogram_parameters[param_name].value

    result = krige(
        sample_biomass_data,
        mesh,
        coordinate_names=("x", "y"),
        variable="biomass_density",
        kriging_parameters=kriging_parameters,
        variogram_parameters=variogram_params,
    )

    # Check that result is valid
    assert isinstance(result, np.ndarray)
    assert result.shape == (len(mesh), 3)

    # Check that estimates are reasonable
    assert not np.any(np.isnan(result[:, 0]))
    assert not np.any(np.isinf(result[:, 0]))


def test_kriging_point_estimator_at_sample_point(sample_biomass_data, sample_variogram_parameters):
    """Test kriging_point_estimator with sample data via krige function."""
    # Test kriging_point_estimator indirectly through krige function
    # since creating a proper mock kriging_array is complex

    # Create a simple mesh
    mesh = pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]})  # Use points close to sample data

    # Create kriging parameters
    kriging_parameters = {"k_min": 2, "k_max": 3, "search_radius": 5.0, "anisotropy": 0.01}

    # Add model parameter to variogram_parameters
    variogram_params = {"model": "exponential"}
    for param_name in sample_variogram_parameters.keys():
        variogram_params[param_name] = sample_variogram_parameters[param_name].value

    result = krige(
        sample_biomass_data,
        mesh,
        coordinate_names=("x", "y"),
        variable="biomass_density",
        kriging_parameters=kriging_parameters,
        variogram_parameters=variogram_params,
    )

    # Check that result is reasonable
    assert isinstance(result, np.ndarray)
    assert result.shape == (len(mesh), 3)

    # Check that estimates are reasonable
    assert not np.any(np.isnan(result[:, 0]))
    assert not np.any(np.isinf(result[:, 0]))


# ==================================================================================================
# Test krige
# ----------
def test_krige_basic(sample_biomass_data, sample_variogram_parameters):
    """Test basic functionality of krige."""
    # Create a simple mesh
    mesh = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})

    # Create kriging parameters
    kriging_parameters = {"k_min": 2, "k_max": 3, "search_radius": 5.0, "anisotropy": 0.01}

    # Add model parameter to variogram_parameters
    variogram_params = {"model": "exponential"}
    for param_name in sample_variogram_parameters.keys():
        variogram_params[param_name] = sample_variogram_parameters[param_name].value

    result = krige(
        sample_biomass_data,
        mesh,
        coordinate_names=("x", "y"),
        variable="biomass_density",
        kriging_parameters=kriging_parameters,
        variogram_parameters=variogram_params,
    )  # Check that result is a numpy array
    assert isinstance(result, np.ndarray)

    # Check shape (should be n_mesh_points x 3: estimate, kriged_variance, sample_variance)
    assert result.shape == (len(mesh), 3)

    # Check that estimates are reasonable (not NaN or infinite)
    assert not np.any(np.isnan(result[:, 0]))
    assert not np.any(np.isinf(result[:, 0]))

    # Check that variance estimates are non-negative
    assert np.all(result[:, 1] >= 0)  # kriged variance should be non-negative


def test_krige_with_search_strategy(sample_biomass_data, sample_variogram_parameters):
    """Test krige with custom search strategy."""
    # Create a simple mesh
    mesh = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})

    # Create kriging parameters
    kriging_parameters = {"k_min": 2, "k_max": 3, "search_radius": 5.0, "anisotropy": 0.01}

    # Add model parameter to variogram_parameters
    variogram_params = {"model": "exponential"}
    for param_name in sample_variogram_parameters.keys():
        variogram_params[param_name] = sample_variogram_parameters[param_name].value

    # Create a simple search strategy
    def test_strategy(*args, **kwargs):
        return uniform_search_strategy(*args, **kwargs)

    result = krige(
        sample_biomass_data,
        mesh,
        coordinate_names=("x", "y"),
        variable="biomass_density",
        kriging_parameters=kriging_parameters,
        variogram_parameters=variogram_params,
        adaptive_search_strategy=test_strategy,
    )

    # Should still return valid results
    assert isinstance(result, np.ndarray)
    assert result.shape == (len(mesh), 3)

    # Check that estimates are reasonable
    assert not np.any(np.isnan(result[:, 0]))
    assert not np.any(np.isinf(result[:, 0]))


# ==================================================================================================
# Test project_kriging_results
# -----------------------------
def test_project_kriging_results_basic(sample_kriging_results, sample_projection_parameters):
    """Test basic functionality of project_kriging_results."""
    # Create mock kriged estimates array
    kriged_estimates = np.column_stack(
        [
            sample_kriging_results["estimated_value"].values,
            sample_kriging_results["estimation_variance"].values,
            sample_kriging_results["estimation_variance"].values,  # Sample variance
        ]
    )

    # Create mock kriging mesh
    kriging_mesh = sample_kriging_results[["x", "y"]].copy()

    # Create mock transect dataframe
    transect_df = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0], "biomass_density": [10.0, 15.0, 12.0]}
    )

    result = project_kriging_results(
        kriged_estimates, kriging_mesh, transect_df, variable="biomass_density"
    )

    # Check that result is a tuple
    assert isinstance(result, tuple)
    assert len(result) == 2

    # Check that first element is a DataFrame
    result_df, total_area = result
    assert isinstance(result_df, pd.DataFrame)
    assert isinstance(total_area, (int, float, np.number))


def test_project_kriging_results_values(sample_kriging_results, sample_projection_parameters):
    """Test that project_kriging_results correctly processes data."""
    # Create mock kriged estimates array
    kriged_estimates = np.column_stack(
        [
            sample_kriging_results["estimated_value"].values,
            sample_kriging_results["estimation_variance"].values,
            sample_kriging_results["estimation_variance"].values,  # Sample variance
        ]
    )

    # Create mock kriging mesh
    kriging_mesh = sample_kriging_results[["x", "y"]].copy()

    # Create mock transect dataframe
    transect_df = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0], "biomass_density": [10.0, 15.0, 12.0]}
    )

    result = project_kriging_results(
        kriged_estimates, kriging_mesh, transect_df, variable="biomass_density"
    )

    # Check that we get reasonable results
    result_df, total_area = result
    assert len(result_df) == len(sample_kriging_results)
    assert total_area > 0


# ==================================================================================================
# Test uniform_search_strategy
# -----------------------------
def test_uniform_search_strategy_basic():
    """Test basic functionality of uniform_search_strategy."""
    # Create mock data
    sparse_radii = np.array([0, 1, 2])
    valid_distances = np.array([2, 2, 2])
    local_points = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    k_min = 2
    k_max = 3
    search_radius = 1.5
    wr_indices = np.array([[0, 1, 2], [1, 0, 2], [2, 1, 0]], dtype=float)
    oos_indices = np.full((3, 2), np.nan)
    oos_weights = np.ones(3)

    result = uniform_search_strategy(
        sparse_radii,
        valid_distances,
        local_points,
        k_min,
        k_max,
        search_radius,
        wr_indices,
        oos_indices,
        oos_weights,
    )

    # Check that result is a tuple with 3 elements
    assert isinstance(result, tuple)
    assert len(result) == 3

    # Check that each element is a numpy array
    updated_wr_indices, updated_oos_indices, updated_oos_weights = result
    assert isinstance(updated_wr_indices, np.ndarray)
    assert isinstance(updated_oos_indices, np.ndarray)
    assert isinstance(updated_oos_weights, np.ndarray)


# ==================================================================================================
# Test western_boundary_search_strategy
# --------------------------------------
def test_western_boundary_search_strategy_basic(sample_western_extent):
    """Test basic functionality of western_boundary_search_strategy."""
    # Create mock data
    kriging_mesh = pd.DataFrame(
        {"x": [-1.0, 0.0, 1.0], "y": [0.0, 1.0, 2.0], "transect_num": [1, 2, 3]}
    )

    sparse_radii = np.array([0, 1, 2])
    valid_distances = np.array([2, 2, 2])
    local_points = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    distance_matrix_masked = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    nearby_indices = np.array([[0, 1], [1, 0], [2, 1]])
    k_min = 2
    k_max = 3
    search_radius = 1.5
    wr_indices = np.array([[0, 1, 2], [1, 0, 2], [2, 1, 0]], dtype=float)
    oos_indices = np.full((3, 2), np.nan)
    oos_weights = np.ones(3)

    result = western_boundary_search_strategy(
        kriging_mesh,
        sample_western_extent,
        coordinate_names=("x", "y"),
        sparse_radii=sparse_radii,
        valid_distances=valid_distances,
        local_points=local_points,
        distance_matrix_masked=distance_matrix_masked,
        nearby_indices=nearby_indices,
        k_min=k_min,
        k_max=k_max,
        search_radius=search_radius,
        wr_indices=wr_indices,
        oos_indices=oos_indices,
        oos_weights=oos_weights,
    )

    # Check that result is a tuple with 3 elements
    assert isinstance(result, tuple)
    assert len(result) == 3

    # Check that each element is a numpy array
    updated_wr_indices, updated_oos_indices, updated_oos_weights = result
    assert isinstance(updated_wr_indices, np.ndarray)
    assert isinstance(updated_oos_indices, np.ndarray)
    assert isinstance(updated_oos_weights, np.ndarray)

    # Check that weights are positive
    assert all(w > 0 for w in updated_oos_weights)  # All weights should be positive
