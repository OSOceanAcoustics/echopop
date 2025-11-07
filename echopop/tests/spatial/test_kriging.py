import numpy as np
import pandas as pd

from echopop.geostatistics.kriging import (
    Kriging,
    adaptive_search_radius,
    count_within_radius,
    krige,
    kriging_lambda,
    ordinary_kriging_matrix,
    parse_stacked_kriging_array,
    project_kriging_results,
    search_radius_mask,
    uniform_search_strategy,
)


# ==================================================================================================
# Test search_radius_mask
# -----------------------
def test_search_radius_mask_basic():
    """Test basic search radius mask functionality."""
    distance_matrix = np.array(
        [[0.0, 1.0, 2.0, 3.0], [1.0, 0.0, 1.0, 2.0], [2.0, 1.0, 0.0, 1.0], [3.0, 2.0, 1.0, 0.0]]
    )
    search_radius = 1.5

    result = search_radius_mask(distance_matrix, search_radius)

    # Check output type and shape
    assert isinstance(result, np.ndarray)
    assert result.shape == distance_matrix.shape

    # Check that distances within radius are preserved and beyond radius are NaN
    expected = distance_matrix.copy()
    expected[distance_matrix > search_radius] = np.nan
    np.testing.assert_array_equal(result, expected)


def test_search_radius_mask_zero_radius():
    """Test search radius mask with zero radius."""
    distance_matrix = np.array([[0.0, 1.0], [1.0, 0.0]])
    search_radius = 0.0

    result = search_radius_mask(distance_matrix, search_radius)

    # Only diagonal elements should be preserved, off-diagonal should be NaN
    expected = distance_matrix.copy()
    expected[distance_matrix > search_radius] = np.nan
    np.testing.assert_array_equal(result, expected)


# ==================================================================================================
# Test count_within_radius
# ------------------------
def test_count_within_radius_basic():
    """Test basic count within radius functionality."""
    # Create a masked distance matrix with NaN for out-of-range distances
    distance_matrix_masked = np.array(
        [
            [0.0, 1.0, np.nan, np.nan],
            [1.0, 0.0, 1.0, np.nan],
            [np.nan, 1.0, 0.0, 1.0],
            [np.nan, np.nan, 1.0, 0.0],
        ]
    )

    result = count_within_radius(distance_matrix_masked)

    # Check output type
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.int64 or result.dtype == np.int32

    # Check expected counts per row
    expected = np.array([2, 3, 3, 2])
    np.testing.assert_array_equal(result, expected)


def test_count_within_radius_empty():
    """Test count within radius with no valid points."""
    distance_matrix_masked = np.array([[np.nan, np.nan], [np.nan, np.nan]])

    result = count_within_radius(distance_matrix_masked)

    expected = np.array([0, 0])
    np.testing.assert_array_equal(result, expected)


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
# Test adaptive_search_radius
# ---------------------------
def test_adaptive_search_radius_basic():
    """Test adaptive search radius functionality."""
    distance_matrix = np.array(
        [[0.0, 1.0, 2.0, 3.0], [1.0, 0.0, 1.5, 2.5], [2.0, 1.5, 0.0, 1.0], [3.0, 2.5, 1.0, 0.0]]
    )
    k_min = 2
    k_max = 3
    search_radius = 2.0

    result = adaptive_search_radius(distance_matrix, k_min, k_max, search_radius)

    # Check that result is a tuple with 4 elements
    assert isinstance(result, tuple)
    assert len(result) == 4

    # Check that all elements are numpy arrays
    for element in result:
        assert isinstance(element, np.ndarray)


def test_adaptive_search_radius_with_strategy(
    sample_sparse_distance_matrix,
    sample_kriging_settings,
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
    sample_kriging_settings,
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
# ----------------------------
def test_ordinary_kriging_matrix_basic():
    """Test ordinary kriging matrix construction."""
    local_distance_matrix = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.5], [2.0, 1.5, 0.0]])
    variogram_parameters = {
        "model": "exponential",
        "nugget": 0.1,
        "sill": 1.0,
        "correlation_range": 0.5,
    }

    result = ordinary_kriging_matrix(local_distance_matrix, variogram_parameters)

    # Check output is numpy array
    assert isinstance(result, np.ndarray)

    # Check that it's a square matrix with one additional row/column for Lagrange multiplier
    expected_size = local_distance_matrix.shape[0] + 1
    assert result.shape == (expected_size, expected_size)


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
# -------------------
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
# --------------------------------
def test_parse_stacked_kriging_array_basic():
    """Test parsing of stacked kriging array."""
    stacked_array = np.array(
        [[1.0, 2.0, 3.0, 4.0, 5.0], [1.5, 2.5, 3.5, 4.5, 5.5], [2.0, 3.0, 4.0, 5.0, 6.0]]
    )
    k_min = 2
    k_max = 4

    result = parse_stacked_kriging_array(stacked_array, k_min, k_max)

    # Check that result is a tuple with 4 elements
    assert isinstance(result, tuple)
    assert len(result) == 4

    # Check that all elements are numpy arrays
    for element in result:
        assert isinstance(element, np.ndarray)


# ==================================================================================================
# Test kriging_point_estimator
# ----------------------------
def test_kriging_point_estimator_basic(sample_biomass_data, sample_variogram_parameters):
    """Test basic functionality of kriging_point_estimator via krige function."""
    # Test kriging_point_estimator indirectly through krige function
    # since creating a proper mock kriging_array is complex

    # Create a simple mesh
    mesh = pd.DataFrame({"x": [1.0, 2.0], "y": [1.0, 2.0]})

    # Create kriging parameters
    kriging_parameters = {"k_min": 2, "k_max": 3, "search_radius": 5.0, "aspect_ratio": 0.01}

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


# ==================================================================================================
# Test krige
# ----------
def test_krige_basic(sample_transect_df):
    """Test basic kriging functionality."""
    # Create kriging mesh
    kriging_mesh = pd.DataFrame(
        {
            "x": np.linspace(0, 2, 5),
            "y": np.linspace(0, 2, 5),
            "area": np.ones(5) * 1.0,  # Required area column
        }
    )

    coordinate_names = ("x", "y")
    variable = "biomass_density"
    kriging_parameters = {
        "k_max": 5,
        "k_min": 2,
        "search_radius": 2.0,
        "aspect_ratio": 0.001,  # Missing parameter
    }
    variogram_parameters = {
        "model": "exponential",
        "nugget": 0.1,
        "sill": 1.0,
        "correlation_range": 0.5,
    }

    result = krige(
        sample_transect_df,
        kriging_mesh,
        coordinate_names,
        variable,
        kriging_parameters,
        variogram_parameters,
    )

    # Check output type
    assert isinstance(result, np.ndarray)

    # Check that we get one estimate per mesh point
    assert len(result) == len(kriging_mesh)


# ==================================================================================================
# Test project_kriging_results
# ----------------------------
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
# Test Kriging class
# ------------------
def test_kriging_class_initialization():
    """Test Kriging class initialization."""
    mesh = pd.DataFrame(
        {
            "x": [0, 1, 2],
            "y": [0, 1, 2],
            "area": [1.0, 1.0, 1.0],  # Required area column
        }
    )
    kriging_params = {
        "k_max": 5,
        "k_min": 2,
        "search_radius": 2.0,
        "aspect_ratio": 0.001,  # Required parameter
    }
    variogram_params = {
        "model": "exponential",
        "nugget": 0.1,
        "sill": 1.0,
        "correlation_range": 0.5,
    }
    coordinate_names = ("x", "y")

    krg = Kriging(
        mesh=mesh,
        kriging_params=kriging_params,
        variogram_params=variogram_params,
        coordinate_names=coordinate_names,
    )

    # Check that attributes are set correctly
    assert np.allclose(krg.mesh, mesh)
    assert krg.kriging_params == kriging_params
    assert krg.variogram_params == variogram_params
    assert krg.coordinate_names == coordinate_names


def test_kriging_class_register_search_strategy():
    """Test registering custom search strategy."""
    mesh = pd.DataFrame({"x": [0, 1], "y": [0, 1], "area": [1.5, 1.5]})  # Required area column
    kriging_params = {
        "k_max": 5,
        "k_min": 2,
        "search_radius": 2.0,  # Required parameter
        "aspect_ratio": 0.001,  # Required parameter
    }
    variogram_params = {
        "model": "exponential",
        "nugget": 0.1,
        "sill": 1.0,
        "correlation_range": 0.5,
    }

    krg = Kriging(mesh, kriging_params, variogram_params, ("x", "y"))

    # Define a dummy search strategy
    def dummy_strategy(**kwargs):
        return np.array([0]), np.array([0]), np.array([1.0])

    # Register the strategy
    krg.register_search_strategy("dummy", dummy_strategy)

    # Check that it was registered
    strategies = krg.list_search_strategies()
    assert "dummy" in strategies


def test_kriging_class_crop_mesh():
    """Test mesh cropping functionality."""
    mesh = pd.DataFrame(
        {
            "x": np.linspace(0, 5, 20),
            "y": np.linspace(0, 5, 20),
            "area": np.ones(20) * 2.5,  # Required area column
        }
    )
    kriging_params = {
        "k_max": 5,
        "k_min": 2,
        "search_radius": 3.0,  # Required parameter
        "aspect_ratio": 0.001,  # Required parameter
    }
    variogram_params = {
        "model": "exponential",
        "nugget": 0.1,
        "sill": 1.0,
        "correlation_range": 0.5,
    }

    krg = Kriging(mesh, kriging_params, variogram_params, ("x", "y"))

    # Define a simple crop function
    def simple_crop(mesh, **kwargs):
        return mesh[mesh["x"] < 3.0]

    # Apply cropping
    krg.crop_mesh(simple_crop, coordinate_names=("x", "y"))

    # Check that mesh was cropped
    assert len(krg.mesh) == len(mesh)
    assert len(krg.mesh_cropped) < len(mesh)
    assert krg.mesh_cropped["x"].max() < 3.0


def test_kriging_class_krige_method(sample_transect_df):
    """Test Kriging class krige method."""
    mesh = pd.DataFrame(
        {
            "x": np.linspace(0, 2, 5),
            "y": np.linspace(0, 2, 5),
            "area": np.ones(5) * 1.0,  # Required area column
        }
    )
    kriging_params = {
        "k_max": 5,
        "k_min": 2,
        "search_radius": 2.0,
        "aspect_ratio": 0.001,  # Required parameter
    }
    variogram_params = {
        "model": "exponential",
        "nugget": 0.1,
        "sill": 1.0,
        "correlation_range": 0.5,
    }

    krg = Kriging(mesh, kriging_params, variogram_params, ("x", "y"))

    result = krg.krige(
        transects=sample_transect_df, variable="biomass_density", default_mesh_cell_area=1.0
    )

    # Check output type
    assert isinstance(result, pd.DataFrame)

    # Check that kriged values are included
    assert "biomass_density" in result.columns


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
