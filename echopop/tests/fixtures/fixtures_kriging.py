import numpy as np
import pandas as pd
import pytest
from lmfit import Parameters


# ==================================================================================================
# Fixtures for kriging testing
# ----------------------------
@pytest.fixture
def sample_distance_matrix():
    """Create a sample distance matrix for kriging tests."""
    return np.array(
        [
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [1.0, 0.0, 1.0, 2.0, 3.0],
            [2.0, 1.0, 0.0, 1.0, 2.0],
            [3.0, 2.0, 1.0, 0.0, 1.0],
            [4.0, 3.0, 2.0, 1.0, 0.0],
        ]
    )


@pytest.fixture
def sample_mesh_data():
    """Create sample mesh data for kriging tests."""
    return pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            "y": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            "transect_num": [1, 1, 2, 2, 3, 3],
        }
    )


@pytest.fixture
def sample_kriging_settings():
    """Create sample kriging settings dictionary."""
    return {"kriging_parameters": {"kmin": 3, "kmax": 8, "search_radius": 5.0}, "verbose": False}


@pytest.fixture
def sample_variogram_parameters():
    """Create sample variogram parameters for kriging."""
    variogram_parameters = Parameters()
    variogram_parameters.add_many(
        ("nugget", 0.1, True, 0.0, None),
        ("sill", 1.0, True, 0.0, None),
        ("correlation_range", 2.0, True, 0.0, None),
        ("hole_effect_range", 0.0, True, 0.0, None),
        ("decay_power", 1.5, True, 0.0, None),
    )
    return variogram_parameters


@pytest.fixture
def sample_local_distance_matrix():
    """Create a sample local distance matrix for kriging matrix tests."""
    return np.array([[0.0, 1.2, 2.1], [1.2, 0.0, 1.5], [2.1, 1.5, 0.0]])


@pytest.fixture
def sample_lagged_semivariogram():
    """Create sample lagged semivariogram for kriging tests."""
    return np.array([0.3, 0.7, 0.9, 1.0])


@pytest.fixture
def sample_kriging_covariance_matrix():
    """Create sample kriging covariance matrix."""
    return np.array(
        [[1.0, 0.7, 0.3, 1.0], [0.7, 1.0, 0.5, 1.0], [0.3, 0.5, 1.0, 1.0], [1.0, 1.0, 1.0, 0.0]]
    )


@pytest.fixture
def sample_sparse_distance_matrix():
    """Create a distance matrix with sparse areas (for testing adaptive search)."""
    # Create a matrix where some areas have very few nearby points
    matrix = np.full((10, 20), 10.0)  # Fill with large distances

    # Add some closer points for the first few rows
    matrix[0, :3] = [0.0, 1.0, 2.0]
    matrix[1, :3] = [1.0, 0.0, 1.0]
    matrix[2, :3] = [2.0, 1.0, 0.0]

    # Make some rows very sparse (only 1-2 nearby points)
    matrix[8, :2] = [0.0, 1.0]
    matrix[9, :1] = [0.0]

    return matrix


@pytest.fixture
def sample_survey_extents():
    """Create sample survey extents for testing."""
    return pd.DataFrame(
        {
            "transect_num": [1, 2, 3, 4, 5],
            "x": [-2.5, -2.3, -2.7, -2.4, -2.6],
            "y": [-2.0, -1.0, 0.0, 1.0, 2.0],
        }
    )


@pytest.fixture
def sample_biomass_data():
    """Create sample biomass data for kriging tests."""
    return pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0, 4.0],
            "y": [0.0, 1.0, 2.0, 3.0, 4.0],
            "biomass_density": [10.5, 15.2, 12.8, 18.1, 14.3],
            "transect_num": [1, 1, 2, 2, 3],
        }
    )


@pytest.fixture
def sample_anisotropy_value():
    """Create sample anisotropy value for kriging tests."""
    return 0.01


@pytest.fixture
def sample_stacked_kriging_array():
    """Create sample stacked kriging array for parsing tests."""
    return np.array(
        [
            [1.0, 2.0, 3.0, 0.5, 12.5],
            [2.0, 3.0, 4.0, 0.3, 15.2],
            [3.0, 4.0, 5.0, 0.8, 18.1],
            [4.0, 5.0, 6.0, 0.2, 14.3],
        ]
    )


@pytest.fixture
def sample_kriging_results():
    """Create sample kriging results for projection tests."""
    return pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": [0.0, 1.0, 2.0, 3.0],
            "estimated_value": [10.5, 15.2, 12.8, 18.1],
            "estimation_variance": [0.5, 0.3, 0.8, 0.2],
        }
    )


@pytest.fixture
def sample_projection_parameters():
    """Create sample projection parameters."""
    return {
        "delta_longitude": 0.1,
        "delta_latitude": 0.05,
        "longitude_offset": -124.5,
        "latitude_offset": 46.2,
    }
