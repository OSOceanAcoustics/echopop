import numpy as np
import pandas as pd
import pytest


# ==================================================================================================
# Fixtures for spatial coordinate testing
# ---------------------------------------
@pytest.fixture
def sample_coordinates_df():
    """Create a sample DataFrame with longitude and latitude coordinates."""
    return pd.DataFrame(
        {
            "longitude": [-124.5, -124.3, -124.1, -123.9],
            "latitude": [46.2, 46.4, 46.6, 46.8],
            "value": [10.5, 15.2, 12.8, 18.1],
        }
    )


@pytest.fixture
def sample_coordinates_array():
    """Create a sample numpy array with x, y coordinates."""
    return np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])


@pytest.fixture
def sample_reference_df():
    """Create a reference DataFrame for coordinate standardization."""
    return pd.DataFrame(
        {"longitude": [-124.6, -124.4, -124.2, -124.0], "latitude": [46.0, 46.2, 46.4, 46.6]}
    )


@pytest.fixture
def sample_lag_matrix():
    """Create a sample lag matrix for testing."""
    return np.array([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]])


@pytest.fixture
def sample_mask_matrix():
    """Create a sample boolean mask matrix."""
    return np.array(
        [
            [False, True, True, True],
            [False, False, True, True],
            [False, False, False, True],
            [False, False, False, False],
        ]
    )


@pytest.fixture
def sample_azimuth_matrix():
    """Create a sample azimuth matrix."""
    return np.array(
        [
            [0.0, 45.0, 90.0, 135.0],
            [225.0, 0.0, 45.0, 90.0],
            [270.0, 225.0, 0.0, 45.0],
            [315.0, 270.0, 225.0, 0.0],
        ]
    )


@pytest.fixture
def sample_mesh_with_area():
    """Create a sample mesh DataFrame with required area column."""
    return pd.DataFrame(
        {
            "longitude": np.linspace(-124.6, -123.0, 10),
            "latitude": np.linspace(46.0, 46.8, 10),
            "x": np.linspace(-0.5, 1.0, 10),
            "y": np.linspace(-0.2, 1.2, 10),
            "area": np.ones(10) * 1.5,  # Required area column
            "mesh_id": range(10),
        }
    )


@pytest.fixture
def sample_mesh_with_fraction():
    """Create a sample mesh DataFrame with fraction column instead of area."""
    return pd.DataFrame(
        {
            "longitude": np.linspace(-124.6, -123.0, 8),
            "latitude": np.linspace(46.0, 46.8, 8),
            "x": np.linspace(-0.5, 1.0, 8),
            "y": np.linspace(-0.2, 1.2, 8),
            "fraction": np.ones(8) * 0.25,  # Fraction column for validation
            "mesh_id": range(8),
        }
    )
