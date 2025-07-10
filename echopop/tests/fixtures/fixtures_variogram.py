import numpy as np
import pandas as pd
import pytest


# ==================================================================================================
# Fixtures for variogram testing
# ------------------------------
@pytest.fixture
def sample_estimates():
    """Create sample field estimates for variogram analysis."""
    return np.array([10.5, 15.2, 12.8, 18.1, 14.3, 16.7, 11.9, 13.4])


@pytest.fixture
def sample_lag_counts():
    """Create sample lag counts."""
    return np.array([20, 18, 15, 12, 8, 5])


@pytest.fixture
def sample_lag_estimates():
    """Create sample lag estimates."""
    return np.array([280.0, 250.0, 210.0, 180.0, 120.0, 75.0])


@pytest.fixture
def sample_lag_estimates_squared():
    """Create sample squared lag estimates."""
    return np.array([3920.0, 3600.0, 3150.0, 2700.0, 1800.0, 1125.0])


@pytest.fixture
def sample_lag_deviations():
    """Create sample lag deviations."""
    return np.array([25.5, 22.8, 18.9, 15.2, 10.4, 6.7])


@pytest.fixture
def sample_head_index():
    """Create sample head index matrix."""
    return np.array(
        [
            [3, 2, 1, 1, 0, 0],
            [2, 3, 2, 1, 1, 0],
            [1, 2, 3, 2, 1, 0],
            [1, 1, 2, 3, 2, 1],
            [0, 1, 1, 2, 3, 2],
            [0, 0, 1, 1, 2, 3],
            [0, 0, 0, 1, 1, 2],
            [0, 0, 0, 0, 1, 1],
        ]
    )


@pytest.fixture
def sample_transect_df():
    """Create sample transect DataFrame for empirical variogram testing."""
    return pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            "y": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
            "biomass_density": [10.5, 15.2, 12.8, 18.1, 14.3, 16.7, 11.9, 13.4],
        }
    )


@pytest.fixture
def sample_variogram_lag_matrix():
    """Create sample lag matrix for variogram analysis."""
    return np.array(
        [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 0, 1, 2, 3, 4, 5, 6],
            [2, 1, 0, 1, 2, 3, 4, 5],
            [3, 2, 1, 0, 1, 2, 3, 4],
            [4, 3, 2, 1, 0, 1, 2, 3],
            [5, 4, 3, 2, 1, 0, 1, 2],
            [6, 5, 4, 3, 2, 1, 0, 1],
            [7, 6, 5, 4, 3, 2, 1, 0],
        ]
    )


@pytest.fixture
def sample_variogram_mask_matrix():
    """Create sample mask matrix for variogram analysis."""
    return np.array(
        [
            [False, True, True, True, True, True, True, True],
            [False, False, True, True, True, True, True, True],
            [False, False, False, True, True, True, True, True],
            [False, False, False, False, True, True, True, True],
            [False, False, False, False, False, True, True, True],
            [False, False, False, False, False, False, True, True],
            [False, False, False, False, False, False, False, True],
            [False, False, False, False, False, False, False, False],
        ]
    )


@pytest.fixture
def sample_variogram_azimuth_matrix():
    """Create sample azimuth matrix for variogram analysis."""
    return np.array(
        [
            [0.0, 30.0, 45.0, 60.0, 90.0, 120.0, 135.0, 150.0],
            [210.0, 0.0, 30.0, 45.0, 60.0, 90.0, 120.0, 135.0],
            [225.0, 210.0, 0.0, 30.0, 45.0, 60.0, 90.0, 120.0],
            [240.0, 225.0, 210.0, 0.0, 30.0, 45.0, 60.0, 90.0],
            [270.0, 240.0, 225.0, 210.0, 0.0, 30.0, 45.0, 60.0],
            [300.0, 270.0, 240.0, 225.0, 210.0, 0.0, 30.0, 45.0],
            [315.0, 300.0, 270.0, 240.0, 225.0, 210.0, 0.0, 30.0],
            [330.0, 315.0, 300.0, 270.0, 240.0, 225.0, 210.0, 0.0],
        ]
    )
