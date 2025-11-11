import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_western_extent():
    """Create sample western extent data."""
    return pd.DataFrame(
        {
            "transect_num": [1, 2, 3, 4, 5],
            "x": [-2.0, -1.8, -2.2, -1.9, -2.1],
            "y": [0.0, 2.0, 4.0, 6.0, 8.0],
        }
    )


@pytest.fixture
def sample_mesh_coordinates():
    """Create sample mesh coordinates for testing."""
    return pd.DataFrame(
        {
            "x": np.linspace(-2, 2, 10),
            "y": np.linspace(-2, 2, 10),
            "transect_num": np.repeat([1, 2, 3, 4, 5], 2),
        }
    )
