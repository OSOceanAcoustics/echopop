import pandas as pd
import pytest


@pytest.fixture
def basic_nasc_data():
    """Basic NASC DataFrame for testing"""
    return pd.DataFrame(
        {
            "distance_s": [0.0, 1.0, 2.0, 3.0, 4.0],
            "distance_e": [1.0, 2.0, 3.0, 4.0, 5.0],
            "transect_spacing": [0.1, 0.1, 0.1, 0.1, 0.1],
        }
    )


@pytest.fixture
def transect_single_row():
    """Single row DataFrame for edge case testing"""
    return pd.DataFrame({"distance_s": [0.0], "distance_e": [1.5], "transect_spacing": [0.1]})


@pytest.fixture
def irregular_spacing_data():
    """DataFrame with irregular spacing for testing threshold logic"""
    return pd.DataFrame(
        {
            "distance_s": [0.0, 1.0, 2.0, 5.0, 6.0],  # Gap between 2.0 and 5.0
            "distance_e": [1.0, 2.0, 3.0, 6.0, 7.0],
            "transect_spacing": [0.1, 0.1, 0.1, 0.1, 0.1],
        }
    )
