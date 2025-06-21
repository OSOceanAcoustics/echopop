import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def simple_bins():
    """Create simple evenly spaced bins for testing."""
    return np.array([1, 2, 3, 4, 5])


@pytest.fixture
def float_bins():
    """Create floating point bins for testing."""
    return np.array([1.5, 2.7, 3.9, 5.1, 6.3])


@pytest.fixture
def linspace_bins():
    """Create bins using np.linspace for testing."""
    return np.linspace(0, 10, 11)


@pytest.fixture
def uneven_bins():
    """Create unevenly spaced bins for testing."""
    return np.array([0, 1, 3, 7, 15, 31])


@pytest.fixture
def age_bins():
    """Create age bins similar to those used in biological data."""
    return np.linspace(1.0, 22.0, 22)


@pytest.fixture
def length_bins():
    """Create length bins similar to those used in biological data."""
    return np.arange(2.0, 82.0, 2.0)


@pytest.fixture
def minimal_bins():
    """Create minimal two-element bins for edge case testing."""
    return np.array([0, 1])


@pytest.fixture
def large_bins():
    """Create a large array of bins for performance testing."""
    return np.linspace(0, 1000, 1001)


@pytest.fixture
def numeric_bins():
    """Create numeric bins for testing."""
    return np.array([10, 15, 20, 25, 30])


@pytest.fixture
def secondary_bins():
    """Create secondary bins for testing."""
    return np.array([1, 2, 3, 4, 5])


@pytest.fixture
def target_dataframe():
    """Create DataFrame with target columns for binning."""
    return pd.DataFrame(
        {
            "numeric_col": [12.5, 18.3, 22.7, 27.1, 15.8, 24.2],
            "secondary_col": [1.2, 2.8, 3.5, 4.1, 1.9, 3.2],
            "value_col": [25.3, 45.7, 68.2, 89.5, 32.1, 72.8],
            "category_col": ["A", "B", "A", "B", "A", "B"],
            "id_col": [1, 2, 3, 4, 5, 6],
        }
    )


@pytest.fixture
def non_target_dataframe():
    """Create DataFrame without target columns."""
    return pd.DataFrame(
        {
            "id_col": [1, 2, 3, 4],
            "coord_x": [45.2, 46.1, 44.8, 45.7],
            "coord_y": [-125.3, -124.8, -125.7, -125.1],
            "depth_col": [100, 150, 120, 180],
        }
    )


@pytest.fixture
def partial_target_dataframe():
    """Create DataFrame with some but not all target columns."""
    return pd.DataFrame(
        {
            "id_col": [1, 1, 2, 2],
            "group_col": [100, 101, 100, 101],
            "secondary_col": [2.5, 3.1, 1.8, 4.2],
            "measurement_col": [125.5, 89.3, 156.7, 203.4],
        }
    )


@pytest.fixture
def mixed_dataframes_dict(target_dataframe, non_target_dataframe, partial_target_dataframe):
    """Create dictionary of mixed DataFrames."""
    return {
        "target_data": target_dataframe,
        "non_target_data": non_target_dataframe,
        "partial_data": partial_target_dataframe,
    }


@pytest.fixture
def empty_dataframe():
    """Create empty DataFrame."""
    return pd.DataFrame()


@pytest.fixture
def single_row_dataframe():
    """Create single-row DataFrame."""
    return pd.DataFrame({"numeric_col": [15.5], "secondary_col": [2.3], "value_col": [35.2]})


@pytest.fixture
def large_dataframe():
    """Create larger DataFrame for performance testing."""
    np.random.seed(42)
    n = 1000
    return pd.DataFrame(
        {
            "numeric_col": np.random.uniform(10, 30, n),
            "secondary_col": np.random.uniform(1, 5, n),
            "value_col": np.random.uniform(20, 100, n),
            "category_col": np.random.choice(["A", "B"], n),
            "id_col": np.arange(n),
        }
    )


@pytest.fixture
def invalid_bins():
    """Create invalid bins for error testing."""
    return np.array([10])  # Single element array


@pytest.fixture
def small_bins():
    """Create small bins for testing."""
    return np.array([10, 20, 30])


@pytest.fixture
def mixed_objects_dict(target_dataframe, non_target_dataframe):
    """Create dictionary with mixed DataFrame and non-DataFrame objects."""
    return {
        "dataframe_1": target_dataframe,
        "dataframe_2": non_target_dataframe,
        "metadata": {"source": "test", "version": 1.0},  # Non-DataFrame
        "config": [1, 2, 3],  # Non-DataFrame
    }


@pytest.fixture
def dataframe_with_na():
    """Create DataFrame with NaN values in target column."""
    return pd.DataFrame(
        {
            "numeric_col": [10.5, np.nan, 20.3, 25.1, np.nan],
            "value_col": [25.0, 30.0, 45.0, 60.0, 35.0],
        }
    )
