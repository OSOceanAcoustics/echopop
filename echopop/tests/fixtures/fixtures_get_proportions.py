import numpy as np
import pandas as pd
import pytest
from scipy import interpolate as interp


@pytest.fixture
def sample_specimen_data():
    """Create sample specimen data for testing."""
    return pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2, 1, 2, 1, 2, 1, 2],
            "length_bin": pd.Categorical(
                [
                    "(10,15]",
                    "(15,20]",
                    "(20,25]",
                    "(25,30]",
                    "(10,15]",
                    "(15,20]",
                    "(20,25]",
                    "(25,30]",
                    "(15,20]",
                    "(20,25]",
                ]
            ),
            "age_bin": pd.Categorical([1, 2, 3, 4, 1, 2, 3, 4, 2, 3]),
            "sex": [
                "male",
                "female",
                "male",
                "female",
                "unsexed",
                "male",
                "female",
                "unsexed",
                "male",
                "female",
            ],
            "length": [12.5, 18.3, 22.7, 27.1, 13.2, 17.8, 23.5, 28.4, 19.1, 21.6],
            "weight": [25.3, 45.7, 68.2, 89.5, 28.1, 42.3, 71.8, 92.1, 48.2, 63.7],
        }
    )


@pytest.fixture
def sample_length_data():
    """Create sample length data for testing."""
    return pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2, 1, 2, 1, 2],
            "length_bin": pd.Categorical(
                [
                    "(10,15]",
                    "(15,20]",
                    "(20,25]",
                    "(25,30]",
                    "(10,15]",
                    "(15,20]",
                    "(20,25]",
                    "(25,30]",
                ]
            ),
            "sex": ["male", "female", "male", "female", "unsexed", "male", "female", "unsexed"],
            "length": [12.8, 17.5, 22.1, 26.9, 14.1, 18.7, 24.3, 27.8],
            "length_count": [15, 23, 18, 12, 8, 19, 14, 9],
        }
    )


@pytest.fixture
def minimal_specimen_data():
    """Create minimal specimen data for edge case testing."""
    return pd.DataFrame(
        {
            "stratum_num": [1, 1],
            "length_bin": pd.Categorical(["(10,15]", "(15,20]"]),
            "length": [12.0, 17.0],
            "sex": ["male", "female"],
        }
    )


@pytest.fixture
def null_specimen_data():
    """Create empty specimen DataFrame."""
    return pd.DataFrame(columns=["stratum_num", "length_bin", "age_bin", "sex", "length", "weight"])


@pytest.fixture
def data_with_missing_columns():
    """Create data missing some expected columns."""
    return pd.DataFrame(
        {
            "stratum_num": [1, 2, 1, 2],
            "length_bin": pd.Categorical(["(10,15]", "(15,20]", "(20,25]", "(25,30]"]),
            "length": [12.5, 18.3, 22.7, 27.1],
            "other_col": ["A", "B", "C", "D"],
        }
    )


@pytest.fixture
def data_all_unsexed():
    """Create data where all fish are unsexed."""
    return pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2],
            "length_bin": pd.Categorical(["(10,15]", "(15,20]", "(20,25]", "(25,30]"]),
            "sex": ["unsexed", "unsexed", "unsexed", "unsexed"],
            "length": [12.5, 18.3, 22.7, 27.1],
        }
    )


@pytest.fixture
def expanded_specimen_dataset():
    """Create larger dataset for performance testing."""
    np.random.seed(42)
    n = 1000

    return pd.DataFrame(
        {
            "stratum_num": np.random.randint(1, 6, n),
            "length_bin": pd.Categorical(
                np.random.choice(["(10,15]", "(15,20]", "(20,25]", "(25,30]", "(30,35]"], n)
            ),
            "age_bin": pd.Categorical(np.random.randint(1, 6, n)),
            "sex": np.random.choice(["male", "female", "unsexed"], n),
            "length": np.random.uniform(10, 35, n),
        }
    )


@pytest.fixture
def data_multiple_groups():
    """Create data with multiple grouping variables."""
    return pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2, 1, 2],
            "length_bin": pd.Categorical(
                ["(10,15]", "(15,20]", "(20,25]", "(25,30]", "(10,15]", "(15,20]"]
            ),
            "age_bin": pd.Categorical([1, 2, 3, 4, 1, 2]),
            "sex": ["male", "female", "male", "female", "unsexed", "male"],
            "species": ["A", "A", "B", "B", "A", "B"],
            "region": ["North", "South", "North", "South", "North", "South"],
            "length": [12.5, 18.3, 22.7, 27.1, 13.2, 17.8],
        }
    )


@pytest.fixture
def aged_dataframe():
    """Create sample aged count data with stratum, length, age, and sex."""
    data = {
        "stratum_num": [1, 1, 1, 2, 2, 2],
        "length_bin": ["(10, 20]", "(20, 30]", "(30, 40]", "(10, 20]", "(20, 30]", "(30, 40]"],
        "age_bin": ["(1, 2]", "(2, 3]", "(3, 4]", "(1, 2]", "(2, 3]", "(3, 4]"],
        "sex": ["female", "male", "unsexed", "female", "male", "unsexed"],
        "species_id": [1, 1, 1, 1, 1, 1],
        "count": [10, 20, 5, 15, 25, 10],
    }
    return pd.DataFrame(data)


@pytest.fixture
def unaged_dataframe():
    """Create sample unaged count data with stratum, length and sex."""
    data = {
        "stratum_num": [1, 1, 1, 2, 2, 2],
        "length_bin": ["(10, 20]", "(20, 30]", "(30, 40]", "(10, 20]", "(20, 30]", "(30, 40]"],
        "sex": ["female", "male", "unsexed", "female", "male", "unsexed"],
        "species_id": [1, 1, 1, 1, 1, 1],
        "count": [30, 40, 10, 35, 45, 15],
    }
    return pd.DataFrame(data)


@pytest.fixture
def length_dataset_with_bins(grouped_length_weight_data):
    """Create dataset with length bins and counts using pandas Interval."""
    from pandas import Interval

    df = grouped_length_weight_data.copy()

    # Add length bins using pandas Interval objects
    bins = []
    for length in df["length"]:
        bin_min = length - 2.5
        bin_max = length + 2.5
        bins.append(Interval(bin_min, bin_max, closed="right"))

    df["length_bin"] = bins

    # Add length counts
    df["length_count"] = [10, 15, 20, 5, 8, 12, 18, 6]

    return df


@pytest.fixture
def length_weight_dataset_with_bins(grouped_length_weight_data):
    """Create length-weight dataset with bins and fitted weights using pandas Interval."""
    from pandas import Interval

    df = grouped_length_weight_data.copy()

    # Add length bins using pandas Interval
    bins = []
    for length in df["length"]:
        bin_min = length - 2.5
        bin_max = length + 2.5
        bins.append(Interval(bin_min, bin_max, closed="right"))

    df["length_bin"] = bins

    # Add fitted weights (simulate a simple length-weight relationship)
    df["weight_fitted"] = df["length"] ** 3 * 0.01

    return df


@pytest.fixture
def real_interpolators():
    """Create real interpolator functions for testing."""
    # Male interpolator: f(x) = x * 3.0
    x_male = np.array([10.0, 20.0, 30.0])
    y_male = np.array([30.0, 60.0, 90.0])
    male_interp = interp.interp1d(
        x_male, y_male, kind="linear", bounds_error=False, fill_value=(y_male[0], y_male[-1])
    )

    # Female interpolator: f(x) = x * 2.0
    x_female = np.array([10.0, 20.0, 30.0])
    y_female = np.array([20.0, 40.0, 60.0])
    female_interp = interp.interp1d(
        x_female,
        y_female,
        kind="linear",
        bounds_error=False,
        fill_value=(y_female[0], y_female[-1]),
    )

    # Global interpolator: f(x) = x * 2.5
    x_global = np.array([10.0, 20.0, 30.0])
    y_global = np.array([25.0, 50.0, 75.0])
    global_interp = interp.interp1d(
        x_global,
        y_global,
        kind="linear",
        bounds_error=False,
        fill_value=(y_global[0], y_global[-1]),
    )

    return {"male": male_interp, "female": female_interp, "_global_": global_interp}


@pytest.fixture
def proportion_dict():
    """Create sample proportion dictionaries for testing."""
    # First group - aged data
    aged_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 2, 2, 2],
            "sex": ["female", "male", "male", "female", "male", "female"],
            "length_bin": ["(10, 20]", "(10, 20]", "(20, 30]", "(10, 20]", "(20, 30]", "(20, 30]"],
            "proportion": [0.3, 0.4, 0.1, 0.05, 0.15, 0.1],
            "proportion_overall": [0.15, 0.2, 0.05, 0.025, 0.075, 0.05],
        }
    )

    # Second group - unaged data
    unaged_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 2, 2, 2],
            "sex": ["female", "male", "male", "female", "male", "female"],
            "length_bin": ["(10, 20]", "(10, 20]", "(20, 30]", "(10, 20]", "(20, 30]", "(20, 30]"],
            "proportion": [0.25, 0.35, 0.2, 0.1, 0.2, 0.1],
            "proportion_overall": [0.125, 0.175, 0.1, 0.05, 0.1, 0.05],
        }
    )

    return {"aged": aged_data, "unaged": unaged_data}


@pytest.fixture
def weight_table():
    """Create sample weight table for testing."""
    return pd.DataFrame(
        {
            "sex": ["female", "male", "all", "female", "male", "all"],
            "length_bin": ["(10, 20]", "(10, 20]", "(10, 20]", "(20, 30]", "(20, 30]", "(20, 30]"],
            "weight_fitted": [0.5, 0.4, 0.45, 1.2, 1.0, 1.1],
        }
    )


@pytest.fixture
def proportion_test_dict():
    """Create sample proportion dictionary for testing stratum_averaged_weight."""
    # First group - aged data
    aged_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 2, 2, 2],
            "sex": ["female", "male", "male", "female", "male", "female"],
            "length_bin": ["(10, 20]", "(10, 20]", "(20, 30]", "(10, 20]", "(20, 30]", "(20, 30]"],
            "proportion": [0.3, 0.4, 0.1, 0.05, 0.15, 0.1],
            "proportion_overall": [0.15, 0.2, 0.05, 0.025, 0.075, 0.05],
        }
    )

    # Second group - unaged data
    unaged_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 2, 2, 2],
            "sex": ["female", "male", "male", "female", "male", "female"],
            "length_bin": ["(10, 20]", "(10, 20]", "(20, 30]", "(10, 20]", "(20, 30]", "(20, 30]"],
            "proportion": [0.25, 0.35, 0.2, 0.1, 0.2, 0.1],
            "proportion_overall": [0.125, 0.175, 0.1, 0.05, 0.1, 0.05],
        }
    )

    return {"aged": aged_data, "unaged": unaged_data}


@pytest.fixture
def test_weight_table():
    """Create sample weight table for testing stratum_averaged_weight."""
    return pd.DataFrame(
        {
            "sex": ["female", "male", "all", "female", "male", "all"],
            "length_bin": ["(10, 20]", "(10, 20]", "(10, 20]", "(20, 30]", "(20, 30]", "(20, 30]"],
            "weight_fitted": [0.5, 0.4, 0.45, 1.2, 1.0, 1.1],
        }
    )
