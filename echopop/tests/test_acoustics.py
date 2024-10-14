import numpy as np
import pandas as pd
import pytest

from ..acoustics import impute_missing_sigma_bs, to_dB, to_linear, ts_length_regression
from .conftest import assert_dataframe_equal


def test_ts_length_regression():

    # -------------------------
    # Mock values
    # ---- length values [ ARRAY input ]
    mock_length_array = np.array([1.0, 2.0, 3.0])
    # ---- x values [ FLOAT input ]
    mock_length_float = np.array([1.0])
    # ---- Slope
    mock_slope = 5.0
    # ---- Intercept
    mock_intercept = -2.0

    # -------------------------
    # Evaluate [ ARRAY ]
    test_results_array = ts_length_regression(mock_length_array, mock_slope, mock_intercept)
    # Evaluate [ FLOAT ]
    test_results_float = ts_length_regression(mock_length_float, mock_slope, mock_intercept)

    # -------------------------
    # Expected outcomes
    # Expected [ ARRAY ]
    expected_array = np.array([-2.0, -0.49485002, 0.38560627])
    # Expected [ FLOAT ]
    expected_float = np.array([-2.0])

    # -------------------------
    # Run tests [ ARRAY ]
    # ---- Type
    assert isinstance(test_results_array, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_array, expected_array)
    # Run tests [ FLOAT ]
    # ---- Type
    assert isinstance(test_results_float, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_float, expected_float)


def test_to_linear():

    # -------------------------
    # Mock values
    # ---- length values [ ARRAY input ]
    mock_db_array = np.array([-80.0, -60.0, -40.0])
    # ---- x values [ FLOAT input ]
    mock_db_float = np.array([-60.0])

    # -------------------------
    # Evaluate [ ARRAY ]
    test_results_array = to_linear(mock_db_array)
    # Evaluate [ FLOAT ]
    test_results_float = to_linear(mock_db_float)

    # -------------------------
    # Expected outcomes
    # Expected [ ARRAY ]
    expected_array = np.array([1e-8, 1e-6, 1e-4])
    # Expected [ FLOAT ]
    expected_float = np.array([1e-6])

    # -------------------------
    # Run tests [ ARRAY ]
    # ---- Type
    assert isinstance(test_results_array, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_array, expected_array)
    # Run tests [ FLOAT ]
    # ---- Type
    assert isinstance(test_results_float, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_float, expected_float)


def test_to_dB():

    # -------------------------
    # Mock values
    # ---- length values [ ARRAY input ]
    mock_linear_array = np.array([1e-8, 1e-6, 1e-4])
    # ---- x values [ FLOAT input ]
    mock_linear_float = np.array([1e-6])

    # -------------------------
    # Evaluate [ ARRAY ]
    test_results_array = to_dB(mock_linear_array)
    # Evaluate [ FLOAT ]
    test_results_float = to_dB(mock_linear_float)

    # -------------------------
    # Expected outcomes
    # Expected [ ARRAY ]
    expected_array = np.array([-80.0, -60.0, -40.0])
    # Expected [ FLOAT ]
    expected_float = np.array([-60.0])

    # -------------------------
    # Run tests [ ARRAY ]
    # ---- Type
    assert isinstance(test_results_array, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_array, expected_array)
    # Run tests [ FLOAT ]
    # ---- Type
    assert isinstance(test_results_float, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_float, expected_float)


@pytest.mark.parametrize(
    "strata, dataframe, expected",
    [
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [1.0, 2.0, 3.0, 4.0, 5.0],
                }
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [1.0, 2.0, 3.0, 4.0, 5.0],
                }
            ),
        ),
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 4, 5],
                    "species_id": np.repeat(94832, 4),
                    "sigma_bs_mean": [1.0, 2.0, 4.0, 5.0],
                }
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [1.0, 2.0, 3.0, 4.0, 5.0],
                }
            ),
        ),
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {
                    "stratum_num": [2, 3, 4, 5],
                    "species_id": np.repeat(94832, 4),
                    "sigma_bs_mean": [2.0, 3.0, 4.0, 5.0],
                }
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [2.0, 2.0, 3.0, 4.0, 5.0],
                }
            ),
        ),
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4],
                    "species_id": np.repeat(94832, 4),
                    "sigma_bs_mean": [1.0, 2.0, 3.0, 4.0],
                }
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [1.0, 2.0, 3.0, 4.0, 4.0],
                }
            ),
        ),
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {
                    "stratum_num": [2, 3, 4],
                    "species_id": np.repeat(94832, 3),
                    "sigma_bs_mean": [2.0, 3.0, 4.0],
                }
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [2.0, 2.0, 3.0, 4.0, 4.0],
                }
            ),
        ),
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3],
                    "species_id": np.repeat(94832, 3),
                    "sigma_bs_mean": [1.0, 2.0, 3.0],
                }
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [1.0, 2.0, 3.0, 3.0, 3.0],
                }
            ),
        ),
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {
                    "stratum_num": [3, 4, 5],
                    "species_id": np.repeat(94832, 3),
                    "sigma_bs_mean": [3.0, 4.0, 5.0],
                }
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [3.0, 3.0, 3.0, 4.0, 5.0],
                }
            ),
        ),
        (
            np.array([1, 2, 3, 4, 5]),
            pd.DataFrame(
                {"stratum_num": [3], "species_id": np.repeat(94832, 1), "sigma_bs_mean": [3.0]}
            ),
            pd.DataFrame(
                {
                    "stratum_num": [1, 2, 3, 4, 5],
                    "species_id": np.repeat(94832, 5),
                    "sigma_bs_mean": [3.0, 3.0, 3.0, 3.0, 3.0],
                }
            ),
        ),
    ],
    ids=[
        "No missing strata (valid)",
        "Missing stratum in middle: 3 (valid)",
        "Missing stratum at top: 1 (valid)",
        "Missing stratum at bottom: 5 (valid)",
        "Missing stratum at top and bottom: (1,5) (valid)",
        "Missing strata at bottom: (4,5) (valid)",
        "Missing strata at top: (1,2) (valid)",
        "Only 1 stratum: 3 (valid)",
    ],
)
def test_impute_missing_sigma_bs(strata, dataframe, expected):

    # -------------------------
    # Expected outcomes
    # ---- Types [~ALL]
    expected_dtypes = {
        "stratum_num": np.integer,
        "species_id": np.integer,
        "sigma_bs_mean": np.floating,
    }

    # -------------------------
    # COMPUTE
    results = impute_missing_sigma_bs(strata, dataframe)

    # -------------------------
    # ASSERT
    # ---- Is the output a DataFrame?
    assert isinstance(results, pd.DataFrame)
    # ---- Are the expected and results DataFrames equal and of the expected datatypes?
    assert_dataframe_equal(results, expected_dtypes, expected)
