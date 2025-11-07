import warnings

import numpy as np
import pandas as pd
import pytest

from echopop.survey.biology import fit_length_weight_regression, quantize_length_data


def test_fit_length_weight_regression_basic(length_weight_data):
    """Test basic functionality of length-weight regression."""
    result = fit_length_weight_regression(length_weight_data)

    assert isinstance(result, pd.Series)
    assert list(result.index) == ["slope", "intercept"]
    assert len(result) == 2
    assert not result.isna().any()


def test_fit_length_weight_regression_coefficients(length_weight_data):
    """Test that regression coefficients are reasonable for biological data."""
    result = fit_length_weight_regression(length_weight_data)

    # For fish, slope should typically be between 2.5 and 3.5 (close to 3 for isometric growth)
    assert 2.0 < result["slope"] < 4.0
    # Intercept should be negative (log scale)
    assert result["intercept"] < 0


def test_fit_length_weight_regression_with_groupby(grouped_length_weight_data):
    """Test function works correctly with pandas groupby."""
    result = grouped_length_weight_data.groupby("sex").apply(
        fit_length_weight_regression, include_groups=False
    )

    assert isinstance(result, pd.DataFrame)
    assert result.index.names == ["sex"]
    assert list(result.columns) == ["slope", "intercept"]
    assert len(result) == 2  # Two sex groups
    assert not result.isna().any().any()


def test_fit_length_weight_regression_multiple_groups(grouped_length_weight_data):
    """Test with multiple grouping variables."""
    result = grouped_length_weight_data.groupby(["sex", "stratum"]).apply(
        fit_length_weight_regression, include_groups=False
    )

    assert isinstance(result, pd.DataFrame)
    assert result.index.names == ["sex", "stratum"]
    assert len(result) == 2


def test_fit_length_weight_regression_missing_values(data_with_missing_values):
    """Test that missing values are properly handled."""
    result = fit_length_weight_regression(data_with_missing_values)

    assert isinstance(result, pd.Series)
    assert not result.isna().any()
    # Should use only the 2 complete cases (rows 0 and 3)


def test_fit_length_weight_regression_minimal_data(minimal_data):
    """Test with minimal two-point dataset."""
    result = fit_length_weight_regression(minimal_data)

    assert isinstance(result, pd.Series)
    assert len(result) == 2
    assert not result.isna().any()


def test_fit_length_weight_regression_single_row_error(single_row_data):
    """Test error handling with insufficient data."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        fit_length_weight_regression(single_row_data)
        # Just check that we got at least one warning about poorly conditioned polyfit
        assert any("poorly conditioned" in str(warning.message) for warning in w)


def test_fit_length_weight_regression_empty_after_dropna():
    """Test behavior when all data is dropped due to missing values."""
    data = pd.DataFrame({"length": [np.nan, np.nan], "weight": [np.nan, np.nan]})

    with pytest.raises(TypeError):
        fit_length_weight_regression(data)


def test_fit_length_weight_regression_missing_columns():
    """Test error handling for missing required columns."""
    data = pd.DataFrame({"size": [10.0, 20.0], "mass": [5.0, 35.0]})

    with pytest.raises(KeyError):
        fit_length_weight_regression(data)


def test_fit_length_weight_regression_zero_negative_values(zero_negative_data):
    """Test error handling for zero and negative values."""
    # log10 of zero/negative values should raise an error
    with pytest.raises((ValueError, RuntimeWarning)):
        fit_length_weight_regression(zero_negative_data)


def test_fit_length_weight_regression_realistic_data(realistic_fish_data):
    """Test with realistic fish data."""
    result = fit_length_weight_regression(realistic_fish_data)

    assert isinstance(result, pd.Series)
    # Realistic slope for fish should be close to 3
    assert 2.5 < result["slope"] < 3.5
    assert result["intercept"] < 0


def test_fit_length_weight_regression_large_dataset(large_dataset):
    """Test performance and accuracy with larger dataset."""
    result = fit_length_weight_regression(large_dataset)

    assert isinstance(result, pd.Series)
    assert not result.isna().any()
    # With synthetic data where weight = 0.01 * length^3, slope should be close to 3
    assert 2.8 < result["slope"] < 3.2


def test_fit_length_weight_regression_preserves_input():
    """Test that input DataFrame is not modified."""
    data = pd.DataFrame(
        {"length": [10.0, np.nan, 20.0], "weight": [5.0, 15.0, 35.0], "other_col": ["A", "B", "C"]}
    )
    data_copy = data.copy()

    # Actually use the result to make the test meaningful
    result = fit_length_weight_regression(data)

    # Verify the function worked
    assert isinstance(result, pd.Series)

    # Verify input wasn't modified
    pd.testing.assert_frame_equal(data, data_copy)


def test_fit_length_weight_regression_deterministic():
    """Test that function produces deterministic results."""
    data = pd.DataFrame({"length": [10.0, 15.0, 20.0, 25.0], "weight": [5.0, 15.0, 35.0, 70.0]})

    result1 = fit_length_weight_regression(data)
    result2 = fit_length_weight_regression(data)

    pd.testing.assert_series_equal(result1, result2)


def test_fit_length_weight_regression_column_types():
    """Test that function works with different numeric column types."""
    data = pd.DataFrame(
        {
            "length": pd.array([10, 15, 20, 25], dtype="Int64"),
            "weight": pd.array([5.0, 15.0, 35.0, 70.0], dtype="float32"),
        }
    )

    result = fit_length_weight_regression(data)

    assert isinstance(result, pd.Series)
    assert not result.isna().any()


# ==============================================================================
# TESTS FOR utils.quantize_length_data
# ==============================================================================


def test_quantize_length_data_no_count(specimen_df_no_count, group_columns):
    """Test quantization when no length_count column exists."""
    result = quantize_length_data(specimen_df_no_count, group_columns)

    assert "length_count" in result.columns
    assert isinstance(result.index, pd.MultiIndex)

    # Check that it correctly counts occurrences
    # Length 20.5 appears twice in stratum 1, haul 101, sex M
    mask = (
        (result.index.get_level_values("stratum_ks") == 1)
        & (result.index.get_level_values("haul_num") == 101)
        & (result.index.get_level_values("sex") == "M")
        & (result.index.get_level_values("length") == 20.5)
    )

    assert result[mask]["length_count"].iloc[0] == 2


def test_quantize_length_data_with_count(length_df, group_columns):
    """Test quantization when length_count column already exists."""
    result = quantize_length_data(length_df, group_columns)

    assert "length_count" in result.columns
    assert isinstance(result.index, pd.MultiIndex)

    # Should sum existing counts, not count rows
    total_original = length_df["length_count"].sum()
    total_result = result["length_count"].sum()
    assert total_original == total_result


def test_quantize_length_data_empty():
    """Test quantization with empty DataFrame."""
    empty_df = pd.DataFrame(
        {"stratum_ks": [], "haul_num": [], "sex": [], "length": [], "length_count": []}
    )
    group_columns = ["stratum_ks", "haul_num", "sex"]
    result = quantize_length_data(empty_df, group_columns)

    assert len(result) == 0
    assert "length_count" in result.columns


def test_quantize_length_data_single_row(single_row_df, group_columns):
    """Test quantization with single row."""
    result = quantize_length_data(single_row_df, group_columns)

    assert len(result) == 1
    assert result["length_count"].iloc[0] == 1


def test_quantize_length_data_different_groups():
    """Test quantization with different grouping columns."""
    df = pd.DataFrame(
        {"stratum": [1, 1, 2, 2], "sex": ["M", "M", "F", "F"], "length": [20.0, 20.0, 22.0, 24.0]}
    )

    result = quantize_length_data(df, ["stratum", "sex"])

    # Should have 3 unique combinations
    assert len(result) == 3

    # Stratum 1, sex M, length 20.0 should have count 2
    mask = (
        (result.index.get_level_values("stratum") == 1)
        & (result.index.get_level_values("sex") == "M")
        & (result.index.get_level_values("length") == 20.0)
    )
    assert result[mask]["length_count"].iloc[0] == 2
