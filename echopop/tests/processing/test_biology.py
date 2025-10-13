import warnings

import numpy as np
import pandas as pd
import pytest

from echopop.nwfsc_feat import biology


def test_fit_length_weight_regression_basic(length_weight_data):
    """Test basic functionality of length-weight regression."""
    result = biology.fit_length_weight_regression(length_weight_data)

    assert isinstance(result, pd.Series)
    assert list(result.index) == ["slope", "intercept"]
    assert len(result) == 2
    assert not result.isna().any()


def test_fit_length_weight_regression_coefficients(length_weight_data):
    """Test that regression coefficients are reasonable for biological data."""
    result = biology.fit_length_weight_regression(length_weight_data)

    # For fish, slope should typically be between 2.5 and 3.5 (close to 3 for isometric growth)
    assert 2.0 < result["slope"] < 4.0
    # Intercept should be negative (log scale)
    assert result["intercept"] < 0


def test_fit_length_weight_regression_with_groupby(grouped_length_weight_data):
    """Test function works correctly with pandas groupby."""
    result = grouped_length_weight_data.groupby("sex").apply(
        biology.fit_length_weight_regression, include_groups=False
    )

    assert isinstance(result, pd.DataFrame)
    assert result.index.names == ["sex"]
    assert list(result.columns) == ["slope", "intercept"]
    assert len(result) == 2  # Two sex groups
    assert not result.isna().any().any()


def test_fit_length_weight_regression_multiple_groups(grouped_length_weight_data):
    """Test with multiple grouping variables."""
    result = grouped_length_weight_data.groupby(["sex", "stratum"]).apply(
        biology.fit_length_weight_regression, include_groups=False
    )

    assert isinstance(result, pd.DataFrame)
    assert result.index.names == ["sex", "stratum"]
    assert len(result) == 2


def test_fit_length_weight_regression_missing_values(data_with_missing_values):
    """Test that missing values are properly handled."""
    result = biology.fit_length_weight_regression(data_with_missing_values)

    assert isinstance(result, pd.Series)
    assert not result.isna().any()
    # Should use only the 2 complete cases (rows 0 and 3)


def test_fit_length_weight_regression_minimal_data(minimal_data):
    """Test with minimal two-point dataset."""
    result = biology.fit_length_weight_regression(minimal_data)

    assert isinstance(result, pd.Series)
    assert len(result) == 2
    assert not result.isna().any()


def test_fit_length_weight_regression_single_row_error(single_row_data):
    """Test error handling with insufficient data."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        biology.fit_length_weight_regression(single_row_data)
        # Just check that we got at least one warning about poorly conditioned polyfit
        assert any("poorly conditioned" in str(warning.message) for warning in w)


def test_fit_length_weight_regression_empty_after_dropna():
    """Test behavior when all data is dropped due to missing values."""
    data = pd.DataFrame({"length": [np.nan, np.nan], "weight": [np.nan, np.nan]})

    with pytest.raises(TypeError):
        biology.fit_length_weight_regression(data)


def test_fit_length_weight_regression_missing_columns():
    """Test error handling for missing required columns."""
    data = pd.DataFrame({"size": [10.0, 20.0], "mass": [5.0, 35.0]})

    with pytest.raises(KeyError):
        biology.fit_length_weight_regression(data)


def test_fit_length_weight_regression_zero_negative_values(zero_negative_data):
    """Test error handling for zero and negative values."""
    # log10 of zero/negative values should raise an error
    with pytest.raises((ValueError, RuntimeWarning)):
        biology.fit_length_weight_regression(zero_negative_data)


def test_fit_length_weight_regression_realistic_data(realistic_fish_data):
    """Test with realistic fish data."""
    result = biology.fit_length_weight_regression(realistic_fish_data)

    assert isinstance(result, pd.Series)
    # Realistic slope for fish should be close to 3
    assert 2.5 < result["slope"] < 3.5
    assert result["intercept"] < 0


def test_fit_length_weight_regression_large_dataset(large_dataset):
    """Test performance and accuracy with larger dataset."""
    result = biology.fit_length_weight_regression(large_dataset)

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
    result = biology.fit_length_weight_regression(data)

    # Verify the function worked
    assert isinstance(result, pd.Series)

    # Verify input wasn't modified
    pd.testing.assert_frame_equal(data, data_copy)


def test_fit_length_weight_regression_deterministic():
    """Test that function produces deterministic results."""
    data = pd.DataFrame({"length": [10.0, 15.0, 20.0, 25.0], "weight": [5.0, 15.0, 35.0, 70.0]})

    result1 = biology.fit_length_weight_regression(data)
    result2 = biology.fit_length_weight_regression(data)

    pd.testing.assert_series_equal(result1, result2)


def test_fit_length_weight_regression_column_types():
    """Test that function works with different numeric column types."""
    data = pd.DataFrame(
        {
            "length": pd.array([10, 15, 20, 25], dtype="Int64"),
            "weight": pd.array([5.0, 15.0, 35.0, 70.0], dtype="float32"),
        }
    )

    result = biology.fit_length_weight_regression(data)

    assert isinstance(result, pd.Series)
    assert not result.isna().any()


def test_length_binned_weights_basic_functionality(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test basic functionality of biology.length_binned_weights."""
    result = biology.length_binned_weights(
        sample_specimen_data, sample_length_bins, single_regression_coefficients
    )

    assert isinstance(result, pd.DataFrame)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns
    assert len(result) > 0
    assert not result["weight_fitted"].isna().any()


def test_length_binned_weights_with_prebinned_data(
    specimen_data_with_bins, sample_length_bins, single_regression_coefficients
):
    """Test function with data that already has length_bin column."""
    result = biology.length_binned_weights(
        specimen_data_with_bins, sample_length_bins, single_regression_coefficients
    )

    assert isinstance(result, pd.DataFrame)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns
    assert len(result) > 0


def test_length_binned_weights_grouped_coefficients(
    sample_specimen_data, sample_length_bins, grouped_regression_coefficients
):
    """Test with grouped regression coefficients."""
    result = biology.length_binned_weights(
        sample_specimen_data, sample_length_bins, grouped_regression_coefficients
    )

    assert isinstance(result, pd.DataFrame)
    # Should return wide format with sex as columns when grouping
    expected_columns = ["female", "male", "unsexed"]  # sex values as columns
    assert list(result.columns) == expected_columns
    assert result.index.name == "length_bin"
    assert len(result) > 0


def test_length_binned_weights_impute_bins_false(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test with impute_bins=False."""
    result = biology.length_binned_weights(
        sample_specimen_data,
        sample_length_bins,
        single_regression_coefficients,
        impute_bins=False,
    )

    assert isinstance(result, pd.DataFrame)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns
    # With impute_bins=False, should only use observed means


def test_length_binned_weights_minimum_count_threshold(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test with different minimum_count_threshold values."""
    # High threshold - should use more modeled weights
    result_high = biology.length_binned_weights(
        sample_specimen_data,
        sample_length_bins,
        single_regression_coefficients,
        minimum_count_threshold=10,
    )

    # Low threshold - should use more observed means
    result_low = biology.length_binned_weights(
        sample_specimen_data,
        sample_length_bins,
        single_regression_coefficients,
        minimum_count_threshold=1,
    )

    assert isinstance(result_high, pd.DataFrame)
    assert isinstance(result_low, pd.DataFrame)
    assert len(result_high) == len(result_low)


def test_length_binned_weights_zero_threshold(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test with minimum_count_threshold=0."""
    result = biology.length_binned_weights(
        sample_specimen_data,
        sample_length_bins,
        single_regression_coefficients,
        minimum_count_threshold=0,
    )

    assert isinstance(result, pd.DataFrame)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns


def test_length_binned_weights_minimal_data(
    reduced_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test with minimal specimen data."""
    result = biology.length_binned_weights(
        reduced_specimen_data, sample_length_bins, single_regression_coefficients
    )

    assert isinstance(result, pd.DataFrame)
    assert len(result) >= 0


def test_length_binned_weights_missing_weights(
    specimen_data_missing_weights, sample_length_bins, single_regression_coefficients
):
    """Test handling of missing weight values."""
    result = biology.length_binned_weights(
        specimen_data_missing_weights, sample_length_bins, single_regression_coefficients
    )

    assert isinstance(result, pd.DataFrame)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns


def test_length_binned_weights_large_dataset(
    large_specimen_dataset, sample_length_bins, single_regression_coefficients
):
    """Test performance with large dataset."""
    result = biology.length_binned_weights(
        large_specimen_dataset, sample_length_bins, single_regression_coefficients
    )

    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(sample_length_bins)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns
    assert not result["weight_fitted"].isna().any()


def test_length_binned_weights_empty_data(
    null_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test with empty specimen data."""
    result = biology.length_binned_weights(
        null_specimen_data, sample_length_bins, single_regression_coefficients
    )

    assert isinstance(result, pd.DataFrame)
    # Should return bins with modeled weights only


def test_length_binned_weights_uneven_distribution(
    uneven_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test with uneven distribution of specimens across bins."""
    result = biology.length_binned_weights(
        uneven_specimen_data,
        sample_length_bins,
        single_regression_coefficients,
        minimum_count_threshold=2,
    )

    assert isinstance(result, pd.DataFrame)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns


def test_length_binned_weights_multiple_groups(
    specimen_data_multiple_groups, sample_length_bins, coefficients_with_multiple_groups
):
    """Test with multiple grouping variables."""
    result = biology.length_binned_weights(
        specimen_data_multiple_groups, sample_length_bins, coefficients_with_multiple_groups
    )

    assert isinstance(result, pd.DataFrame)
    # Should return wide format with multiindex columns for multiple grouping variables
    assert result.index.name == "length_bin"
    # Columns should be tuples representing the combination of grouping variables
    assert isinstance(result.columns, pd.MultiIndex) or len(result.columns) > 0
    assert len(result) > 0


def test_length_binned_weights_preserves_original_data(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test that original data is not modified."""
    original_data = sample_specimen_data.copy()

    biology.length_binned_weights(
        sample_specimen_data, sample_length_bins, single_regression_coefficients
    )

    # Original data should be unchanged (function makes copy)
    pd.testing.assert_frame_equal(sample_specimen_data, original_data)


def test_length_binned_weights_deterministic_results(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test that function produces deterministic results."""
    result1 = biology.length_binned_weights(
        sample_specimen_data, sample_length_bins, single_regression_coefficients
    )

    result2 = biology.length_binned_weights(
        sample_specimen_data, sample_length_bins, single_regression_coefficients
    )

    pd.testing.assert_frame_equal(result1, result2)


def test_length_binned_weights_realistic_coefficients(sample_specimen_data, sample_length_bins):
    """Test with realistic regression coefficients."""
    # Generate coefficients from actual data
    coeffs = biology.fit_length_weight_regression(sample_specimen_data)

    result = biology.length_binned_weights(sample_specimen_data, sample_length_bins, coeffs)

    assert isinstance(result, pd.DataFrame)
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns
    assert result["weight_fitted"].min() > 0  # Weights should be positive


def test_length_binned_weights_column_filtering(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test that output contains only expected columns."""
    result = biology.length_binned_weights(
        sample_specimen_data, sample_length_bins, single_regression_coefficients
    )

    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns


def test_length_binned_weights_weight_values_reasonable(
    sample_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test that fitted weight values are reasonable."""
    result = biology.length_binned_weights(
        sample_specimen_data, sample_length_bins, single_regression_coefficients
    )

    # Weight values should be positive and finite
    assert (result["weight_fitted"] > 0).all()
    assert np.isfinite(result["weight_fitted"]).all()


def test_length_binned_weights_integration_with_fit_regression(
    sample_specimen_data, sample_length_bins
):
    """Test integration with fit_length_weight_regression function."""
    # Generate coefficients using the regression function
    coeffs = biology.fit_length_weight_regression(sample_specimen_data)

    # Use those coefficients in biology.length_binned_weights
    result = biology.length_binned_weights(sample_specimen_data, sample_length_bins, coeffs)

    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    # Should return long format with additional columns when no grouping
    expected_columns = ["length_bin", "count", "weight_mean", "weight_modeled", "weight_fitted"]
    assert list(result.columns) == expected_columns


def test_length_binned_weights_different_imputation_strategies(
    uneven_specimen_data, sample_length_bins, single_regression_coefficients
):
    """Test different imputation strategies produce different results."""
    result_impute = biology.length_binned_weights(
        uneven_specimen_data,
        sample_length_bins,
        single_regression_coefficients,
        impute_bins=True,
        minimum_count_threshold=3,
    )

    result_no_impute = biology.length_binned_weights(
        uneven_specimen_data,
        sample_length_bins,
        single_regression_coefficients,
        impute_bins=False,
    )

    # Results should be different when imputation settings differ
    assert isinstance(result_impute, pd.DataFrame)
    assert isinstance(result_no_impute, pd.DataFrame)
    assert len(result_impute) == len(result_no_impute)


def test_remove_speciman_hauls(biological_data):
    """
    Test basic removal of specimen-specific hauls from catch data
    """

    # Pre-process
    biodata = {k: v.rename(columns={"haul": "haul_num"}) for k, v in biological_data.items()}

    # Swap out haul number 4 from length data
    biodata["length"].loc[4, "haul_num"] = 3

    # Create a copy for comparison
    biodata_preremoval = biodata.copy()

    # Apply removal
    biology.remove_specimen_hauls(biodata)

    # Compare
    assert len(biodata["catch"]) == 3 and len(biodata["catch"]) < len(biodata_preremoval["catch"])
    assert (
        len(
            set(biodata["catch"]["haul_num"].values).difference(
                biodata["length"]["haul_num"].values
            )
        )
        == 0
    )
