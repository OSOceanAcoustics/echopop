import numpy as np
import pandas as pd
import pytest

import echopop.nwfsc_feat.utils as echoutils

def test_create_centered_bins_dataframe_output(simple_bins):
    """Test basic functionality with DataFrame output."""
    result = echoutils.binned_distribution(simple_bins)
    
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["bin", "interval"]
    assert len(result) == len(simple_bins)
    np.testing.assert_array_equal(result["bin"].values, simple_bins)


def test_create_centered_bins_array_output(simple_bins):
    """Test functionality with array tuple output."""
    bins_out, centered_bins = echoutils.binned_distribution(simple_bins, return_dataframe=False)
    
    assert isinstance(bins_out, np.ndarray)
    assert isinstance(centered_bins, np.ndarray)
    np.testing.assert_array_equal(bins_out, simple_bins)
    assert len(centered_bins) == len(simple_bins) + 1


def test_binwidth_calculation(simple_bins):
    """Test that binwidth is calculated correctly."""
    bins_out, centered_bins = echoutils.binned_distribution(simple_bins, return_dataframe=False)
    
    # For evenly spaced bins [1,2,3,4,5], diff = [1,1,1,1], binwidth = 0.5
    expected_binwidth = 0.5
    actual_binwidth = centered_bins[1] - simple_bins[0]
    
    assert np.isclose(actual_binwidth, expected_binwidth)


def test_centered_bins_extend_range(simple_bins):
    """Test that centered bins properly extend the original range."""
    bins_out, centered_bins = echoutils.binned_distribution(simple_bins, return_dataframe=False)
    
    # First centered bin should be less than first original bin
    assert centered_bins[0] < simple_bins[0]
    # Last centered bin should be greater than last original bin
    assert centered_bins[-1] > simple_bins[-1]


def test_float_bins(float_bins):
    """Test functionality with floating point bins."""
    result = echoutils.binned_distribution(float_bins)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(float_bins)
    np.testing.assert_array_equal(result["bin"].values, float_bins)


def test_uneven_bins(uneven_bins):
    """Test functionality with unevenly spaced bins."""
    result = echoutils.binned_distribution(uneven_bins)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(uneven_bins)
    # Should still create valid intervals
    assert all(pd.notna(result["interval"]))


def test_minimal_bins(minimal_bins):
    """Test edge case with minimal two-element array."""
    result = echoutils.binned_distribution(minimal_bins)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 2
    np.testing.assert_array_equal(result["bin"].values, minimal_bins)


def test_age_bins_realistic(age_bins):
    """Test with realistic age bins used in biological data."""
    result = echoutils.binned_distribution(age_bins)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(age_bins)
    assert all(pd.notna(result["interval"]))


def test_length_bins_realistic(length_bins):
    """Test with realistic length bins used in biological data."""
    result = echoutils.binned_distribution(length_bins)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(length_bins)
    assert all(pd.notna(result["interval"]))


def test_interval_coverage(simple_bins):
    """Test that all original bins fall within created intervals."""
    result = echoutils.binned_distribution(simple_bins)
    
    for i, (bin_val, interval) in result.iterrows():
        assert interval.left <= bin_val <= interval.right

def test_input_type_conversion():
    """Test that function accepts list input and converts to numpy array."""
    bins_list = [1, 2, 3, 4, 5]
    result = echoutils.binned_distribution(bins_list)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(bins_list)


def test_large_array_performance(large_bins):
    """Test performance with large arrays."""
    result = echoutils.binned_distribution(large_bins)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(large_bins)
    assert all(pd.notna(result["interval"]))


def test_return_dataframe_parameter():
    """Test that return_dataframe parameter works correctly."""
    bins = np.array([1, 2, 3, 4, 5])
    
    # Test DataFrame return
    df_result = echoutils.binned_distribution(bins, return_dataframe=True)
    assert isinstance(df_result, pd.DataFrame)
    
    # Test tuple return
    tuple_result = echoutils.binned_distribution(bins, return_dataframe=False)
    assert isinstance(tuple_result, tuple)
    assert len(tuple_result) == 2


def test_binwidth_consistency_across_methods():
    """Test that binwidth calculation is consistent."""
    bins = np.array([1, 2, 3, 4, 5])
    
    # Get results from both methods
    df_result = echoutils.binned_distribution(bins, return_dataframe=True)
    bins_out, centered_bins = echoutils.binned_distribution(bins, return_dataframe=False)
    
    # Extract intervals from DataFrame and compare with centered_bins
    first_interval = df_result["interval"].iloc[0]
    expected_left = centered_bins[0]
    
    assert np.isclose(first_interval.left, expected_left)

def test_binify_single_dataframe_inplace(target_dataframe, numeric_bin_distribution):
    """Test binify with single DataFrame inplace operation."""
    original_shape = target_dataframe.shape
    result = echoutils.binify(target_dataframe, numeric_bin_distribution, "numeric_col", inplace=True)
    
    assert result is None
    assert "numeric_col_bin" in target_dataframe.columns
    assert target_dataframe.shape == (original_shape[0], original_shape[1] + 1)
    assert target_dataframe["numeric_col_bin"].notna().all()


def test_binify_single_dataframe_copy(target_dataframe, numeric_bin_distribution):
    """Test binify with single DataFrame returning copy."""
    original_columns = list(target_dataframe.columns)
    result = echoutils.binify(target_dataframe, numeric_bin_distribution, "numeric_col", inplace=False)
    
    assert isinstance(result, pd.DataFrame)
    assert "numeric_col_bin" in result.columns
    assert "numeric_col_bin" not in target_dataframe.columns
    assert list(target_dataframe.columns) == original_columns
    assert result["numeric_col_bin"].notna().all()


def test_binify_dictionary_inplace(mixed_dataframes_dict, numeric_bin_distribution):
    """Test binify with dictionary of DataFrames inplace."""
    result = echoutils.binify(mixed_dataframes_dict, numeric_bin_distribution, "numeric_col", inplace=True)
    
    assert result is None
    assert "numeric_col_bin" in mixed_dataframes_dict["target_data"].columns
    assert "numeric_col_bin" not in mixed_dataframes_dict["non_target_data"].columns
    assert "numeric_col_bin" not in mixed_dataframes_dict["partial_data"].columns


def test_binify_dictionary_copy(mixed_dataframes_dict, numeric_bin_distribution):
    """Test binify with dictionary of DataFrames returning copy."""
    original_target_cols = list(mixed_dataframes_dict["target_data"].columns)
    result = echoutils.binify(mixed_dataframes_dict, numeric_bin_distribution, "numeric_col", inplace=False)
    
    assert isinstance(result, dict)
    assert "numeric_col_bin" in result["target_data"].columns
    assert "numeric_col_bin" not in result["non_target_data"].columns
    assert "numeric_col_bin" not in result["partial_data"].columns
    
    # Original data unchanged
    assert list(mixed_dataframes_dict["target_data"].columns) == original_target_cols
    assert "numeric_col_bin" not in mixed_dataframes_dict["target_data"].columns


def test_binify_secondary_column(mixed_dataframes_dict, secondary_bin_distribution):
    """Test binify with secondary column instead of primary."""
    echoutils.binify(mixed_dataframes_dict, secondary_bin_distribution, "secondary_col", inplace=True)
    
    assert "secondary_col_bin" in mixed_dataframes_dict["target_data"].columns
    assert "secondary_col_bin" in mixed_dataframes_dict["partial_data"].columns
    assert "secondary_col_bin" not in mixed_dataframes_dict["non_target_data"].columns


def test_binify_missing_column_single_df(non_target_dataframe, numeric_bin_distribution):
    """Test binify with single DataFrame missing target column."""
    original_columns = list(non_target_dataframe.columns)
    result = echoutils.binify(non_target_dataframe, numeric_bin_distribution, "numeric_col", inplace=False)
    
    # Should return copy of original data
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == original_columns
    assert "numeric_col_bin" not in result.columns


def test_binify_missing_column_single_df_inplace(non_target_dataframe, numeric_bin_distribution):
    """Test binify inplace with single DataFrame missing target column."""
    original_columns = list(non_target_dataframe.columns)
    result = echoutils.binify(non_target_dataframe, numeric_bin_distribution, "numeric_col", inplace=True)
    
    assert result is None
    assert list(non_target_dataframe.columns) == original_columns
    assert "numeric_col_bin" not in non_target_dataframe.columns


def test_binify_empty_dataframe(empty_dataframe, numeric_bin_distribution):
    """Test binify with empty DataFrame."""
    result = echoutils.binify(empty_dataframe, numeric_bin_distribution, "numeric_col", inplace=False)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0


def test_binify_single_row(single_row_dataframe, numeric_bin_distribution):
    """Test binify with single-row DataFrame."""
    result = echoutils.binify(single_row_dataframe, numeric_bin_distribution, "numeric_col", inplace=False)
    
    assert isinstance(result, pd.DataFrame)
    assert "numeric_col_bin" in result.columns
    assert len(result) == 1


def test_binify_large_dataset(large_dataframe, numeric_bin_distribution):
    """Test binify performance with large dataset."""
    result = echoutils.binify(large_dataframe, numeric_bin_distribution, "numeric_col", inplace=False)
    
    assert isinstance(result, pd.DataFrame)
    assert "numeric_col_bin" in result.columns
    assert len(result) == len(large_dataframe)


def test_binify_mixed_objects_dict(mixed_objects_dict, numeric_bin_distribution):
    """Test binify with dictionary containing non-DataFrame objects."""
    result = echoutils.binify(mixed_objects_dict, numeric_bin_distribution, "numeric_col", inplace=False)
    
    assert isinstance(result, dict)
    assert "numeric_col_bin" in result["dataframe_1"].columns
    assert "numeric_col_bin" not in result["dataframe_2"].columns
    assert result["metadata"] == {"source": "test", "version": 1.0}
    assert result["config"] == [1, 2, 3]


def test_binify_invalid_bin_distribution(target_dataframe, invalid_bin_distribution):
    """Test error handling for bin distribution without interval column."""
    with pytest.raises(KeyError, match="interval"):
        echoutils.binify(target_dataframe, invalid_bin_distribution, "numeric_col")


def test_binify_non_interval_distribution(target_dataframe, non_interval_bin_distribution):
    """Test error handling for non-interval bin distribution."""
    with pytest.raises(ValueError, match="must contain pandas Interval objects"):
        echoutils.binify(target_dataframe, non_interval_bin_distribution, "numeric_col")


def test_binify_invalid_data_type(numeric_bin_distribution):
    """Test error handling for invalid data type."""
    with pytest.raises(TypeError, match="data must be DataFrame or dict of DataFrames"):
        echoutils.binify("invalid_data", numeric_bin_distribution, "numeric_col")


def test_binify_bin_column_name_format(target_dataframe, numeric_bin_distribution):
    """Test that bin column name is formatted correctly."""
    echoutils.binify(target_dataframe, numeric_bin_distribution, "numeric_col", inplace=True)
    assert "numeric_col_bin" in target_dataframe.columns
    
    # Test with different column name
    echoutils.binify(target_dataframe, numeric_bin_distribution, "secondary_col", inplace=True)
    assert "secondary_col_bin" in target_dataframe.columns


def test_binify_preserves_original_data(target_dataframe, numeric_bin_distribution):
    """Test that original data structure is preserved."""
    original_index = target_dataframe.index.copy()
    original_dtypes = target_dataframe.dtypes.copy()
    
    echoutils.binify(target_dataframe, numeric_bin_distribution, "numeric_col", inplace=True)
    
    # Check that original columns and index are preserved
    pd.testing.assert_index_equal(target_dataframe.index, original_index)
    for col in original_dtypes.index:
        assert target_dataframe[col].dtype == original_dtypes[col]


def test_binify_deterministic_results(target_dataframe, numeric_bin_distribution):
    """Test that binify produces deterministic results."""
    result1 = echoutils.binify(target_dataframe.copy(), numeric_bin_distribution, "numeric_col", inplace=False)
    result2 = echoutils.binify(target_dataframe.copy(), numeric_bin_distribution, "numeric_col", inplace=False)
    
    pd.testing.assert_frame_equal(result1, result2)


def test_binify_handles_na_values(dataframe_with_na, numeric_bin_distribution):
    """Test binify behavior with NaN values in target column."""
    result = echoutils.binify(dataframe_with_na, numeric_bin_distribution, "numeric_col", inplace=False)
    
    assert "numeric_col_bin" in result.columns
    assert result["numeric_col_bin"].isna().sum() == 2  # Two NaN values
    assert result["numeric_col_bin"].notna().sum() == 3  # Three valid values


def test_binify_interval_assignment(target_dataframe, numeric_bin_distribution):
    """Test that intervals are assigned correctly."""
    result = echoutils.binify(target_dataframe, numeric_bin_distribution, "numeric_col", inplace=False)
    
    # Check that binned values are within expected intervals
    for i, numeric_val in enumerate(result["numeric_col"]):
        if pd.notna(numeric_val) and pd.notna(result["numeric_col_bin"].iloc[i]):
            interval = result["numeric_col_bin"].iloc[i]
            assert interval.left < numeric_val <= interval.right


def test_binify_multiple_columns_sequentially(target_dataframe, numeric_bin_distribution, secondary_bin_distribution):
    """Test binify with multiple columns applied sequentially."""
    echoutils.binify(target_dataframe, numeric_bin_distribution, "numeric_col", inplace=True)
    echoutils.binify(target_dataframe, secondary_bin_distribution, "secondary_col", inplace=True)
    
    assert "numeric_col_bin" in target_dataframe.columns
    assert "secondary_col_bin" in target_dataframe.columns
    assert target_dataframe["numeric_col_bin"].notna().all()
    assert target_dataframe["secondary_col_bin"].notna().all()


def test_binify_different_bin_distributions(target_dataframe):
    """Test binify with different bin distribution sizes."""
    # Create small bin distribution
    small_bins = np.array([10, 20, 30])
    small_binwidth = np.mean(np.diff(small_bins) / 2.0)
    small_centered = np.concatenate([[small_bins[0] - small_binwidth], small_bins + small_binwidth])
    small_intervals = pd.cut(small_bins, small_centered)
    small_dist = pd.DataFrame({"bin": small_bins, "interval": small_intervals})
    
    result = echoutils.binify(target_dataframe, small_dist, "numeric_col", inplace=False)
    
    assert "numeric_col_bin" in result.columns
    assert result["numeric_col_bin"].notna().any()