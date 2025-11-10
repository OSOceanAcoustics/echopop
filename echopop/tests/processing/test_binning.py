import numpy as np
import pandas as pd
import pytest

import echopop.utils as echoutils


def test_create_centered_bins_dataframe_output(simple_bins):
    """Test basic functionality with DataFrame output."""
    result = echoutils.binned_distribution(simple_bins)

    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["bin", "interval"]
    assert len(result) == len(simple_bins)
    np.testing.assert_array_equal(result["bin"].values, simple_bins)


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


def test_binify_single_dataframe_inplace(target_dataframe, numeric_bins):
    """Test binify with single DataFrame inplace operation."""
    original_shape = target_dataframe.shape
    result = echoutils.binify(target_dataframe, numeric_bins, "numeric_col")

    assert result is None
    assert "numeric_col_bin" in target_dataframe.columns
    assert target_dataframe.shape == (original_shape[0], original_shape[1] + 1)
    assert target_dataframe["numeric_col_bin"].notna().all()


def test_binify_dictionary_inplace(mixed_dataframes_dict, numeric_bins):
    """Test binify with dictionary of DataFrames inplace."""
    result = echoutils.binify(mixed_dataframes_dict, numeric_bins, "numeric_col")

    assert result is None
    assert "numeric_col_bin" in mixed_dataframes_dict["target_data"].columns
    assert "numeric_col_bin" not in mixed_dataframes_dict["non_target_data"].columns
    assert "numeric_col_bin" not in mixed_dataframes_dict["partial_data"].columns


def test_binify_secondary_column(mixed_dataframes_dict, secondary_bins):
    """Test binify with secondary column instead of primary."""
    echoutils.binify(mixed_dataframes_dict, secondary_bins, "secondary_col")

    assert "secondary_col_bin" in mixed_dataframes_dict["target_data"].columns
    assert "secondary_col_bin" in mixed_dataframes_dict["partial_data"].columns
    assert "secondary_col_bin" not in mixed_dataframes_dict["non_target_data"].columns


def test_binify_missing_column_single_df(non_target_dataframe, numeric_bins):
    """Test binify with single DataFrame missing target column."""
    original_columns = list(non_target_dataframe.columns)
    original_data = non_target_dataframe.copy()

    result = echoutils.binify(non_target_dataframe, numeric_bins, "numeric_col")

    # Should modify dataframe in place but not add the bin column since target column is missing
    assert result is None
    assert list(non_target_dataframe.columns) == original_columns
    assert "numeric_col_bin" not in non_target_dataframe.columns
    # Data should remain unchanged
    pd.testing.assert_frame_equal(non_target_dataframe, original_data)


def test_binify_empty_dataframe(empty_dataframe, numeric_bins):
    """Test binify with empty DataFrame."""
    original_data = empty_dataframe.copy()
    result = echoutils.binify(empty_dataframe, numeric_bins, "numeric_col")

    assert result is None
    assert len(empty_dataframe) == 0
    # Should remain unchanged since no target column exists
    pd.testing.assert_frame_equal(empty_dataframe, original_data)


def test_binify_single_row(single_row_dataframe, numeric_bins):
    """Test binify with single-row DataFrame."""
    result = echoutils.binify(single_row_dataframe, numeric_bins, "numeric_col")

    assert result is None
    assert "numeric_col_bin" in single_row_dataframe.columns
    assert len(single_row_dataframe) == 1


def test_binify_large_dataset(large_dataframe, numeric_bins):
    """Test binify performance with large dataset."""
    result = echoutils.binify(large_dataframe, numeric_bins, "numeric_col")

    assert result is None
    assert "numeric_col_bin" in large_dataframe.columns
    assert len(large_dataframe) == 1000  # Original size


def test_binify_mixed_objects_dict(mixed_objects_dict, numeric_bins):
    """Test binify with dictionary containing non-DataFrame objects."""
    result = echoutils.binify(mixed_objects_dict, numeric_bins, "numeric_col")

    assert result is None
    assert "numeric_col_bin" in mixed_objects_dict["dataframe_1"].columns
    assert "numeric_col_bin" not in mixed_objects_dict["dataframe_2"].columns
    # Non-DataFrame objects should remain unchanged
    assert mixed_objects_dict["metadata"] == {"source": "test", "version": 1.0}
    assert mixed_objects_dict["config"] == [1, 2, 3]


def test_binify_invalid_bins(target_dataframe):
    """Test error handling for invalid bins array."""
    invalid_bins = np.array([10])  # Single element
    with pytest.raises(Exception):  # Will fail in binned_distribution
        echoutils.binify(target_dataframe, invalid_bins, "numeric_col")


def test_binify_invalid_data_type(numeric_bins):
    """Test error handling for invalid data type."""
    with pytest.raises(TypeError, match="data must be DataFrame or dict of DataFrames"):
        echoutils.binify("invalid_data", numeric_bins, "numeric_col")


def test_binify_bin_column_name_format(target_dataframe, numeric_bins, secondary_bins):
    """Test that bin column name is formatted correctly."""
    echoutils.binify(target_dataframe, numeric_bins, "numeric_col")
    assert "numeric_col_bin" in target_dataframe.columns

    # Test with different column name
    echoutils.binify(target_dataframe, secondary_bins, "secondary_col")
    assert "secondary_col_bin" in target_dataframe.columns


def test_binify_preserves_original_data(target_dataframe, numeric_bins):
    """Test that original data structure is preserved."""
    original_index = target_dataframe.index.copy()
    original_dtypes = target_dataframe.dtypes.copy()

    echoutils.binify(target_dataframe, numeric_bins, "numeric_col")

    # Check that original columns and index are preserved
    pd.testing.assert_index_equal(target_dataframe.index, original_index)
    for col in original_dtypes.index:
        assert target_dataframe[col].dtype == original_dtypes[col]


def test_binify_deterministic_results(target_dataframe, numeric_bins):
    """Test that binify produces deterministic results."""
    df_copy1 = target_dataframe.copy()
    df_copy2 = target_dataframe.copy()

    echoutils.binify(df_copy1, numeric_bins, "numeric_col")
    echoutils.binify(df_copy2, numeric_bins, "numeric_col")

    pd.testing.assert_frame_equal(df_copy1, df_copy2)


def test_binify_handles_na_values(dataframe_with_na, numeric_bins):
    """Test binify behavior with NaN values in target column."""
    result = echoutils.binify(dataframe_with_na, numeric_bins, "numeric_col")

    assert result is None
    assert "numeric_col_bin" in dataframe_with_na.columns
    assert dataframe_with_na["numeric_col_bin"].isna().sum() == 2  # Two NaN values
    assert dataframe_with_na["numeric_col_bin"].notna().sum() == 3  # Three valid values


def test_binify_interval_assignment(target_dataframe, numeric_bins):
    """Test that intervals are assigned correctly."""
    original_data = target_dataframe.copy()
    echoutils.binify(target_dataframe, numeric_bins, "numeric_col")

    # Check that binned values are within expected intervals
    for i, numeric_val in enumerate(original_data["numeric_col"]):
        if pd.notna(numeric_val) and pd.notna(target_dataframe["numeric_col_bin"].iloc[i]):
            interval = target_dataframe["numeric_col_bin"].iloc[i]
            assert interval.left < numeric_val <= interval.right


def test_binify_multiple_columns_sequentially(target_dataframe, numeric_bins, secondary_bins):
    """Test binify with multiple columns applied sequentially."""
    echoutils.binify(target_dataframe, numeric_bins, "numeric_col")
    echoutils.binify(target_dataframe, secondary_bins, "secondary_col")

    assert "numeric_col_bin" in target_dataframe.columns
    assert "secondary_col_bin" in target_dataframe.columns
    assert target_dataframe["numeric_col_bin"].notna().all()
    assert target_dataframe["secondary_col_bin"].notna().all()


def test_binify_different_bin_sizes(target_dataframe, small_bins):
    """Test binify with different bin array sizes."""
    result = echoutils.binify(target_dataframe, small_bins, "numeric_col")

    assert result is None
    assert "numeric_col_bin" in target_dataframe.columns
    assert target_dataframe["numeric_col_bin"].notna().any()


def test_binify_with_linspace_bins(target_dataframe):
    """Test binify with numpy linspace bins."""
    age_bins = np.linspace(start=1.0, stop=22.0, num=22)
    # Add age column to test data
    target_dataframe["age"] = [5, 8, 12, 15, 18, 20]

    result = echoutils.binify(target_dataframe, age_bins, "age")

    assert result is None
    assert "age_bin" in target_dataframe.columns
    assert target_dataframe["age_bin"].notna().all()
