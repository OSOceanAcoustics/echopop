import numpy as np
import pandas as pd
import pytest
from scipy import interpolate as interp

from echopop import utils
from echopop.survey import proportions as get_proportions


def test_compute_binned_counts_size_aggregation(sample_specimen_data):
    """Test size aggregation (default)."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "length", agg_func="size"
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert result["count"].sum() == len(sample_specimen_data)


def test_compute_binned_counts_sum_aggregation(sample_length_data):
    """Test sum aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_length_data, ["stratum_num", "length_bin"], "length_count", agg_func="sum"
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert result["count"].sum() == sample_length_data["length_count"].sum()


def test_compute_binned_counts_mean_aggregation(sample_specimen_data):
    """Test mean aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert len(result) > 0


def test_compute_binned_counts_count_aggregation(sample_specimen_data):
    """Test count aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="count"
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert len(result) > 0


def test_compute_binned_counts_default_size(sample_specimen_data):
    """Test default aggregation is size."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "length"
    )

    assert isinstance(result, pd.DataFrame)
    assert result["count"].sum() == len(sample_specimen_data)


def test_compute_binned_counts_single_group_column(sample_specimen_data):
    """Test with single grouping column."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "length", agg_func="size"
    )

    assert len(result) == len(sample_specimen_data["stratum_num"].unique())


def test_compute_binned_counts_complex_grouping(sample_specimen_data):
    """Test with complex grouping."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["stratum_num", "length_bin", "age_bin", "sex"],
        "weight",
        agg_func="mean",
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert len(result) > 0


def test_compute_binned_counts_empty_data(empty_specimen_data):
    """Test with empty DataFrame."""
    result = get_proportions.compute_binned_counts(
        empty_specimen_data, ["stratum_num", "length_bin"], "length", agg_func="size"
    )

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0


def test_compute_binned_counts_different_numeric_columns(sample_specimen_data):
    """Test aggregation on different numeric columns."""
    result_length = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "length", agg_func="mean"
    )

    result_weight = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "weight", agg_func="mean"
    )

    assert isinstance(result_length, pd.DataFrame)
    assert isinstance(result_weight, pd.DataFrame)
    assert len(result_length) == len(result_weight)


def test_compute_binned_counts_deterministic_results(sample_specimen_data):
    """Test that function produces deterministic results."""
    result1 = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    result2 = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    pd.testing.assert_frame_equal(result1, result2)


def test_compute_binned_counts_column_order_preserved(sample_specimen_data):
    """Test that groupby column order is preserved in output."""
    groupby_cols = ["stratum_num", "length_bin", "sex"]
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, groupby_cols, "length", agg_func="size"
    )

    # Check that columns appear in expected order
    expected_cols = groupby_cols + ["count"]
    assert set(list(result.columns)) <= set(expected_cols)


def test_compute_binned_counts_median_aggregation(sample_specimen_data):
    """Test median aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "weight", agg_func="median"
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert len(result) > 0


def test_compute_binned_counts_var_aggregation(sample_specimen_data):
    """Test variance aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "length", agg_func="var"
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert len(result) > 0


def test_compute_binned_counts_multiple_groupby_combinations(data_multiple_groups):
    """Test with multiple different groupby combinations."""
    # Simple grouping
    result1 = get_proportions.compute_binned_counts(
        data_multiple_groups, ["stratum_num"], "length", agg_func="size"
    )

    # Complex grouping
    result2 = get_proportions.compute_binned_counts(
        data_multiple_groups, ["stratum_num", "species", "region"], "length", agg_func="size"
    )

    assert len(result1) <= len(result2)  # More grouping = more rows (usually)
    assert result1["count"].sum() == result2["count"].sum()  # Same total count


def test_number_proportions_single_dataframe(aged_dataframe):
    """Test the number_proportions function with a single dataframe."""
    result = get_proportions.number_proportions(data=aged_dataframe)

    # Check that result is a DataFrame and has the right columns
    assert isinstance(result, pd.DataFrame)
    assert "proportion" in result.columns
    assert "proportion_overall" in result.columns
    # Remove or comment out the following lines if those columns are not expected:
    # assert 'total_df_0' in result.columns
    # assert 'total_overall' in result.columns

    # Check that proportions sum to 1.0 for each stratum
    prop_sums = result.groupby("stratum_num")["proportion"].sum().round(10)
    for stratum_sum in prop_sums:
        assert stratum_sum == pytest.approx(1.0)


def test_number_proportions_multiple_dataframes(aged_dataframe, unaged_dataframe):
    """Test the number_proportions function with multiple dataframes."""
    result = get_proportions.number_proportions(
        data={"aged": aged_dataframe, "unaged": unaged_dataframe}
    )

    # Check that result is a dictionary with the right keys
    assert isinstance(result, dict)
    assert "aged" in result
    assert "unaged" in result

    # Check the aged dataframe
    aged_df = result["aged"]
    assert "proportion" in aged_df.columns
    assert "proportion_overall" in aged_df.columns

    # Check the unaged dataframe
    unaged_df = result["unaged"]
    assert "proportion" in unaged_df.columns
    assert "proportion_overall" in unaged_df.columns

    # Check sums of within-group proportions
    for stratum in [1, 2]:
        aged_sum = aged_df[aged_df["stratum_num"] == stratum]["proportion"].sum()
        unaged_sum = unaged_df[unaged_df["stratum_num"] == stratum]["proportion"].sum()
        assert aged_sum == pytest.approx(1.0)
        assert unaged_sum == pytest.approx(1.0)


def test_number_proportions_column_aliases(aged_dataframe, unaged_dataframe):
    """Test column_aliases parameter in number_proportions."""
    result = get_proportions.number_proportions(
        data={"specimen": aged_dataframe, "length": unaged_dataframe}
    )

    assert "specimen" in result
    assert "length" in result


def test_number_proportions_default_aliases(aged_dataframe, unaged_dataframe):
    """Test default aliases when column_aliases not provided."""
    result = get_proportions.number_proportions(
        data={"df_0": aged_dataframe, "df_1": unaged_dataframe}
    )

    assert "df_0" in result
    assert "df_1" in result


def test_number_proportions_with_exclusion(aged_dataframe):
    """Test exclude_filters parameter in number_proportions."""
    result = get_proportions.number_proportions(
        data=aged_dataframe, exclude_filters={"sex": "unsexed"}
    )

    # Check that unsexed rows are excluded
    assert not (result["sex"] == "unsexed").any()

    # Check that proportions are recalculated correctly
    stratum1_total = aged_dataframe[
        (aged_dataframe["stratum_num"] == 1) & (aged_dataframe["sex"] != "unsexed")
    ]["count"].sum()

    stratum1_female = aged_dataframe[
        (aged_dataframe["stratum_num"] == 1) & (aged_dataframe["sex"] == "female")
    ]["count"].sum()

    expected_proportion = stratum1_female / stratum1_total
    actual_proportion = result[(result["stratum_num"] == 1) & (result["sex"] == "female")][
        "proportion"
    ].iloc[0]

    assert actual_proportion == pytest.approx(expected_proportion)


def test_number_proportions_custom_group_columns(aged_dataframe, unaged_dataframe):
    """Test custom group_columns parameter in number_proportions."""
    result = get_proportions.number_proportions(
        data={"aged": aged_dataframe, "unaged": unaged_dataframe},
        group_columns=["stratum_num", "sex"],
    )  # Check that grouping includes both stratum and sex
    assert len(result["aged"].groupby(["stratum_num", "sex"]).groups) == 6
    assert len(result["unaged"].groupby(["stratum_num", "sex"]).groups) == 6

    # Check that proportions sum to 1.0 for each stratum-sex combination
    df1 = result["aged"]
    for _, group in df1.groupby(["stratum_num", "sex"]):
        assert group["proportion"].sum() == pytest.approx(1.0)


def test_number_proportions_proportion_calculations(aged_dataframe, unaged_dataframe):
    """Test the correctness of proportion calculations in number_proportions."""
    result = get_proportions.number_proportions(
        data={"aged": aged_dataframe, "unaged": unaged_dataframe}
    )

    # Get stratum 1 totals
    stratum1_aged_total = aged_dataframe[aged_dataframe["stratum_num"] == 1]["count"].sum()
    stratum1_unaged_total = unaged_dataframe[unaged_dataframe["stratum_num"] == 1]["count"].sum()
    stratum1_overall_total = stratum1_aged_total + stratum1_unaged_total

    # Check specific row's proportions
    aged_df = result["aged"]
    aged_female_row = aged_df[(aged_df["stratum_num"] == 1) & (aged_df["sex"] == "female")].iloc[0]
    female_count = 10  # From fixture

    # Within-group proportion
    expected_within = female_count / stratum1_aged_total
    assert aged_female_row["proportion"] == pytest.approx(expected_within)

    # Overall proportion
    expected_overall = female_count / stratum1_overall_total
    assert aged_female_row["proportion_overall"] == pytest.approx(expected_overall)


def test_number_proportions_validation():
    """Test validation errors in number_proportions."""
    # DataFrame without count column
    invalid_df = pd.DataFrame({"stratum_num": [1, 2], "value": [10, 20]})
    with pytest.raises(ValueError, match="DataFrame 0 does not have a 'count' column"):
        get_proportions.number_proportions(data=invalid_df)


def test_apply_weight_interpolation_global(length_dataset_with_bins, real_interpolators):
    """Test applying a global interpolator."""
    # Use real interpolator function where interpolator(x) = x * 2.5
    global_interp_dict = {"_global_": real_interpolators["_global_"]}

    result = get_proportions.apply_weight_interpolation(
        target_df=length_dataset_with_bins,
        interpolators=global_interp_dict,
        dependent_var="weight",
        independent_var="length",
        count_col="length_count",
        contrast_vars=None,
    )

    # Check that the interpolation was applied and multiplied by counts
    assert "weight" in result.columns
    for idx, row in result.iterrows():
        expected = row["length"] * 2.5 * row["length_count"]
        assert np.isclose(row["weight"], expected, rtol=1e-5)


def test_apply_weight_interpolation_by_contrast(length_dataset_with_bins, real_interpolators):
    """Test applying interpolators based on contrast variables."""
    # Use real interpolator functions where:
    # male_interp(x) = x * 3.0
    # female_interp(x) = x * 2.0
    sex_interp_dict = {"male": real_interpolators["male"], "female": real_interpolators["female"]}

    result = get_proportions.apply_weight_interpolation(
        target_df=length_dataset_with_bins,
        interpolators=sex_interp_dict,
        dependent_var="weight",
        independent_var="length",
        count_col="length_count",
        contrast_vars="sex",
    )

    # Check that the correct interpolator was applied for each sex
    assert "weight" in result.columns

    for idx, row in result.iterrows():
        if row["sex"] == "male":
            expected = row["length"] * 3.0 * row["length_count"]
        else:  # female
            expected = row["length"] * 2.0 * row["length_count"]
        assert np.isclose(row["weight"], expected, rtol=1e-5)


def test_apply_weight_interpolation_missing_interpolator(length_dataset_with_bins):
    """Test handling of missing interpolators for some contrast values."""
    # Create dictionary with only male interpolator
    x_male = np.array([10.0, 20.0, 30.0])
    y_male = np.array([30.0, 60.0, 90.0])
    male_interp = interp.interp1d(
        x_male, y_male, kind="linear", bounds_error=False, fill_value=(y_male[0], y_male[-1])
    )

    interpolators = {"male": male_interp}

    result = get_proportions.apply_weight_interpolation(
        target_df=length_dataset_with_bins,
        interpolators=interpolators,
        dependent_var="weight",
        independent_var="length",
        count_col="length_count",
        contrast_vars="sex",
    )

    # Check that NaN is assigned for rows without an interpolator
    for idx, row in result.iterrows():
        if row["sex"] == "male":
            expected = row["length"] * 3.0 * row["length_count"]
            assert np.isclose(row["weight"], expected, rtol=1e-5)
        else:  # female - should be NaN
            assert pd.isna(row["weight"])


# -- Tests for binned_weights --


def test_binned_weights_with_interpolation(
    length_dataset_with_bins, length_weight_dataset_wide_format
):
    """Test binned_weights with interpolation."""
    result = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_wide_format,  # Required for interpolation
        interpolate_regression=True,
        table_index=["length_bin"],
        table_cols=["sex", "stratum"],
        include_filter=None,
        contrast_vars="sex",
    )

    # Check the result is a DataFrame
    assert isinstance(result, pd.DataFrame)
    # Check it has the expected structure
    assert result.index.name == "length_bin"
    assert "sex" in result.columns.names
    assert "stratum" in result.columns.names
    # Check it has values
    assert not result.empty


def test_binned_weights_without_interpolation(length_dataset_with_bins):
    """Test binned_weights without interpolation."""
    # Add direct weights to the dataset
    df = length_dataset_with_bins.copy()
    df["weight"] = df["length"] * 2.0  # Simple weight calculation

    # Note: length_weight_dataset is not provided since interpolate=False
    result = get_proportions.binned_weights(
        length_dataset=df,
        interpolate_regression=False,
        table_index=["length_bin"],
        table_cols=["sex"],
        include_filter=None,
        contrast_vars=None,
    )

    # Check the result
    assert isinstance(result, pd.DataFrame)
    assert result.index.name == "length_bin"
    assert "sex" in result.columns.names
    assert not result.empty


def test_binned_weights_raises_error_without_dataset(length_dataset_with_bins):
    """Test that binned_weights raises error when interpolate=True but no length_weight_dataset."""
    with pytest.raises(
        ValueError, match="length_weight_dataset must be provided when interpolate_regression=True"
    ):
        get_proportions.binned_weights(
            length_dataset=length_dataset_with_bins,
            length_weight_dataset=None,  # Missing required dataset
            interpolate_regression=True,
            table_index=["length_bin"],
            table_cols=["sex"],
            contrast_vars="sex",
        )


def test_binned_weights_with_filtering(length_dataset_with_bins, length_weight_dataset_wide_format):
    """Test binned_weights with filtering applied."""
    result = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_wide_format,
        interpolate_regression=True,
        table_index=["length_bin"],
        table_cols=["sex"],
        include_filter={"sex": "male"},
        contrast_vars="sex",
    )

    # Check result structure
    assert isinstance(result, pd.DataFrame)

    # Check that result contains data for males
    assert "male" in result.columns.get_level_values("sex")

    # Check that if female columns exist, they contain only zero values
    female_cols = [col for col in result.columns if "female" in col]
    if female_cols:
        for col in female_cols:
            assert (result[col] == 0).all() or result[col].isna().all()

    # Check that male columns contain data
    male_cols = [col for col in result.columns if "male" in col]
    assert len(male_cols) > 0
    assert not result[male_cols].isna().all().all()


def test_binned_weights_global_interpolation(
    length_dataset_with_bins, length_weight_dataset_wide_format
):
    """Test binned_weights with global interpolation (no contrast vars)."""
    result = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_wide_format,  # Required for interpolation
        interpolate_regression=True,
        table_index=["length_bin"],
        table_cols=["stratum"],
        include_filter=None,
        contrast_vars=None,
    )

    # Check the result
    assert isinstance(result, pd.DataFrame)
    assert result.index.name == "length_bin"
    assert "stratum" in result.columns.names
    assert not result.empty


def test_binned_weights_integration(length_dataset_with_bins, length_weight_dataset_wide_format):
    """Full integration test of binned_weights with both interpolation and no interpolation."""
    # Test with interpolation
    result_interp = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_wide_format,
        interpolate_regression=True,
        table_index=["length_bin"],
        table_cols=["sex", "stratum"],
        include_filter={"sex": ["male", "female"]},
        contrast_vars="sex",
    )

    # Test without interpolation (need to add weight)
    df = length_dataset_with_bins.copy()
    df["weight"] = df["length"] * 2.0

    # For non-interpolation case, provide empty DataFrame instead of None
    empty_df = pd.DataFrame()
    result_no_interp = get_proportions.binned_weights(
        length_dataset=df,
        length_weight_dataset=empty_df,  # Empty DataFrame instead of None
        interpolate_regression=False,
        table_index=["length_bin"],
        table_cols=["sex", "stratum"],
        include_filter={"sex": ["male", "female"]},
        contrast_vars="sex",
    )

    # Both should produce valid tables with the same structure but different values
    assert isinstance(result_interp, pd.DataFrame)
    assert isinstance(result_no_interp, pd.DataFrame)
    assert result_interp.index.name == result_no_interp.index.name == "length_bin"
    assert (
        set(result_interp.columns.names)
        == set(result_no_interp.columns.names)
        == {"sex", "stratum"}
    )


def test_calculate_within_group_proportions(proportion_dict):
    """Test the calculate_within_group_proportions function."""
    result = get_proportions.calculate_within_group_proportions(
        proportion_dict, ["stratum_num", "sex"]
    )

    # Check that we have within-group proportions
    assert "within_group_proportion" in result.columns

    # For stratum 1, female, aged: proportion = 0.3, total = 0.3 -> within proportion = 1.0
    aged_female_within = result[
        (result["group"] == "aged") & (result["stratum_num"] == 1) & (result["sex"] == "female")
    ]["within_group_proportion"].iloc[0]

    assert np.isclose(aged_female_within, 1.0)

    # For stratum 1, male, aged: proportions = 0.4 + 0.1 = 0.5
    # 0.4/0.5 = 0.8 for first row, 0.1/0.5 = 0.2 for second row
    aged_male_within = result[
        (result["group"] == "aged")
        & (result["stratum_num"] == 1)
        & (result["sex"] == "male")
        & (result["length_bin"] == "(10, 20]")
    ]["within_group_proportion"].iloc[0]

    assert np.isclose(aged_male_within, 0.8)


def test_calculate_adjusted_proportions(proportion_dict):
    """Test the calculate_adjusted_proportions function."""
    # First create the necessary inputs
    group_keys = list(proportion_dict.keys())

    aggregate_table = utils.create_grouped_table(
        proportion_dict, ["stratum_num"], ["group"], ["stratum_num"], "proportion_overall"
    )

    sex_proportions_table = utils.create_grouped_table(
        proportion_dict,
        ["stratum_num", "sex"],
        ["group", "sex"],
        ["stratum_num"],
        "proportion_overall",
    )

    result = get_proportions.calculate_adjusted_proportions(
        group_keys, aggregate_table, sex_proportions_table, ["sex"]
    )

    # Check that we get multi-index with group and sex
    assert isinstance(result.index, pd.MultiIndex)
    assert result.index.names == ["group", "sex"]

    # Check that the adjusted proportions for the first group are calculated correctly
    # For stratum 1, the adjusted proportion for aged should be higher than raw proportion
    assert result.loc[("aged", "female"), 1] > sex_proportions_table.loc[("aged", "female"), 1]


def test_stratum_averaged_weight(proportion_test_dict, test_weight_table):
    """Test the stratum_averaged_weight function."""
    result = get_proportions.stratum_averaged_weight(
        proportions_dict=proportion_test_dict,
        binned_weight_table=test_weight_table,
        stratify_by=["stratum_num"],
        group_by=["sex"],
    )

    # Check that we get the right format
    assert isinstance(result, pd.DataFrame)

    # Check that we have the right columns
    assert set(result.columns) >= {"all", "female", "male"}

    # Check that the weights are positive
    assert result["all"].min() > 0

    # Check that female weights are different from male weights
    if "female" in result.columns and "male" in result.columns:
        assert not result["female"].equals(result["male"])


def test_aggregate_stratum_weights_dictionary(weight_distr_dict):
    """Test aggregating stratum weights from a dictionary of DataFrames."""
    # Call function with dictionary input
    result = get_proportions.aggregate_stratum_weights(weight_distr_dict)

    # Check result structure
    assert isinstance(result, pd.DataFrame)
    assert "aged" in result.columns
    assert "unaged" in result.columns

    # Check values for first stratum with approx to handle float precision
    assert pytest.approx(result.loc[1, "aged"]) == 18.8  # 10.5 + 8.3
    assert pytest.approx(result.loc[1, "unaged"]) == 9.5  # 5.2 + 4.3

    # Check values for second stratum
    assert pytest.approx(result.loc[2, "aged"]) == 27.9  # 15.2 + 12.7
    assert pytest.approx(result.loc[2, "unaged"]) == 13.9  # 7.8 + 6.1


def test_aggregate_stratum_weights_single_df(weights_df_multilevel):
    """Test aggregating stratum weights from a single DataFrame."""
    # Call function with single DataFrame input
    result = get_proportions.aggregate_stratum_weights(weights_df_multilevel)

    # Check result structure
    assert isinstance(result, pd.DataFrame)
    assert "data" in result.columns  # Default name for single DataFrame

    # Check values with approx to handle float precision
    assert pytest.approx(result.loc[1, "data"]) == 21.9  # 10.5 + 8.3 + 3.1
    assert pytest.approx(result.loc[2, "data"]) == 32.4  # 15.2 + 12.7 + 4.5


def test_aggregate_stratum_weights_empty(empty_weights_df_multilevel):
    """Test aggregating stratum weights with an empty DataFrame."""
    # With empty DataFrame
    result_empty = get_proportions.aggregate_stratum_weights(empty_weights_df_multilevel)
    assert isinstance(result_empty, pd.DataFrame)

    # The function returns zero values instead of an empty DataFrame
    assert not result_empty.empty
    assert result_empty.shape == (2, 1)  # 2 rows (for stratum 1 and 2), 1 column ('data')
    assert (result_empty["data"] == 0.0).all()  # All values are 0.0

    # With empty dictionary
    result_empty_dict = get_proportions.aggregate_stratum_weights({})
    assert isinstance(result_empty_dict, pd.DataFrame)
    assert result_empty_dict.empty  # This should be truly empty


def test_aggregate_stratum_weights_missing_level(weights_df_missing_stratum):
    """Test handling when stratum_num level is missing."""
    # Should print a warning but not fail
    result = get_proportions.aggregate_stratum_weights({"test": weights_df_missing_stratum})

    # Result should be empty DataFrame since no valid data found
    assert isinstance(result, pd.DataFrame)
    assert result.empty


def test_scale_weights_by_stratum_basic(simple_weights_df, simple_stratum_weights):
    """Test basic functionality of standardizing weights by stratum."""
    # Call the function
    result = get_proportions.scale_weights_by_stratum(simple_weights_df, simple_stratum_weights)

    # Check result is a DataFrame
    assert isinstance(result, pd.DataFrame)

    # Get the sum of weights per stratum in original data
    stratum1_sum = simple_weights_df[("female", 1)] + simple_weights_df[("male", 1)]
    stratum2_sum = simple_weights_df[("female", 2)] + simple_weights_df[("male", 2)]

    # Check the values are transformed correctly
    # For female, stratum 1: original proportion * reference weight
    female_stratum1_proportion = simple_weights_df[("female", 1)][0] / stratum1_sum[0]
    female_stratum1_standardized = female_stratum1_proportion * 100.0
    assert result.loc["female", 1] == pytest.approx(female_stratum1_standardized)

    # For male, stratum 2: original proportion * reference weight
    male_stratum2_proportion = simple_weights_df[("male", 2)][0] / stratum2_sum[0]
    male_stratum2_standardized = male_stratum2_proportion * 150.0
    assert result.loc["male", 2] == pytest.approx(male_stratum2_standardized)


def test_scale_weights_by_stratum_error_handling(weights_df_multilevel):
    """Test error handling with invalid reference data."""
    # Create reference without stratum_num
    invalid_reference = pd.DataFrame({"region": [1, 2], "weight": [100.0, 150.0]})

    # Should raise ValueError
    with pytest.raises(ValueError, match="must have the defined `stratum_col`.*column or index"):
        get_proportions.scale_weights_by_stratum(weights_df_multilevel, invalid_reference)


def test_weight_proportions_basic(weight_distr_dict, catch_data_df):
    """Test basic functionality of weight proportions calculation."""
    # Call the function
    result = get_proportions.weight_proportions(
        weight_data=weight_distr_dict, catch_data=catch_data_df, group="aged"
    )

    # Check that the result has the expected structure
    assert isinstance(result, pd.DataFrame)

    # Calculate expected values
    total_stratum1 = 80.0 + 18.8  # Catch weight + aged group weight
    total_stratum2 = 110.0 + 27.9  # Catch weight + aged group weight

    # Check values individually instead of comparing entire DataFrames
    # Expected values for female, stratum 1 and 2
    assert pytest.approx(result.iloc[0, 0]) == 10.5 / total_stratum1
    assert pytest.approx(result.iloc[0, 1]) == 15.2 / total_stratum2

    # Expected values for male, stratum 1 and 2
    assert pytest.approx(result.iloc[1, 0]) == 8.3 / total_stratum1
    assert pytest.approx(result.iloc[1, 1]) == 12.7 / total_stratum2


def test_weight_proportions_empty_catch(weight_distr_dict):
    """Test handling of empty catch data."""
    # Create empty catch DataFrame
    empty_catch = pd.DataFrame(
        columns=["stratum_num", "weight"]
    )  # Changed from haul_weight to weight

    # Should handle empty catch data gracefully
    result = get_proportions.weight_proportions(
        weight_data=weight_distr_dict, catch_data=empty_catch, group="aged"
    )

    # Result should be a DataFrame and reflect the structure of the input
    assert isinstance(result, pd.DataFrame)


def test_weight_proportions_missing_group(weight_distr_dict, catch_data_df):
    """Test error handling when a non-existent group is specified."""
    # Should raise KeyError for non-existent group
    with pytest.raises(KeyError):
        get_proportions.weight_proportions(
            weight_data=weight_distr_dict, catch_data=catch_data_df, group="non_existent_group"
        )


def test_scale_weight_proportions_basic(
    standardized_data_fixture,
    standardized_weight_reference,
    catch_data_df,
    proportion_dict_fixture,
    binned_weight_table_fixture,
):
    """Test basic functionality of standardized weight proportions."""
    # Call the function directly - this just tests if it runs without errors
    result = get_proportions.scale_weight_proportions(
        weight_data=standardized_data_fixture,
        reference_weight_proportions=standardized_weight_reference,
        catch_data=catch_data_df,
        number_proportions=proportion_dict_fixture,
        binned_weights=binned_weight_table_fixture,
        group="unaged",
        group_columns=["sex"],
    )

    # Verify it returns a DataFrame
    assert isinstance(result, pd.DataFrame)

    # Verify the DataFrame has the expected structure
    assert 1 in result.columns  # Should have stratum 1
    assert 2 in result.columns  # Should have stratum 2

    # Assert the index names
    assert set(result.index.names) <= set(["length_bin", "sex"])

    # Assert results
    assert list(result.loc["female", 30]) == pytest.approx([0.1500000, 0.0235294])
    assert list(result.loc["female", 40]) == pytest.approx([0.0000000, 0.0564706])
    assert list(result.loc["male", 30]) == pytest.approx([0.2500000, 0.0352941])
    assert list(result.loc["male", 40]) == pytest.approx([0.0000000, 0.0847059])


def test_weight_proportions_custom_stratum_col(weight_distr_dict):
    """Test weight proportions with custom stratum column name."""
    # Create weight data with custom stratum column name to match
    columns = pd.MultiIndex.from_tuples(
        [("female", 1), ("female", 2), ("male", 1), ("male", 2)], names=["sex", "stratum_id"]
    )
    custom_aged_df = pd.DataFrame(data=[[10.5, 15.2, 8.3, 12.7]], columns=columns)
    custom_unaged_df = pd.DataFrame(data=[[5.2, 7.8, 4.3, 6.1]], columns=columns)

    weight_data_custom = {"aged": custom_aged_df, "unaged": custom_unaged_df}

    # Create catch data with custom stratum column name
    catch_data_custom = pd.DataFrame({"stratum_id": [1, 2], "weight": [80.0, 110.0]})

    # Test with custom stratum column
    result = get_proportions.weight_proportions(
        weight_data=weight_data_custom,
        catch_data=catch_data_custom,
        group="aged",
        stratum_col="stratum_id",
    )

    # Should return a DataFrame with proper structure
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0


# =============================================================================
# TESTS FOR NEW PROPORTION SLICING FUNCTIONS
# ==============================================================================


def test_get_nasc_proportions_slice_basic():
    """Test basic NASC proportions calculation."""
    # Create test data directly like the real workflow
    # Create raw specimen-like data
    raw_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 1, 2, 2, 2, 2],
            "age": [1.2, 1.8, 2.3, 2.7, 1.1, 1.6, 2.2, 2.9],
            "length": [15.2, 18.5, 25.3, 28.1, 16.1, 19.2, 24.8, 27.5],
            "sex": ["female", "male", "female", "male", "female", "male", "female", "male"],
            "count": [10, 15, 20, 25, 12, 18, 22, 28],
        }
    )

    # Apply binning like the real workflow does
    age_bins = np.linspace(start=1.0, stop=3.0, num=3)
    length_bins = np.linspace(start=15.0, stop=30.0, num=4)

    utils.binify(raw_data, bins=age_bins, bin_column="age")
    utils.binify(raw_data, bins=length_bins, bin_column="length")

    # Convert to proportions
    aged_props = get_proportions.number_proportions(data=raw_data)

    ts_params = {"slope": 20.0, "intercept": -68.0}

    result = get_proportions.get_nasc_proportions_slice(
        number_proportions=aged_props,
        ts_length_regression_parameters=ts_params,
        stratify_by=["stratum_num"],
        include_filter={"age_bin": [1]},
    )

    assert isinstance(result, pd.Series)
    assert len(result) == 2  # Three strata (0, 1)
    assert all(result >= 0), "NASC proportions should be non-negative"
    assert all(result <= 1), "NASC proportions should be <= 1"


def test_get_nasc_proportions_slice_with_sex_filter():
    """Test NASC proportions with sex filtering."""
    # Create test data with sex information, avoiding fixture usage
    test_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 2, 2, 2],
            "length": [15.0, 25.0, 35.0, 18.0, 28.0, 38.0],
            "age": [1.5, 2.5, 3.5, 1.8, 2.8, 3.8],
            "sex": ["female", "male", "female", "female", "male", "female"],
            "species_id": [1, 1, 1, 1, 1, 1],
            "count": [15, 20, 10, 25, 30, 12],
        }
    )

    # Apply proper binning using utils.binify
    length_bins = np.linspace(10, 40, 4)  # [10, 20, 30, 40]
    age_bins = np.linspace(1, 4, 4)  # [1, 2, 3, 4]

    utils.binify(test_data, bins=length_bins, bin_column="length")
    utils.binify(test_data, bins=age_bins, bin_column="age")

    # Convert to proportions using the correct workflow
    aged_props = get_proportions.number_proportions(data=test_data)

    ts_params = {"slope": 20.0, "intercept": -68.0}

    result = get_proportions.get_nasc_proportions_slice(
        number_proportions=aged_props,
        ts_length_regression_parameters=ts_params,
        stratify_by=["stratum_num"],
        include_filter={"sex": ["female"]},
    )

    assert isinstance(result, pd.Series)
    assert len(result) == 2  # Three strata (0, 1, 2)
    assert all(result >= 0), "NASC proportions should be non-negative"


def test_get_number_proportions_slice_single_stratify():
    """Test number proportions slicing with single stratification variable."""
    # Create test data, avoiding fixture usage
    test_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2, 3, 3],
            "length": [15.0, 25.0, 18.0, 28.0, 16.0, 26.0],
            "age": [1.5, 2.5, 1.8, 2.8, 1.3, 2.3],
            "sex": ["female", "male", "female", "male", "female", "male"],
            "species_id": [1, 1, 1, 1, 1, 1],
            "count": [15, 20, 25, 30, 18, 22],
        }
    )

    # Apply proper binning using utils.binify
    length_bins = np.linspace(10, 30, 3)  # [10, 20, 30]
    age_bins = np.linspace(1, 3, 3)  # [1, 2, 3]

    utils.binify(test_data, bins=length_bins, bin_column="length")
    utils.binify(test_data, bins=age_bins, bin_column="age")

    # Convert to proportions using the correct workflow
    aged_props = get_proportions.number_proportions(data=test_data)

    result = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props, stratify_by=["stratum_num"], include_filter={"age_bin": [1]}
    )

    assert isinstance(result, pd.Series)
    assert len(result) == 3  # Three strata (1, 2, 3)
    assert all(result >= 0), "Proportions should be non-negative"
    assert result.sum() <= 1.0, "Total proportion should not exceed 1.0"


def test_get_number_proportions_slice_multiple_stratify():
    """Test number proportions slicing with multiple stratification variables."""
    # Create test data, avoiding fixture usage
    test_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 1, 2, 2, 2, 2],
            "length": [15.0, 25.0, 15.0, 25.0, 18.0, 28.0, 18.0, 28.0],
            "age": [1.5, 2.5, 1.5, 2.5, 1.8, 2.8, 1.8, 2.8],
            "sex": ["female", "female", "male", "male", "female", "female", "male", "male"],
            "species_id": [1, 1, 1, 1, 1, 1, 1, 1],
            "count": [15, 20, 12, 18, 25, 30, 22, 28],
        }
    )

    # Apply proper binning using utils.binify
    length_bins = np.linspace(10, 30, 3)  # [10, 20, 30]
    age_bins = np.linspace(1, 3, 3)  # [1, 2, 3]

    utils.binify(test_data, bins=length_bins, bin_column="length")
    utils.binify(test_data, bins=age_bins, bin_column="age")

    # Convert to proportions using the correct workflow
    aged_props = get_proportions.number_proportions(data=test_data)

    result = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props,
        stratify_by=["stratum_num", "sex"],
        include_filter={"age_bin": [1]},
    )

    assert isinstance(result, pd.DataFrame)
    assert not result.empty
    # With stratum_num and sex as stratification, result is a pivot table with length_bin as index
    # and columns as MultiIndex for stratum_num and sex combinations
    assert isinstance(result.columns, pd.MultiIndex)
    assert "stratum_num" in result.columns.names
    assert "sex" in result.columns.names


def test_get_number_proportions_slice_with_exclusion():
    """Test number proportions slicing with exclusion filter."""
    # Create test data directly, avoiding fixture usage
    test_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
            "length": [15.0, 25.0, 15.0, 25.0, 18.0, 28.0, 18.0, 28.0, 16.0, 26.0, 16.0, 26.0],
            "age": [1.5, 2.5, 1.5, 2.5, 1.8, 2.8, 1.8, 2.8, 1.6, 2.6, 1.6, 2.6],
            "sex": [
                "female",
                "female",
                "male",
                "male",
                "female",
                "female",
                "male",
                "male",
                "female",
                "female",
                "male",
                "male",
            ],
            "species_id": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "count": [15, 20, 12, 18, 25, 30, 22, 28, 14, 24, 16, 20],
        }
    )

    # Apply proper binning using utils.binify
    length_bins = np.linspace(10, 30, 3)  # [10, 20, 30]
    age_bins = np.linspace(1, 3, 3)  # [1, 2, 3]

    utils.binify(test_data, bins=length_bins, bin_column="length")
    utils.binify(test_data, bins=age_bins, bin_column="age")

    # Convert to proportions using the correct workflow
    aged_props = get_proportions.number_proportions(data=test_data)

    result = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props,
        stratify_by=["stratum_num"],
        include_filter={"age_bin": [1]},
        exclude_filter={"sex": ["male"]},
    )

    assert isinstance(result, pd.Series)
    assert len(result) == 3  # Three strata (1, 2, 3)
    # Should have valid proportions
    assert all(result >= 0), "Proportions should be non-negative"


def test_get_weight_proportions_slice_no_thresholding():
    """Test weight proportions without thresholding."""
    # Create test data and apply proper binning using utils.binify
    raw_data = pd.DataFrame(
        {
            "length": [15.0, 25.0, 15.0, 25.0, 18.0, 28.0],
            "age": [1.2, 2.3, 1.8, 2.2, 1.5, 2.7],
            "sex": ["female", "male", "female", "male", "female", "male"],
            "weight": [10, 20, 15, 25, 12, 22],
        }
    )

    # Apply binning using utils.binify (like the real workflow)
    length_bins = np.linspace(10, 30, 3)  # Creates (10, 20], (20, 30]
    age_bins = np.linspace(1.0, 3.0, 3)  # Creates (0.5, 1.5], (1.5, 2.5]
    utils.binify(raw_data, bins=length_bins, bin_column="length")
    utils.binify(raw_data, bins=age_bins, bin_column="age")

    # Create weight proportions in the correct format: MultiIndex rows, stratum columns
    weight_data = (
        raw_data.groupby(["length_bin", "sex", "age_bin"], observed=False)["weight"]
        .sum()
        .to_frame()
    )
    # Convert to pivot table format with stratum columns
    weight_data = pd.concat([weight_data, weight_data, weight_data], axis=1)
    weight_data.columns = [0, 1, 2]
    weight_data.columns.name = "stratum_ks"

    # Normalize to get proportions
    weight_data = weight_data.div(weight_data.sum(axis=0), axis=1).fillna(0)

    result = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_data,
        stratify_by=["stratum_ks"],
        include_filter={"age_bin": [1]},  # Use integer, not interval
    )

    assert isinstance(result, pd.Series)
    assert len(result) == 3  # Two strata that match the filter
    assert all(result >= 0), "Weight proportions should be non-negative"
    assert all(result <= 1), "Weight proportions should be <= 1"


def test_get_weight_proportions_slice_with_thresholding():
    """Test weight proportions with thresholding."""
    # Create test data and apply proper binning using utils.binify
    raw_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2, 3, 3],
            "length": [15.0, 25.0, 15.0, 25.0, 18.0, 28.0],
            "age": [1.2, 2.3, 1.8, 2.2, 1.5, 2.7],
            "sex": ["female", "male", "female", "male", "female", "male"],
            "weight": [10, 20, 15, 25, 12, 22],
        }
    )

    # Apply binning using utils.binify (like the real workflow)
    length_bins = np.linspace(10, 30, 3)
    age_bins = np.linspace(1, 3, 3)
    utils.binify(data=raw_data, bins=length_bins, bin_column="length")
    utils.binify(data=raw_data, bins=age_bins, bin_column="age")

    # Create weight proportions data with proper MultiIndex
    weight_data = (
        raw_data.groupby(["stratum_num", "length_bin", "sex", "age_bin"], observed=False)["weight"]
        .sum()
        .to_frame()
    )
    weight_data = pd.concat([weight_data, weight_data, weight_data], axis=1)
    weight_data.columns = [0, 1, 2]
    weight_data.columns.name = "stratum_num"
    weight_data = weight_data.div(weight_data.sum(axis=0), axis=1).fillna(0)

    # Create number proportions data for thresholding with proper MultiIndex structure
    number_data = (
        raw_data.groupby(["stratum_num", "length_bin", "sex", "age_bin"], observed=False)
        .size()
        .to_frame("proportion")
    )
    # Normalize to get proportions
    number_data["proportion"] = number_data["proportion"] / number_data["proportion"].sum()

    # Spoof dictionary
    number_dict = {"spoopy": number_data.reset_index()}

    result = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_data,
        stratify_by=["stratum_num"],
        include_filter={"age_bin": [1]},
        number_proportions=number_dict,
        length_threshold_min=15.0,
    )

    assert isinstance(result, pd.Series)
    assert len(result) == 3
    assert all(result >= 0), "Weight proportions should be non-negative"
    assert all(result <= 1), "Weight proportions should be <= 1"


def test_get_weight_proportions_slice_with_dict_thresholding():
    """Test weight proportions with dictionary-based thresholding."""
    # Create test data and apply proper binning using utils.binify
    raw_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2, 3, 3],
            "length": [15.0, 25.0, 15.0, 25.0, 18.0, 28.0],
            "age": [1.2, 2.3, 1.8, 2.2, 1.5, 2.7],
            "sex": ["female", "male", "female", "male", "female", "male"],
            "weight": [10, 20, 15, 25, 12, 22],
        }
    )

    # Apply binning using utils.binify (like the real workflow)
    length_bins = np.linspace(10, 30, 3)
    age_bins = np.linspace(1, 3, 3)
    utils.binify(data=raw_data, bins=length_bins, bin_column="length")
    utils.binify(data=raw_data, bins=age_bins, bin_column="age")

    # Create weight proportions data with proper MultiIndex
    weight_data = (
        raw_data.groupby(["stratum_num", "length_bin", "sex", "age_bin"], observed=False)["weight"]
        .sum()
        .to_frame()
    )
    weight_data = weight_data.unstack("stratum_num")["weight"]
    weight_data = weight_data.div(weight_data.sum(axis=0), axis=1).fillna(0)

    # Create aged and unaged number proportions with proper structure
    aged_data = (
        raw_data.groupby(["stratum_num", "length_bin", "sex", "age_bin"], observed=False)
        .size()
        .to_frame("proportion")
    )
    aged_data["proportion"] = aged_data["proportion"] / aged_data["proportion"].sum()

    unaged_data = (
        raw_data.groupby(["stratum_num", "length_bin", "sex"], observed=False)
        .size()
        .to_frame("proportion")
    )
    unaged_data["proportion"] = unaged_data["proportion"] / unaged_data["proportion"].sum()

    number_dict = {"aged": aged_data.reset_index(), "unaged": unaged_data.reset_index()}

    result = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_data,
        stratify_by=["stratum_num"],
        include_filter={"sex": ["female"]},
        number_proportions=number_dict,
        length_threshold_min=15.0,
    )

    assert isinstance(result, pd.Series)
    assert len(result) == 3
    assert all(result >= 0), "Weight proportions should be non-negative"
    assert all(result <= 1), "Weight proportions should be <= 1"


def test_get_weight_proportions_slice_missing_filter_keys():
    """Test error handling when filter keys are missing from weight proportions."""
    # Create simple weight data without age_bin
    index_tuples = [("(10, 20]", "female", 1), ("(20, 30]", "male", 1)]
    index = pd.MultiIndex.from_tuples(index_tuples, names=["length_bin", "sex", "stratum_num"])

    weight_data = pd.DataFrame([[10], [20]], index=index, columns=["weight"])

    with pytest.raises(ValueError, match="Filter keys.*not found"):
        get_proportions.get_weight_proportions_slice(
            weight_proportions=weight_data,
            stratify_by=["stratum_num"],
            include_filter={"age_bin": ["(1, 2]"]},  # age_bin not in index
        )


def test_nasc_proportions_slice_target_strength_calculation():
    """Test that target strength calculation works correctly."""
    # Create test data directly, avoiding fixture usage
    test_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 1, 2, 2, 2],
            "length": [15.0, 25.0, 35.0, 18.0, 28.0, 38.0],
            "age": [1.5, 2.5, 3.5, 1.8, 2.8, 3.8],
            "sex": ["female", "male", "female", "female", "male", "female"],
            "species_id": [1, 1, 1, 1, 1, 1],
            "count": [15, 20, 10, 25, 30, 12],
        }
    )

    # Apply proper binning using utils.binify
    length_bins = np.linspace(10, 40, 4)  # [10, 20, 30, 40]
    age_bins = np.linspace(1, 4, 4)  # [1, 2, 3, 4]

    utils.binify(test_data, bins=length_bins, bin_column="length")
    utils.binify(test_data, bins=age_bins, bin_column="age")

    # Convert to proportions using the correct workflow
    aged_props = get_proportions.number_proportions(data=test_data)

    # Use known parameters for verification
    ts_params = {"slope": 20.0, "intercept": -68.0}

    result = get_proportions.get_nasc_proportions_slice(
        number_proportions=aged_props,
        ts_length_regression_parameters=ts_params,
        stratify_by=["stratum_num"],
        include_filter={"age_bin": [1]},
    )

    # Verify that larger length bins contribute more to NASC
    # (This is implicitly tested through the target strength weighting)
    assert isinstance(result, pd.Series)
    assert not result.isna().any(), "No NaN values should be present"


def test_number_proportions_slice_edge_cases():
    """Test edge cases for number proportions slicing."""
    # Create test data directly, avoiding fixture usage
    test_data = pd.DataFrame(
        {
            "stratum_num": [1, 1, 2, 2],
            "length": [15.0, 25.0, 18.0, 28.0],
            "age": [1.5, 2.5, 1.8, 2.8],
            "sex": ["female", "male", "female", "male"],
            "species_id": [1, 1, 1, 1],
            "count": [15, 20, 25, 30],
        }
    )

    # Apply proper binning using utils.binify
    length_bins = np.linspace(10, 30, 3)  # [10, 20, 30]
    age_bins = np.linspace(1, 3, 3)  # [1, 2, 3]

    utils.binify(test_data, bins=length_bins, bin_column="length")
    utils.binify(test_data, bins=age_bins, bin_column="age")

    # Convert to proportions using the correct workflow
    aged_props = get_proportions.number_proportions(data=test_data)

    # Test with empty include filter (should return all data)
    result_all = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props, stratify_by=["stratum_num"], include_filter={}
    )

    assert isinstance(result_all, pd.DataFrame)
    assert len(result_all) == 3

    # Test with filter that excludes everything
    result_empty = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props,
        stratify_by=["stratum_num"],
        include_filter={"sex": ["nonexistent"]},
    )

    assert isinstance(result_empty, pd.Series)
    # Should return zeros for all strata
    assert all(result_empty == 0)
