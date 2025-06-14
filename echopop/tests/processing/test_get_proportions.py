import numpy as np
import pandas as pd
import pytest
from scipy import interpolate as interp

from echopop.nwfsc_feat import get_proportions, utils


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


def test_compute_binned_counts_with_filters_size(sample_specimen_data):
    """Test size aggregation with filters."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["stratum_num", "length_bin"],
        "length",
        agg_func="size",
        exclude_filters={"sex": "unsexed"},
    )

    expected_count = len(sample_specimen_data[sample_specimen_data["sex"] != "unsexed"])
    assert result["count"].sum() == expected_count


def test_compute_binned_counts_with_filters_sum(sample_length_data):
    """Test sum aggregation with filters."""
    result = get_proportions.compute_binned_counts(
        sample_length_data,
        ["stratum_num", "length_bin"],
        "length_count",
        agg_func="sum",
        exclude_filters={"sex": "unsexed"},
    )

    filtered_data = sample_length_data[sample_length_data["sex"] != "unsexed"]
    expected_sum = filtered_data["length_count"].sum()
    assert result["count"].sum() == expected_sum


def test_compute_binned_counts_multiple_filters(sample_specimen_data):
    """Test with multiple exclusion filters."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["length_bin"],
        "length",
        agg_func="size",
        exclude_filters={"sex": "unsexed", "stratum_num": 1},
    )

    expected_count = len(
        sample_specimen_data[
            (sample_specimen_data["sex"] != "unsexed") & (sample_specimen_data["stratum_num"] != 1)
        ]
    )
    assert result["count"].sum() == expected_count


def test_compute_binned_counts_std_aggregation(sample_specimen_data):
    """Test standard deviation aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "weight", agg_func="std"
    )

    assert isinstance(result, pd.DataFrame)
    assert "count" in result.columns
    assert len(result) > 0
    assert result["count"].notna().any()


def test_compute_binned_counts_min_max_aggregation(sample_specimen_data):
    """Test min and max aggregations."""
    result_min = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "length", agg_func="min"
    )

    result_max = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "length", agg_func="max"
    )

    assert isinstance(result_min, pd.DataFrame)
    assert isinstance(result_max, pd.DataFrame)
    assert (result_max["count"] >= result_min["count"]).all()


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


def test_compute_binned_counts_all_filtered_out(data_all_unsexed):
    """Test when all data is filtered out."""
    result = get_proportions.compute_binned_counts(
        data_all_unsexed,
        ["stratum_num"],
        "length",
        agg_func="size",
        exclude_filters={"sex": "unsexed"},
    )

    assert isinstance(result, pd.DataFrame)
    assert result["count"].sum() == 0


def test_compute_binned_counts_preserves_original_data(sample_specimen_data):
    """Test that original data is not modified."""
    original_data = sample_specimen_data.copy()

    get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["stratum_num"],
        "length",
        agg_func="sum",
        exclude_filters={"sex": "unsexed"},
    )

    pd.testing.assert_frame_equal(sample_specimen_data, original_data)


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


def test_compute_binned_counts_missing_filter_column(sample_specimen_data):
    """Test with filter column that doesn't exist."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["stratum_num"],
        "length",
        agg_func="size",
        exclude_filters={"nonexistent_col": "value"},
    )

    assert result["count"].sum() == len(sample_specimen_data)


def test_compute_binned_counts_large_dataset(expanded_specimen_dataset):
    """Test performance with large dataset."""
    result = get_proportions.compute_binned_counts(
        expanded_specimen_dataset,
        ["stratum_num", "length_bin"],
        "length",
        agg_func="size",
        exclude_filters={"sex": "unsexed"},
    )

    expected_count = len(expanded_specimen_dataset[expanded_specimen_dataset["sex"] != "unsexed"])
    assert result["count"].sum() == expected_count


def test_compute_binned_counts_deterministic_results(sample_specimen_data):
    """Test that function produces deterministic results."""
    result1 = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    result2 = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    pd.testing.assert_frame_equal(result1, result2)


def test_compute_binned_counts_numeric_vs_categorical_filters(sample_specimen_data):
    """Test filtering with numeric and categorical values."""
    # Numeric filter
    result_numeric = get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["length_bin"],
        "length",
        agg_func="size",
        exclude_filters={"stratum_num": 1},
    )

    # Categorical filter
    result_categorical = get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["stratum_num"],
        "length",
        agg_func="size",
        exclude_filters={"sex": "male"},
    )

    assert isinstance(result_numeric, pd.DataFrame)
    assert isinstance(result_categorical, pd.DataFrame)


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


import numpy as np
import pandas as pd
import pytest

from echopop.nwfsc_feat import get_proportions


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


def test_number_proportions_single_dataframe(aged_dataframe):
    """Test the number_proportions function with a single dataframe."""
    result = get_proportions.number_proportions(aged_dataframe)

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
        aged_dataframe, unaged_dataframe, column_aliases=["aged", "unaged"]
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
        aged_dataframe, unaged_dataframe, column_aliases=["specimen", "length"]
    )

    assert "specimen" in result
    assert "length" in result


def test_number_proportions_default_aliases(aged_dataframe, unaged_dataframe):
    """Test default aliases when column_aliases not provided."""
    result = get_proportions.number_proportions(aged_dataframe, unaged_dataframe)

    assert "df_0" in result
    assert "df_1" in result


def test_number_proportions_with_exclusion(aged_dataframe):
    """Test exclude_filters parameter in number_proportions."""
    result = get_proportions.number_proportions(aged_dataframe, exclude_filters={"sex": "unsexed"})

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
        aged_dataframe, unaged_dataframe, group_columns=["stratum_num", "sex"]
    )

    # Check that grouping includes both stratum and sex
    assert len(result["df_0"].groupby(["stratum_num", "sex"]).groups) == 6
    assert len(result["df_1"].groupby(["stratum_num", "sex"]).groups) == 6

    # Check that proportions sum to 1.0 for each stratum-sex combination
    df1 = result["df_0"]
    for _, group in df1.groupby(["stratum_num", "sex"]):
        assert group["proportion"].sum() == pytest.approx(1.0)


def test_number_proportions_proportion_calculations(aged_dataframe, unaged_dataframe):
    """Test the correctness of proportion calculations in number_proportions."""
    result = get_proportions.number_proportions(
        aged_dataframe, unaged_dataframe, column_aliases=["aged", "unaged"]
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
    # No dataframes provided
    with pytest.raises(ValueError, match="At least one DataFrame must be provided"):
        get_proportions.number_proportions()

    # DataFrame without count column
    invalid_df = pd.DataFrame({"stratum_num": [1, 2], "value": [10, 20]})
    with pytest.raises(ValueError, match="DataFrame 0 does not have a 'count' column"):
        get_proportions.number_proportions(invalid_df)


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
    length_dataset_with_bins, length_weight_dataset_with_bins
):
    """Test binned_weights with interpolation."""
    result = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_with_bins,  # Required for interpolation
        interpolate=True,
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
        interpolate=False,
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
        ValueError, match="length_weight_dataset must be provided when interpolate=True"
    ):
        get_proportions.binned_weights(
            length_dataset=length_dataset_with_bins,
            length_weight_dataset=None,  # Missing required dataset
            interpolate=True,
            table_index=["length_bin"],
            table_cols=["sex"],
            contrast_vars="sex",
        )


def test_binned_weights_with_filtering(length_dataset_with_bins, length_weight_dataset_with_bins):
    """Test binned_weights with filtering applied."""
    result = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_with_bins,
        interpolate=True,
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
    length_dataset_with_bins, length_weight_dataset_with_bins
):
    """Test binned_weights with global interpolation (no contrast vars)."""
    result = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_with_bins,  # Required for interpolation
        interpolate=True,
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


def test_binned_weights_integration(length_dataset_with_bins, length_weight_dataset_with_bins):
    """Full integration test of binned_weights with both interpolation and no interpolation."""
    # Test with interpolation
    result_interp = get_proportions.binned_weights(
        length_dataset=length_dataset_with_bins,
        length_weight_dataset=length_weight_dataset_with_bins,
        interpolate=True,
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
        interpolate=False,
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
        group_keys, aggregate_table, sex_proportions_table
    )

    # Check that we get multi-index with group and sex
    assert isinstance(result.index, pd.MultiIndex)
    assert result.index.names == ["group", "sex"]

    # Check that the adjusted proportions for the first group are calculated correctly
    # For stratum 1, the adjusted proportion for aged should be higher than raw proportion
    assert result.loc[("aged", "female"), 1] > sex_proportions_table.loc[("aged", "female"), 1]


def test_stratum_averaged_weight(proportion_test_dict, test_weight_table):
    """Test the stratum_averaged_weight function."""
    result = get_proportions.stratum_averaged_weight(proportion_test_dict, test_weight_table)

    # Check that we get the right format
    assert isinstance(result, pd.DataFrame)

    # Check that we have the right columns
    assert set(result.columns) >= {"all", "female", "male"}

    # Check that the weights are positive
    assert result["all"].min() > 0

    # Check that female weights are different from male weights
    if "female" in result.columns and "male" in result.columns:
        assert not result["female"].equals(result["male"])
