import numpy as np
import pandas as pd
import xarray as xr
import pytest
from scipy import interpolate as interp

from echopop import utils
from echopop.survey import proportions as get_proportions

def test_compute_binned_counts_size_aggregation(sample_specimen_data):
    """Test size aggregation (default)."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "length", agg_func="size"
    )

    assert isinstance(result, xr.DataArray)
    assert set(result.coords) == set(["variable", "stratum_num", "length_bin"])
    assert result.sum() == len(sample_specimen_data)


def test_compute_binned_counts_sum_aggregation(sample_length_data):
    """Test sum aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_length_data, ["stratum_num", "length_bin"], "length_count", agg_func="sum"
    )

    assert isinstance(result, xr.DataArray)
    assert set(result.coords) == set(["variable", "stratum_num", "length_bin"])
    assert result.sum() == sample_length_data["length_count"].sum()


def test_compute_binned_counts_mean_aggregation(sample_specimen_data):
    """Test mean aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    assert isinstance(result, xr.DataArray)
    assert set(result.coords) == set(["variable", "stratum_num", "length_bin"])
    assert len(result) > 0


def test_compute_binned_counts_count_aggregation(sample_specimen_data):
    """Test count aggregation."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="count"
    )

    assert isinstance(result, xr.DataArray)
    assert set(result.coords) == set(["variable", "stratum_num", "length_bin"])
    assert result.sum() == len(sample_specimen_data)


def test_compute_binned_counts_default_size(sample_specimen_data):
    """Test default aggregation is size."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "length"
    )
    
    assert isinstance(result, xr.DataArray)
    assert set(result.coords) == set(["variable", "stratum_num", "length_bin"])
    assert result.sum() == len(sample_specimen_data)
    

def test_compute_binned_counts_single_group_column(sample_specimen_data):
    """Test with single grouping column."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "length", agg_func="size"
    )
    
    assert set(result.coords) == set(["variable", "stratum_num"])
    assert len(result.sel(variable="count")) == len(sample_specimen_data["stratum_num"].unique())

def test_compute_binned_counts_complex_grouping(sample_specimen_data):
    """Test with complex grouping."""
    result = get_proportions.compute_binned_counts(
        sample_specimen_data,
        ["stratum_num", "length_bin", "age_bin", "sex"],
        "weight",
        agg_func="mean",
    )

    assert isinstance(result, xr.DataArray)
    assert set(result.coords) == set(["variable", "stratum_num", "length_bin", "age_bin", "sex"])
    assert result.size == sample_specimen_data.nunique().filter(set(result.coords)).product()


def test_compute_binned_counts_empty_data(empty_specimen_data):
    """Test with empty DataFrame."""
    result = get_proportions.compute_binned_counts(
        empty_specimen_data, ["stratum_num", "length_bin"], "length", agg_func="size"
    )

    assert isinstance(result, xr.DataArray)
    assert result.shape == (1, 0, 0)

def test_compute_binned_counts_different_numeric_columns(sample_specimen_data):
    """Test aggregation on different numeric columns."""
    result_length = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "length", agg_func="mean"
    )

    result_weight = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num"], "weight", agg_func="mean"
    )

    assert isinstance(result_length, xr.DataArray)
    assert isinstance(result_weight, xr.DataArray)
    assert result_length.size == result_weight.size


def test_compute_binned_counts_deterministic_results(sample_specimen_data):
    """Test that function produces deterministic results."""
    result1 = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    result2 = get_proportions.compute_binned_counts(
        sample_specimen_data, ["stratum_num", "length_bin"], "weight", agg_func="mean"
    )

    assert result1.equals(result2)
    

def test_compute_binned_counts_column_order_preserved(sample_specimen_data):
    """Test that groupby column order is preserved in output."""
    groupby_cols = ["stratum_num", "length_bin", "sex"]
    result = get_proportions.compute_binned_counts(
        sample_specimen_data, groupby_cols, "length", agg_func="size"
    )

    # Check that columns appear in expected order
    expected_cols = groupby_cols + ["variable"]
    assert set(list(result.coords)) <= set(expected_cols)


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

    assert result1.size <= result2.size
    assert result1.sum() == result2.sum()  # Same total count


def test_number_proportions_single_dataarray(aged_dataarray):
    """Test the number_proportions function with a single DataArray."""
    result = get_proportions.number_proportions(data=aged_dataarray, group_columns=["stratum_num"])
    
    # Check that result is a Dataset and has the right coords
    assert isinstance(result, xr.Dataset)
    assert set(result.data_vars) == {"proportion_overall", "proportion", "count"}
    assert set(result.coords) == (set(aged_dataarray.coords) - {"variable"})
    
    # Check that proportions sum to 1.0 for each stratum
    assert all(result["proportion"].sum(dim=["length_bin", "age_bin", "sex"]) == 1)
    assert all(result["proportion_overall"].sum(dim=["length_bin", "age_bin", "sex"]) == 1)
    

def test_number_proportions_multiple_dataarrays(aged_dataarray, unaged_dataarray):
    """Test the number_proportions function with multiple DataArrays."""
    result = get_proportions.number_proportions(
        data=xr.Dataset({"aged": aged_dataarray, "unaged": unaged_dataarray}),
        group_columns=["stratum_num"]
    )

    # Check that result is a dictionary with the right keys
    assert isinstance(result, dict)
    assert "aged" in result
    assert "unaged" in result

    # Check the aged DataArray
    aged = result["aged"]
    assert isinstance(aged, xr.Dataset)
    assert set(aged.coords) == ({*aged_dataarray.coords} - {"variable"})    
    assert set(aged.data_vars)  == {"count", "proportion", "proportion_overall"}   

    # Check the unaged dataframe
    unaged = result["unaged"]
    assert isinstance(unaged, xr.Dataset)
    assert set(unaged.coords) == ({*unaged_dataarray.coords} - {"variable"})    
    assert set(unaged.data_vars)  == {"count", "proportion", "proportion_overall"}  

    # Check sums of within-group proportions
    assert all(aged["proportion"].sum(dim=["length_bin", "age_bin", "sex"]) == 1)
    assert all(unaged["proportion"].sum(dim=["length_bin", "sex"]) == 1)
    
    # Check sums across group proportions
    assert all(
        aged["proportion_overall"].sum(dim=["length_bin", "age_bin", "sex"]) + 
        unaged["proportion_overall"].sum(dim=["length_bin", "sex"])
        == 1
    )

def test_number_proportions_column_aliases(aged_dataarray, unaged_dataarray):
    """Test column_aliases parameter in number_proportions."""
    result = get_proportions.number_proportions(
        data=xr.Dataset({"specimen": aged_dataarray, "length": unaged_dataarray}),
        group_columns=["stratum_num"]
    )

    assert "specimen" in result
    assert "length" in result


def test_number_proportions_with_exclusion(aged_dataarray):
    """Test exclude_filters parameter in number_proportions."""
    result = get_proportions.number_proportions(
        data=aged_dataarray, exclude_filters={"sex": "unsexed"}, group_columns=["stratum_num"]
    )

    # Check that unsexed rows are excluded
    assert {*result["sex"].values} == {"female", "male"}

    # Check that proportions are recalculated correctly
    stratum1_total = aged_dataarray.sel(stratum_num=1, sex=["female", "male"]).sum()
    stratum1_female = aged_dataarray.sel(stratum_num=1, sex=["female"]).sum()
    expected_proportion = stratum1_female / stratum1_total
    actual_proportion = result.sel(stratum_num=1, sex="female")["proportion"]
    assert actual_proportion.sum() == expected_proportion


def test_number_proportions_custom_group_columns(aged_dataarray, unaged_dataarray):
    """Test custom group_columns parameter in number_proportions."""
    result = get_proportions.number_proportions(
        data=xr.Dataset({"aged": aged_dataarray, "unaged": unaged_dataarray}),
        group_columns=["stratum_num", "sex"],
    )  
    
    # Check that grouping includes both stratum and sex
    assert len(result["aged"].groupby(["stratum_num", "sex"]).groups) == 6
    assert len(result["unaged"].groupby(["stratum_num", "sex"]).groups) == 6

    # Check that proportions sum to 1.0 for each stratum-sex combination
    assert (result["aged"]["proportion"].sum(dim=["length_bin", "age_bin"]).values == 1).all()
    assert (result["unaged"]["proportion"].sum(dim=["length_bin"]).values == 1).all()

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
        length_data=length_dataset_with_bins,
        length_weight_data=length_weight_dataset_wide_format,  # Required for interpolation
        interpolate_regression=True,
        group_columns=["stratum", "sex"]
    )

    # Check the result is a DataFrame
    assert isinstance(result, xr.DataArray)
    # Check it has the expected structure
    assert {*result.coords} == {"length_bin", "stratum", "sex"}
    assert result.shape == tuple(length_dataset_with_bins.nunique().loc[[*result.coords]])


def test_binned_weights_without_interpolation(length_dataset_with_bins):
    """Test binned_weights without interpolation."""
    
    # Add direct weights to the dataset
    df = length_dataset_with_bins.copy()
    df["weight"] = df["length"] * 2.0  # Simple weight calculation

    # Note: length_weight_dataset is not provided since interpolate=False
    result = get_proportions.binned_weights(
        length_data=df,
        interpolate_regression=False,
        group_columns=["sex"]
    )

    # Check the result
    assert isinstance(result, xr.DataArray)
    assert {*result.coords} == {"length_bin", "sex"}


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
        number_proportions=proportion_test_dict,
        length_weight_data=test_weight_table,
        group_columns=["stratum_num"],
    )

    # Check that we get the right format
    assert isinstance(result, xr.DataArray)

    # Check that we have the right coordinates
    assert {*result.coords} == {"sex", "stratum_num"}
    assert (result["stratum_num"] == np.array([1, 2])).all()
    assert (result["sex"] == np.array(["all", "female", "male"])).all()
    
    # Check that the weights are positive
    assert result.min() > 0

    # Check that female weights are different from male weights
    assert (result.sel(sex="female") != result.sel(sex="male")).all()

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
    result = get_proportions.scale_weights_by_stratum(
        weight_data=simple_weights_df,
        catch_data=simple_stratum_weights, 
        group_columns=["stratum_num"]
    )

    # Check result is a DataFrame
    assert isinstance(result,xr.DataArray)

    # Get the sum of weights per stratum in original data
    stratum_sums = simple_weights_df.sum(dim=["sex"])

    # Check the values are transformed correctly
    # For female, stratum 1: original proportion * reference weight
    female_stratum1_proportion = (
        simple_weights_df.sel(stratum_num=1, sex="female") / stratum_sums.sel(stratum_num=1)
    )
    female_stratum1_standardized = female_stratum1_proportion * 100.0
    assert result.sel(stratum_num=1, sex="female") == female_stratum1_standardized 

    # For male, stratum 2: original proportion * reference weight
    male_stratum2_proportion = (
        simple_weights_df.sel(stratum_num=2, sex="male") / stratum_sums.sel(stratum_num=2)
    )
    male_stratum2_standardized = male_stratum2_proportion * 150.0
    assert result.sel(stratum_num=2, sex="male") == male_stratum2_standardized 


def test_weight_proportions_basic(weight_distr_dict, catch_data_df):
    """Test basic functionality of weight proportions calculation."""
    
    # Call the function
    result = get_proportions.weight_proportions(
        weight_data=weight_distr_dict["aged"], 
        catch_data=catch_data_df, 
        group_columns=["stratum_num"]
    )["proportion_overall"]

    # Check that the result has the expected structure
    assert isinstance(result, xr.DataArray)

    # Calculate expected values
    total_stratum1 = 80.0 + 18.8  # Catch weight + aged group weight
    total_stratum2 = 110.0 + 27.9  # Catch weight + aged group weight

    # Check values individually instead of comparing entire DataFrames
    # Expected values for female, stratum 1 and 2
    assert (result.sel(stratum_num=1).values == np.array([10.5, 8.3]) / total_stratum1).all()
    assert (result.sel(stratum_num=2).values == np.array([15.2, 12.7]) / total_stratum2).all()
    

def test_scale_weight_proportions_basic(
    simple_weights_df, 
    simple_stratum_weights,
    weight_distr_dict, 
    catch_data_df,
    aged_dataarray,
    unaged_dataarray
):
    """Test basic functionality of standardized weight proportions."""

    # Scale weights by stratum first
    strata_weights = get_proportions.scale_weights_by_stratum(
        weight_data=simple_weights_df,
        catch_data=simple_stratum_weights, 
        group_columns=["stratum_num"]
    )
    
    # Tabulate weight proportions
    weight_props = get_proportions.weight_proportions(
        weight_data=weight_distr_dict["aged"],
        catch_data=catch_data_df,
        group_columns=["stratum_num"]
    )
    
    # Get number proportions
    number_props = get_proportions.number_proportions(
        data=xr.Dataset({"aged": aged_dataarray, "unaged": unaged_dataarray}),
        group_columns=["stratum_num", "sex"],
    )  
    
    # Get the length-binned weights
    lb_weights = pd.DataFrame({
        "length_bin": number_props["unaged"]["length_bin"].values,
        "sex": np.repeat("all", 3),
        "weight_fitted": [1, 2, 3]
    }).set_index(
        ["length_bin", "sex"]
    ).to_xarray().to_dataarray().squeeze("variable").reset_coords("variable", drop=True)
    
    # Call the function directly - this just tests if it runs without errors
    result = get_proportions.scale_weight_proportions(
        scaled_weight_data=strata_weights,
        reference_weight_proportions=weight_props,
        catch_data=catch_data_df,
        number_proportions=number_props["unaged"],
        binned_weights=lb_weights,
        group_columns=["stratum_num"]
    )

    # Verify it returns a DataFrame
    assert isinstance(result, xr.Dataset)
    
    # Convert to Array for ease of comparisons
    result_da = result["proportion_overall"]

    # Verify the DataFrame has the expected structure
    assert result_da.shape == (2, 3, 2)

    # Assert results
    assert (
        result_da.sel(sex="female", length_bin="(20, 30]") == 
        pytest.approx([0.15074511, 0.14485936])
    ).all()
    assert (
        result_da.sel(sex="female", length_bin="(30, 40]") == 
        pytest.approx([0.22611767, 0.21728903])
    ).all()
    assert(
        result_da.sel(sex="male", length_bin="(20, 30]") == 
        pytest.approx([0.11916042, 0.1210338 ])        
    ).all()
    assert(
        result_da.sel(sex="male", length_bin="(30, 40]") == 
        pytest.approx([0.17874063, 0.18155071])        
    ).all()


# =============================================================================
# TESTS FOR NEW PROPORTION SLICING FUNCTIONS
# ==============================================================================


def test_get_nasc_proportions_slice():
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
    data_cnv = (
        raw_data.set_index(["length_bin", "age_bin", "stratum_num", "sex"])["count"].to_xarray()
    )
    data_cnv = data_cnv.assign_coords(variable="count")
    aged_props = get_proportions.number_proportions(data=data_cnv, group_columns=["stratum_num"])

    ts_params = {"slope": 20.0, "intercept": -68.0}

    result = get_proportions.get_nasc_proportions_slice(
        number_proportions=aged_props,
        group_columns=["stratum_num"],
        ts_length_regression_parameters=ts_params,
        include_filter={"age_bin": [1]},
    )

    assert isinstance(result, xr.DataArray)
    assert result.shape == (2,)
    assert all(result >= 0), "NASC proportions should be non-negative"
    assert all(result <= 1), "NASC proportions should be <= 1"
    
    # Compare to non-filtered
    result_full = get_proportions.get_nasc_proportions_slice(
        number_proportions=aged_props,
        group_columns=["stratum_num"],
        ts_length_regression_parameters=ts_params,
    )
    
    assert isinstance(result_full, xr.DataArray)
    assert result_full.shape == (2,)
    assert all(result_full ==  1)
    
    # Use a different contrast for filtering
    result_female = get_proportions.get_nasc_proportions_slice(
        number_proportions=aged_props,
        ts_length_regression_parameters=ts_params,
        group_columns=["stratum_num"],
        include_filter={"sex": ["female"]},
    )

    assert isinstance(result_female, xr.DataArray)
    assert result_female.shape == (2,)
    assert all(result_female >= 0), "NASC proportions should be non-negative"
    assert all(result_female <= 1), "NASC proportions should be <= 1"
    
    # Test inverse with exclusion filter
    result_male_excl = get_proportions.get_nasc_proportions_slice(
        number_proportions=aged_props,
        ts_length_regression_parameters=ts_params,
        group_columns=["stratum_num"],
        exclude_filter={"sex": ["male"]},
    )
    
    assert isinstance(result_male_excl, xr.DataArray)
    assert result_male_excl.shape == (2,)
    assert all(result_male_excl >= 0), "NASC proportions should be non-negative"
    assert all(result_male_excl <= 1), "NASC proportions should be <= 1"
    assert all(result_male_excl == result_female)
    
    # Try different grouping
    sexed_props = get_proportions.number_proportions(data=data_cnv, group_columns=["sex"])
    result_sex = get_proportions.get_nasc_proportions_slice(
        number_proportions=sexed_props,
        ts_length_regression_parameters=ts_params,
        group_columns=["sex"],
    )
    
    assert isinstance(result_sex, xr.DataArray)
    assert result_sex.shape == (2,)
    assert all(result_sex ==  1)


def test_get_number_proportions_slice():
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

    # Convert to proportions
    data_cnv = (
        test_data.set_index(["length_bin", "age_bin", "stratum_num", "sex"])["count"].to_xarray()
    )
    data_cnv = data_cnv.assign_coords(variable="count")
    aged_props = get_proportions.number_proportions(data=data_cnv, group_columns=["stratum_num"])

    # Convert to proportions using the correct workflow
    result = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props, 
        group_columns=["stratum_num"], 
        include_filter={"age_bin": [1]}
    )

    assert isinstance(result, xr.DataArray)
    assert result.size == 3  # Three strata (1, 2, 3)
    assert result.shape == (3,)
    assert all(result >= 0), "Proportions should be non-negative"
    assert result.sum() <= 1.0, "Total proportion should not exceed 1.0"
    
    # Compare to non-filtered
    result_full = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props, 
        group_columns=["stratum_num"]
    )
    
    assert isinstance(result_full, xr.DataArray)
    assert result_full.shape == (3,)
    assert all(result_full ==  1)
    
    # Use a different contrast for filtering
    result_female = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props, 
        group_columns=["stratum_num"], 
        include_filter={"sex": ["female"]}
    )
    
    assert isinstance(result_female, xr.DataArray)
    assert result_female.size == 3
    assert all(result_female >= 0), "NASC proportions should be non-negative"
    assert all(result_female <= 1), "NASC proportions should be <= 1"
    
    # Test inverse with exclusion filter
    result_male_excl = get_proportions.get_number_proportions_slice(
        number_proportions=aged_props, 
        group_columns=["stratum_num"], 
        exclude_filter={"sex": ["male"]}
    )
    
    assert isinstance(result_male_excl, xr.DataArray)
    assert result_male_excl.shape == (3,)
    assert all(result_male_excl >= 0), "NASC proportions should be non-negative"
    assert all(result_male_excl <= 1), "NASC proportions should be <= 1"
    assert all(result_male_excl == result_female)
    
    # Try different grouping
    sexed_props = get_proportions.number_proportions(data=data_cnv, group_columns=["sex"])
    result_sex = get_proportions.get_number_proportions_slice(
        number_proportions=sexed_props, 
        group_columns=["sex"],
    )

    assert isinstance(result_sex, xr.DataArray)
    assert result_sex.shape == (2,)
    assert all(result_sex ==  1)

def test_get_weight_proportions_slice():
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
    weight_data.columns.name = "stratum_num"

    # Normalize to get proportions
    weight_data = weight_data.div(weight_data.sum(axis=0), axis=1).fillna(0)
    # ---- Convert to DataArray
    weight_cnv = weight_data.stack().to_xarray().rename("proportion_overall")
    
    # Apply binning using utils.binify (like the real workflow)
    utils.binify(test_data, bins=length_bins, bin_column="length")
    utils.binify(test_data, bins=age_bins, bin_column="age")
    
    # Convert to proportions
    data_cnv = (
        test_data.set_index(["length_bin", "age_bin", "stratum_num", "sex"])["count"].to_xarray()
    )
    data_cnv = data_cnv.assign_coords(variable="count")
    aged_props = {
        "aged": get_proportions.number_proportions(data=data_cnv, group_columns=["stratum_num"])  
    }  

    # Convert to proportions using the correct workflow
    result = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["stratum_num"],
        include_filter={"age_bin": [1]},
    )

    assert isinstance(result, xr.DataArray)
    assert result.size == 3  # Three strata (1, 2, 3)
    assert result.shape == (3,)
    assert all(result >= 0), "Proportions should be non-negative"
    assert result.sum() <= 1.0, "Total proportion should not exceed 1.0"
    
    # Compare to non-filtered
    result_full = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["stratum_num"]
    )
    
    assert isinstance(result_full, xr.DataArray)
    assert result_full.shape == (3,)
    assert all(result_full ==  1)
    
    # Use a different contrast for filtering
    result_female = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["stratum_num"], 
        include_filter={"sex": ["female"]}
    )
    
    assert isinstance(result_female, xr.DataArray)
    assert result_female.size == 3
    assert all(result_female >= 0), "NASC proportions should be non-negative"
    assert all(result_female <= 1), "NASC proportions should be <= 1"
    
    # Test inverse with exclusion filter
    result_male_excl = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["stratum_num"], 
        exclude_filter={"sex": ["male"]}
    )
    
    assert isinstance(result_male_excl, xr.DataArray)
    assert result_male_excl.shape == (3,)
    assert all(result_male_excl >= 0), "NASC proportions should be non-negative"
    assert all(result_male_excl <= 1), "NASC proportions should be <= 1"
    assert all(result_male_excl == result_female)
    
    # Try different grouping
    result_sex = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["sex"]
    )

    assert isinstance(result_sex, xr.DataArray)
    assert result_sex.shape == (2,)
    assert all(result_sex ==  1)
    
    # Try thresholding based solely on weight
    result_thresh = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["stratum_num"],
        number_proportions=aged_props
    )
    
    assert isinstance(result_thresh, xr.DataArray)
    assert result_thresh.shape == (3,)
    assert all(result_thresh ==  1)
    
    # Try thresholding
    result_thresh = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["stratum_num"],
        weight_proportion_threshold = 1
    )
    
    assert isinstance(result_thresh, xr.DataArray)
    assert result_thresh.shape == (3,)
    assert all(result_thresh ==  0)
    
    # Set minimum lengths
    result_lenthresh = get_proportions.get_weight_proportions_slice(
        weight_proportions=weight_cnv,
        group_columns=["stratum_num"],
        number_proportions=aged_props,
        include_filter={"age_bin": [1]},
        length_threshold_min=20
    )
    
    assert isinstance(result_lenthresh, xr.DataArray)
    assert result_lenthresh.shape == (3,)
    assert all(result_lenthresh >= 0), "NASC proportions should be non-negative"
    assert all(result_lenthresh <= 1), "NASC proportions should be <= 1"
    assert all((result_lenthresh != result) == [True, True, False])
