import numpy as np
import pandas as pd

from echopop import utils


def test_apply_filters_include():
    """Test include filter functionality."""
    df = pd.DataFrame(
        {"species": ["cod", "hake", "pollock", "cod", "hake"], "length": [10, 15, 20, 25, 30]}
    )

    result = utils.apply_filters(df, include_filter={"species": "cod"})
    assert len(result) == 2
    assert all(result["species"] == "cod")

    result = utils.apply_filters(df, include_filter={"species": ["cod", "hake"]})
    assert len(result) == 4
    assert all(result["species"].isin(["cod", "hake"]))


def test_apply_filters_exclude():
    """Test exclude filter functionality."""
    df = pd.DataFrame(
        {"species": ["cod", "hake", "pollock", "cod", "hake"], "length": [10, 15, 20, 25, 30]}
    )

    result = utils.apply_filters(df, exclude_filter={"species": "cod"})
    assert len(result) == 3
    assert not any(result["species"] == "cod")

    result = utils.apply_filters(df, exclude_filter={"species": ["cod", "hake"]})
    assert len(result) == 1
    assert all(result["species"] == "pollock")


def test_apply_filters_combined():
    """Test combining include and exclude filters."""
    df = pd.DataFrame(
        {"species": ["cod", "hake", "pollock", "cod", "hake"], "length": [10, 15, 20, 25, 30]}
    )

    result = utils.apply_filters(
        df, include_filter={"species": ["cod", "hake"]}, exclude_filter={"length": 30}
    )
    assert len(result) == 3
    assert all(result["species"].isin(["cod", "hake"]))
    assert all(result["length"] != 30)


def test_group_interpolator_creator_global(length_weight_dataset_with_bins):
    """Test creating a global interpolator with no contrast variables."""
    interpolators = utils.group_interpolator_creator(
        grouped_data=length_weight_dataset_with_bins,
        independent_var="length",
        dependent_var="weight_fitted",
        contrast_vars=None,
    )

    assert len(interpolators) == 1
    assert "_global_" in interpolators
    assert callable(interpolators["_global_"])

    # Test the interpolator function
    global_interp = interpolators["_global_"]
    assert np.isclose(global_interp(15.0), 15**3 * 0.01, rtol=1e-5)


def test_group_interpolator_creator_by_sex(length_weight_dataset_with_bins):
    """Test creating interpolators grouped by sex."""
    interpolators = utils.group_interpolator_creator(
        grouped_data=length_weight_dataset_with_bins,
        independent_var="length",
        dependent_var="weight_fitted",
        contrast_vars="sex",
    )

    assert len(interpolators) == 2
    assert "male" in interpolators
    assert "female" in interpolators
    assert callable(interpolators["male"])
    assert callable(interpolators["female"])

    # Test interpolator functions
    male_df = length_weight_dataset_with_bins[length_weight_dataset_with_bins["sex"] == "male"]
    male_length = male_df["length"].iloc[0]
    male_weight = male_df["weight_fitted"].iloc[0]

    assert np.isclose(interpolators["male"](male_length), male_weight, rtol=1e-5)


def test_group_interpolator_creator_multiple_groups(length_weight_dataset_with_bins):
    """Test creating interpolators with multiple grouping variables."""
    interpolators = utils.group_interpolator_creator(
        grouped_data=length_weight_dataset_with_bins,
        independent_var="length",
        dependent_var="weight_fitted",
        contrast_vars=["sex", "stratum"],
    )

    assert len(interpolators) == 2  # male-A, female-B combinations in the fixture
    assert ("male", "A") in interpolators or ("male", "stratum_A") in interpolators
    assert ("female", "B") in interpolators or ("female", "stratum_B") in interpolators
    assert callable(interpolators[list(interpolators.keys())[0]])


def test_group_interpolator_creator_insufficient_data():
    """Test handling of groups with insufficient data for interpolation."""
    df = pd.DataFrame(
        {"length": [10.0, 15.0], "weight_fitted": [5.0, 10.0], "sex": ["male", "female"]}
    )

    interpolators = utils.group_interpolator_creator(
        grouped_data=df,
        independent_var="length",
        dependent_var="weight_fitted",
        contrast_vars="sex",
    )

    # Function returns empty dictionary when all groups have insufficient data
    assert isinstance(interpolators, dict)
    assert len(interpolators) == 0


def test_create_grouped_series(proportion_dict):
    """Test the create_grouped_series function."""
    result = utils.create_grouped_series(
        proportion_dict, ["stratum_num", "sex"], "proportion_overall"
    )

    # Check shape and column types
    assert result.shape[0] == 8  # 2 strata * 2 sexes * 2 groups
    assert set(result.columns) == {"stratum_num", "sex", "proportion_overall", "group"}

    # Check that we have both groups in the result
    assert set(result["group"].unique()) == {"aged", "unaged"}

    # Check aggregation values for specific combinations
    female_aged_strat1 = result[
        (result["group"] == "aged") & (result["sex"] == "female") & (result["stratum_num"] == 1)
    ]["proportion_overall"].iloc[0]
    assert female_aged_strat1 == 0.15


def test_create_pivot_table(proportion_dict):
    """Test the create_pivot_table function."""
    grouped = utils.create_grouped_series(
        proportion_dict, ["stratum_num", "sex"], "proportion_overall"
    )

    result = utils.create_pivot_table(
        grouped, ["group", "sex"], ["stratum_num"], "proportion_overall"
    )

    # Check that the pivot table has the right shape and structure
    assert isinstance(result.index, pd.MultiIndex)
    assert result.index.names == ["group", "sex"]
    assert result.columns.names == ["stratum_num"]

    # Check specific values
    assert result.loc[("aged", "female"), 1] == 0.15
    assert result.loc[("unaged", "male"), 2] == 0.1


def test_create_grouped_table(proportion_dict):
    """Test the create_grouped_table function."""
    result = utils.create_grouped_table(
        proportion_dict, ["stratum_num", "sex"], ["group"], ["stratum_num"], "proportion_overall"
    )

    # Check that we get a pivot table with the right structure
    assert result.index.name == "group"
    assert result.columns.name == "stratum_num"

    # Check the values (sum of proportions by group and stratum)
    assert result.loc["aged", 1] == 0.4  # 0.15 + 0.2 + 0.05
    assert result.loc["unaged", 2] == 0.2  # 0.05 + 0.1 + 0.05
