import pandas as pd
import pytest

from echopop.nwfsc_feat import utils


def test_get_transect_basic_functionality(basic_nasc_data):
    """Test basic functionality of utils.compute_interval_distance"""
    utils.compute_interval_distance(basic_nasc_data)

    # Check that distance_interval column was added
    assert "distance_interval" in basic_nasc_data.columns

    # Check that original columns are preserved
    assert "distance_s" in basic_nasc_data.columns
    assert "distance_e" in basic_nasc_data.columns
    assert "transect_spacing" in basic_nasc_data.columns

    # Check expected interval values (forward differences)
    expected_intervals = [1.0, 1.0, 1.0, 1.0, 1.0]  # Last one calculated as distance_e - distance_s
    pd.testing.assert_series_equal(
        basic_nasc_data["distance_interval"],
        pd.Series(expected_intervals, name="distance_interval"),
        check_dtype=False,
    )


def test_get_transect_single_row(transect_single_row):
    """Test utils.compute_interval_distance with single row DataFrame"""
    utils.compute_interval_distance(transect_single_row)

    # Check that distance_interval column was added
    assert "distance_interval" in transect_single_row.columns

    # For single row, interval should be distance_e - distance_s
    expected_interval = 1.5
    assert transect_single_row["distance_interval"].iloc[0] == expected_interval


def test_get_transect_preserves_original(basic_nasc_data):
    """Test that original DataFrame is not modified"""
    original_columns = basic_nasc_data.columns.tolist()
    utils.compute_interval_distance(basic_nasc_data)

    # Original DataFrame should not have distance_interval column
    assert set(basic_nasc_data.columns) - set(original_columns) == {"distance_interval"}

    # Result should have the new column
    assert "distance_interval" in basic_nasc_data.columns


def test_get_transect_irregular_spacing(irregular_spacing_data):
    """Test utils.compute_interval_distance with irregular spacing"""
    utils.compute_interval_distance(irregular_spacing_data)

    # Check that distance_interval column was added
    assert "distance_interval" in irregular_spacing_data.columns

    # Check that all values are non-negative
    assert (irregular_spacing_data["distance_interval"] >= 0).all()

    # Check that last interval uses distance_e - distance_s
    last_interval = irregular_spacing_data["distance_interval"].iloc[-1]
    expected_last = (
        irregular_spacing_data["distance_e"].iloc[-1]
        - irregular_spacing_data["distance_s"].iloc[-1]
    )
    assert last_interval == expected_last


def test_get_transect_empty_dataframe():
    """Test utils.compute_interval_distance with empty DataFrame"""
    empty_df = pd.DataFrame(columns=["distance_s", "distance_e", "transect_spacing"])

    with pytest.raises(IndexError):
        utils.compute_interval_distance(empty_df)


def test_get_transect_missing_columns():
    """Test utils.compute_interval_distance with missing required columns"""
    incomplete_df = pd.DataFrame({"distance_s": [1, 2, 3]})

    with pytest.raises(KeyError):
        utils.compute_interval_distance(incomplete_df)
