import pandas as pd
import pytest

from echopop.workflows.nwfsc_feat.functions import filter_transect_intervals


def test_filter_transect_intervals_basic():
    """Test basic filtering functionality."""
    # Create test dataframes
    nasc_df = pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "distance_s": [0, 1, 0, 1],
            "distance_e": [1, 2, 1, 2],
            "nasc": [10, 20, 30, 40],
        }
    )

    filter_df = pd.DataFrame(
        {
            "transect_num": [1, 2],
            "log_start": [0.5, 0.5],
            "log_end": [1.5, 1.5],
            "region_id": ["A", "B"],
        }
    )

    # Apply filter
    result = filter_transect_intervals(nasc_df, filter_df)

    # Check results
    assert len(result) == 4  # All rows should match
    assert "nasc" in result.columns
    assert "log_start" not in result.columns  # Should only contain original columns


def test_filter_transect_intervals_with_subset_filter():
    """Test filtering with subset filter."""
    nasc_df = pd.DataFrame(
        {
            "transect_num": [1, 1, 1, 2, 2, 3, 3, 3],
            "distance_s": [0, 1, 2, 0, 1, 0, 1, 2],
            "distance_e": [1, 2, 3, 1, 2, 1, 2, 3],
            "nasc": [10, 20, 30, 40, 50, 60, 70, 80],
        }
    )

    filter_df = pd.DataFrame(
        {
            "transect_num": [1, 2],
            "log_start": [0.5, 0.5],
            "log_end": [1.5, 1.5],
            "region_id": ["A", "B"],
        }
    )

    # Apply filter with subset
    result = filter_transect_intervals(nasc_df, filter_df, "region_id == 'A'")

    # Check results - should include 1 fewer row
    assert len(result) == 7
    assert all(result["transect_num"] == [1, 1, 2, 2, 3, 3, 3])


def test_filter_transect_intervals_no_overlap():
    """Test when there's no overlap between intervals."""
    nasc_df = pd.DataFrame(
        {"transect_num": [1, 1], "distance_s": [0, 1], "distance_e": [1, 2], "nasc": [10, 20]}
    )

    filter_df = pd.DataFrame(
        {
            "transect_num": [1],
            "log_start": [3],  # No overlap with any interval
            "log_end": [4],
            "region_id": ["A"],
        }
    )

    # Apply filter
    result = filter_transect_intervals(nasc_df, filter_df)

    # Check results - should be empty
    assert len(result) == 0


def test_filter_transect_intervals_column_rename():
    """Test automatic renaming of 'transect' to 'transect_num'."""
    nasc_df = pd.DataFrame(
        {"transect_num": [1, 2], "distance_s": [0, 0], "distance_e": [1, 1], "nasc": [10, 20]}
    )

    # Use 'transect' instead of 'transect_num'
    filter_df = pd.DataFrame(
        {
            "transect": [1, 2],  # Should be automatically renamed
            "log_start": [0, 0],
            "log_end": [1, 1],
        }
    )

    result = filter_transect_intervals(nasc_df, filter_df)

    # Should match both rows
    assert len(result) == 2
    assert set(result["transect_num"]) == {1, 2}


def test_filter_transect_intervals_invalid_column():
    """Test error handling with invalid column in subset filter."""
    nasc_df = pd.DataFrame({"transect_num": [1, 2], "distance_s": [0, 0], "distance_e": [1, 1]})

    filter_df = pd.DataFrame({"transect_num": [1, 2], "log_start": [0, 0], "log_end": [1, 1]})

    # Use non-existent column in filter
    with pytest.raises(pd.errors.UndefinedVariableError):
        filter_transect_intervals(nasc_df, filter_df, "nonexistent_column == 'A'")
