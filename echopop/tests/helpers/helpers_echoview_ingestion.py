import tempfile
from pathlib import Path

import pandas as pd


# Helper function to create temporary CSV
def create_temp_csv(data):
    """Create a temporary CSV file with provided data."""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w+") as f:
        data.to_csv(f.name, index=False)
        temp_path = f.name

    return Path(temp_path)


# Mock reader function to substitute for read_echoview_nasc
def mock_read_echoview_nasc(filepath, transect_num):
    """Mock function that simulates read_echoview_nasc for testing."""
    # Return different DataFrames based on transect number
    if transect_num == 1.0:
        return pd.DataFrame({"data": [1, 2, 3], "transect_num": 1.0})
    elif transect_num == 2.0:
        return pd.DataFrame({"data": [4, 5, 6], "transect_num": 2.0})
    elif transect_num == 3.0:
        return pd.DataFrame({"data": [7, 8, 9], "transect_num": 3.0})
    else:
        return pd.DataFrame({"data": [0], "transect_num": transect_num})


# Direct test implementation with dependency injection
def echoview_nasc_to_df_test(filtered_df: pd.DataFrame, reader_func) -> list[pd.DataFrame]:
    """Test version with injectable reader function."""
    return [reader_func(row["file_path"], row["transect_num"]) for _, row in filtered_df.iterrows()]


def mock_reader_with_error(filepath, transect_num):
    """Mock reader that raises an error for transect 2.0"""
    if transect_num == 2.0:
        raise ValueError("Test error reading file")
    return pd.DataFrame({"data": [1, 2, 3], "transect_num": transect_num})
