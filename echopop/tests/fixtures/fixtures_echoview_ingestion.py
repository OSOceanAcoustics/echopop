import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


# ==================================================================================================
# Mock Paths
# ----------
@pytest.fixture
def mock_export_paths():
    """Create mock Path objects for testing."""
    paths = {
        "analysis": [
            Path("T1_survey_data_(analysis).csv"),
            Path("T2_survey_data_(analysis).csv"),
            Path("T3_survey_data_(analysis).csv"),
        ],
        "cells": [
            Path("T1_survey_data_(cells).csv"),
            Path("T2_survey_data_(cells).csv"),
            Path("T3_survey_data_(cells).csv"),
        ],
        "intervals": [
            Path("T1_survey_data_(intervals).csv"),
            Path("T2_survey_data_(intervals).csv"),
            Path("T3_survey_data_(intervals).csv"),
        ],
        "layers": [
            Path("T1_survey_data_(layers).csv"),
            Path("T2_survey_data_(layers).csv"),
            # T3 missing layers file
        ],
    }

    # Convert lists to generators to mimic Path.glob behavior
    return {key: (path for path in path_list) for key, path_list in paths.items()}


# ==================================================================================================
# Retrieving transect numbers and validating filesets
# ---------------------------------------------------
@pytest.fixture
def mock_transect_df():
    """Create a mock transect DataFrame for testing validate_transect_exports."""
    data = [
        # Complete transect T1
        {
            "file_type": "analysis",
            "file_path": Path("/mock/Transect_T1_analysis.csv"),
            "transect_num": 1.0,
        },
        {
            "file_type": "cells",
            "file_path": Path("/mock/Transect_T1_cells.csv"),
            "transect_num": 1.0,
        },
        {
            "file_type": "intervals",
            "file_path": Path("/mock/Transect_T1_intervals.csv"),
            "transect_num": 1.0,
        },
        {
            "file_type": "layers",
            "file_path": Path("/mock/Transect_T1_layers.csv"),
            "transect_num": 1.0,
        },
        # Complete transect T2
        {
            "file_type": "analysis",
            "file_path": Path("/mock/Transect_T2_analysis.csv"),
            "transect_num": 2.0,
        },
        {
            "file_type": "cells",
            "file_path": Path("/mock/Transect_T2_cells.csv"),
            "transect_num": 2.0,
        },
        {
            "file_type": "intervals",
            "file_path": Path("/mock/Transect_T2_intervals.csv"),
            "transect_num": 2.0,
        },
        {
            "file_type": "layers",
            "file_path": Path("/mock/Transect_T2_layers.csv"),
            "transect_num": 2.0,
        },
        # Incomplete transect T3 (missing layers)
        {
            "file_type": "analysis",
            "file_path": Path("/mock/Transect_T3_analysis.csv"),
            "transect_num": 3.0,
        },
        {
            "file_type": "cells",
            "file_path": Path("/mock/Transect_T3_cells.csv"),
            "transect_num": 3.0,
        },
        {
            "file_type": "intervals",
            "file_path": Path("/mock/Transect_T3_intervals.csv"),
            "transect_num": 3.0,
        },
    ]
    return pd.DataFrame(data)


# ==================================================================================================
# Fileset reading generator
# -------------------------
@pytest.fixture
def mock_filtered_df():
    """Create a mock filtered DataFrame for testing."""
    return pd.DataFrame(
        {
            "file_type": ["intervals", "intervals", "intervals"],
            "file_path": [
                Path("/mock/T1_survey_data_(intervals).csv"),
                Path("/mock/T2_survey_data_(intervals).csv"),
                Path("/mock/T3_survey_data_(intervals).csv"),
            ],
            "transect_num": [1.0, 2.0, 3.0],
        }
    )


@pytest.fixture
def mock_empty_df():
    """Create an empty DataFrame with the correct column structure."""
    return pd.DataFrame(columns=["file_type", "file_path", "transect_num"])


# ==================================================================================================
# Coordinate imputation/repair
# ----------------------------
@pytest.fixture
def bad_coords_at_start():
    return pd.DataFrame(
        {
            "latitude": [999.0, 999.0, 45.2, 45.3, 45.4],
            "longitude": [-125.1, -125.2, -125.3, -125.4, -125.5],
        }
    )


@pytest.fixture
def bad_coords_in_middle():
    return pd.DataFrame(
        {
            "latitude": [45.1, 45.2, 999.0, 999.0, 45.5],
            "longitude": [-125.1, -125.2, -125.3, -125.4, -125.5],
        }
    )


@pytest.fixture
def bad_coords_at_end():
    return pd.DataFrame(
        {
            "latitude": [45.1, 45.2, 45.3, 999.0, 999.0],
            "longitude": [-125.1, -125.2, -125.3, -125.4, -125.5],
        }
    )


@pytest.fixture
def multiple_bad_coord_groups():
    return pd.DataFrame(
        {
            "latitude": [999.0, 45.1, 45.2, 999.0, 999.0, 45.5, 999.0],
            "longitude": [-125.1, -125.2, -125.3, -125.4, -125.5, -125.6, -125.7],
        }
    )


@pytest.fixture
def no_bad_coords():
    return pd.DataFrame(
        {
            "latitude": [45.1, 45.2, 45.3, 45.4, 45.5],
            "longitude": [-125.1, -125.2, -125.3, -125.4, -125.5],
        }
    )


# ==================================================================================================
# Fix region column strings
# -------------------------
@pytest.fixture
def test_cells_df():
    """Create test DataFrame for clean_echoview_cells_df tests."""
    return pd.DataFrame(
        {
            "region_class": ['"Fish"', ' "Plankton" ', "  Krill  ", np.nan],
            "region_name": ['"Region1"', " Area51 ", "  Zone99  ", np.nan],
            "other_column": [1, 2, 3, 4],
        }
    )


# ==================================================================================================
# Sorting and reindexing the export data rows
# -------------------------------------------
@pytest.fixture
def test_export_df():
    """Create a test DataFrame for Echoview export sorting."""
    return pd.DataFrame(
        {"interval": [3, 1, 2], "region_name": ["B", "A", "C"], "data_value": [10, 20, 30]}
    )


# ==================================================================================================
# Transect spacing imputation/calculation
# ---------------------------------------
@pytest.fixture
def test_transect_data():
    """Create test transect data for spacing calculations."""
    # Create data for 4 transects
    return pd.DataFrame(
        {
            "transect_num": [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0],
            "latitude": [34.5, 34.6, 34.7, 34.8, 34.9, 35.0, 35.1, 35.2],
            "longitude": [-121.1, -121.2, -121.3, -121.4, -121.5, -121.6, -121.7, -121.8],
            "process_id": [148293] * 8,
            "interval": range(1489, 1497),
        }
    )


# ==================================================================================================
# Merge Echoview exports into a consolidated file
# -----------------------------------------------
@pytest.fixture
def sample_intervals_df():
    """Create a sample intervals DataFrame."""
    return pd.DataFrame(
        {
            "transect_num": [1.0, 1.0, 2.0],
            "interval": [1, 2, 3],
            "process_id": [101, 101, 102],
            "latitude": [34.5, 34.6, 34.7],
            "longitude": [-121.1, -121.2, -121.3],
            "vessel_log_start": [100.1, 100.2, 200.1],
        }
    )


@pytest.fixture
def sample_cells_df():
    """Create a sample cells DataFrame."""
    return pd.DataFrame(
        {
            "transect_num": [1.0, 1.0, 2.0],
            "interval": [1, 2, 3],
            "process_id": [101, 101, 102],
            "region_class": ["fish", "plankton", "fish"],
            "region_name": ["A", "B", "C"],
            "nasc": [10.5, 15.2, 20.7],
        }
    )


@pytest.fixture
def sample_layers_df():
    """Create a sample layers DataFrame."""
    return pd.DataFrame(
        {
            "transect_num": [1.0, 1.0, 2.0],
            "interval": [1, 2, 3],
            "process_id": [101, 101, 102],
            "layer_name": ["surface", "middle", "bottom"],
            "min_depth": [0.0, 10.5, 0.0],
            "max_depth": [10.5, 50.2, 25.8],
        }
    )


# ==================================================================================================
# Read in, clean, and merge Echoview exports (database format)
# ------------------------------------------------------------
@pytest.fixture
def mock_nasc_directory():
    """Create a temporary directory with mock Echoview export files."""
    # Create a temporary directory
    temp_dir = tempfile.mkdtemp()

    # Create sample data for different file types
    intervals_data = pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-02"],
            "lat_s": [45.1, 45.2, 46.1],
            "lon_s": [-125.1, -125.2, -126.1],
            "interval": [1, 2, 1],
            "process_id": [101, 101, 102],
            "vl_start": [10.0, 10.5, 20.0],
            "vl_end": [10.5, 11.0, 20.5],
        }
    )

    cells_data = pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-02"],
            "lat_s": [45.1, 45.2, 46.1],
            "lon_s": [-125.1, -125.2, -126.1],
            "interval": [1, 2, 1],
            "process_id": [101, 101, 102],
            "region_class": ['"Fish"', '"Plankton"', '"Fish"'],
            "region_name": ['"RegionA"', '"RegionB"', '"RegionA"'],
            "prc_nasc": [125.5, 130.2, 140.8],
        }
    )

    layers_data = pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-02"],
            "lat_s": [45.1, 45.2, 46.1],
            "lon_s": [-125.1, -125.2, -126.1],
            "interval": [1, 2, 1],
            "process_id": [101, 101, 102],
            "layer_name": ["Surface", "Bottom", "Surface"],
            "layer_depth_min": [0.0, 50.0, 0.0],
            "layer_depth_max": [50.0, 100.0, 50.0],
        }
    )

    analysis_data = pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-02"],
            "lat_s": [45.1, 45.2, 46.1],
            "lon_s": [-125.1, -125.2, -126.1],
            "interval": [1, 2, 1],
            "process_id": [101, 101, 102],
        }
    )

    # Create files for two transects
    for t in [1, 2]:
        # Save each file type
        intervals_file = os.path.join(temp_dir, f"T{t}_data_(intervals).csv")
        cells_file = os.path.join(temp_dir, f"T{t}_data_(cells).csv")
        layers_file = os.path.join(temp_dir, f"T{t}_data_(layers).csv")
        analysis_file = os.path.join(temp_dir, f"T{t}_data_(analysis).csv")

        # Filter data by transect (using process_id as a proxy)
        transect_filter = intervals_data["process_id"] == (100 + t)

        # Save to files
        intervals_data[transect_filter].to_csv(intervals_file, index=False)
        cells_data[transect_filter].to_csv(cells_file, index=False)
        layers_data[transect_filter].to_csv(layers_file, index=False)
        analysis_data[transect_filter].to_csv(analysis_file, index=False)

    yield Path(temp_dir)

    # Clean up
    shutil.rmtree(temp_dir)


# ==================================================================================================
# Read transect-region-haul key file
# ----------------------------------
@pytest.fixture
def sample_mapping_data():
    """Create sample data for transect-region-haul mapping."""
    return pd.DataFrame(
        {
            "transect_num": [1, 2, 3, 4],
            "region_id": ["A1", "B2", "C3", "D4"],
            "haul_num": [101, 102, 103, 104],
            "extra_col": ["extra1", "extra2", "extra3", "extra4"],
        }
    )


@pytest.fixture
def csv_mapping_file(sample_mapping_data):
    """Create a temporary CSV file with sample mapping data."""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as temp_file:
        temp_filename = temp_file.name

    sample_mapping_data.to_csv(temp_filename, index=False)

    yield Path(temp_filename)

    # Clean up
    os.unlink(temp_filename)


@pytest.fixture
def excel_mapping_file(sample_mapping_data):
    """Create a temporary Excel file with sample mapping data."""
    with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as temp_file:
        temp_filename = temp_file.name

    # Create Excel file with the data
    with pd.ExcelWriter(temp_filename) as writer:
        sample_mapping_data.to_excel(writer, sheet_name="MappingSheet", index=False)

    yield Path(temp_filename)

    # Clean up
    os.unlink(temp_filename)


# ==================================================================================================
# Reading Echoview data files
# ---------------------------
@pytest.fixture
def sample_echoview_data():
    """Create sample Echoview export data for testing."""
    # Create data with original Echoview column names
    return pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-01"],
            "exclude_below_line_depth_mean": [200.5, 198.2, 201.3],
            "lat_s": [45.1, 45.2, 45.3],
            "lon_s": [-125.1, -125.2, -125.3],
            "prc_nasc": [120.5, 135.8, 118.2],
            "time_s": ["10:15:30", "10:16:00", "10:16:30"],
            "vl_end": [10.5, 11.0, 11.5],
            "vl_start": [10.0, 10.5, 11.0],
            "other_column": ["a", "b", "c"],  # Column not in mapping
        }
    )


@pytest.fixture
def echoview_temp_csv(sample_echoview_data):
    """Create a temporary CSV file with sample Echoview data."""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w+") as f:
        sample_echoview_data.to_csv(f.name, index=False)
        temp_path = f.name

    yield Path(temp_path)

    # Clean up the temporary file
    os.unlink(temp_path)


@pytest.fixture
def empty_echoview_data():
    """Create empty dataframe with Echoview headers."""
    return pd.DataFrame(
        columns=[
            "date_s",
            "exclude_below_line_depth_mean",
            "lat_s",
            "lon_s",
            "prc_nasc",
            "time_s",
            "vl_end",
            "vl_start",
        ]
    )


@pytest.fixture
def missing_columns_data():
    """Create dataframe with some columns missing."""
    return pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01"],
            # Missing exclude_below_line_depth_mean
            "lat_s": [45.1, 45.2],
            "lon_s": [-125.1, -125.2],
            # Missing prc_nasc
            "time_s": ["10:15:30", "10:16:00"],
            # Missing vl_end
            "vl_start": [10.0, 10.5],
        }
    )


@pytest.fixture
def extreme_values_data():
    """Create dataframe with extreme/edge values."""
    return pd.DataFrame(
        {
            "date_s": ["9999-12-31", "0001-01-01", "2023-01-01"],
            "exclude_below_line_depth_mean": [np.inf, -np.inf, np.nan],
            "lat_s": [90.0, -90.0, 0.0],  # Max, min, and zero latitude
            "lon_s": [180.0, -180.0, 0.0],  # Max, min, and zero longitude
            "prc_nasc": [1e10, -1e10, 0.0],  # Very large, very negative, zero
            "time_s": ["23:59:59", "00:00:00", "12:00:00"],
            "vl_end": [1e6, -1e6, 0.0],
            "vl_start": [1e6, -1e6, 0.0],
        }
    )


@pytest.fixture
def duplicate_columns_data():
    """Create dataframe with duplicate column names (after reading)."""
    # This will simulate a CSV with duplicate columns after reading
    df = pd.DataFrame()
    df["date_s"] = ["2023-01-01", "2023-01-01"]
    df["lat_s"] = [45.1, 45.2]
    df["lat_s.1"] = [45.3, 45.4]  # Duplicate column that pandas would create
    df["lon_s"] = [-125.1, -125.2]
    return df


@pytest.fixture
def bad_coords_middle():
    """Create dataframe with bad coordinates in middle."""
    return pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-01", "2023-01-01", "2023-01-01"],
            "lat_s": [45.1, 45.2, 999.0, 999.0, 45.5],
            "lon_s": [-125.1, -125.2, -125.3, -125.4, -125.5],
        }
    )


@pytest.fixture
def bad_coords_end():
    """Create dataframe with bad coordinates at end."""
    return pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-01", "2023-01-01", "2023-01-01"],
            "lat_s": [45.1, 45.2, 45.3, 999.0, 999.0],
            "lon_s": [-125.1, -125.2, -125.3, -125.4, -125.5],
        }
    )


@pytest.fixture
def bad_coords_mixed():
    """Create dataframe with mixed patterns of bad coordinates."""
    return pd.DataFrame(
        {
            "date_s": [
                "2023-01-01",
                "2023-01-01",
                "2023-01-01",
                "2023-01-01",
                "2023-01-01",
                "2023-01-01",
                "2023-01-01",
            ],
            "lat_s": [999.0, 45.1, 45.2, 999.0, 999.0, 45.5, 999.0],
            "lon_s": [-125.1, -125.2, -125.3, -125.4, -125.5, -125.6, -125.7],
        }
    )


@pytest.fixture
def bad_coords_all():
    """Create dataframe with all bad coordinates."""
    return pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-01"],
            "lat_s": [999.0, 999.0, 999.0],
            "lon_s": [999.0, 999.0, 999.0],
            "prc_nasc": [100, 200, 300],
        }
    )


@pytest.fixture
def bad_coords_start():
    """Create dataframe with bad coordinates at start."""
    return pd.DataFrame(
        {
            "date_s": ["2023-01-01", "2023-01-01", "2023-01-01", "2023-01-01", "2023-01-01"],
            "lat_s": [999.0, 999.0, 45.2, 45.3, 45.4],
            "lon_s": [-125.1, -125.2, -125.3, -125.4, -125.5],
        }
    )


# ==================================================================================================
# Export layer intervals from Echoview exports
# --------------------------------------------
@pytest.fixture
def sample_transect_data():
    """Create sample transect data for testing."""
    return pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "interval": [1, 2, 1, 2],
            "max_depth": [100, 120, 90, 110],
            "layer_depth_min": [10, 20, 15, 25],
            "layer_depth_max": [50, 60, 45, 55],
        }
    )


# ==================================================================================================
# Consolidated XLSX file
# ----------------------
@pytest.fixture
def sample_excel_path(tmp_path):
    """Create a simple Excel file for testing"""
    df = pd.DataFrame(
        {
            "transect": [1, 2, 3],
            "region id": ["A", "B", "C"],
            "latitude": [45.0, 46.0, 47.0],
            "longitude": [-125.0, -126.0, -127.0],
        }
    )

    # Create file in the pytest-provided temp directory
    file_path = tmp_path / "test_data.xlsx"
    df.to_excel(file_path, sheet_name="Sheet1", index=False)

    # Return path as string
    return str(file_path)


# ==================================================================================================
# Transect-region-haul key map generator
# --------------------------------------
@pytest.fixture
def sample_region_data():
    """Create sample region data for testing."""
    return pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2, 3],
            "haul_num": [1, 1, 2, 2, 3],
            "region_id": [1, 2, 1, 2, 1],
            "region_class": ["Age-1 Hake", "Hake", "Age-1 Hake Mix", "Age-1 Hake", "Hake Mix"],
            "region_name": ["R1", "R2", "R3", "R4", "R5"],
        }
    )


# ==================================================================================================
# NASC consolidation
# ------------------
@pytest.fixture
def sample_merged_data():
    """Create sample merged Echoview data."""
    return pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "interval": [1, 2, 1, 2],
            "region_id": [1, 2, 1, 2],
            "region_class": ["Hake", "Age-1 Hake", "Hake Mix", "Age-1 Hake"],
            "nasc": [100.0, 150.0, 200.0, 250.0],
            "layer_depth_min": [10, 20, 15, 25],
            "layer_depth_max": [50, 60, 45, 55],
            "max_depth": [100, 120, 90, 110],
        }
    )


@pytest.fixture
def sample_interval_data():
    """Create sample interval data."""
    return pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "interval": [1, 2, 1, 2],
            "distance_s": [0, 1, 0, 1],
            "distance_e": [1, 2, 1, 2],
            "latitude": [45.0, 45.1, 45.2, 45.3],
            "longitude": [-124.0, -124.1, -124.2, -124.3],
            "transect_spacing": [1.0, 1.0, 1.0, 1.0],
        }
    )


@pytest.fixture
def sample_haul_key():
    """Create sample haul key data."""
    return pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "interval": [1, 2, 1, 2],
            "region_id": [1, 2, 1, 2],
            "haul_num": [10, 11, 12, 13],
        }
    )


@pytest.fixture
def sample_sv_csv_file():
    """Create a temporary SV CSV file for testing."""
    df = pd.DataFrame(
        {
            "Process_ID": [7217, 7217, 7217, 7217],
            "Interval": [6, 6, 7, 7],
            "Layer": [2, 3, 2, 3],
            "Sv_mean": [-83.381331, -75.2, -82.1, -999],
            "NASC": [1.979553, 2.5, 1.8, 0],
            "Thickness_mean": [6.0088, 10.0048, 6.0088, 10.0048],
            "Depth_mean": [16.994, 25.0008, 16.994, 25.0008],
            "Lat_M": [34.64083725, 34.64083725, 34.64064364, 34.64064364],
            "Lon_M": [-120.6977052, -120.6977052, -120.7081345, -120.7081345],
            "Frequency": [18, 18, 18, 18],
            "VL_start": [2.504577, 2.504577, 3.004596, 3.004596],
            "VL_end": [3.001603, 3.001603, 3.503897, 3.503897],
        }
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        df.to_csv(f.name, index=False)
        filename = f.name
    yield Path(filename)
    os.unlink(filename)


@pytest.fixture
def empty_sv_csv_file():
    """Create an empty SV CSV file for testing."""
    df = pd.DataFrame(
        columns=["Process_ID", "Interval", "Layer", "Sv_mean", "NASC", "Thickness_mean"]
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        df.to_csv(f.name, index=False)
        filename = f.name
    yield Path(filename)
    os.unlink(filename)


@pytest.fixture
def sv_csv_with_coordinates():
    """Create SV CSV file with latitude/longitude columns for coordinate testing."""
    df = pd.DataFrame(
        {
            "Process_ID": [7217, 7217],
            "Interval": [6, 6],
            "Layer": [2, 3],
            "Sv_mean": [-83.381331, -75.2],
            "NASC": [1.979553, 2.5],
            "Thickness_mean": [6.0088, 10.0048],
            "latitude": [34.64083725, 34.64083725],
            "longitude": [-120.6977052, -120.6977052],
            "Frequency": [18, 18],
        }
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        df.to_csv(f.name, index=False)
        filename = f.name
    yield Path(filename)
    os.unlink(filename)


@pytest.fixture
def sv_directory_with_files(tmp_path):
    """Create a directory with multiple SV CSV files for testing."""
    sv_dir = tmp_path / "sv_data"
    sv_dir.mkdir()

    # Create first CSV file
    df1 = pd.DataFrame(
        {
            "Process_ID": [7217, 7217],
            "Interval": [6, 6],
            "Layer": [2, 3],
            "Sv_mean": [-83.381331, -75.2],
            "NASC": [1.979553, 2.5],
            "Thickness_mean": [6.0088, 10.0048],
            "Lat_M": [34.641, 34.641],
            "Lon_M": [-120.698, -120.698],
            "Frequency": [18, 18],
            "VL_start": [100.0, 100.0],
            "VL_end": [200.0, 200.0],
        }
    )

    # Create second CSV file
    df2 = pd.DataFrame(
        {
            "Process_ID": [7218, 7218],
            "Interval": [7, 7],
            "Layer": [2, 3],
            "Sv_mean": [-82.1, -78.5],
            "NASC": [1.8, 2.2],
            "Thickness_mean": [6.0088, 10.0048],
            "Lat_M": [34.642, 34.642],
            "Lon_M": [-120.699, -120.699],
            "Frequency": [18, 18],
            "VL_start": [200.0, 200.0],
            "VL_end": [300.0, 300.0],
        }
    )

    file1 = sv_dir / "x01.csv"
    file2 = sv_dir / "x02.csv"

    df1.to_csv(file1, index=False)
    df2.to_csv(file2, index=False)

    return sv_dir
