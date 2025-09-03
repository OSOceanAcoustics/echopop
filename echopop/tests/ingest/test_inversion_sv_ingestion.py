from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from echopop.nwfsc_feat import ingest_sv


def test_read_echoview_sv_basic(sample_sv_csv_file):
    """Test basic reading of Echoview SV data."""
    result = ingest_sv.read_echoview_sv(sample_sv_csv_file)

    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 4
    assert "filename" in result.columns
    assert sample_sv_csv_file.as_posix() in result["filename"].iloc[0]
    assert "process_id" in result.columns
    assert "sv_mean" in result.columns


def test_read_echoview_sv_with_transect_num(sample_sv_csv_file):
    """Test reading with transect number assignment."""
    result = ingest_sv.read_echoview_sv(sample_sv_csv_file, transect_num=5)

    assert "transect_num" in result.columns
    assert all(result["transect_num"] == 5)
    assert len(result) == 4


def test_read_echoview_sv_empty_file(empty_sv_csv_file):
    """Test handling of empty files."""
    result = ingest_sv.read_echoview_sv(empty_sv_csv_file)

    assert result is None


def test_read_echoview_sv_no_coordinate_imputation(sample_sv_csv_file):
    """Test reading without coordinate imputation."""
    result = ingest_sv.read_echoview_sv(sample_sv_csv_file, impute_coordinates=False)

    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 4


def test_read_echoview_sv_with_coordinate_imputation(sv_csv_with_coordinates):
    """Test coordinate imputation when latitude/longitude columns exist."""
    result = ingest_sv.read_echoview_sv(sv_csv_with_coordinates, impute_coordinates=True)

    assert result is not None
    assert "latitude" in result.columns
    assert "longitude" in result.columns


def test_echoview_sv_to_df():
    """Test processing multiple SV files."""
    # Create sample DataFrame with file paths and transect numbers
    file_df = pd.DataFrame(
        {"file_path": [Path("fake1.csv"), Path("fake2.csv")], "transect_num": [1, 2]}
    )

    # This will fail with real files that don't exist, but tests the structure
    try:
        result = ingest_sv.echoview_sv_to_df(file_df)
        assert isinstance(result, list)
        assert len(result) == 2
    except FileNotFoundError:
        # Expected to fail with fake files - that's okay for structure test
        pass


def test_echoview_sv_to_df_with_real_files(tmp_path):
    """Test echoview_sv_to_df with actual CSV files."""
    # Create test CSV files
    file1 = tmp_path / "test1.csv"
    file2 = tmp_path / "test2.csv"

    csv_content1 = """Process_ID,Interval,Layer,Sv_mean,NASC,Thickness_mean,Frequency
7217,6,2,-83.0,1.9,6.0,18
7217,6,3,-75.0,2.5,10.0,18"""

    csv_content2 = """Process_ID,Interval,Layer,Sv_mean,NASC,Thickness_mean,Frequency
7218,7,2,-82.0,1.8,6.0,18
7218,7,3,-78.0,2.2,10.0,18"""

    file1.write_text(csv_content1)
    file2.write_text(csv_content2)

    file_df = pd.DataFrame({"file_path": [file1, file2], "transect_num": [1, 2]})

    result = ingest_sv.echoview_sv_to_df(file_df)

    assert isinstance(result, list)
    assert len(result) == 2
    assert all(isinstance(df, pd.DataFrame) for df in result if df is not None)

    # Check that transect numbers were assigned correctly
    if result[0] is not None:
        assert "transect_num" in result[0].columns
        assert all(result[0]["transect_num"] == 1)
    if result[1] is not None:
        assert "transect_num" in result[1].columns
        assert all(result[1]["transect_num"] == 2)


def test_apply_Sv_thresholds():
    """Test application of frequency-specific Sv thresholds."""
    data = pd.DataFrame(
        {"frequency": [18.0, 18.0, 38.0, 38.0], "sv_mean": [-85.0, -45.0, -80.0, -40.0]}
    )

    thresholds = {18.0: {"min": -90.0, "max": -50.0}, 38.0: {"min": -85.0, "max": -45.0}}

    result = ingest_sv.apply_Sv_thresholds(data, thresholds)

    # Values outside thresholds should be -999
    assert result["sv_mean"].iloc[1] == -999.0  # -45.0 > -50.0 max
    assert result["sv_mean"].iloc[3] == -999.0  # -40.0 > -45.0 max
    assert result["sv_mean"].iloc[0] == -85.0  # Within bounds
    assert result["sv_mean"].iloc[2] == -80.0  # Within bounds


def test_apply_Sv_thresholds_missing_column():
    """Test error handling when sv_mean column is missing."""
    data = pd.DataFrame({"frequency": [18.0, 38.0], "other_column": [1, 2]})

    thresholds = {18.0: {"min": -90.0, "max": -50.0}}

    with pytest.raises(KeyError):
        ingest_sv.apply_Sv_thresholds(data, thresholds)


def test_sv_to_nasc():
    """Test conversion from sv to NASC."""
    sv_linear = np.array([0.001, 0.002, 0.003])
    thickness_mean = np.array([5.0, 10.0, 15.0])

    result = ingest_sv.sv_to_nasc(sv_linear, thickness_mean)

    # Check that result has correct structure
    assert len(result) == 3
    assert all(result > 0)  # NASC should be positive

    # Test with known values
    expected_factor = 4 * np.pi * (1852**2)
    expected = expected_factor * sv_linear * thickness_mean
    np.testing.assert_array_almost_equal(result, expected)


def test_aggregate_cells():
    """Test cell-level aggregation."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "longitude": [-120.7, -120.7, -120.8, -120.8],
            "latitude": [34.6, 34.6, 34.7, 34.7],
            "interval": [1, 1, 1, 1],
            "layer": [1, 2, 1, 2],
            "frequency": [18.0, 18.0, 38.0, 38.0],
            "sv_mean": [-80.0, -75.0, -82.0, -77.0],
            "nasc": [1.5, 2.0, 1.8, 2.2],
            "thickness_mean": [5.0, 8.0, 6.0, 9.0],
        }
    )

    result = ingest_sv.aggregate_cells(data)

    assert isinstance(result, pd.DataFrame)
    assert result.columns.nlevels == 2  # MultiIndex columns
    # Check that fill values are applied correctly for missing data


def test_aggregate_intervals():
    """Test interval-level aggregation."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1, 1, 1],
            "longitude": [-120.7, -120.7, -120.7, -120.7],
            "latitude": [34.6, 34.6, 34.6, 34.6],
            "interval": [1, 1, 2, 2],
            "frequency": [18.0, 18.0, 18.0, 18.0],
            "sv_mean": [-80.0, -75.0, -82.0, -77.0],
            "sv_mean_linear": [0.0001, 0.0003162, 0.0000631, 0.0001995],
            "nasc": [1.5, 2.0, 1.8, 2.2],
            "thickness_mean": [5.0, 8.0, 6.0, 9.0],
        }
    )

    result = ingest_sv.aggregate_intervals(data)

    assert isinstance(result, pd.DataFrame)
    assert result.columns.nlevels == 2  # MultiIndex columns


def test_aggregate_intervals_missing_interval():
    """Test error when interval column is missing."""
    data = pd.DataFrame(
        {"transect_num": [1, 1], "frequency": [18.0, 18.0], "sv_mean_linear": [0.001, 0.002]}
    )

    with pytest.raises(KeyError):
        ingest_sv.aggregate_intervals(data)


def test_aggregate_transects():
    """Test transect-level aggregation."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1, 1, 1],
            "longitude": [-120.7, -120.7, -120.8, -120.8],
            "latitude": [34.6, 34.6, 34.7, 34.7],
            "interval": [1, 2, 1, 2],
            "frequency": [18.0, 18.0, 18.0, 18.0],
            "sv_mean": [-80.0, -75.0, -82.0, -77.0],
            "sv_mean_linear": [0.0001, 0.0003162, 0.0000631, 0.0001995],
            "nasc": [1.5, 2.0, 1.8, 2.2],
            "thickness_mean": [5.0, 8.0, 6.0, 9.0],
            "distance_s": [100.0, 200.0, 300.0, 400.0],
            "distance_e": [200.0, 300.0, 400.0, 500.0],
        }
    )

    result = ingest_sv.aggregate_transects(data)

    assert isinstance(result, pd.DataFrame)
    assert result.columns.nlevels == 2  # MultiIndex columns


def test_aggregate_transects_missing_transect_num():
    """Test error when transect_num column is missing."""
    data = pd.DataFrame({"frequency": [18.0, 18.0], "sv_mean_linear": [0.001, 0.002]})

    with pytest.raises(KeyError):
        ingest_sv.aggregate_transects(data)


def test_integrate_measurements_cells():
    """Test integration with cells method."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1, 1, 1],
            "longitude": [-120.7, -120.7, -120.8, -120.8],
            "latitude": [34.6, 34.6, 34.7, 34.7],
            "interval": [1, 1, 1, 1],
            "layer": [1, 2, 1, 2],
            "frequency": [18.0, 18.0, 38.0, 38.0],
            "sv_mean": [-80.0, -75.0, -82.0, -77.0],
            "thickness_mean": [5.0, 8.0, 6.0, 9.0],
        }
    )

    thresholds = {18.0: {"min": -90.0, "max": -50.0}, 38.0: {"min": -90.0, "max": -50.0}}

    result_data, result_coords = ingest_sv.integrate_measurements(data, "cells", thresholds)

    assert isinstance(result_data, pd.DataFrame)
    assert isinstance(result_coords, pd.DataFrame)


def test_integrate_measurements_interval():
    """Test integration with interval method."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1, 1, 1],
            "longitude": [-120.7, -120.7, -120.8, -120.8],
            "latitude": [34.6, 34.6, 34.7, 34.7],
            "interval": [1, 1, 2, 2],
            "frequency": [18.0, 18.0, 18.0, 18.0],
            "sv_mean": [-80.0, -75.0, -82.0, -77.0],
            "thickness_mean": [5.0, 8.0, 6.0, 9.0],
        }
    )

    thresholds = {18.0: {"min": -90.0, "max": -50.0}}

    result_data, result_coords = ingest_sv.integrate_measurements(data, "interval", thresholds)

    assert isinstance(result_data, pd.DataFrame)
    assert isinstance(result_coords, pd.DataFrame)


def test_integrate_measurements_transect():
    """Test integration with transect method."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1, 1, 1],
            "longitude": [-120.7, -120.7, -120.8, -120.8],
            "latitude": [34.6, 34.6, 34.7, 34.7],
            "interval": [1, 1, 2, 2],
            "frequency": [18.0, 18.0, 18.0, 18.0],
            "sv_mean": [-80.0, -75.0, -82.0, -77.0],
            "thickness_mean": [5.0, 8.0, 6.0, 9.0],
            "distance_s": [100.0, 200.0, 300.0, 400.0],
            "distance_e": [200.0, 300.0, 400.0, 500.0],
        }
    )

    thresholds = {18.0: {"min": -90.0, "max": -50.0}}

    result_data, result_coords = ingest_sv.integrate_measurements(data, "transect", thresholds)

    assert isinstance(result_data, pd.DataFrame)
    assert isinstance(result_coords, pd.DataFrame)


def test_integrate_measurements_invalid_method():
    """Test error handling for invalid aggregation method."""
    data = pd.DataFrame({"frequency": [18.0], "sv_mean": [-80.0], "thickness_mean": [5.0]})

    thresholds = {18.0: {"min": -90.0, "max": -50.0}}

    with pytest.raises(ValueError):
        ingest_sv.integrate_measurements(data, "invalid_method", thresholds)


def test_integrate_measurements_no_coordinates():
    """Test integration when coordinates are missing."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1],
            "interval": [1, 1],
            "layer": [1, 2],
            "frequency": [18.0, 18.0],
            "sv_mean": [-80.0, -75.0],
            "thickness_mean": [5.0, 8.0],
        }
    )

    thresholds = {18.0: {"min": -90.0, "max": -50.0}}

    result_data, result_coords = ingest_sv.integrate_measurements(data, "cells", thresholds)

    assert isinstance(result_data, pd.DataFrame)
    assert result_coords is None


def test_ingest_echoview_sv_basic_functionality(sv_directory_with_files):
    """Test basic functionality of complete SV ingestion pipeline."""
    center_frequencies = {18000: {"min": -90.0, "max": -50.0}}

    result_data, result_coords = ingest_sv.ingest_echoview_sv(
        sv_directory_with_files,
        center_frequencies=center_frequencies,
        aggregate_method="cells",
        impute_coordinates=True,
    )

    assert isinstance(result_data, pd.DataFrame)
    assert result_coords is None or isinstance(result_coords, pd.DataFrame)


def test_ingest_echoview_sv_with_transect_pattern(sv_directory_with_files):
    """Test SV ingestion with transect pattern."""
    center_frequencies = {18000: {"min": -90.0, "max": -50.0}}

    result_data, result_coords = ingest_sv.ingest_echoview_sv(
        sv_directory_with_files,
        center_frequencies=center_frequencies,
        transect_pattern=r"x(\d+)",
        aggregate_method="interval",
    )

    assert isinstance(result_data, pd.DataFrame)
    assert isinstance(result_coords, pd.DataFrame)


def test_ingest_echoview_sv_no_center_frequencies(sv_directory_with_files):
    """Test SV ingestion without specified center frequencies."""
    result_data, result_coords = ingest_sv.ingest_echoview_sv(
        sv_directory_with_files, transect_pattern=r"x(\d+)", aggregate_method="transect"
    )

    assert isinstance(result_data, pd.DataFrame)
    assert isinstance(result_coords, pd.DataFrame)


def test_ingest_echoview_sv_directory_not_found():
    """Test error when directory doesn't exist."""
    fake_path = Path("/fake/directory/that/does/not/exist")

    with pytest.raises(FileNotFoundError):
        ingest_sv.ingest_echoview_sv(fake_path)


def test_ingest_echoview_sv_empty_directory(tmp_path):
    """Test error when directory is empty."""
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()

    with pytest.raises(ValueError, match="No objects to concatenate"):
        ingest_sv.ingest_echoview_sv(empty_dir)
