from pathlib import Path

import pandas as pd
import pytest

from echopop.ingest import (
    join_geostrata_by_latitude,
    join_strata_by_haul,
    load_geostrata,
    load_strata,
)


# Stratification tests
def test_load_strata_basic(strata_excel_file, strata_sheet_map):
    """Test basic loading of stratification data."""
    result = load_strata(strata_excel_file, strata_sheet_map)

    assert isinstance(result, dict)
    assert set(result.keys()) == set(strata_sheet_map.keys())

    for df in result.values():
        assert isinstance(df, pd.DataFrame)
        assert not df.empty


def test_load_strata_with_column_map(strata_excel_file, strata_sheet_map, strata_column_map):
    """Test loading with column name mapping."""
    result = load_strata(strata_excel_file, strata_sheet_map, strata_column_map)

    for df in result.values():
        assert "nasc_proportion" in df.columns
        assert "haul_num" in df.columns
        assert "stratum_num" in df.columns

        assert "fraction_hake" not in df.columns
        assert "haul" not in df.columns
        assert "stratum" not in df.columns


def test_load_strata_file_not_found(strata_sheet_map):
    """Test error handling for non-existent file."""
    non_existent_file = Path("non_existent_file.xlsx")

    with pytest.raises(FileNotFoundError):
        load_strata(non_existent_file, strata_sheet_map)


# Geographic stratification tests
def test_load_geostrata_basic(geostrata_excel_file, geostrata_sheet_map):
    """Test basic loading of geographic stratification data."""
    result = load_geostrata(geostrata_excel_file, geostrata_sheet_map)

    assert isinstance(result, dict)
    assert set(result.keys()) == set(geostrata_sheet_map.keys())

    for df in result.values():
        assert isinstance(df, pd.DataFrame)
        assert not df.empty


def test_load_geostrata_with_column_map(
    geostrata_excel_file, geostrata_sheet_map, geostrata_column_map
):
    """Test loading with column name mapping."""
    result = load_geostrata(geostrata_excel_file, geostrata_sheet_map, geostrata_column_map)

    for df in result.values():
        assert "northlimit_latitude" in df.columns
        assert "stratum_num" in df.columns

        assert "latitude (upper limit)" not in df.columns
        assert "stratum" not in df.columns


def test_load_geostrata_latitude_intervals(
    geostrata_excel_file, geostrata_sheet_map, geostrata_column_map
):
    """Test that latitude intervals are created."""
    result = load_geostrata(geostrata_excel_file, geostrata_sheet_map, geostrata_column_map)

    for df in result.values():
        assert "latitude_interval" in df.columns
        assert isinstance(df["latitude_interval"].dtype, pd.CategoricalDtype)


# Join stratification tests
def test_join_strata_by_haul_dataframe(biological_data, strata_data):
    """Test joining stratum data to a single DataFrame."""
    # Get length data and rename haul to haul_num to match the expected join column
    df = biological_data["length"].rename(columns={"haul": "haul_num"})
    strata_df = strata_data["inpfc"].rename(columns={"haul": "haul_num", "stratum": "stratum_num"})

    result = join_strata_by_haul(df, strata_df)

    assert "stratum_num" in result.columns

    for idx, row in result.iterrows():
        if row["haul_num"] in strata_df["haul_num"].values:
            assert not pd.isna(row["stratum_num"])


def test_join_strata_by_haul_dictionary(biological_data, strata_data):
    """Test joining stratum data to a dictionary of DataFrames."""
    strata_df = strata_data["inpfc"].rename(columns={"haul": "haul_num", "stratum": "stratum_num"})

    result = join_strata_by_haul(biological_data, strata_df)

    assert isinstance(result, dict)
    assert set(result.keys()) == set(biological_data.keys())

    for key, df in result.items():
        if "haul_num" in df.columns:
            assert "stratum_num" in df.columns


def test_join_geostrata_by_latitude_dataframe(mesh_data, geostrata_data):
    """Test joining geostrata to a single DataFrame."""
    geostrata_df = geostrata_data["inpfc"].rename(
        columns={"latitude (upper limit)": "northlimit_latitude", "stratum": "stratum_num"}
    )

    df = mesh_data.rename(columns={"centroid_latitude": "latitude"})

    result = join_geostrata_by_latitude(df, geostrata_df)

    assert "stratum_num" in result.columns
    assert not result["stratum_num"].isna().any()
