import pandas as pd

from echopop.ingest import (
    load_kriging_variogram_params,
    load_mesh_data,
)


# Mesh data tests
def test_load_mesh_data_basic(mesh_excel_file):
    """Test basic loading of mesh data."""
    result = load_mesh_data(mesh_excel_file, "krigedgrid2_5nm_forChu")

    assert isinstance(result, pd.DataFrame)
    assert not result.empty

    for col in ["centroid_latitude", "centroid_longitude", "fraction_cell_in_polygon"]:
        assert col in result.columns


def test_load_mesh_data_with_column_map(mesh_excel_file, mesh_column_map):
    """Test loading mesh data with column mapping."""
    result = load_mesh_data(mesh_excel_file, "krigedgrid2_5nm_forChu", mesh_column_map)

    assert "latitude" in result.columns
    assert "longitude" in result.columns
    assert "fraction" in result.columns

    assert "centroid_latitude" not in result.columns
    assert "centroid_longitude" not in result.columns
    assert "fraction_cell_in_polygon" not in result.columns


# Kriging parameter tests
def test_load_kriging_variogram_params_basic(kriging_params_excel_file):
    """Test basic loading of kriging parameters."""
    kriging_params, variogram_params = load_kriging_variogram_params(
        kriging_params_excel_file, "Sheet1"
    )

    assert isinstance(kriging_params, dict)
    assert isinstance(variogram_params, dict)

    assert "ratio" in kriging_params
    assert "srad" in kriging_params

    assert "lscl" in variogram_params
    assert "nugt" in variogram_params
    assert "powr" in variogram_params
    assert "hole" in variogram_params
    assert "res" in variogram_params


def test_load_kriging_variogram_params_with_column_map(
    kriging_params_excel_file, kriging_column_map
):
    """Test loading with column mapping."""
    kriging_params, variogram_params = load_kriging_variogram_params(
        kriging_params_excel_file, "Sheet1", kriging_column_map
    )

    assert "anisotropy" in kriging_params
    assert "search_radius" in kriging_params

    assert "correlation_range" in variogram_params
    assert "nugget" in variogram_params
    assert "decay_power" in variogram_params
    assert "hole_effect_range" in variogram_params
    assert "lag_resolution" in variogram_params
