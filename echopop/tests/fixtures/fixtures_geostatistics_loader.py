import pandas as pd
import pytest


@pytest.fixture
def mesh_data():
    """Create sample mesh grid data for testing."""
    return pd.DataFrame(
        {
            "centroid_latitude": [45.1, 45.2, 45.3, 45.4],
            "centroid_longitude": [-125.1, -125.2, -125.3, -125.4],
            "fraction_cell_in_polygon": [1.0, 0.8, 0.95, 0.7],
        }
    )


@pytest.fixture
def mesh_excel_file(mesh_data, tmp_path):
    """Create a temporary Excel file with mesh grid data."""
    file_path = tmp_path / "test_mesh.xlsx"

    with pd.ExcelWriter(file_path) as writer:
        mesh_data.to_excel(writer, sheet_name="krigedgrid2_5nm_forChu", index=False)

    return file_path


@pytest.fixture
def mesh_column_map():
    """Create column mapping for mesh grid data."""
    return {
        "centroid_latitude": "latitude",
        "centroid_longitude": "longitude",
        "fraction_cell_in_polygon": "fraction",
    }


@pytest.fixture
def kriging_params_data():
    """Create sample kriging and variogram parameters data."""
    params = pd.DataFrame(columns=["Parameter", "Value"])
    params.loc[0] = ["krig.ratio", 1.5]
    params.loc[1] = ["krig.srad", 25.0]
    params.loc[2] = ["vario.lscl", 15.0]
    params.loc[3] = ["vario.nugt", 0.1]
    params.loc[4] = ["vario.powr", 1.5]
    params.loc[5] = ["vario.hole", 30.0]
    params.loc[6] = ["vario.res", 2.0]

    return params


@pytest.fixture
def kriging_params_excel_file(kriging_params_data, tmp_path):
    """Create a temporary Excel file with kriging parameters data."""
    file_path = tmp_path / "test_kriging_params.xlsx"

    # Transpose data for the correct format
    kriging_params_data.to_excel(file_path, sheet_name="Sheet1", index=False, header=False)

    return file_path


@pytest.fixture
def kriging_column_map():
    """Create column mapping for kriging parameters."""
    return {
        "hole": "hole_effect_range",
        "lscl": "correlation_range",
        "nugt": "nugget",
        "powr": "decay_power",
        "ratio": "anisotropy",
        "res": "lag_resolution",
        "srad": "search_radius",
    }
