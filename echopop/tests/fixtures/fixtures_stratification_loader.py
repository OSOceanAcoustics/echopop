import pandas as pd
import pytest


@pytest.fixture
def strata_data():
    """Create sample stratification data for testing."""
    inpfc_data = pd.DataFrame(
        {
            "haul": [1, 2, 3, 4, 5],
            "stratum": [10, 20, 30, 10, 20],
            "fraction_hake": [0.8, 0.7, 0.9, 0.75, 0.85],
        }
    )

    ks_data = pd.DataFrame(
        {
            "haul": [1, 2, 3, 4, 5],
            "stratum": [101, 102, 103, 101, 102],
            "fraction_hake": [0.75, 0.65, 0.85, 0.7, 0.8],
        }
    )

    return {"inpfc": inpfc_data, "ks": ks_data}


@pytest.fixture
def strata_excel_file(strata_data, tmp_path):
    """Create a temporary Excel file with stratification data."""
    file_path = tmp_path / "test_strata.xlsx"

    with pd.ExcelWriter(file_path) as writer:
        strata_data["inpfc"].to_excel(writer, sheet_name="INPFC", index=False)
        strata_data["ks"].to_excel(writer, sheet_name="Base KS", index=False)

    return file_path


@pytest.fixture
def strata_sheet_map():
    """Create sheet mapping for stratification data."""
    return {"inpfc": "INPFC", "ks": "Base KS"}


@pytest.fixture
def strata_column_map():
    """Create column mapping for stratification data."""
    return {"fraction_hake": "nasc_proportion", "haul": "haul_num", "stratum": "stratum_num"}


@pytest.fixture
def geostrata_data():
    """Create sample geographic stratification data for testing."""
    inpfc_data = pd.DataFrame(
        {"latitude (upper limit)": [40.0, 42.5, 45.0, 47.5, 50.0], "stratum": [10, 20, 30, 40, 50]}
    )

    ks_data = pd.DataFrame(
        {
            "latitude (upper limit)": [41.0, 43.5, 46.0, 48.5, 51.0],
            "stratum": [101, 102, 103, 104, 105],
        }
    )

    return {"inpfc": inpfc_data, "ks": ks_data}


@pytest.fixture
def geostrata_excel_file(geostrata_data, tmp_path):
    """Create a temporary Excel file with geographic stratification data."""
    file_path = tmp_path / "test_geostrata.xlsx"

    with pd.ExcelWriter(file_path) as writer:
        geostrata_data["inpfc"].to_excel(writer, sheet_name="INPFC", index=False)
        geostrata_data["ks"].to_excel(writer, sheet_name="stratification1", index=False)

    return file_path


@pytest.fixture
def geostrata_sheet_map():
    """Create sheet mapping for geographic stratification data."""
    return {"inpfc": "INPFC", "ks": "stratification1"}


@pytest.fixture
def geostrata_column_map():
    """Create column mapping for geographic stratification data."""
    return {"latitude (upper limit)": "northlimit_latitude", "stratum": "stratum_num"}
