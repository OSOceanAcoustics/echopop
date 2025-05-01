import pytest
import tempfile
import os
from echopop.ingest import read_csv_file

# import copy
# from pathlib import Path

# import numpy as np
# import pytest
# import yaml

# from .. import Survey
# from ..core_tmp import LAYER_NAME_MAP
# from ..utils.load import load_configuration
# from .conftest import assert_dictionary_structure_equal

# ==================================================================================================
# FIXTURES
# --------
@pytest.fixture
def sample_csv_content():
    """
    Provide sample CSV content with mixed-case column names for testing.
    
    Returns
    -------
        str: Sample CSV data with 3 columns and 2 rows
    """
    return "Column1,COLUMN2,CoLuMn3\n1,2,3\n4,5,6"

@pytest.fixture
def sample_csv_file(sample_csv_content):
    """
    Create a temporary CSV file for testing.
    
    Parameters
    ----------
        sample_csv_content: Content to write to the CSV file
        
    Yields
    ------
        str: Path to the temporary CSV file
        
    Notes
    -----
        File is automatically deleted after the test
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(sample_csv_content)  # Now properly indented under 'with'
        filename = f.name
    yield filename
    os.unlink(filename)

# ==================================================================================================
# READ_CSV_FILE
# -------------
def test_read_csv_file(sample_csv_file):
    """
    Test that read_csv_file correctly reads a CSV and converts column names to lowercase.
    
    This test verifies:
    1. All column names are converted to lowercase
    2. The expected column names are present
    3. Data is read correctly with proper dimensions and values
    
    Parameters
    ----------
        sample_csv_file: Path to test CSV file
    """
    df = read_csv_file(sample_csv_file)
    
    # Check that all column names are lowercase
    assert all(col == col.lower() for col in df.columns), "Column names should be lowercase"
    assert list(df.columns) == ['column1', 'column2', 'column3'], "Incorrect column names"
    
    # Check data was read correctly
    assert df.shape == (2, 3), "DataFrame should have 2 rows and 3 columns"
    assert df.iloc[0, 0] == 1, "First value should be 1"
    assert df.iloc[1, 2] == 6, "Last value should be 6"

def test_read_csv_file_nonexistent():
    """
    Test that read_csv_file raises an appropriate error for nonexistent files.
    
    This test verifies that the function properly handles the case of a missing file
    by raising a FileNotFoundError.
    """
    with pytest.raises(FileNotFoundError):
        read_csv_file("nonexistent_file.csv")

# def test_load_configuration(test_path, tmp_path):

#     # Read in the initialization and file configuration
#     # ---- Initialization
#     init_config = yaml.safe_load(Path(test_path["CONFIG"] / "config_init.yml").read_text())
#     # ---- Files
#     files_config = yaml.safe_load(Path(test_path["CONFIG"] / "config_survey.yml").read_text())

#     # Swap out test data root path
#     files_config["data_root_dir"] = str(test_path["INPUT"])

#     # Write a new temp yaml file with correct data path
#     temp_config_survey_path = tmp_path / "config_survey_local.yaml"
#     with open(temp_config_survey_path, "w") as yf:
#         yaml.safe_dump(files_config, yf)

#     # Use class method
#     config = load_configuration(
#         init_config_path=Path(test_path["CONFIG"] / "config_init.yml"),
#         survey_year_config_path=temp_config_survey_path,
#     )

#     # Check parsed values (to be completed!)
#     assert (
#         config["stratified_survey_mean_parameters"]["strata_transect_proportion"]
#         == init_config["stratified_survey_mean_parameters"]["strata_transect_proportion"]
#     )


# def test_init(mock_survey):
#     objS = mock_survey
#     assert isinstance(objS, Survey)


# @pytest.mark.skip(reason="Function has since been updated!")
# def test_load_survey_data(mock_survey, test_path):

#     # Pull in configuration values
#     mock_survey.config = load_configuration(
#         Path(test_path["CONFIG"] / "config_init.yml"),
#         Path(test_path["CONFIG"] / "config_survey.yml"),
#     )

#     # Initialize data attributes
#     mock_survey.acoustics = copy.deepcopy(LAYER_NAME_MAP["NASC"]["data_tree"])
#     mock_survey.biology = copy.deepcopy(LAYER_NAME_MAP["biological"]["data_tree"])
#     mock_survey.spatial = copy.deepcopy(LAYER_NAME_MAP["stratification"]["data_tree"])
#     mock_survey.statistics = copy.deepcopy(LAYER_NAME_MAP["kriging"]["data_tree"])

#     # Load in data using the `load_survey_data` method
#     mock_survey.load_survey_data()

#     # -----------------
#     # Evaluate results
#     # Evaluate results
#     # -----------------
#     # Dictionary structure
#     # !!! TODO: based on the original data structure -- will need to be updated once the core data
#     # structure is also updated
#     # Dictionary structure
#     # !!! TODO: based on the original data structure -- will need to be updated once the core data
#     # structure is also updated
#     # ---- Check attributes
#     assert set(["acoustics", "biology", "spatial", "statistics"]) <= set(dir(mock_survey))
#     # ---- Check sub-directory keys
#     assert_dictionary_structure_equal(mock_survey.acoustics, LAYER_NAME_MAP["NASC"]["data_tree"])
#     assert_dictionary_structure_equal(
#         mock_survey.biology, LAYER_NAME_MAP["biological"]["data_tree"]
#     )
#     assert_dictionary_structure_equal(
#         mock_survey.spatial, LAYER_NAME_MAP["stratification"]["data_tree"]
#     )
#     assert_dictionary_structure_equal(
#         mock_survey.statistics, LAYER_NAME_MAP["kriging"]["data_tree"]
#     )
#     # ++++ acoustics
#     assert mock_survey.acoustics["nasc"]["nasc_df"].shape == tuple([1, 10])
#     # ++++ biology
#     assert mock_survey.biology["catch_df"].shape == tuple([2, 7])
#     assert mock_survey.biology["distributions"]["age_bins_arr"].shape == tuple(
#         [
#             0,
#         ]
#     )
#     assert mock_survey.biology["distributions"]["length_bins_arr"].shape == tuple(
#         [
#             0,
#         ]
#     )
#     assert mock_survey.biology["haul_to_transect_df"].shape == tuple([2, 5])
#     assert mock_survey.biology["length_df"].shape == tuple([2, 10])
#     assert mock_survey.biology["specimen_df"].shape == tuple([2, 11])
#     # ++++ spatial
#     assert mock_survey.spatial["strata_df"].shape == tuple([1, 3])
#     assert mock_survey.spatial["geo_strata_df"].shape == tuple([1, 2])
#     assert mock_survey.spatial["inpfc_strata_df"].shape == tuple([1, 2])
#     # ++++ statistics
#     assert mock_survey.statistics["kriging"]["mesh_df"].shape == tuple([19843, 3])
#     assert mock_survey.statistics["kriging"]["isobath_200m_df"].shape == tuple([147, 2])
#     assert len(mock_survey.statistics["kriging"]["model_config"]) == 39
#     assert len(mock_survey.statistics["variogram"]["model_config"]) == 13
#     # Test merged outputs
#     assert set(mock_survey.biology["haul_to_transect_df"].columns) <= set(
#         mock_survey.biology["catch_df"].columns
#     )
#     assert set(mock_survey.biology["haul_to_transect_df"].columns) <= set(
#         mock_survey.biology["length_df"].columns
#     )
#     assert set(mock_survey.biology["haul_to_transect_df"].columns) <= set(
#         mock_survey.biology["specimen_df"].columns
#     )
#     # Test biological data (sex definition)
#     assert np.all(
#         (mock_survey.biology["length_df"].sex == "female")
#         & (mock_survey.biology["length_df"].group == "sexed")
#     )
#     assert np.all(
#         (mock_survey.biology["specimen_df"].sex == ["male", "female"])
#         & (mock_survey.biology["specimen_df"].group == "sexed")
#     )
