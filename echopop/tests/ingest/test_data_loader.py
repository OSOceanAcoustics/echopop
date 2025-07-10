import pandas as pd
import pytest

from echopop.ingest import read_csv_file, read_xlsx_file


# ==================================================================================================
# Test *.csv reader
# -----------------
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
    assert list(df.columns) == ["column1", "column2", "column3"], "Incorrect column names"

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


# ==================================================================================================
# Test *.xlsx reader
# ------------------
def test_read_xlsx_file_basic(sample_excel_file):
    """Test basic functionality of reading an Excel file."""
    # Call the function
    result = read_xlsx_file(sample_excel_file, "Sheet1")

    # Check that result is a DataFrame
    assert isinstance(result, pd.DataFrame)

    # Check data content
    assert len(result) == 3
    assert result["column1"].tolist() == [1, 2, 3]
    assert result["column2"].tolist() == ["a", "b", "c"]
    assert result["mixedcase_column"].tolist() == [1.1, 2.2, 3.3]


def test_read_xlsx_file_column_names(sample_excel_file):
    """Test that column names are converted to lowercase."""
    # Call the function
    result = read_xlsx_file(sample_excel_file, "Sheet1")

    # Check that all column names are lowercase
    assert all(col == col.lower() for col in result.columns)

    # Check specific column names
    assert "column1" in result.columns
    assert "column2" in result.columns
    assert "mixedcase_column" in result.columns


def test_read_xlsx_file_different_sheet(sample_excel_file):
    """Test reading a different sheet from the Excel file."""
    # Call the function with Sheet2
    result = read_xlsx_file(sample_excel_file, "Sheet2")

    # Check data content from Sheet2
    assert len(result) == 3
    assert "sheet2_col1" in result.columns
    assert "sheet2_col2" in result.columns
    assert result["sheet2_col1"].tolist() == [4, 5, 6]
    assert result["sheet2_col2"].tolist() == ["d", "e", "f"]


def test_read_xlsx_file_nonexistent_file():
    """Test behavior when file doesn't exist."""
    nonexistent_file = "nonexistent_file.xlsx"

    # Check that attempting to read a non-existent file raises FileNotFoundError
    with pytest.raises(FileNotFoundError):
        read_xlsx_file(nonexistent_file, "Sheet1")


def test_read_xlsx_file_nonexistent_sheet(sample_excel_file):
    """Test behavior when sheet doesn't exist."""
    # Check that attempting to read a non-existent sheet raises ValueError
    with pytest.raises(ValueError):
        read_xlsx_file(sample_excel_file, "NonexistentSheet")


# ==================================================================================================
# Test load_isobath_data
# ----------------------
def test_load_isobath_data_basic(sample_isobath_file):
    """Test basic functionality of loading isobath data."""
    from echopop.nwfsc_feat.load_data import load_isobath_data

    df = load_isobath_data(sample_isobath_file, "isobath_data")

    # Check that all column names are lowercase
    assert all(col == col.lower() for col in df.columns), "Column names should be lowercase"
    assert list(df.columns) == [
        "longitude",
        "latitude",
        "depth_200m",
        "other_data",
    ], "Incorrect column names"

    # Check data was read correctly
    assert df.shape == (4, 4), "DataFrame should have 4 rows and 4 columns"
    assert df.iloc[0, 0] == -124.5, "First longitude should be -124.5"
    assert df.iloc[0, 1] == 46.2, "First latitude should be 46.2"


def test_load_isobath_data_with_column_mapping(sample_isobath_file):
    """Test load_isobath_data with column name mapping."""
    from echopop.nwfsc_feat.load_data import load_isobath_data

    column_map = {"depth_200m": "depth", "other_data": "value"}
    df = load_isobath_data(sample_isobath_file, "isobath_data", column_name_map=column_map)

    # Check that columns were renamed correctly
    assert "depth" in df.columns, "Column should be renamed to 'depth'"
    assert "value" in df.columns, "Column should be renamed to 'value'"
    assert "depth_200m" not in df.columns, "Original column name should be replaced"
    assert "other_data" not in df.columns, "Original column name should be replaced"


def test_load_isobath_data_nonexistent_file():
    """Test that load_isobath_data raises FileNotFoundError for nonexistent files."""
    from echopop.nwfsc_feat.load_data import load_isobath_data

    with pytest.raises(FileNotFoundError):
        load_isobath_data("nonexistent_file.xlsx", "Sheet1")


def test_load_isobath_data_nonexistent_sheet(sample_isobath_file):
    """Test that load_isobath_data raises error for nonexistent sheet."""
    from echopop.nwfsc_feat.load_data import load_isobath_data

    with pytest.raises(ValueError):
        load_isobath_data(sample_isobath_file, "NonexistentSheet")


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
