import copy
from pathlib import Path

import numpy as np
import yaml

from echopop import Survey
from echopop.core import LAYER_NAME_MAP
from echopop.tests.conftest import assert_dictionary_structure_equal
from echopop.utils.data_file_validation import load_configuration


def test_load_configuration(test_path, tmp_path):
    init_params = yaml.safe_load(Path(test_path["CONFIG"] / "config_init.yml").read_text())
    survey_params = yaml.safe_load(Path(test_path["CONFIG"] / "config_survey.yml").read_text())

    # Swap out test data root path
    survey_params["data_root_dir"] = str(test_path["INPUT"])

    # Write a new temp yaml file with correct data path
    temp_config_survey_path = tmp_path / "config_survey_local.yaml"
    with open(temp_config_survey_path, "w") as yf:
        yaml.safe_dump(survey_params, yf)

    # Use class method
    config = load_configuration(
        init_config_path=Path(test_path["CONFIG"] / "config_init.yml"),
        survey_year_config_path=temp_config_survey_path,
    )

    # Check parsed values (to be completed!)
    assert (
        config["stratified_survey_mean_parameters"]["strata_transect_proportion"]
        == init_params["stratified_survey_mean_parameters"]["strata_transect_proportion"]
    )


def test_init(mock_survey):
    objS = mock_survey
    assert isinstance(objS, Survey)


def test_load_survey_data(mock_survey, test_path):

    # Pull in configuration values
    mock_survey.config = load_configuration(
        Path(test_path["CONFIG"] / "config_init.yml"),
        Path(test_path["CONFIG"] / "config_survey.yml"),
    )

    # Initialize data attributes
    mock_survey.acoustics = copy.deepcopy(LAYER_NAME_MAP["NASC"]["data_tree"])
    mock_survey.biology = copy.deepcopy(LAYER_NAME_MAP["biological"]["data_tree"])
    mock_survey.spatial = copy.deepcopy(LAYER_NAME_MAP["stratification"]["data_tree"])
    mock_survey.statistics = copy.deepcopy(LAYER_NAME_MAP["kriging"]["data_tree"])

    # Load in data using the `load_survey_data` method
    mock_survey.load_survey_data()

    # -----------------
    # Evaluate results
    # -----------------
    # Dictionary structure
    # !!! TODO: based on the original data structure -- will need to be updated once the core data
    # structure is also updated
    # ---- Check attributes
    assert set(["acoustics", "biology", "spatial", "statistics"]) <= set(dir(mock_survey))
    # ---- Check sub-directory keys
    assert_dictionary_structure_equal(mock_survey.acoustics, LAYER_NAME_MAP["NASC"]["data_tree"])
    assert_dictionary_structure_equal(
        mock_survey.biology, LAYER_NAME_MAP["biological"]["data_tree"]
    )
    assert_dictionary_structure_equal(
        mock_survey.spatial, LAYER_NAME_MAP["stratification"]["data_tree"]
    )
    assert_dictionary_structure_equal(
        mock_survey.statistics, LAYER_NAME_MAP["kriging"]["data_tree"]
    )
    # Data structure
    # ++++ acoustics
    assert mock_survey.acoustics["nasc"]["nasc_df"].shape == tuple([1, 10])
    # ++++ biology
    assert mock_survey.biology["catch_df"].shape == tuple([2, 7])
    assert mock_survey.biology["distributions"]["age_bins_arr"].shape == tuple(
        [
            0,
        ]
    )
    assert mock_survey.biology["distributions"]["length_bins_arr"].shape == tuple(
        [
            0,
        ]
    )
    assert mock_survey.biology["haul_to_transect_df"].shape == tuple([2, 5])
    assert mock_survey.biology["length_df"].shape == tuple([2, 10])
    assert mock_survey.biology["specimen_df"].shape == tuple([2, 11])
    # ++++ spatial
    assert mock_survey.spatial["strata_df"].shape == tuple([1, 3])
    assert mock_survey.spatial["geo_strata_df"].shape == tuple([1, 2])
    assert mock_survey.spatial["inpfc_strata_df"].shape == tuple([1, 2])
    # ++++ statistics
    assert mock_survey.statistics["kriging"]["mesh_df"].shape == tuple([19843, 3])
    assert mock_survey.statistics["kriging"]["isobath_200m_df"].shape == tuple([147, 2])
    assert len(mock_survey.statistics["kriging"]["model_config"]) == 39
    assert len(mock_survey.statistics["variogram"]["model_config"]) == 13
    # Test merged outputs
    assert set(mock_survey.biology["haul_to_transect_df"].columns) <= set(
        mock_survey.biology["catch_df"].columns
    )
    assert set(mock_survey.biology["haul_to_transect_df"].columns) <= set(
        mock_survey.biology["length_df"].columns
    )
    assert set(mock_survey.biology["haul_to_transect_df"].columns) <= set(
        mock_survey.biology["specimen_df"].columns
    )
    # Test biological data (sex definition)
    assert np.all(
        (mock_survey.biology["length_df"].sex == "female")
        & (mock_survey.biology["length_df"].group == "sexed")
    )
    assert np.all(
        (mock_survey.biology["specimen_df"].sex == ["male", "female"])
        & (mock_survey.biology["specimen_df"].group == "sexed")
    )


def test_biometric_distributions(mock_survey, test_path):

    # Pull in configuration values
    mock_survey.config = load_configuration(
        Path(test_path["CONFIG"] / "config_init.yml"),
        Path(test_path["CONFIG"] / "config_survey.yml"),
    )

    # Initialize data attributes
    mock_survey.acoustics = copy.deepcopy(LAYER_NAME_MAP["NASC"]["data_tree"])
    mock_survey.biology = copy.deepcopy(LAYER_NAME_MAP["biological"]["data_tree"])
    mock_survey.spatial = copy.deepcopy(LAYER_NAME_MAP["stratification"]["data_tree"])
    mock_survey.statistics = copy.deepcopy(LAYER_NAME_MAP["kriging"]["data_tree"])

    # Load in data using the `load_survey_data` method
    mock_survey.load_survey_data()

    # Generate length and age distributions
    mock_survey.biometric_distributions()

    # -----------------
    # Evaluate results
    # -----------------
    # Data structure
    assert mock_survey.biology["distributions"]["age"]["age_interval_arr"].shape == tuple(
        [
            23,
        ]
    )
    assert mock_survey.biology["distributions"]["age"]["age_bins_arr"].shape == tuple(
        [
            22,
        ]
    )
    assert mock_survey.biology["distributions"]["length"]["length_interval_arr"].shape == tuple(
        [
            41,
        ]
    )
    assert mock_survey.biology["distributions"]["length"]["length_bins_arr"].shape == tuple(
        [
            40,
        ]
    )
    # Data equality
    assert np.all(
        mock_survey.biology["distributions"]["age"]["age_interval_arr"]
        == np.linspace(0.5, 22.5, 23)
    )
    assert np.all(
        mock_survey.biology["distributions"]["age"]["age_bins_arr"] == np.linspace(1, 22, 22)
    )
    assert np.all(
        mock_survey.biology["distributions"]["length"]["length_interval_arr"]
        == np.linspace(1, 81, 41)
    )
    assert np.all(
        mock_survey.biology["distributions"]["length"]["length_bins_arr"] == np.linspace(2, 80, 40)
    )
