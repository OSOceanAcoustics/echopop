import yaml
from pathlib import Path
from EchoPro.utils.data_file_validation import load_configuration , validate_data_columns

def test_load_configuration(test_path, tmp_path):
    init_params = yaml.safe_load(
        Path(test_path["CONFIG"] / "config_init.yml").read_text()
    )
    survey_params = yaml.safe_load(
        Path(test_path["CONFIG"] / "config_survey.yml").read_text()
    )

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
