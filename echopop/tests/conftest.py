import pytest
from pathlib import Path
from echopop import Survey
from _pytest.assertion.util import assertrepr_compare

def pytest_assertrepr_compare( config , op , left , right ):
    # hack the verbosity so we always show full diffs on assertion failures,
    # even if we're otherwise not fully verbose
    config.option.verbose = 2
    return assertrepr_compare( config , op , left , right) 

# Set up path to test_data folder
HERE = Path(__file__).parent.absolute()
TEST_DATA_ROOT = HERE.parent / "test_data"

@pytest.fixture(scope="session")
def test_path():
    return {
        "ROOT": TEST_DATA_ROOT,
        "CONFIG": TEST_DATA_ROOT / "config_files",
        "INPUT": TEST_DATA_ROOT / "input_files",  # this doesn't exist yet
    }


@pytest.fixture(scope="session")
def mock_survey(test_path) -> Survey:
    return Survey(
        init_config_path=Path(test_path["CONFIG"] / "config_init.yml"),
        survey_year_config_path=Path(test_path["CONFIG"] / "config_survey.yml"),
    )


# ============ below from previous version, remove after revamping is complete ============
@pytest.fixture(scope="session")
def config_base_path() -> Path:
    """
    Defines the base directory path for the
    configuration files.
    
    Returns
    -------
    pathlib.Path
        The base directory path for the configuration files
    """
    return HERE / "../config_files"


@pytest.fixture(scope="session")
def reports_base_path() -> Path:
    """
    Defines the base directory path were all reports
    generated should be saved.
    Returns
    -------
    pathlib.Path
        The base directory path for the reports
    """
    return HERE / "tests/reports/echopop_python_output"


@pytest.fixture(scope="session")
def matlab_output_base_path() -> Path:
    """
    Defines the base directory path for the
    Matlab output files.
    Returns
    -------
    pathlib.Path
        The base directory path for the Matlab output files
    """
    return Path("<YOUR-MATLAB-OUTPUTS-BASEPATH>")
