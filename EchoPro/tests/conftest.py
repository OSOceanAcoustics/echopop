import pytest
from pathlib import Path


# Set up path to test_data folder
HERE = Path(__file__).parent.absolute()
TEST_DATA_ROOT = HERE.parent / "test_data"

@pytest.fixture(scope="session")
def test_path():
    return {
        "ROOT": TEST_DATA_ROOT,
        "CONFIG": TEST_DATA_ROOT / "config_files",
        "INPUT": TEST_DATA_ROOT / "input_files",
    }


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
    return HERE / "tests/reports/EchoPro_python_output"


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
