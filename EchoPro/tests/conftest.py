import pytest
from pathlib import Path


# Set up path to test_data folder
HERE = Path(__file__).parent.absolute()
TEST_DATA_FOLDER = HERE / "../test_data"


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
    return pathlib.Path("<YOUR-MATLAB-OUTPUTS-BASEPATH>")
