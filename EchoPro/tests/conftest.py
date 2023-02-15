"""``pytest`` configuration."""

import pytest
import pathlib


@pytest.fixture(scope="session")
def config_base_path() -> pathlib.Path:
    """
    Defines the base directory path for the
    configuration files.
    Returns
    -------
    pathlib.Path
        The base directory path for the configuration files
    """
    return pathlib.Path("../../../example_notebooks")


@pytest.fixture(scope="session")
def reports_base_path() -> pathlib.Path:
    """
    Defines the base directory path were all reports
    generated should be saved.
    Returns
    -------
    pathlib.Path
        The base directory path for the reports
    """
    return pathlib.Path("../EchoPro/tests/reports/EchoPro_python_output")


@pytest.fixture(scope="session")
def matlab_output_base_path() -> pathlib.Path:
    """
    Defines the base directory path for the
    Matlab output files.
    Returns
    -------
    pathlib.Path
        The base directory path for the Matlab output files
    """
    return pathlib.Path("/Users/brandonreyes/UW_work/EchoPro_work/UW_EchoProMatlab_Repackaged/outputs/EchoPro_matlab_output_brandon_age_22_end_bin")
