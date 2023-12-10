import pathlib
from typing import List, Optional

import pandas as pd
import pytest
import numpy as np

import EchoPro
from EchoPro.computation import SemiVariogram as SV


@pytest.fixture(scope="module")
def generate_reports(config_base_path: pathlib.Path,
                     reports_base_path: pathlib.Path) -> None:
    """
    The purpose of this function is to generate all reports
    and write them into the ``reports_base_path`` directory.
    Once this is done, all tests constructed in this script
    can be run using the generated reports.

    Parameters
    ----------
    config_base_path: pathlib.Path
        The base directory path for the configuration files
    reports_base_path: pathlib.Path
        The base directory path for the reports
    """

    # initialize Survey object
    survey_2019 = EchoPro.Survey(
        init_file_path=config_base_path / 'initialization_config.yml',
        survey_year_file_path=config_base_path / 'survey_year_2019_config.yml',
        source=3,
        exclude_age1=True
    )

    # load all data
    survey_2019.load_survey_data()

    # compute all transect variables
    survey_2019.compute_transect_results()

    # initialize Kriging mesh object and apply appropriate coordinate transformations
    krig_mesh = survey_2019.get_kriging_mesh()
    krig_mesh.apply_coordinate_transformation(coord_type='transect')
    krig_mesh.apply_coordinate_transformation(coord_type='mesh')

    kriging_params = dict(
        # kriging parameters
        k_max=10,
        k_min=3,
        R=0.0226287,
        ratio=0.001,

        # parameters for semi-variogram model
        s_v_params={'nugget': 0.0, 'sill': 0.95279, 'ls': 0.0075429,
                    'exp_pow': 1.5, 'ls_hole_eff': 0.0},

        # grab appropriate semi-variogram model
        s_v_model=SV.generalized_exp_bessel
    )

    # initialize kriging routine
    krig = survey_2019.get_kriging(kriging_params)

    # run Kriging routine
    krig.run_biomass_kriging(krig_mesh)

    # compute additional variables at each Kriging mesh point
    krig.compute_kriging_variables()

    # compute Transect and Kriging based variables at each length and age bin
    survey_2019.compute_length_age_variables(data="all")

    # generate and write reports
    survey_2019.create_and_write_reports(output_path=reports_base_path)


def _compare_truth_produced_files(truth_base_path: pathlib.Path, produced_base_path: pathlib.Path,
                                  file_names_truth: List[str], file_names_produced: List[str],
                                  sheet_names_truth: List[List[str]], sheet_names_produced: List[List[str]],
                                  skiprows_truth: List[Optional[List[int]]],
                                  skiprows_produced: List[Optional[List[int]]],
                                  usecols: List[str], index_col_truth: List[Optional[int]],
                                  index_col_produced: List[Optional[int]],
                                  sortby: List[Optional[List[List[str]]]]) -> None:
    """
    This function reads in the true data (i.e. Matlab EchoPro generated files) and
    compares it against the produced (i.e. Python generated) data. This function
    is utilized by all tests within this script.

    Parameters
    ----------
    truth_base_path: pathlib.Path
        The base directory path for the Matlab output files
    produced_base_path: pathlib.Path
        The base directory path for the Python generated reports
    file_names_truth: list of str
        The names of the files for the true data
    file_names_produced: list of str
        The names of the files for the produced data
    sheet_names_truth: list of lists of str
        The sheet names for each true data file
    sheet_names_produced: list of lists of str
        The sheet names for each produced data file
    skiprows_truth: list of lists of int
        The rows to skip when reading in file for true data
    skiprows_produced: list of lists of int
        The rows to skip when reading in file for produced data
    usecols: list of str
        The columns that should be read in the data
    index_col_truth: list of int
        The column that should be used as the index for the true data
    index_col_produced: list of int
        The column that should be used as the index for the produced data
    sortby: list of lists of list of str
        The column(s) that should be used to sort both the true and produced
        data, where the first list element corresponds to the true data
        and the second corresponds to the produced data

    Notes
    -----
    If the true and produced data do not match, then an assert statement will
    be triggered.
    """

    for file_ind in range(len(file_names_truth)):
        # obtain file_path pointing to the known data
        file_path_truth = truth_base_path / file_names_truth[file_ind]

        # obtain file_path pointing to the produced data
        file_path_produced = produced_base_path / file_names_produced[file_ind]

        for sheet_ind in range(len(sheet_names_truth[file_ind])):
            # gather known solution data produced by the Matlab version of EchoPro
            df_truth = pd.read_excel(file_path_truth, index_col=index_col_truth[file_ind],
                                     sheet_name=sheet_names_truth[file_ind][sheet_ind],
                                     skiprows=skiprows_truth[file_ind], usecols=usecols[file_ind])

            # gather produced solution
            df_produced = pd.read_excel(file_path_produced, index_col=index_col_produced[file_ind],
                                        sheet_name=sheet_names_produced[file_ind][sheet_ind],
                                        skiprows=skiprows_produced[file_ind], usecols=usecols[file_ind])

            if sortby[file_ind]:
                # sort produced and truth results, so they can be compared
                df_truth = df_truth.sort_values(by=sortby[file_ind][0])
                df_produced = df_produced.sort_values(by=sortby[file_ind][1])

            # compare known and produced values for abundance
            assert np.all(np.isclose(df_truth.to_numpy(), df_produced.to_numpy()))


@pytest.mark.usefixtures("generate_reports")
def test_length_age_reports(config_base_path: pathlib.Path, matlab_output_base_path: pathlib.Path,
                            reports_base_path: pathlib.Path):
    """
    Ensures that the length age reports for Kriging and Transect based data
    produced match the Matlab generated output.

    Parameters
    ----------
    config_base_path: pathlib.Path
        The base directory path for the configuration files
    matlab_output_base_path: pathlib.Path
        The base directory path for the Matlab output files
    reports_base_path: pathlib.Path
        The base directory path for the reports
    """

    # specify file names
    file_names_truth = ["un-kriged_len_age_abundance_table.xlsx", "kriged_len_age_abundance_table.xlsx",
                        "un-kriged_len_age_biomass_table.xlsx", "kriged_len_age_biomass_table.xlsx"]
    file_names_produced = ["transect_based_len_age_abundance.xlsx", "kriging_based_len_age_abundance.xlsx",
                           "transect_based_len_age_biomass.xlsx", "kriging_based_len_age_biomass.xlsx"]

    # specify sheet names
    sheet_names_truth = [["Sheet3", "Sheet1", "Sheet2"]]*4
    sheet_names_produced = [["all genders", "male", "female"]]*4

    # specify skiprows
    skiprows_truth = [[0, 42, 43, 44, 45, 46]]*4
    skiprows_produced = [None] * 4

    # specify usecols
    usecols = ["A:X"]*2 + ["A:W"]*2

    # specify sortby columns
    sortby = [None] * 4

    # specify the index column
    index_col_truth = [0] * 4
    index_col_produced = [0] * 4

    _compare_truth_produced_files(matlab_output_base_path, reports_base_path, file_names_truth,
                                  file_names_produced, sheet_names_truth, sheet_names_produced,
                                  skiprows_truth, skiprows_produced, usecols, index_col_truth,
                                  index_col_produced, sortby)


@pytest.mark.usefixtures("generate_reports")
def test_biomass_ages_reports(config_base_path: pathlib.Path, matlab_output_base_path: pathlib.Path,
                              reports_base_path: pathlib.Path):
    """
    Ensures that the biomass at different age bins reports for Kriging and Transect based data
    produced match the Matlab generated output.

    Parameters
    ----------
    config_base_path: pathlib.Path
        The base directory path for the configuration files
    matlab_output_base_path: pathlib.Path
        The base directory path for the Matlab output files
    reports_base_path: pathlib.Path
        The base directory path for the reports
    """

    # specify file names
    file_names_truth = ["EchoPro_un-kriged_aged_output-2019_0.xlsx", "EchoPro_kriged_aged_output-2019_0.xlsx",
                        "EchoPro_un-kriged_aged_output-2019_1.xlsx", "EchoPro_kriged_aged_output-2019_1.xlsx"]
    file_names_produced = ["transect_based_aged_output_all.xlsx", "kriging_based_aged_output_all.xlsx",
                           "transect_based_aged_output_non_zero.xlsx", "kriging_based_aged_output_non_zero.xlsx"]

    # specify sheet names
    sheet_names_truth = [["Sheet1", "Sheet2", "Sheet3"]] * 4
    sheet_names_produced = [["all genders", "male", "female"]] * 4

    # specify skiprows
    skiprows_truth = [[0]] * 4
    skiprows_produced = [None] * 4

    # specify usecols
    usecols = ["A:AA"] * 4

    # specify sortby columns
    sortby = [[["Lat", "Lon"], ["latitude", "longitude"]],
              [["Lat", "Lon"], ["centroid_latitude", "centroid_longitude"]]] * 2

    # specify the index column
    index_col_truth = [0, None] * 2
    index_col_produced = [0, 0] * 2

    _compare_truth_produced_files(matlab_output_base_path, reports_base_path, file_names_truth,
                                  file_names_produced, sheet_names_truth, sheet_names_produced,
                                  skiprows_truth, skiprows_produced, usecols, index_col_truth,
                                  index_col_produced, sortby)


@pytest.mark.usefixtures("generate_reports")
def test_core_variables_reports(config_base_path: pathlib.Path, matlab_output_base_path: pathlib.Path,
                                reports_base_path: pathlib.Path):
    """
    Ensures that the core variable reports for Kriging and Transect based data
    produced match the Matlab generated output.

    Parameters
    ----------
    config_base_path: pathlib.Path
        The base directory path for the configuration files
    matlab_output_base_path: pathlib.Path
        The base directory path for the Matlab output files
    reports_base_path: pathlib.Path
        The base directory path for the reports
    """

    # specify file names
    file_names_truth = ["EchoPro_un-kriged_output-26-Jan-2023_0.xlsx", "EchoPro_kriged_output-26-Jan-2023_0.xlsx",
                        "EchoPro_un-kriged_output-26-Jan-2023_1.xlsx", "EchoPro_kriged_output-26-Jan-2023_1.xlsx"]
    file_names_produced = ["transect_based_core_output_all.xlsx", "kriging_based_core_output_all.xlsx",
                           "transect_based_core_output_non_zero.xlsx", "kriging_based_core_output_non_zero.xlsx"]

    # specify sheet names
    sheet_names_truth = [["Sheet1"]] * 4
    sheet_names_produced = [["Sheet1"]] * 4

    # specify skiprows
    skiprows_truth = [None] * 4
    skiprows_produced = [None] * 4

    # specify usecols
    usecols = ["A:L"] * 4

    # specify sortby columns
    sortby = [[["Lat", "Lon"], ["latitude", "longitude"]],
              [["Lat", "Lon"], ["centroid_latitude", "centroid_longitude"]]] * 2

    # specify the index column
    index_col_truth = [None] * 4
    index_col_produced = [None] * 4

    _compare_truth_produced_files(matlab_output_base_path, reports_base_path, file_names_truth,
                                  file_names_produced, sheet_names_truth, sheet_names_produced,
                                  skiprows_truth, skiprows_produced, usecols, index_col_truth,
                                  index_col_produced, sortby)


@pytest.mark.usefixtures("generate_reports")
def test_kriging_input_report(config_base_path: pathlib.Path, matlab_output_base_path: pathlib.Path,
                              reports_base_path: pathlib.Path):
    """
    Ensures that the Kriging input report produced matches the Matlab generated output.

    Parameters
    ----------
    config_base_path: pathlib.Path
        The base directory path for the configuration files
    matlab_output_base_path: pathlib.Path
        The base directory path for the Matlab output files
    reports_base_path: pathlib.Path
        The base directory path for the reports
    """

    # specify file names
    file_names_truth = ["kriging_input.xlsx"]
    file_names_produced = ["kriging_input.xlsx"]

    # specify sheet names
    sheet_names_truth = [["Sheet1"]]
    sheet_names_produced = [["Sheet1"]]

    # specify skiprows
    skiprows_truth = [None]
    skiprows_produced = [None]

    # specify usecols
    usecols = ["A:E"]

    # specify sortby columns
    sortby = [[["Lat", "Lon"], ["latitude", "longitude"]]]

    # specify the index column
    index_col_truth = [None]
    index_col_produced = [None]

    _compare_truth_produced_files(matlab_output_base_path, reports_base_path, file_names_truth,
                                  file_names_produced, sheet_names_truth, sheet_names_produced,
                                  skiprows_truth, skiprows_produced, usecols, index_col_truth,
                                  index_col_produced, sortby)


@pytest.mark.usefixtures("generate_reports")
def test_len_haul_count_reports(config_base_path: pathlib.Path, matlab_output_base_path: pathlib.Path,
                                reports_base_path: pathlib.Path):
    """
    Ensures that the length count at each haul reports
    produced match the Matlab generated output.

    Parameters
    ----------
    config_base_path: pathlib.Path
        The base directory path for the configuration files
    matlab_output_base_path: pathlib.Path
        The base directory path for the Matlab output files
    reports_base_path: pathlib.Path
        The base directory path for the reports
    """

    # specify file names
    file_names_truth = ["aged_len_haul_counts_table.xlsx", "total_len_haul_counts_table.xlsx"]
    file_names_produced = ["specimen_length_counts_haul.xlsx", "total_length_counts_haul.xlsx"]

    # specify sheet names
    sheet_names_truth = [["Sheet3", "Sheet1", "Sheet2"]] * 2
    sheet_names_produced = [["all genders", "male", "female"]] * 2

    # specify skiprows
    skiprows_truth = [[0, 42, 43, 44, 45, 46]]*2
    skiprows_produced = [[41]] * 2

    # specify usecols
    usecols = ["A:BX"] * 2

    # specify sortby columns
    sortby = [None] * 2

    # specify the index column
    index_col_truth = [None] * 2
    index_col_produced = [None] * 2

    _compare_truth_produced_files(matlab_output_base_path, reports_base_path, file_names_truth,
                                  file_names_produced, sheet_names_truth, sheet_names_produced,
                                  skiprows_truth, skiprows_produced, usecols, index_col_truth,
                                  index_col_produced, sortby)
