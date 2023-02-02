import os
import pandas as pd
import numpy as np
import EchoPro
from EchoPro.computation import SemiVariogram as SV
import pathlib


def test_generate_reports():

    config_base_path = pathlib.Path("../../../example_notebooks")
    reports_path = "../EchoPro/tests/reports/EchoPro_python_output"

    # TODO: formalize this test

    # change working directory so no initialization files need to be modified
    # TODO: this may not be necessary in the future
    os.chdir(config_base_path)

    # initialize Survey object
    survey_2019 = EchoPro.Survey(init_file_path='../config_files/initialization_config.yml',
                                 survey_year_file_path='../config_files/survey_year_2019_config.yml',
                                 source=3,
                                 exclude_age1=True)

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

    # initialize Reports object
    reports = survey_2019.get_reports()

    # generate and write reports
    reports.create_and_write_reports(output_path=reports_path)


def _compare_truth_produced_files(truth_base_path, produced_base_path, file_names_truth,
                                  file_names_produced, sheet_names_truth,
                                  sheet_names_produced, skiprows_truth, skiprows_produced,
                                  usecols, sortby = None):

    # TODO: document

    for file_ind in range(len(file_names_truth)):

        # obtain file_path pointing to the known data
        file_path_truth = truth_base_path / file_names_truth[file_ind]

        # obtain file_path pointing to the produced data
        file_path_produced = produced_base_path / file_names_produced[file_ind]

        for sheet_ind in range(len(sheet_names_truth[file_ind])):

            # gather known solution data produced by the Matlab version of EchoPro
            df_truth = pd.read_excel(file_path_truth, index_col=0,
                                     sheet_name=sheet_names_truth[file_ind][sheet_ind],
                                     skiprows=skiprows_truth[file_ind], usecols=usecols[file_ind])

            # gather produced solution
            df_produced = pd.read_excel(file_path_produced, index_col=0,
                                        sheet_name=sheet_names_produced[file_ind][sheet_ind],
                                        skiprows=skiprows_produced[file_ind], usecols=usecols[file_ind])

            if sortby:
                # sort produced and truth results, so they can be compared
                df_truth = df_truth.sort_values(by=sortby[file_ind][0])
                df_produced = df_produced.sort_values(by=sortby[file_ind][1])

            # compare known and produced values for abundance
            assert np.all(np.isclose(df_truth.to_numpy(), df_produced.to_numpy()))


def test_length_age_reports():

    # TODO: document

    matlab_output_base_path = pathlib.Path(
        "/Users/brandonreyes/UW_work/EchoPro_work/UW_EchoProMatlab_Repackaged/outputs/EchoPro_matlab_output_brandon_age_22_end_bin")

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

    produced_base_path = pathlib.Path("./EchoPro_python_output")

    _compare_truth_produced_files(matlab_output_base_path, produced_base_path, file_names_truth,
                                  file_names_produced, sheet_names_truth, sheet_names_produced,
                                  skiprows_truth, skiprows_produced, usecols)


# def test_biomass_ages_reports():
#     # TODO: document
#
#     matlab_output_base_path = pathlib.Path(
#         "/Users/brandonreyes/UW_work/EchoPro_work/UW_EchoProMatlab_Repackaged/outputs/EchoPro_matlab_output_brandon_age_22_end_bin")
#
#     # specify file names
#     file_names_truth = ["EchoPro_un-kriged_aged_output-2019_0.xlsx", "EchoPro_kriged_aged_output-2019_0.xlsx",
#                         "EchoPro_un-kriged_aged_output-2019_1.xlsx", "EchoPro_kriged_aged_output-2019_1.xlsx"]
#     file_names_produced = ["transect_based_aged_output_all.xlsx", "kriging_based_aged_output_all.xlsx",
#                            "transect_based_aged_output_non_zero.xlsx", "kriging_based_aged_output_non_zero.xlsx"]
#
#     # specify sheet names
#     sheet_names_truth = [["Sheet1", "Sheet2", "Sheet3"]] * 4
#     sheet_names_produced = [["all genders", "male", "female"]] * 4
#
#     # specify skiprows
#     skiprows_truth = [[0]] * 4
#     skiprows_produced = [None] * 4
#
#     # specify usecols
#     usecols = ["A:AA"] * 4
#
#     # specify sortby columns
#     sortby = [[["Lat", "Lon"], ["centroid_latitude", "centroid_longitude"]]] * 4
#
#     produced_base_path = pathlib.Path("./EchoPro_python_output")
#
#     _compare_truth_produced_files(matlab_output_base_path, produced_base_path, file_names_truth,
#                                   file_names_produced, sheet_names_truth, sheet_names_produced,
#                                   skiprows_truth, skiprows_produced, usecols, sortby)