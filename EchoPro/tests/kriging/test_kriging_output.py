import pytest

import os
import pandas as pd
import numpy as np
import EchoPro
from EchoPro.computation import SemiVariogram as SV


def test_biomass_age_output(config_base_path, matlab_output_base_path):

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

    # set variables to improve readability
    krig_results = survey_2019.bio_calc.kriging_results_gdf
    krig_results_male = survey_2019.bio_calc.kriging_results_male_gdf
    krig_results_female = survey_2019.bio_calc.kriging_results_female_gdf

    # obtain file_path pointing to the known output and associated sheet names
    file_path = matlab_output_base_path / "EchoPro_kriged_aged_output-2019_0.xlsx"
    sheet_name = "Sheet1"
    sheet_name_male = "Sheet2"
    sheet_name_female = "Sheet3"

    # gather known solution data produced by the Matlab version of EchoPro
    df_known = pd.read_excel(file_path, sheet_name=sheet_name, skiprows=1, usecols="A:Z").drop(columns=["stratum"])
    male_df_known = pd.read_excel(file_path, sheet_name=sheet_name_male, skiprows=1, usecols="A:Z").drop(
        columns=["stratum"])
    female_df_known = pd.read_excel(file_path, sheet_name=sheet_name_female, skiprows=1, usecols="A:Z").drop(
        columns=["stratum"])

    # sort known solution dfs by latitude and longitude
    df_known.sort_values(by=["Lat", "Lon"], inplace=True)
    male_df_known.sort_values(by=["Lat", "Lon"], inplace=True)
    female_df_known.sort_values(by=["Lat", "Lon"], inplace=True)

    # define columns to compare in the produced solution
    produced_wanted_columns = ["centroid_latitude", "centroid_longitude", "biomass_adult"] + ["biomass_age_bin_" + str(i + 1) for i in range(22)]

    df_produced = krig_results.sort_values(by=["centroid_latitude",
                                               "centroid_longitude"])[produced_wanted_columns]
    male_df_produced = krig_results_male.sort_values(by=["centroid_latitude",
                                                         "centroid_longitude"])[produced_wanted_columns]
    female_df_produced = krig_results_female.sort_values(by=["centroid_latitude",
                                                             "centroid_longitude"])[produced_wanted_columns]

    # compare known and produced values
    assert np.all(np.isclose(df_known.to_numpy(), df_produced.to_numpy()))

    assert np.all(np.isclose(male_df_known.to_numpy(), male_df_produced.to_numpy()))

    assert np.all(np.isclose(female_df_known.to_numpy(), female_df_produced.to_numpy()))


def test_core_output(config_base_path, matlab_output_base_path):

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

    # set variables to improve readability
    krig_results = survey_2019.bio_calc.kriging_results_gdf
    krig_results_male = survey_2019.bio_calc.kriging_results_male_gdf
    krig_results_female = survey_2019.bio_calc.kriging_results_female_gdf

    # obtain file_path pointing to the known output and associated sheet name
    file_path = matlab_output_base_path / "EchoPro_kriged_output-26-Jan-2023_0.xlsx"
    sheet_name = "Sheet1"

    # gather known solution data produced by the Matlab version of EchoPro
    df_known = pd.read_excel(file_path, sheet_name=sheet_name, usecols="A:L")

    # sort known solution dfs by latitude and longitude
    df_known.sort_values(by=["Lat", "Lon"], inplace=True)

    # define columns grab from the produced results
    wanted_columns = ["centroid_latitude", "centroid_longitude", "stratum_num", "NASC", "biomass_adult",
                      "abundance_adult", "sig_b", "biomass_adult_cell_CV"]
    gender_wanted_columns = ["centroid_latitude", "centroid_longitude", "biomass_adult", "abundance_adult"]

    # sort produced results and format them, so they can be compared to known results
    df_produced = krig_results.sort_values(by=["centroid_latitude",
                                               "centroid_longitude"])[wanted_columns]
    male_df_produced = krig_results_male.sort_values(by=["centroid_latitude",
                                                         "centroid_longitude"])[gender_wanted_columns].drop(
        columns=["centroid_latitude", "centroid_longitude"]).rename(
        columns={"biomass_adult": "biomass_male_adult", "abundance_adult": "abundance_male_adult"})
    female_df_produced = krig_results_female.sort_values(by=["centroid_latitude",
                                                             "centroid_longitude"])[gender_wanted_columns].drop(
        columns=["centroid_latitude", "centroid_longitude"]).rename(
        columns={"biomass_adult": "biomass_female_adult", "abundance_adult": "abundance_female_adult"})

    # put together results to match the known solution DataFrame
    final_df_produced = pd.concat([df_produced, male_df_produced, female_df_produced], axis=1)

    # produced column names in the same order as the known solution
    ordered_columns = ["centroid_latitude", "centroid_longitude", "stratum_num", "NASC", "abundance_male_adult",
                       "abundance_female_adult", "abundance_adult", "biomass_male_adult", "biomass_female_adult",
                       "biomass_adult", "sig_b", "biomass_adult_cell_CV"]

    # compare known and produced values
    assert np.all(np.isclose(df_known.to_numpy(), final_df_produced[ordered_columns].to_numpy()))
