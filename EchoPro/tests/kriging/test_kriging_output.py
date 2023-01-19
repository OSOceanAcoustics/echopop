import pytest

import os
import pandas as pd
import numpy as np
import EchoPro
from EchoPro.computation import SemiVariogram as SV


def test_biomass_age_output():

    # change working directory so no initialization files need to be modified
    os.chdir("../../../example_notebooks")

    # initialize Survey object
    survey_2019 = EchoPro.Survey(init_file_path='../config_files/initialization_config.yml',
                                 survey_year_file_path='../config_files/survey_year_2019_config.yml',
                                 source=3,
                                 exclude_age1=True)

    # load all data
    survey_2019.load_survey_data()

    # compute all transect variables
    survey_2019.compute_biomass_density()

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

    # TODO: replace this with input in the future
    file_path = "/Users/brandonreyes/UW_work/EchoPro_work/UW_EchoProMatlab_Repackaged/outputs/Runs-brandon/EchoPro_kriged_aged_output-2019_0.xlsx"
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

    assert np.all(np.isclose(df_known.to_numpy(), df_produced.to_numpy()))

    assert np.all(np.isclose(male_df_known.to_numpy(), male_df_produced.to_numpy()))

    assert np.all(np.isclose(female_df_known.to_numpy(), female_df_produced.to_numpy()))




