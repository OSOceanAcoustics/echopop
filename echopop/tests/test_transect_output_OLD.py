import pytest

pytestmark = pytest.mark.skip("Disable all existing tests to revamp testing mechanism.")


import pandas as pd
import numpy as np
import echopop


def test_biomass_age_output(config_base_path, matlab_output_base_path):
    # TODO: formalize this test

    # initialize Survey object
    survey_2019 = echopop.Survey(
        init_file_path=config_base_path / 'initialization_config.yml',
        survey_year_file_path=config_base_path / 'survey_year_2019_config.yml',
        source=3,
        exclude_age1=True
    )

    # load all data
    survey_2019.load_survey_data()

    # compute all transect variables
    survey_2019.compute_transect_results()

    # set variables to improve readability
    transect_results = survey_2019.bio_calc.transect_results_gdf
    transect_results_male = survey_2019.bio_calc.transect_results_male_gdf
    transect_results_female = survey_2019.bio_calc.transect_results_female_gdf

    # obtain file_path pointing to the known output and associated sheet names
    file_path = matlab_output_base_path / "EchoPro_un-kriged_aged_output-2019_0.xlsx"
    sheet_name = "Sheet1"
    sheet_name_male = "Sheet2"
    sheet_name_female = "Sheet3"

    # gather known solution data produced by the Matlab version of EchoPro
    def _read_transform_biomass_age(sheet_name):
        df = (
            pd.read_excel(
                file_path,
                sheet_name=sheet_name,
                skiprows=1,
                usecols="A:AA"
            )
            .drop(columns=["stratum", "Transect"])
            # sort known solution dfs by latitude and longitude
            .sort_values(by=["Lat", "Lon"])
        )
        return df

    df_known = _read_transform_biomass_age(sheet_name)
    male_df_known = _read_transform_biomass_age(sheet_name_male)
    female_df_known = _read_transform_biomass_age(sheet_name_female)

    # define columns to compare in the produced solution
    produced_wanted_columns = ["latitude", "longitude", "biomass_adult"] + [
        "biomass_age_bin_" + str(i + 1) for i in range(22)]

    df_produced = transect_results.sort_values(
        by=["latitude", "longitude"]
    )[produced_wanted_columns]
    male_df_produced = transect_results_male.sort_values(
        by=["latitude", "longitude"]
    )[produced_wanted_columns]
    female_df_produced = transect_results_female.sort_values(
        by=["latitude", "longitude"]
    )[produced_wanted_columns]

    # compare known and produced values
    assert np.all(np.isclose(df_known.to_numpy(), df_produced.to_numpy()))

    assert np.all(np.isclose(male_df_known.to_numpy(), male_df_produced.to_numpy()))

    assert np.all(np.isclose(female_df_known.to_numpy(), female_df_produced.to_numpy()))


def test_core_output(config_base_path, matlab_output_base_path):
    # TODO: formalize this test

    # initialize Survey object
    survey_2019 = echopop.Survey(
        init_file_path=config_base_path / 'initialization_config.yml',
        survey_year_file_path=config_base_path / 'survey_year_2019_config.yml',
        source=3,
        exclude_age1=True
    )

    # load all data
    survey_2019.load_survey_data()

    # compute all transect variables
    survey_2019.compute_transect_results()

    # set variables to improve readability
    transect_results = survey_2019.bio_calc.transect_results_gdf
    transect_results_male = survey_2019.bio_calc.transect_results_male_gdf
    transect_results_female = survey_2019.bio_calc.transect_results_female_gdf

    # obtain file_path pointing to the known output and associated sheet name
    file_path = matlab_output_base_path / "EchoPro_un-kriged_output-26-Jan-2023_0.xlsx"
    sheet_name = "Sheet1"

    # gather known solution data produced by the Matlab version of EchoPro
    df_known = pd.read_excel(file_path, sheet_name=sheet_name, usecols="E:F,J:U")

    # sort known solution dfs by latitude and longitude
    df_known.sort_values(by=["Lat", "Lon"], inplace=True)

    # define columns grab from the produced results
    wanted_columns = ["latitude", "longitude", "stratum_num", "biomass_adult", "biomass_density_adult",
                      "numerical_density_adult", "abundance_adult", "transect_spacing", "interval"]
    gender_wanted_columns = ["latitude", "longitude", "biomass_adult", "biomass_density",
                             "abundance_adult", "numerical_density"]

    # sort produced results and format them, so they can be compared to known results
    df_produced = transect_results.sort_values(by=["latitude",
                                                   "longitude"])[wanted_columns]
    male_df_produced = transect_results_male.sort_values(by=["latitude",
                                                             "longitude"])[gender_wanted_columns].drop(
        columns=["latitude", "longitude"]).rename(
        columns={"biomass_adult": "biomass_male_adult", "biomass_density": "biomass_density_male",
                 "abundance_adult": "abundance_male_adult",
                 "numerical_density": "numerical_density_male"})
    female_df_produced = transect_results_female.sort_values(by=["latitude",
                                                                 "longitude"])[gender_wanted_columns].drop(
        columns=["latitude", "longitude"]).rename(
        columns={"biomass_adult": "biomass_female_adult", "biomass_density": "biomass_density_female",
                 "abundance_adult": "abundance_female_adult",
                 "numerical_density": "numerical_density_female"})

    # put together results to match the known solution DataFrame
    final_df_produced = pd.concat([df_produced, male_df_produced, female_df_produced], axis=1)

    # produced column names in the same order as the known solution
    ordered_columns = ["latitude", "longitude", "abundance_male_adult",
                       "abundance_female_adult", "abundance_adult", "biomass_male_adult", "biomass_female_adult",
                       "biomass_adult", "numerical_density_male", "numerical_density_female",
                       "numerical_density_adult", "biomass_density_male", "biomass_density_female",
                       "biomass_density_adult"]

    # compare known and produced values
    assert np.all(np.isclose(df_known.to_numpy(), final_df_produced[ordered_columns].to_numpy()))
