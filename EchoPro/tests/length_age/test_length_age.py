import pytest

import os
import pandas as pd
import numpy as np
import EchoPro


def test_transect_based_length_age():

    # TODO: formalize this test

    # change working directory so no initialization files need to be modified
    # TODO: this may not be necessary in the future
    os.chdir("../../../example_notebooks")

    # initialize Survey object
    survey_2019 = EchoPro.Survey(init_file_path='../config_files/initialization_config.yml',
                                 survey_year_file_path='../config_files/survey_year_2019_config.yml',
                                 source=3,
                                 exclude_age1=True)

    # load all data
    survey_2019.load_survey_data()

    # compute all transect variables
    survey_2019.compute_transect_results()

    # compute abundance and biomass for each length and age bin
    survey_2019.compute_length_age_variables(data="transect")

    # set variables to improve readability
    abundance_df = survey_2019.bio_calc.transect_len_age_abundance
    male_abundance_df = survey_2019.bio_calc.transect_len_age_abundance_male
    female_abundance_df = survey_2019.bio_calc.transect_len_age_abundance_female
    biomass_df = survey_2019.bio_calc.transect_len_age_biomass
    male_biomass_df = survey_2019.bio_calc.transect_len_age_biomass_male
    female_biomass_df = survey_2019.bio_calc.transect_len_age_biomass_female

    # obtain file_path pointing to the known abundance data
    # TODO: replace this with input in the future
    file_path_abundance = "/Users/brandonreyes/UW_work/EchoPro_work/UW_EchoProMatlab_Repackaged/outputs/EchoPro_matlab_output_brandon_age_22_end_bin/un-kriged_len_age_abundance_table.xlsx"
    sheet_name = "Sheet3"
    sheet_name_male = "Sheet1"
    sheet_name_female = "Sheet2"

    # gather known solution data produced by the Matlab version of EchoPro
    df_known = pd.read_excel(file_path_abundance, index_col=0, sheet_name=sheet_name, skiprows=[0, 42, 43, 44, 45, 46], usecols="A:X")
    male_df_known = pd.read_excel(file_path_abundance, index_col=0, sheet_name=sheet_name_male,
                                  skiprows=[0, 42, 43, 44, 45, 46], usecols="A:X")
    female_df_known = pd.read_excel(file_path_abundance, index_col=0, sheet_name=sheet_name_female,
                                    skiprows=[0, 42, 43, 44, 45, 46], usecols="A:X")

    # compare known and produced values for abundance
    assert np.all(np.isclose(df_known.to_numpy(), abundance_df.to_numpy()))

    assert np.all(np.isclose(male_df_known.to_numpy(), male_abundance_df.to_numpy()))

    assert np.all(np.isclose(female_df_known.to_numpy(), female_abundance_df.to_numpy()))

    # obtain file_path pointing to the known biomass data
    # TODO: replace this with input in the future
    file_path_biomass = "/Users/brandonreyes/UW_work/EchoPro_work/UW_EchoProMatlab_Repackaged/outputs/EchoPro_matlab_output_brandon_age_22_end_bin/un-kriged_len_age_biomass_table.xlsx"

    # gather known solution data produced by the Matlab version of EchoPro
    df_known = pd.read_excel(file_path_biomass, index_col=0, sheet_name=sheet_name,
                             skiprows=[0, 42, 43, 44, 45, 46], usecols="A:W")
    male_df_known = pd.read_excel(file_path_biomass, index_col=0, sheet_name=sheet_name_male,
                                  skiprows=[0, 42, 43, 44, 45, 46], usecols="A:W")
    female_df_known = pd.read_excel(file_path_biomass, index_col=0, sheet_name=sheet_name_female,
                                    skiprows=[0, 42, 43, 44, 45, 46], usecols="A:W")

    # compare known and produced values for biomass
    assert np.all(np.isclose(df_known.to_numpy(), 1e-9*biomass_df.to_numpy()))

    assert np.all(np.isclose(male_df_known.to_numpy(), 1e-9*male_biomass_df.to_numpy()))

    assert np.all(np.isclose(female_df_known.to_numpy(), 1e-9*female_biomass_df.to_numpy()))


def test_kriging_based_length_age():

    # TODO: construct this test once Kriging results are available

    pass
