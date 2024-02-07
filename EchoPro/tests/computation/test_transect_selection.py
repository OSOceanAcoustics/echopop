import pytest

pytestmark = pytest.mark.skip("Disable all existing tests to revamp testing mechanism.")



import math

import numpy as np
import EchoPro
from EchoPro.computation import SemiVariogram as SV


# 12/1/2023: The test typically succeeds for a few iterations, then fails
@pytest.mark.skip(reason="Fails due to https://github.com/uw-echospace/EchoPro/issues/95")
@pytest.mark.parametrize(["removal_percentage"], [
    [38.0], [50.0], [60.0]
])
def test_transects_selection_success(removal_percentage, config_base_path):
    # initialize Survey object
    survey_2019 = EchoPro.Survey(
        init_file_path=config_base_path / 'initialization_config.yml',
        survey_year_file_path=config_base_path / 'survey_year_2019_config.yml',
        source=3,
        exclude_age1=True
    )

    # load all data
    survey_2019.load_survey_data()

    # obtain all unique transects in nasc_df
    unique_transects = survey_2019.nasc_df.index.unique().values

    # determine the number of transects that should be selected
    num_sel_transects = math.floor(len(unique_transects) * (1.0 - removal_percentage / 100.0))

    # initialize the random number generator object and fix the seed
    rng = np.random.default_rng(seed=1234)

    num_iterations = 5
    for iteration in range(num_iterations):
        # randomly select transects without replacement
        selected_transects = list(rng.choice(unique_transects, num_sel_transects, replace=False))

        print("Iteration #", iteration, len(selected_transects))

        # compute all transect variables
        survey_2019.compute_transect_results(selected_transects=selected_transects)

        nasc_fraction_adult_df_set = set(survey_2019.bio_calc.nasc_fraction_adult_df.index)
        nasc_stratum_num_set = set(survey_2019.bio_calc.nasc_df.stratum_num.unique())

        # This test is intended mainly to verify that survey_2019.compute_transect_results
        # runs successfully. Not sure what other assertion test would be better
        assert nasc_stratum_num_set <= nasc_fraction_adult_df_set
        # assert nasc_stratum_num_set <= nasc_fraction_adult_df_set | {0}


def test_transect_selection_output(config_base_path):
    # initialize Survey object
    survey_2019 = EchoPro.Survey(
        init_file_path=config_base_path / 'initialization_config.yml',
        survey_year_file_path=config_base_path / 'survey_year_2019_config.yml',
        source=3,
        exclude_age1=True
    )

    # load all data
    survey_2019.load_survey_data()

    # hard code a selection of transects to test output
    selected_transects = [1, 5, 6, 9, 15, 18, 22, 27, 28, 30, 35, 38, 41, 42, 46, 47, 52,
                          55, 57, 59, 60, 61, 63, 65, 66, 67, 68, 69, 70, 75, 77, 78, 79,
                          92, 94, 97, 98, 100, 105, 109, 117, 119, 121, 135, 140]

    # compute all transect variables
    survey_2019.compute_transect_results(selected_transects=selected_transects)

    # compare known and produced values for transect results
    assert np.isclose(survey_2019.bio_calc.transect_results_gdf.biomass_adult.sum(), 1589234887.734642)

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

    # compare known and produced values for kriging results
    assert np.isclose(survey_2019.bio_calc.kriging_results_gdf.biomass_adult.sum(), 1422611645.1972983)


def test_all_transects_selected_output(config_base_path):
    # The biomass produced should be the same as the case where no transects are selected

    # initialize Survey object
    survey_2019 = EchoPro.Survey(
        init_file_path=config_base_path / 'initialization_config.yml',
        survey_year_file_path=config_base_path / 'survey_year_2019_config.yml',
        source=3,
        exclude_age1=True
    )

    # load all data
    survey_2019.load_survey_data()

    # hard code a selection of transects to test output
    selected_transects = list(survey_2019.nasc_df.index.unique())

    # compute all transect variables
    survey_2019.compute_transect_results(selected_transects=selected_transects)

    # compare known and produced values for transect results
    assert np.isclose(survey_2019.bio_calc.transect_results_gdf.biomass_adult.sum(), 1643221106.9632864)
