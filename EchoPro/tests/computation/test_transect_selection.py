import pytest

import os
import numpy as np
import EchoPro
from EchoPro.computation import SemiVariogram as SV


def test_transect_selection_output(config_base_path):

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
