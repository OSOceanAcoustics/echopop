import copy

import numpy as np
import pandas as pd

from echopop.computation.biology import index_transect_age_sex_proportions
from echopop.computation.spatial import (
    calculate_start_end_coordinates,
    calculate_transect_distance,
    correct_transect_intervals,
)
from echopop.tests.conftest import assert_dataframe_equal


def test_index_transect_age_sex_proportions(mock_survey):

    # Initialize various attributes
    mock_survey.acoustics["sigma_bs"] = {}
    mock_survey.statistics["length_weight"] = {}
    mock_survey.biology["weight"] = {}
    mock_survey.biology["population"] = {}

    # Create mock data for `age_proportions_df`
    mock_survey.biology["weight"]["proportions"] = {}
    mock_survey.biology["weight"]["proportions"]["age_proportions_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 2).astype(np.int64),
            "age": np.tile([1, 2], 2).astype(np.int64),
            "count_age_proportion_all": np.repeat(0.5, 4),
            "count_age_proportion_adult": [0.0, 1.0, 0.0, 1.0],
        }
    )

    # Create mock data for `age_weight_proportions_df`
    mock_survey.biology["weight"]["proportions"]["age_weight_proportions_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 2).astype(np.int64),
            "age": np.tile([1, 2], 2).astype(np.int64),
            "weight_age_proportion_all": [0.50, 0.50, 0.50, 0.50],
            "weight_age_proportion_adult": [0.0, 1.0, 0.0, 1.0],
        }
    )

    # Create mock data for `sex_age_weight_proportions_df`
    mock_survey.biology["weight"]["proportions"]["sex_age_weight_proportions_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 6).astype(np.int64),
            "age": np.tile([1, 1, 1, 2, 2, 2], 2).astype(np.int64),
            "sex": np.tile(["all", "female", "male"], 4),
            "weight_sex_proportion_all": [
                0.5,
                0.6,
                0.4,
                0.5,
                0.4,
                0.6,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
            ],
            "weight_sex_proportion_adult": np.tile([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2),
        }
    )

    # Create mock data for 'length_weight_df'
    mock_survey.statistics["length_weight"]["length_weight_df"] = pd.DataFrame(
        {
            "length_bin": pd.cut(np.repeat([12, 18], 3), np.linspace(9, 21, 3)),
            "sex": np.repeat(["all", "female", "male"], 2),
            "n_length": [4, 2, 2, 4, 2, 2],
            "mean_weight": [2.5, 3.5, 1.5, 7.5, 6.5, 8.5],
            "n_weight": [4, 2, 2, 4, 2, 2],
            "rate": [2.63, 1.36, 3.90, 2.63, 1.36, 3.90],
            "initial": [-2.49, -0.93, -4.06, -2.49, -0.93, -4.06],
            "weight_fitted": [2.21, 3.46, 1.41, 6.43, 6.02, 6.87],
            "weight_modeled": [2.21, 3.46, 1.41, 6.43, 6.02, 6.87],
        }
    )

    # Create mock data for `weight_strata_df`
    mock_survey.biology["weight"]["weight_strata_df"] = pd.DataFrame(
        {
            "stratum_num": [0, 1],
            "proportion_female": [0.592593, 0.407407],
            "proportion_male": [0.407407, 0.592593],
            "proportion_station_1": [0.925926, 0.925926],
            "proportion_station_2": [0.074074, 0.074074],
            "average_weight_female": [4.719110, 2.707892],
            "average_weight_male": [6.640487, 6.299942],
            "average_weight_total": [3.066481, 2.603519],
        }
    )

    # Create mock data for `strata_mean` (sigma_bs)
    mock_survey.acoustics["sigma_bs"]["strata_mean"] = pd.DataFrame(
        {"stratum_num": [0, 1], "species_id": np.repeat(8675309, 2), "sigma_bs_mean": 1.630277e-8}
    )

    # Create mock data for `nasc_df`
    mock_survey.acoustics["nasc"]["nasc_df"] = pd.DataFrame(
        {
            "transect_num": [1, 2, 3, 4],
            "stratum_num": [0, 0, 1, 1],
            "vessel_log_start": [0.0, 10.1, 20.1, 30.1],
            "vessel_log_end": [10.0, 20.0, 30.0, 40.0],
            "latitude": [20.0, 30.0, 40.0, 50.0],
            "longitude": [-180.0, -120.0, -170.0, -110.0],
            "transect_spacing": np.repeat(1.0, 4),
            "NASC_no_age1": [0.0, 1e1, 1e2, 1e3],
            "haul_num": [1, 1, 2, 2],
            "NASC_all_ages": [1e1, 1e2, 1e2, 1e3],
        }
    )

    # Create mock data for `strata_df`
    mock_survey.spatial["strata_df"] = pd.DataFrame(
        {"stratum_num": [0, 1], "haul_num": [1, 2], "fraction_hake": [1.000, 0.500]}
    )

    # Bundle the mocked data into their respective inputs for `index_transect_age_sex_proportions`
    test_acoustics_dict = copy.deepcopy(mock_survey.acoustics)
    test_biology_dict = copy.deepcopy(mock_survey.biology)
    test_info_strata = mock_survey.spatial["strata_df"].copy()

    # Evaluate object for later comparison
    eval_nasc_fraction_total_df = index_transect_age_sex_proportions(
        test_acoustics_dict, test_biology_dict, test_info_strata
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dtypes
    expected_dtypes = {
        "latitude": np.floating,
        "longitude": np.floating,
        "transect_num": np.integer,
        "stratum_num": np.integer,
        "haul_num": np.integer,
        "interval": np.floating,
        "interval_area": np.floating,
        "NASC_all_ages": np.floating,
        "NASC_no_age1": np.floating,
        "fraction_hake": np.floating,
        "species_id": np.integer,
        "sigma_bs_mean": np.floating,
        "proportion_female": np.floating,
        "proportion_male": np.floating,
        "proportion_station_1": np.floating,
        "proportion_station_2": np.floating,
        "average_weight_female": np.floating,
        "average_weight_male": np.floating,
        "average_weight_total": np.floating,
        "age": np.integer,
        "count_age_proportion_all": np.floating,
        "count_age_proportion_adult": np.floating,
        "weight_age_proportion_all": np.floating,
        "weight_age_proportion_adult": np.floating,
    }
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            "latitude": np.repeat([20.0, 30.0, 40.0, 50.0], 2),
            "longitude": np.repeat([-180.0, -120.0, -170.0, -110.0], 2),
            "transect_num": np.repeat([1, 2, 3, 4], 2).astype(np.int64),
            "stratum_num": np.repeat([0, 1], 4).astype(np.int64),
            "haul_num": np.repeat([1, 2], 4).astype(np.int64),
            "interval": np.repeat([10.0, 10.0, 10.0, 9.9], 2),
            "interval_area": np.repeat([10.0, 10.0, 10.0, 9.9], 2),
            "NASC_all_ages": np.repeat([1e1, 1e2, 1e2, 1e3], 2),
            "NASC_no_age1": np.repeat([0.0, 1e1, 1e2, 1e3], 2),
            "fraction_hake": np.repeat([1.0, 0.5], 4),
            "species_id": np.repeat(8675309, 8).astype(np.int64),
            "sigma_bs_mean": np.repeat(1.630277e-8, 8),
            "proportion_female": np.repeat([0.592593, 0.407407], 4),
            "proportion_male": np.repeat([0.407407, 0.592593], 4),
            "proportion_station_1": np.repeat(0.925926, 8),
            "proportion_station_2": np.repeat(0.074074, 8),
            "average_weight_female": np.repeat([4.719110, 2.707892], 4),
            "average_weight_male": np.repeat([6.640487, 6.299942], 4),
            "average_weight_total": np.repeat([3.066481, 2.603519], 4),
            "age": np.tile([1, 2], 4).astype(np.int64),
            "count_age_proportion_all": np.repeat(0.5, 8),
            "count_age_proportion_adult": np.tile([0.0, 1.0], 4),
            "weight_age_proportion_all": np.repeat(0.5, 8),
            "weight_age_proportion_adult": np.tile([0.0, 1.0], 4),
        },
    )

    # ----------------------------------
    # Run tests: `index_transect_age_sex_proportions`
    # ----------------------------------
    assert_dataframe_equal(eval_nasc_fraction_total_df, expected_dtypes, expected_output)


def test_correct_transect_intervals():

    # Create mock data for `nasc_df`
    test_nasc_dataframe = pd.DataFrame(
        {
            "transect_num": [1, 2, 3, 4],
            "stratum_num": [0, 0, 1, 1],
            "vessel_log_start": [0.0, 10.1, 20.1, 30.1],
            "vessel_log_end": [10.0, 20.0, 30.0, 40.0],
            "latitude": [20.0, 30.0, 40.0, 50.0],
            "longitude": [-180.0, -120.0, -170.0, -110.0],
            "transect_spacing": np.repeat(1.0, 4),
            "NASC_no_age1": [0.0, 1e1, 1e2, 1e3],
            "haul_num": [1, 1, 2, 2],
            "NASC_all_ages": [1e1, 1e2, 1e2, 1e3],
        },
    )

    # Evaluate object for later comparison
    eval_nasc_interval = correct_transect_intervals(test_nasc_dataframe)

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dtypes
    expected_dtypes = {
        "latitude": np.floating,
        "longitude": np.floating,
        "transect_num": np.integer,
        "stratum_num": np.integer,
        "haul_num": np.integer,
        "interval": np.floating,
        "interval_area": np.floating,
        "NASC_all_ages": np.floating,
        "NASC_no_age1": np.floating,
    }
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            "latitude": [20.0, 30.0, 40.0, 50.0],
            "longitude": [-180.0, -120.0, -170.0, -110.0],
            "transect_num": [1, 2, 3, 4],
            "stratum_num": [0, 0, 1, 1],
            "haul_num": [1, 1, 2, 2],
            "interval": [10.0, 10.0, 10.0, 9.9],
            "interval_area": [10.0, 10.0, 10.0, 9.9],
            "NASC_all_ages": [1e1, 1e2, 1e2, 1e3],
            "NASC_no_age1": [0.0, 1e1, 1e2, 1e3],
        },
    )

    # ----------------------------------
    # Run tests: `correct_transect_intervals`
    # ----------------------------------
    assert_dataframe_equal(eval_nasc_interval, expected_dtypes, expected_output)


def test_calculate_start_end_coordinates():

    # Create mock data for `nasc_df`
    test_nasc_dataframe = pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "stratum_num": [0, 0, 1, 1],
            "vessel_log_start": [0.0, 10.1, 20.1, 30.1],
            "vessel_log_end": [10.0, 20.0, 30.0, 40.0],
            "latitude": [20.0, 30.0, 40.0, 50.0],
            "longitude": [-180.0, -120.0, -170.0, -110.0],
            "transect_spacing": np.repeat(1.0, 4),
            "NASC_no_age1": [0.0, 1e1, 1e2, 1e3],
            "haul_num": [1, 1, 2, 2],
            "NASC_all_ages": [1e1, 1e2, 1e2, 1e3],
        },
    )

    # Evaluate for later comparison
    eval_test_nasc_df = calculate_start_end_coordinates(test_nasc_dataframe, "transect_num")

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dtypes
    expected_dtypes = {
        "transect_num": np.integer,
        "minimum_longitude": np.floating,
        "maximum_longitude": np.floating,
        "center_latitude": np.floating,
    }
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            "transect_num": [1, 2],
            "minimum_longitude": [-180.0, -170.0],
            "maximum_longitude": [-120.0, -110.0],
            "center_latitude": [25.0, 45.0],
        },
    )

    # ----------------------------------
    # Run tests: `calculate_start_end_coordinates`
    # ----------------------------------
    assert_dataframe_equal(eval_test_nasc_df, expected_dtypes, expected_output)


def test_calculate_transect_distance():

    # Create mock data for `nasc_df`
    test_nasc_dataframe = pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2],
            "stratum_num": [0, 0, 1, 1],
            "vessel_log_start": [0.0, 10.1, 20.1, 30.1],
            "vessel_log_end": [10.0, 20.0, 30.0, 40.0],
            "latitude": [20.0, 30.0, 40.0, 50.0],
            "longitude": [-180.0, -120.0, -170.0, -110.0],
            "transect_spacing": np.repeat(2.0, 4),
            "NASC_no_age1": [0.0, 1e1, 1e2, 1e3],
            "haul_num": [1, 1, 2, 2],
            "NASC_all_ages": [1e1, 1e2, 1e2, 1e3],
        },
    )
    # Evaluate for later comparison
    eval_test_nasc_df = calculate_transect_distance(test_nasc_dataframe, "transect_num")

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dtypes
    expected_dtypes = {
        "transect_num": np.integer,
        "minimum_longitude": np.floating,
        "maximum_longitude": np.floating,
        "center_latitude": np.floating,
        "transect_distance": np.floating,
        "transect_spacing": np.floating,
        "transect_area": np.floating,
    }
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            "transect_num": [1, 2],
            "minimum_longitude": [-180.0, -170.0],
            "maximum_longitude": [-120.0, -110.0],
            "center_latitude": [25.0, 45.0],
            "transect_distance": [3241.273891, 2493.203304],
            "transect_spacing": [2.0, 2.0],
            "transect_area": [6482.547781, 4986.406609],
        },
    )

    # ----------------------------------
    # Run tests: `calculate_transect_distance`
    # ----------------------------------
    assert_dataframe_equal(eval_test_nasc_df, expected_dtypes, expected_output)
