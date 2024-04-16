import numpy as np
import pandas as pd

from echopop.tests.conftest import assert_dataframe_equal


def test_fit_binned_length_weight_relationship(mock_survey):

    # Initialize mock_survey for `length_weight`
    mock_survey.statistics["length_weight"] = {}

    # Re-parameterize `specimen_df` with dummy data
    mock_survey.biology["specimen_df"] = pd.DataFrame(
        {
            "stratum_num": [0, 0, 1, 1, 2, 2, 3, 3],
            "haul_num": [1, 1, 2, 2, 3, 3, 4, 4],
            "sex": np.tile(["male", "female"], 4),
            "group": np.repeat("sexed", 8),
            "species_id": np.repeat([8675309], 8),
            "length": [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            "weight": [4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0],
        }
    )

    # Re-parameterize `length_bins` with dummy data
    mock_survey.biology["distributions"]["length"]["length_bins_arr"] = [2.0, 5.0, 8.0, 11.0]

    # Re-parameterize `length_interval` with dummy data
    mock_survey.biology["distributions"]["length"]["length_interval_arr"] = [
        0.5,
        3.5,
        6.5,
        9.5,
        12.5,
    ]

    # Evaluate object for later comparison
    mock_survey.fit_binned_length_weight_relationship(species_id=8675309)

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected data types
    expected_dtypes = {
        "regression_parameters": {
            "sex": object,
            "rate": np.floating,
            "initial": np.floating,
        },
        "length_weight_df": {
            "length_bin": pd.CategoricalDtype(),
            "sex": object,
            "mean_length": np.floating,
            "n_length": np.integer,
            "mean_weight": np.floating,
            "n_weight": np.integer,
            "rate": np.floating,
            "initial": np.floating,
            "weight_fitted": np.floating,
            "weight_modeled": np.floating,
        },
    }
    # ---- Expected output
    expected_output = {
        "regression_parameters": pd.DataFrame(
            {
                "sex": ["all", "female", "male"],
                "rate": [2.0, 2.0, 2.0],
                "initial": [4.710277e-16, -2.220446e-16, 1.110223e-15],
            },
        ),
        "length_weight_df": pd.DataFrame(
            {
                "length_bin": pd.cut(
                    np.repeat([1, 4, 7, 10], 3), np.array([0.5, 3.5, 6.5, 9.5, 12.5])
                ),
                "sex": np.tile(["all", "female", "male"], 4),
                "mean_length": [2.5, 3.0, 2.0, 5.0, 5.0, 5.0, 8.0, 8.0, 8.0, 0.0, 0.0, 0.0],
                "n_length": [2, 1, 1, 3, 1, 2, 3, 2, 1, 0, 0, 0],
                "mean_weight": [
                    6.50,
                    9.00,
                    4.00,
                    25.6666667,
                    25.00,
                    26.00,
                    64.6666667,
                    65.00,
                    64.00,
                    0.00,
                    0.00,
                    0.00,
                ],
                "n_weight": [2, 1, 1, 3, 1, 2, 3, 2, 1, 0, 0, 0],
                "rate": np.repeat(2.0, 12),
                "initial": np.tile([4.710277e-16, -2.220446e-16, 1.110223e-15], 4),
                "weight_fitted": [
                    4.0,
                    4.0,
                    4.0,
                    25.0,
                    25.0,
                    25.0,
                    64.0,
                    64.0,
                    64.0,
                    121.0,
                    121.0,
                    121.0,
                ],
                "weight_modeled": [
                    4.0,
                    4.0,
                    4.0,
                    25.0,
                    25.0,
                    25.0,
                    64.0,
                    64.0,
                    64.0,
                    121.0,
                    121.0,
                    121.0,
                ],
            },
        ),
    }
    # ----------------------------------
    # Run tests: `fit_binned_length_weight_relationship`
    # ----------------------------------
    eval_dictionary = mock_survey.statistics["length_weight"]
    assert_dataframe_equal(eval_dictionary, expected_dtypes, expected_output)


def test_strata_sex_weight_proportions(mock_survey):

    # Initialize mock_survey for `weight`
    mock_survey.biology["weight"] = {}

    # Initialize mock_survey for `length_weight`
    mock_survey.statistics["length_weight"] = {}

    # Re-parameterize `specimen_df` with dummy data
    mock_survey.biology["specimen_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4).astype(np.int64),
            "sex": np.tile(["male", "female"], 4),
            "group": np.repeat("sexed", 8),
            "haul_num": np.tile([1, 2], 4),
            "species_id": np.repeat([8675309], 8),
            "length": [12.0, 12.0, 19.0, 19.0, 12.0, 12.0, 19.0, 19.0],
            "weight": [2.0, 3.0, 3.0, 2.0, 2.0, 3.0, 2.0, 3.0],
            "age": [1, 1, 2, 2, 1, 1, 2, 2],
        }
    )

    # Re-parameterize `length_df` with dummy data
    mock_survey.biology["length_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4).astype(np.int64),
            "haul_num": [1, 1, 2, 2, 3, 3, 4, 4],
            "sex": np.tile(["male", "female"], 4),
            "group": np.repeat("sexed", 8),
            "species_id": np.repeat([8675309], 8),
            "length": [12, 12, 19, 19, 12, 12, 19, 19],
            "length_count": [5, 10, 15, 20, 20, 15, 10, 5],
        }
    )

    # Re-parameterize `fitted_weight` with dummy data
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

    # Re-parameterize `length_df` with dummy data
    mock_survey.biology["length_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4).astype(np.int64),
            "sex": np.tile(["male", "female"], 4),
            "group": np.repeat("sexed", 8),
            "species_id": np.repeat([8675309], 8),
            "length": [12, 12, 19, 19, 12, 12, 19, 19],
            "length_count": [5, 10, 15, 20, 20, 15, 10, 5],
        }
    )

    # Re-parameterize `fitted_weight` with dummy data
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

    # Re-parameterize `length_bins` with dummy data
    mock_survey.biology["distributions"]["length"]["length_interval_arr"] = np.linspace(9, 21, 3)

    # Evaluate object for later comparison
    mock_survey.strata_sex_weight_proportions(species_id=8675309)

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected data types
    expected_dtypes = {
        "stratum_num": np.integer,
        "proportion_female": np.floating,
        "proportion_male": np.floating,
        "proportion_station_1": np.floating,
        "proportion_station_2": np.floating,
        "average_weight_female": np.floating,
        "average_weight_male": np.floating,
        "average_weight_total": np.floating,
    }
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            "stratum_num": np.array([0, 1]).astype(int),
            "proportion_female": [0.592593, 0.407407],
            "proportion_male": [0.407407, 0.592593],
            "proportion_station_1": [0.925926, 0.925926],
            "proportion_station_2": [0.074074, 0.074074],
            "average_weight_female": [4.719110, 2.707892],
            "average_weight_male": [6.640487, 6.299942],
            "average_weight_total": [3.066481, 2.603519],
        },
    )

    # ----------------------------------
    # Run tests: `strata_sex_weight_proportions`
    # ----------------------------------
    eval_dataframe = mock_survey.biology["weight"]["weight_strata_df"]
    assert_dataframe_equal(eval_dataframe, expected_dtypes, expected_output)


def test_strata_age_binned_weight_proportions(mock_survey):

    # Initialize mock_survey for `weight`
    mock_survey.biology["weight"] = {}

    # Re-parameterize `specimen_df` with dummy data
    mock_survey.biology["specimen_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4),
            "sex": np.tile(["male", "female"], 4),
            "group": np.repeat("sexed", 8),
            "haul_num": [1, 1, 2, 2, 3, 3, 4, 4],
            "species_id": np.repeat([8675309], 8),
            "length": [12.0, 12.0, 19.0, 19.0, 12.0, 12.0, 19.0, 19.0],
            "weight": [2.0, 3.0, 3.0, 2.0, 2.0, 3.0, 2.0, 3.0],
            "age": [1, 1, 2, 2, 1, 1, 2, 2],
        },
    )

    # Re-parameterize `length_bins` with dummy data
    mock_survey.biology["distributions"]["length"]["length_interval_arr"] = np.linspace(9, 21, 3)

    # Evaluate object for later comparison
    mock_survey.strata_age_binned_weight_proportions(species_id=8675309)

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dtypes
    expected_dtypes = {
        "age_proportions_df": {
            "stratum_num": np.integer,
            "age": np.integer,
            "count_age_proportion_all": np.floating,
            "count_age_proportion_adult": np.floating,
        },
        "age_weight_proportions_df": {
            "stratum_num": np.integer,
            "age": np.integer,
            "weight_age_proportion_all": np.floating,
            "weight_age_proportion_adult": np.floating,
        },
        "sex_age_weight_proportions_df": {
            "stratum_num": np.integer,
            "age": np.integer,
            "sex": object,
            "weight_sex_proportion_all": np.floating,
            "weight_sex_proportion_adult": np.floating,
        },
        "length_sex_age_weight_proportions_df": {
            "stratum_num": np.integer,
            "age": np.integer,
            "length_bin": pd.CategoricalDtype(),
            "sex": object,
            "count": np.floating,
            "weight_total_all": np.floating,
            "weight_total_adult": np.floating,
            "weight_length_sex_proportion_all": np.floating,
            "weight_length_sex_proportion_adult": np.floating,
        },
    }
    # ---- Expected output
    expected_output = {
        "age_proportions_df": pd.DataFrame(
            {
                "stratum_num": np.repeat([0, 1], 2).astype(np.int64),
                "age": np.tile([1, 2], 2).astype(np.int64),
                "count_age_proportion_all": np.repeat(0.5, 4),
                "count_age_proportion_adult": [0.0, 1.0, 0.0, 1.0],
            }
        ),
        "age_weight_proportions_df": pd.DataFrame(
            {
                "stratum_num": np.repeat([0, 1], 2).astype(np.int64),
                "age": np.tile([1, 2], 2).astype(np.int64),
                "weight_age_proportion_all": [0.50, 0.50, 0.50, 0.50],
                "weight_age_proportion_adult": [0.0, 1.0, 0.0, 1.0],
            }
        ),
        "sex_age_weight_proportions_df": pd.DataFrame(
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
        ),
        "length_sex_age_weight_proportions_df": pd.DataFrame(
            {
                "stratum_num": np.repeat([0, 1], 12).astype(np.int64),
                "age": np.tile([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2], 2).astype(np.int64),
                "length_bin": pd.cut(
                    np.tile([12.0, 12.0, 12.0, 18.0, 18.0, 18.0], 4), np.linspace(9, 21, 3)
                ),
                "sex": np.tile(["all", "female", "male"], 8),
                "count": [
                    5.0,
                    3.0,
                    2.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    5.0,
                    2.0,
                    3.0,
                    5.0,
                    3.0,
                    2.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    5.0,
                    3.0,
                    2.0,
                ],
                "weight_total_all": [
                    10.0,
                    5.0,
                    5.0,
                    10.0,
                    5.0,
                    5.0,
                    10.0,
                    5.0,
                    5.0,
                    10.0,
                    5.0,
                    5.0,
                    10.0,
                    6.0,
                    4.0,
                    10.0,
                    6.0,
                    4.0,
                    10.0,
                    6.0,
                    4.0,
                    10.0,
                    6.0,
                    4.0,
                ],
                "weight_total_adult": [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    5.0,
                    2.0,
                    3.0,
                    5.0,
                    2.0,
                    3.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    5.0,
                    3.0,
                    2.0,
                    5.0,
                    3.0,
                    2.0,
                ],
                "weight_length_sex_proportion_all": [
                    0.5,
                    0.6,
                    0.4,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.5,
                    0.4,
                    0.6,
                    0.5,
                    0.5,
                    0.5,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.5,
                    0.5,
                    0.5,
                ],
                "weight_length_sex_proportion_adult": np.tile(
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2
                ),
            },
        ),
    }

    # ----------------------------------
    # Run tests: `strata_age_binned_weight_proportions`
    # ----------------------------------
    eval_dictionary = mock_survey.biology["weight"]["proportions"]
    assert_dataframe_equal(eval_dictionary, expected_dtypes, expected_output)


def test_nasc_to_biomass_conversion(mock_survey):

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
        },
    )

    # Create mock data for `age_weight_proportions_df`
    mock_survey.biology["weight"]["proportions"]["age_weight_proportions_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 2).astype(np.int64),
            "age": np.tile([1, 2], 2).astype(np.int64),
            "weight_age_proportion_all": [0.50, 0.50, 0.50, 0.50],
            "weight_age_proportion_adult": [0.0, 1.0, 0.0, 1.0],
        },
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
        },
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
        },
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
        },
    )

    # Create mock data for `strata_mean` (sigma_bs)
    mock_survey.acoustics["sigma_bs"]["strata_mean"] = pd.DataFrame(
        {
            "stratum_num": [0, 1],
            "species_id": np.repeat(8675309, 2),
            "sigma_bs_mean": 1.630277e-8,
        },
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
        },
    )

    # Create mock data for `strata_df`
    mock_survey.spatial["strata_df"] = pd.DataFrame(
        {
            "stratum_num": [0, 1],
            "haul_num": [1, 2],
            "fraction_hake": [1.000, 0.500],
        },
    )

    # Evaluate object for later comparison
    mock_survey.nasc_to_biomass_conversion(species_id=8675309)

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dtypes
    expected_dtypes = {
        "areal_density": {
            "number_density_df": {
                "transect_num": np.integer,
                "latitude": np.floating,
                "longitude": np.floating,
                "stratum_num": np.integer,
                "sex": object,
                "rho_a": np.floating,
                "age": np.integer,
                "count_age_proportion_all": np.floating,
                "count_age_proportion_adult": np.floating,
                "rho_a_adult": np.floating,
            },
            "biomass_density_df": {
                "transect_num": np.integer,
                "latitude": np.floating,
                "longitude": np.floating,
                "stratum_num": np.integer,
                "sex": object,
                "B_a": np.floating,
                "age": np.integer,
                "count_age_proportion_all": np.floating,
                "count_age_proportion_adult": np.floating,
                "B_a_adult": np.floating,
            },
        },
        "abundance": {
            "abundance_df": {
                "transect_num": np.integer,
                "latitude": np.floating,
                "longitude": np.floating,
                "stratum_num": np.integer,
                "sex": object,
                "NASC_all_ages": np.floating,
                "NASC_no_age1": np.floating,
                "N": np.floating,
                "age": np.integer,
                "count_age_proportion_all": np.floating,
                "count_age_proportion_adult": np.floating,
                "N_adult": np.floating,
            },
        },
        "biomass": {
            "biomass_df": {
                "transect_num": np.integer,
                "latitude": np.floating,
                "longitude": np.floating,
                "stratum_num": np.integer,
                "sex": object,
                "B": np.floating,
                "age": np.integer,
                "count_age_proportion_all": np.floating,
                "count_age_proportion_adult": np.floating,
                "B_adult": np.floating,
            },
            "biomass_age_df": {
                "transect_num": np.integer,
                "latitude": np.floating,
                "longitude": np.floating,
                "stratum_num": np.integer,
                "age": np.integer,
                "sex": object,
                "age_proportion": np.floating,
                "B_age": np.floating,
            },
        },
    }
    # ----- Expected output
    expected_output = {
        "areal_density": {
            "number_density_df": pd.DataFrame(
                {
                    "transect_num": np.repeat([1, 2, 3, 4], 8).astype(np.int64),
                    "latitude": np.repeat([20.0, 30.0, 40.0, 50.0], 8),
                    "longitude": np.repeat([-180.0, -120.0, -170.0, -110.0], 8),
                    "stratum_num": np.repeat([0, 1], 16).astype(np.int64),
                    "sex": np.tile(
                        ["all", "all", "male", "male", "female", "female", "unsexed", "unsexed"], 4
                    ),
                    "rho_a": np.concatenate(
                        [
                            np.repeat(0.0, 8),
                            [
                                4.881224e7,
                                4.881224e7,
                                1.988645e7,
                                1.988645e7,
                                2.892579e7,
                                2.892579e7,
                                0.0,
                                0.0,
                                2.440612e8,
                                2.440612e8,
                                1.446290e8,
                                1.446290e8,
                                9.943224e7,
                                9.943224e7,
                                0.0,
                                0.0,
                                2.440612e9,
                                2.440612e9,
                                1.446290e9,
                                1.446290e9,
                                9.943224e8,
                                9.943224e8,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                    "age": np.tile([1, 2], 16).astype(np.int64),
                    "count_age_proportion_all": np.repeat(0.5, 32),
                    "count_age_proportion_adult": np.tile([0.0, 1.0], 16),
                    "rho_a_adult": np.concatenate(
                        [
                            np.repeat(0.0, 9),
                            [
                                4.881224e7,
                                0.0,
                                1.988645e7,
                                0.0,
                                2.892579e7,
                                0.0,
                                0.0,
                                0.0,
                                2.440612e8,
                                0.0,
                                1.446290e8,
                                0.0,
                                9.943224e7,
                                0.0,
                                0.0,
                                0.0,
                                2.440612e9,
                                0.0,
                                1.446290e9,
                                0.0,
                                9.943224e8,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                }
            ),
            "biomass_density_df": pd.DataFrame(
                {
                    "transect_num": np.repeat([1, 2, 3, 4], 8).astype(np.int64),
                    "latitude": np.repeat([20.0, 30.0, 40.0, 50.0], 8),
                    "longitude": np.repeat([-180.0, -120.0, -170.0, -110.0], 8),
                    "stratum_num": np.repeat([0, 1], 16).astype(np.int64),
                    "sex": np.tile(
                        ["all", "all", "male", "male", "female", "female", "unsexed", "unsexed"], 4
                    ),
                    "B_a": np.concatenate(
                        [
                            np.repeat(0.0, 8),
                            [
                                1.496818e8,
                                1.496818e8,
                                1.320557e8,
                                1.320557e8,
                                1.365040e8,
                                1.365040e8,
                                0.0,
                                0.0,
                                6.354180e8,
                                6.354180e8,
                                9.111540e8,
                                9.111540e8,
                                2.692518e8,
                                2.692518e8,
                                0.0,
                                0.0,
                                6.354180e9,
                                6.354180e9,
                                9.111540e9,
                                9.111540e9,
                                2.692518e9,
                                2.692518e9,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                    "age": np.tile([1, 2], 16).astype(np.int64),
                    "count_age_proportion_all": np.repeat(0.5, 32),
                    "count_age_proportion_adult": np.tile([0.0, 1.0], 16),
                    "B_a_adult": np.concatenate(
                        [
                            np.repeat(0.0, 9),
                            [
                                1.496818e8,
                                0.0,
                                1.320557e8,
                                0.0,
                                1.365040e8,
                                0.0,
                                0.0,
                                0.0,
                                6.354180e8,
                                0.0,
                                9.111540e8,
                                0.0,
                                2.692518e8,
                                0.0,
                                0.0,
                                0.0,
                                6.354180e9,
                                0.0,
                                9.111540e9,
                                0.0,
                                2.692518e9,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                }
            ),
        },
        "abundance": {
            "abundance_df": pd.DataFrame(
                {
                    "transect_num": np.repeat([1, 2, 3, 4], 8).astype(np.int64),
                    "latitude": np.repeat([20.0, 30.0, 40.0, 50.0], 8),
                    "longitude": np.repeat([-180.0, -120.0, -170.0, -110.0], 8),
                    "stratum_num": np.repeat([0, 1], 16).astype(np.int64),
                    "sex": np.tile(
                        ["all", "all", "male", "male", "female", "female", "unsexed", "unsexed"], 4
                    ),
                    "NASC_all_ages": np.concatenate(
                        [np.repeat(1e1, 8), np.repeat(1e2, 16), np.repeat(1e3, 8)]
                    ),
                    "NASC_no_age1": np.concatenate(
                        [np.repeat(0, 8), np.repeat(1e1, 8), np.repeat(1e2, 8), np.repeat(1e3, 8)]
                    ),
                    "N": np.concatenate(
                        [
                            np.repeat(0.0, 8),
                            [
                                4.881224e8,
                                4.881224e8,
                                1.988645e8,
                                1.988645e8,
                                2.892579e8,
                                2.892579e8,
                                0.0,
                                0.0,
                                2.440612e9,
                                2.440612e9,
                                1.44629e9,
                                1.44629e9,
                                9.943224e8,
                                9.943224e8,
                                0.0,
                                0.0,
                                2.416206e10,
                                2.416206e10,
                                1.431827e10,
                                1.431827e10,
                                9.843792e9,
                                9.843792e9,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                    "age": np.tile([1, 2], 16).astype(np.int64),
                    "count_age_proportion_all": np.repeat(0.5, 32),
                    "count_age_proportion_adult": np.tile([0.0, 1.0], 16),
                    "N_adult": np.concatenate(
                        [
                            np.repeat(0.0, 9),
                            [
                                4.881224e8,
                                0.0,
                                1.988645e8,
                                0.0,
                                2.892579e8,
                                0.0,
                                0.0,
                                0.0,
                                2.440612e9,
                                0.0,
                                1.44629e9,
                                0.0,
                                9.943224e8,
                                0.0,
                                0.0,
                                0.0,
                                2.416206e10,
                                0.0,
                                1.431827e10,
                                0.0,
                                9.843792e9,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                }
            ),
        },
        "biomass": {
            "biomass_df": pd.DataFrame(
                {
                    "transect_num": np.repeat([1, 2, 3, 4], 8).astype(np.int64),
                    "latitude": np.repeat([20.0, 30.0, 40.0, 50.0], 8),
                    "longitude": np.repeat([-180.0, -120.0, -170.0, -110.0], 8),
                    "stratum_num": np.repeat([0, 1], 16).astype(np.int64),
                    "sex": np.tile(
                        ["all", "all", "male", "male", "female", "female", "unsexed", "unsexed"], 4
                    ),
                    "B": np.concatenate(
                        [
                            np.repeat(0.0, 8),
                            [
                                1.496818e9,
                                1.496818e9,
                                1.320557e9,
                                1.320557e9,
                                1.365040e9,
                                1.365040e9,
                                0.0,
                                0.0,
                                6.354180e9,
                                6.354180e9,
                                9.111540e9,
                                9.111540e9,
                                2.692518e9,
                                2.692518e9,
                                0.0,
                                0.0,
                                6.290638e10,
                                6.290638e10,
                                9.020425e10,
                                9.020425e10,
                                2.665593e10,
                                2.665593e10,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                    "age": np.tile([1, 2], 16).astype(np.int64),
                    "count_age_proportion_all": np.repeat(0.5, 32),
                    "count_age_proportion_adult": np.tile([0.0, 1.0], 16),
                    "B_adult": np.concatenate(
                        [
                            np.repeat(0.0, 9),
                            [
                                1.496818e9,
                                0.0,
                                1.320557e9,
                                0.0,
                                1.365040e9,
                                0.0,
                                0.0,
                                0.0,
                                6.354180e9,
                                0.0,
                                9.111540e9,
                                0.0,
                                2.692518e9,
                                0.0,
                                0.0,
                                0.0,
                                6.290638e10,
                                0.0,
                                9.020425e10,
                                0.0,
                                2.665593e10,
                                0.0,
                                0.0,
                            ],
                        ]
                    ),
                }
            ),
            "biomass_age_df": pd.DataFrame(
                {
                    "transect_num": np.repeat([1, 2, 3, 4], 6).astype(np.int64),
                    "latitude": np.repeat([20.0, 30.0, 40.0, 50.0], 6),
                    "longitude": np.repeat([-180.0, -120.0, -170.0, -110.0], 6),
                    "stratum_num": np.repeat([0, 1], 12).astype(np.int64),
                    "age": np.tile([1, 2], 12).astype(np.int64),
                    "sex": np.tile(["all", "all", "male", "male", "female", "female"], 4),
                    "age_proportion": np.tile([0.0, 1.0], 12),
                    "B_age": np.concatenate(
                        [
                            np.repeat(0.0, 7),
                            [
                                1.496818e9,
                                0.000,
                                1.320557e9,
                                0.000,
                                1.365040e9,
                                0.000,
                                6.354180e9,
                                0.000,
                                9.111540e9,
                                0.000,
                                2.692518e9,
                                0.000,
                                6.290638e10,
                                0.000,
                                9.020425e10,
                                0.000,
                                2.665593e10,
                            ],
                        ]
                    ),
                }
            ),
        },
    }

    # ----------------------------------
    # Run tests: `test_nasc_to_biomass_conversion`
    # ----------------------------------
    eval_dictionary = mock_survey.biology["population"]
    assert_dataframe_equal(eval_dictionary, expected_dtypes, expected_output)
