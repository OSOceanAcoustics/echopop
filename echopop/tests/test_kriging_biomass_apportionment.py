import numpy as np
import pandas as pd

from echopop.computation.biology import (
    aged_sex_weight_proportions,
    aged_weight_proportions,
    apply_age_bins,
    calculate_aged_biomass,
    calculate_aged_unaged_proportions,
    calculate_unaged_biomass,
    distribute_aged_weight_proportions,
    sum_strata_weight,
    unaged_number_proportions,
    unaged_sex_weight_proportions,
    unaged_weight_proportions,
)


def test_sum_strata_weight(mock_survey):

    # Re-parameterize `specimen_df` with dummy data
    mock_survey.biology["specimen_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1, 2, 4, 5], 4),
            "haul_num": np.repeat([1, 2, 3, 4, 5], 4),
            "species_id": np.repeat([19350], 20),
            "length": np.linspace(10, 100, 20),
            "weight": np.linspace(1, 5, 20),
        }
    )

    # Re-parameterize `length_df` with dummy data
    mock_survey.biology["length_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1, 2, 4, 5], 4),
            "haul_num": np.repeat([1, 2, 3, 4, 5], 4),
            "species_id": np.repeat([19350], 20),
            "length": np.linspace(10, 100, 20),
            "length_count": np.linspace(10, 100, 20),
        }
    )

    # Re-parameterize `catch_df` with dummy data
    mock_survey.biology["catch_df"] = pd.DataFrame(
        {
            "stratum_num": [0, 1, 2, 4, 5],
            "haul_num": [1, 2, 3, 4, 5],
            "haul_weight": [51.4, 0.6, 81.7, 16.2, 12.9],
        }
    )

    # Evaluate object for later comparison
    eval_weight_strata = sum_strata_weight(
        mock_survey.biology["catch_df"], mock_survey.biology["specimen_df"]
    )

    # ----------------------------------
    # Run tests: `sum_strata_weight`
    # ----------------------------------
    # Evaluate shape
    # ---- `object_weight_strata`
    assert eval_weight_strata.shape == tuple([5, 2])

    # Evaluate value equality
    # ---- `object_weight_strata`
    check_values = np.array([56.663158, 9.231579, 93.700000, 31.568421, 31.636842])
    assert np.allclose(check_values, eval_weight_strata.weight_stratum_all)


def test_compute_aged_weight_proportions(mock_survey):

    # Re-parameterize `specimen_df` with dummy data
    mock_survey.biology["specimen_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4),
            "sex": np.tile(["male", "female"], 4),
            "haul_num": np.tile([1, 2], 4),
            "species_id": np.repeat([19350], 8),
            "length": [12.0, 12.0, 19.0, 19.0, 12.0, 12.0, 19.0, 19.0],
            "weight": [2.0, 3.0, 8.0, 7.0, 1.0, 4.0, 9.0, 6.0],
            "age": [1, 1, 2, 2, 1, 1, 2, 2],
        }
    )

    # Length interval
    mock_survey.biology["distributions"]["length"]["length_interval_arr"] = np.linspace(9, 21, 3)

    # Age interval
    mock_survey.biology["distributions"]["age"]["age_interval_arr"] = np.array([0.5, 1.5, 2.5])

    # Evaluate object for later comparison
    eval_props_wgt_len_age_sex = aged_weight_proportions(
        mock_survey.biology["specimen_df"],
        mock_survey.biology["distributions"]["length"]["length_interval_arr"],
        mock_survey.biology["distributions"]["age"]["age_interval_arr"],
    )
    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dimensions of `obj_props_wgt_len_age_sex`
    expected_dimensions = tuple([16, 11])

    # ---- Expected dataframe output
    # ---- sex
    sex_array = (
        np.tile(["female", "female", "female", "female", "male", "male", "male", "male"], 2)
        .astype(object)
        .flatten()
    )
    # ---- total_weight_sex_all
    total_weight_sex_all_array = np.repeat(10, 16).astype(np.float64).flatten()
    total_weight_sex_all_array = np.concatenate(
        [total_weight_sex_all_array, np.repeat(20, 8).astype(np.float64).flatten()]
    )

    expected_output = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 8).astype(np.int64),
            "species_id": np.repeat(19350, 16).astype(np.int64),
            "sex": sex_array,
            "length_bin": pd.cut(np.tile([12, 12, 19, 19], 4), np.linspace(9, 21, 3)),
            "age_bin": pd.cut(np.tile([1, 2], 8), np.array([0.5, 1.5, 2.5])),
            "weight_all": [
                3.0,
                0.0,
                0.0,
                7.0,
                2.0,
                0.0,
                0.0,
                8.0,
                4.0,
                0.0,
                0.0,
                6.0,
                1.0,
                0.0,
                0.0,
                9.0,
            ],
            "weight_adult": [
                0.0,
                0.0,
                0.0,
                7.0,
                0.0,
                0.0,
                0.0,
                8.0,
                0.0,
                0.0,
                0.0,
                6.0,
                0.0,
                0.0,
                0.0,
                9.0,
            ],
            "total_weight_sex_all": np.repeat(10, 16).astype(np.float64).flatten(),
            "total_weight_sex_adult": [
                7.0,
                7.0,
                7.0,
                7.0,
                8.0,
                8.0,
                8.0,
                8.0,
                6.0,
                6.0,
                6.0,
                6.0,
                9.0,
                9.0,
                9.0,
                9.0,
            ],
            "proportion_weight_sex_all": [
                0.3,
                0.0,
                0.0,
                0.7,
                0.2,
                0.0,
                0.0,
                0.8,
                0.4,
                0.0,
                0.0,
                0.6,
                0.1,
                0.0,
                0.0,
                0.9,
            ],
            "proportion_weight_sex_adult": [
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
            ],
        }
    )
    expected_output["length_bin"] = pd.IntervalIndex(expected_output["length_bin"])
    expected_output["length_bin"] = pd.Categorical(
        expected_output["length_bin"],
        categories=expected_output["length_bin"].unique(),
        ordered=True,
    )
    expected_output["age_bin"] = pd.IntervalIndex(expected_output["age_bin"])
    expected_output["age_bin"] = pd.Categorical(
        expected_output["age_bin"], categories=expected_output["age_bin"].unique(), ordered=True
    )

    # ----------------------------------
    # Run tests: `compute_index_aged_weight_proportions`
    # ----------------------------------
    # Process the specimen data
    # Check shape
    assert eval_props_wgt_len_age_sex.shape == expected_dimensions

    # Check data value equality
    assert np.all(expected_output == eval_props_wgt_len_age_sex)


def test_calculate_aged_unaged_proportions():

    # Mock data for `specimen_data`
    test_specimen_data = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4),
            "sex": np.tile(["male", "female"], 4),
            "haul_num": np.tile([1, 2], 4),
            "species_id": np.repeat([19350], 8),
            "length": [12.0, 12.0, 19.0, 19.0, 12.0, 12.0, 19.0, 19.0],
            "weight": [2.0, 3.0, 8.0, 7.0, 1.0, 4.0, 9.0, 6.0],
            "age": [1, 1, 2, 2, 1, 1, 2, 2],
        }
    )

    # Mock data for `weight_strata`
    test_weight_strata = pd.DataFrame(
        {
            "stratum_num": [0, 1],
            "weight_stratum_all": [50.0, 100.0],
        }
    )

    # Evaluate for later comparison
    eval_aged_unaged_weight_proportions = calculate_aged_unaged_proportions(
        test_specimen_data, test_weight_strata
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # `eval_aged_sex_proportions`
    # ---- Expected dimensions
    expected_dimensions = tuple([2, 7])
    # ---- Expected dataframe output
    expected_output = pd.DataFrame(
        {
            "stratum_num": np.array([0.0, 1.0]).astype(np.int64),
            "weight": [20.0, 20.0],
            "weight_stratum_all": [50.0, 100.0],
            "proportion_aged_weight_all": [0.40, 0.20],
            "proportion_aged_weight_adult": [0.30, 0.15],
            "proportion_unaged_weight_all": [0.60, 0.80],
            "proportion_unaged_weight_adult": [0.70, 0.85],
        }
    )

    # ----------------------------------
    # Run tests: `compute_index_aged_weight_proportions`
    # ----------------------------------
    # `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_aged_unaged_weight_proportions.shape == expected_dimensions
    # ---- Dataframe equality
    assert np.allclose(eval_aged_unaged_weight_proportions, expected_output)


def test_aged_sex_weight_proportions():

    # Mock data for `proportions_weight_length_age_sex`
    test_proportions_weight_length_age_sex = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 8).astype(np.int64),
            "species_id": np.repeat(19350, 16).astype(np.int64),
            "sex": np.tile(
                ["female", "female", "female", "female", "male", "male", "male", "male"], 2
            )
            .astype(object)
            .flatten(),
            "length_bin": pd.cut(np.tile([12, 12, 19, 19], 4), np.linspace(9, 21, 3)),
            "age_bin": pd.cut(np.tile([1, 2], 8), np.array([0.5, 1.5, 2.5])),
            "weight_all": [
                3.0,
                0.0,
                0.0,
                7.0,
                2.0,
                0.0,
                0.0,
                8.0,
                4.0,
                0.0,
                0.0,
                6.0,
                1.0,
                0.0,
                0.0,
                9.0,
            ],
            "weight_adult": [
                0.0,
                0.0,
                0.0,
                7.0,
                0.0,
                0.0,
                0.0,
                8.0,
                0.0,
                0.0,
                0.0,
                6.0,
                0.0,
                0.0,
                0.0,
                9.0,
            ],
            "total_weight_sex_all": np.repeat(10, 16).astype(np.float64).flatten(),
            "total_weight_sex_adult": [
                7.0,
                7.0,
                7.0,
                7.0,
                8.0,
                8.0,
                8.0,
                8.0,
                6.0,
                6.0,
                6.0,
                6.0,
                9.0,
                9.0,
                9.0,
                9.0,
            ],
            "proportion_weight_sex_all": [
                0.3,
                0.0,
                0.0,
                0.7,
                0.2,
                0.0,
                0.0,
                0.8,
                0.4,
                0.0,
                0.0,
                0.6,
                0.1,
                0.0,
                0.0,
                0.9,
            ],
            "proportion_weight_sex_adult": [
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
            ],
        }
    )
    test_proportions_weight_length_age_sex["length_bin"] = pd.IntervalIndex(
        test_proportions_weight_length_age_sex["length_bin"]
    )
    test_proportions_weight_length_age_sex["length_bin"] = pd.Categorical(
        test_proportions_weight_length_age_sex["length_bin"],
        categories=test_proportions_weight_length_age_sex["length_bin"].unique(),
        ordered=True,
    )
    test_proportions_weight_length_age_sex["age_bin"] = pd.IntervalIndex(
        test_proportions_weight_length_age_sex["age_bin"]
    )
    test_proportions_weight_length_age_sex["age_bin"] = pd.Categorical(
        test_proportions_weight_length_age_sex["age_bin"],
        categories=test_proportions_weight_length_age_sex["age_bin"].unique(),
        ordered=True,
    )

    # Mock data for `aged_proportions`
    test_aged_proportions = pd.DataFrame(
        {
            "stratum_num": [0, 1],
            "weight": [100.0, 200.0],
            "weight_stratum_all": [200.00, 250.0],
            "proportion_aged_weight_all": [0.10, 0.20],
            "proportion_aged_weight_adult": [0.05, 0.10],
            "proportion_unaged_weight_all": [0.90, 0.80],
            "proportion_unaged_weight_adult": [0.45, 0.40],
        }
    )

    # Evaluate for later comparison
    eval_aged_sex_proportions = aged_sex_weight_proportions(
        test_proportions_weight_length_age_sex, test_aged_proportions
    )
    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # `eval_aged_sex_proportions`
    # ---- Expected dimensions
    expected_dimensions = tuple([4, 6])
    # ---- Expected dataframe output
    expected_output = pd.DataFrame(
        {
            "stratum_num": np.repeat([0.0, 1.0], 2).astype(np.int64),
            "sex": np.tile(["female", "male"], 2).astype(object),
            "proportion_aged_weight_all": [0.10, 0.10, 0.20, 0.20],
            "proportion_aged_weight_adult": [0.05, 0.05, 0.10, 0.10],
            "proportion_weight_all": [0.050, 0.050, 0.040, 0.040],
            "proportion_weight_adult": [0.035, 0.040, 0.024, 0.036],
        }
    )

    # ----------------------------------
    # Run tests: `compute_index_aged_weight_proportions`
    # ----------------------------------
    # `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_aged_sex_proportions.shape == expected_dimensions
    # ---- Datatypes
    assert np.all(eval_aged_sex_proportions.dtypes == expected_output.dtypes)
    # ---- Dataframe equality
    assert np.all(eval_aged_sex_proportions == expected_output)


def test_distribute_aged_weight_proportions():

    # Mock data for `proportions_weight_length_age_sex`
    test_proportions_weight_length_age_sex = pd.DataFrame(
        {
            "stratum_num": np.repeat([0.0, 1.0], 2),
            "sex": np.tile(["female", "male"], 2),
            "proportion_weight_sex_all": [0.50, 0.50, 0.25, 0.50],
            "proportion_weight_sex_adult": [0.50, 0.25, 0.00, 0.25],
        }
    )

    # Mock data for `aged_sex_proportions`
    test_aged_sex_proportions = pd.DataFrame(
        {
            "stratum_num": np.repeat([0.0, 1.0], 2).astype(np.int64),
            "sex": np.tile(["female", "male"], 2).astype(object),
            "proportion_aged_weight_all": [0.10, 0.10, 0.20, 0.20],
            "proportion_aged_weight_adult": [0.05, 0.05, 0.10, 0.10],
            "proportion_weight_all": [0.050, 0.050, 0.040, 0.040],
            "proportion_weight_adult": [0.035, 0.040, 0.024, 0.036],
        }
    )

    # Evaluate for later comparison
    eval_distributed_aged_weight_proportions = distribute_aged_weight_proportions(
        test_proportions_weight_length_age_sex, test_aged_sex_proportions
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # `eval_aged_sex_proportions`
    # ---- Expected dimensions
    expected_dimensions = tuple([4, 4])
    # ---- Expected dataframe output
    expected_output = pd.DataFrame(
        {
            "stratum_num": np.repeat([0.0, 1.0], 2),
            "sex": np.tile(["female", "male"], 2),
            "normalized_proportion_weight_all": [0.025, 0.025, 0.010, 0.020],
            "normalized_proportion_weight_adult": [0.0175, 0.0100, 0.00000, 0.0090],
        }
    )

    # ----------------------------------
    # Run tests: `distribute_aged_weight_proportions`
    # ----------------------------------
    # `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_distributed_aged_weight_proportions.shape == expected_dimensions
    # ---- Datatypes
    assert np.all(eval_distributed_aged_weight_proportions.dtypes == expected_output.dtypes)
    # ---- Dataframe equality
    assert np.all(eval_distributed_aged_weight_proportions == expected_output)


def test_calculate_aged_biomass():

    # Mock data for `kriging_biomass_df`
    test_kriging_biomass_df = pd.DataFrame(
        {
            "centroid_latitude": [30.0, 40.0, 50.0],
            "centroid_longitude": [30.0, 31.0, 32.0],
            "fraction_cell_in_polygon": [1.0, 1.0, 1.0],
            "stratum_num": [0.0, 1.0, 2.0],
            "B_a_adult_mean": [1e3, 1e4, 1e5],
            "B_a_adult_prediction_variance": [1e1, 1e2, 1e3],
            "B_a_adult_sample_variance": [1e1, 1e2, 1e3],
            "cell_area_nmi2": [10.0, 10.0, 5.0],
            "B_adult_kriged": [1e4, 1e5, 5e5],
        }
    )

    # Mock data for `specimen_data`
    test_specimen_data = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4),
            "sex": np.tile(["male", "female"], 4),
            "haul_num": np.tile([1, 2], 4),
            "species_id": np.repeat([19350], 8),
            "length": [12.0, 12.0, 19.0, 19.0, 12.0, 12.0, 19.0, 19.0],
            "weight": [2.0, 3.0, 3.0, 2.0, 2.0, 3.0, 2.0, 3.0],
            "age": [1, 1, 2, 2, 1, 1, 2, 2],
        }
    )

    # Mock data for `length_intervals`
    test_length_intervals = np.linspace(9, 21, 3)

    # Mock data for `age_intervals`
    test_age_intervals = np.array([0.5, 1.5, 2.5])

    # Mock data for `aged_proportions`
    test_aged_proportions = pd.DataFrame(
        {
            "stratum_num": np.array([0.0, 1.0]).astype(np.int64),
            "weight": [10.0, 10.0],
            "weight_stratum_all": [100.0, 200.0],
            "proportion_aged_weight_all": [0.10, 0.05],
            "proportion_aged_weight_adult": [0.050, 0.025],
            "proportion_unaged_weight_all": [0.90, 0.95],
            "proportion_unaged_weight_adult": [0.95, 0.975],
        }
    )

    # Evaluate for later comparison
    eval_aged_sex_biomass, eval_aged_biomass = calculate_aged_biomass(
        test_kriging_biomass_df,
        test_specimen_data,
        test_length_intervals,
        test_age_intervals,
        test_aged_proportions,
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # `eval_aged_sex_biomass`
    # ---- Expected dimensions
    expected_dimensions_aged_sex_biomass = tuple([8, 6])
    # ---- Expected dataframe output
    expected_output_aged_sex_biomass = pd.DataFrame(
        {
            "species_id": np.repeat(19350, 8).astype(np.int64),
            "sex": np.repeat(["female", "male"], 4),
            "length_bin": pd.cut(np.tile([12, 12, 19, 19], 2), np.linspace(9, 21, 3)),
            "age_bin": pd.cut(np.tile([1, 2], 4), np.array([0.5, 1.5, 2.5])),
            "biomass_sexed_aged_all": [1800.0, 0.0, 0.0, 1700.0, 1200.0, 0.0, 0.0, 1300.0],
            "biomass_sexed_aged_adult": [0.0, 0.0, 0.0, 1700.0, 0.0, 0.0, 0.0, 1300.0],
        }
    )
    expected_output_aged_sex_biomass["length_bin"] = pd.IntervalIndex(
        expected_output_aged_sex_biomass["length_bin"]
    )
    expected_output_aged_sex_biomass["length_bin"] = pd.Categorical(
        expected_output_aged_sex_biomass["length_bin"],
        categories=expected_output_aged_sex_biomass["length_bin"].unique(),
        ordered=True,
    )
    expected_output_aged_sex_biomass["age_bin"] = pd.IntervalIndex(
        expected_output_aged_sex_biomass["age_bin"]
    )
    expected_output_aged_sex_biomass["age_bin"] = pd.Categorical(
        expected_output_aged_sex_biomass["age_bin"],
        categories=expected_output_aged_sex_biomass["age_bin"].unique(),
        ordered=True,
    )
    # `eval_aged_biomass`
    # ---- Expected dimensions
    expected_dimensions_aged_biomass = tuple([4, 5])
    # ---- Expected dataframe output
    expected_output_aged_biomass = pd.DataFrame(
        {
            "species_id": np.repeat(19350, 4).astype(np.int64),
            "length_bin": pd.cut([12, 12, 19, 19], np.linspace(9, 21, 3)),
            "age_bin": pd.cut(np.tile([1, 2], 2), np.array([0.5, 1.5, 2.5])),
            "biomass_aged_all": [3000.0, 0.0, 0.0, 3000.0],
            "biomass_aged_adult": [0.0, 0.0, 0.0, 3000.0],
        }
    )
    expected_output_aged_biomass["length_bin"] = pd.IntervalIndex(
        expected_output_aged_biomass["length_bin"]
    )
    expected_output_aged_biomass["length_bin"] = pd.Categorical(
        expected_output_aged_biomass["length_bin"],
        categories=expected_output_aged_biomass["length_bin"].unique(),
        ordered=True,
    )
    expected_output_aged_biomass["age_bin"] = pd.IntervalIndex(
        expected_output_aged_biomass["age_bin"]
    )
    expected_output_aged_biomass["age_bin"] = pd.Categorical(
        expected_output_aged_biomass["age_bin"],
        categories=expected_output_aged_biomass["age_bin"].unique(),
        ordered=True,
    )
    # ----------------------------------
    # Run tests: `distribute_aged_weight_proportions`
    # ----------------------------------
    # `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_aged_sex_biomass.shape == expected_dimensions_aged_sex_biomass
    assert eval_aged_biomass.shape == expected_dimensions_aged_biomass
    # ---- Dataframe equality
    assert np.allclose(
        eval_aged_sex_biomass[["biomass_sexed_aged_all", "biomass_sexed_aged_adult"]],
        expected_output_aged_sex_biomass[["biomass_sexed_aged_all", "biomass_sexed_aged_adult"]],
    )
    assert np.allclose(
        eval_aged_biomass[["biomass_aged_all", "biomass_aged_adult"]],
        expected_output_aged_biomass[["biomass_aged_all", "biomass_aged_adult"]],
    )
    # Assess equality between sexed and total biomass estimates
    assert eval_aged_biomass.biomass_aged_all.sum() == 6e3
    assert eval_aged_biomass.biomass_aged_adult.sum() == 3e3
    assert (
        eval_aged_biomass.biomass_aged_all.sum()
        == eval_aged_sex_biomass.biomass_sexed_aged_all.sum()
    )
    assert (
        eval_aged_biomass.biomass_aged_adult.sum()
        == eval_aged_sex_biomass.biomass_sexed_aged_adult.sum()
    )


def test_unaged_number_proportions():

    # Re-parameterize `length_df` with dummy data
    test_length_data = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4),
            "sex": np.tile(["male", "female"], 4),
            "species_id": np.repeat([19350], 8),
            "length": [12, 12, 19, 19, 12, 12, 19, 19],
            "length_count": [5, 10, 15, 20, 20, 15, 10, 5],
        }
    )

    # Length interval
    test_length_intervals = np.linspace(9, 21, 3)

    # Evaluate object for later comparison
    eval_proportions_unaged_length = unaged_number_proportions(
        test_length_data, test_length_intervals
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # `eval_aged_sex_proportions`
    # ---- Expected dimensions
    expected_dimensions = tuple([4, 4])
    # ---- Expected dataframe output
    expected_output = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 2).astype(np.int64),
            "species_id": np.repeat([19350], 4).astype(np.int64),
            "length_bin": pd.IntervalIndex.from_arrays(
                np.tile([9.0, 15.0], 2), np.tile([15.0, 21.0], 2), closed="right"
            ),
            "proportion_number_all": [0.3, 0.7, 0.7, 0.3],
        }
    )
    expected_output["length_bin"] = pd.IntervalIndex(expected_output["length_bin"])
    expected_output["length_bin"] = pd.Categorical(
        expected_output["length_bin"],
        categories=expected_output["length_bin"].unique(),
        ordered=True,
    )
    # ----------------------------------
    # Run tests: `compute_index_unaged_number_proportions`
    # ----------------------------------
    # `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_proportions_unaged_length.shape == expected_dimensions
    # ---- Dataframe equality
    assert np.all(expected_output == eval_proportions_unaged_length)


def test_unaged_weight_proportions():

    # Create mock `proportions_unaged_length`
    test_proportions_unaged_length = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 3),
            "species_id": np.repeat(19880, 6),
            "length_bin": ["small", "medium", "big", "small", "medium", "big"],
            "proportion_number_all": [0.1, 0.85, 0.05, 0.0, 0.75, 0.25],
        }
    )

    # Create mock `length_weight_df`
    test_length_weight_df = pd.DataFrame(
        {
            "length_bin": [
                "small",
                "small",
                "small",
                "medium",
                "medium",
                "medium",
                "big",
                "big",
                "big",
            ],
            "sex": ["all", "male", "female", "all", "male", "female", "all", "male", "female"],
            "weight_modeled": [
                0.5,
                0.242424,
                0.32562252,
                1.0,
                1.141414,
                0.99999999,
                2.0,
                1.942424,
                2.12313131,
            ],
        }
    )

    # Evaluate for later comparison
    eval_proportions_unaged_weight_length = unaged_weight_proportions(
        test_proportions_unaged_length, test_length_weight_df
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # `eval_proportions_unaged_weight_length`
    # ---- Expected dimensions
    expected_dimensions = tuple([6, 4])
    # ---- Expected dataframe output
    expected_output = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 3),
            "species_id": np.repeat([19880], 6),
            "length_bin": ["small", "medium", "big", "small", "medium", "big"],
            "proportion_weight_length": [0.05, 0.85, 0.10, 0.0, 0.60, 0.40],
        }
    )
    # ----------------------------------
    # Run tests: `compute_index_unaged_number_proportions`
    # ----------------------------------
    # `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_proportions_unaged_weight_length.shape == expected_dimensions
    # ---- Datatypes
    assert np.all(eval_proportions_unaged_weight_length.dtypes == expected_output.dtypes)
    # ---- Dataframe equality
    assert eval_proportions_unaged_weight_length.equals(expected_output)


def test_unaged_sex_weight_proportions():

    # Mock data for `length_data`
    test_length_data = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4),
            "sex": np.tile(["male", "female"], 4),
            "species_id": np.repeat([19350], 8),
            "length": [12, 12, 19, 19, 12, 12, 19, 19],
            "length_count": [5, 10, 15, 20, 20, 15, 10, 5],
        }
    )

    # Mock data for `length_intervals`
    test_length_intervals = np.linspace(9, 21, 3)

    # Mock data for `regression_parameters`
    test_regression_parameters = pd.DataFrame(
        {
            "sex": ["female", "male", "all"],
            "rate": [1.0, 3.0, 2.0],
            "initial": [-2.0, -4.0, -3.0],
        }
    )

    # Mock data for `aged_proportions`
    test_aged_proportions = pd.DataFrame(
        {
            "stratum_num": [0, 1],
            "weight": [100.0, 200.0],
            "weight_stratum_all": [200.00, 250.0],
            "proportion_aged_weight_all": [0.10, 0.20],
            "proportion_aged_weight_adult": [0.05, 0.10],
            "proportion_unaged_weight_all": [0.90, 0.80],
            "proportion_unaged_weight_adult": [0.45, 0.40],
        }
    )

    # Evaluate for later comparison
    eval_proportions_unaged_weight_sex = unaged_sex_weight_proportions(
        test_length_data, test_length_intervals, test_regression_parameters, test_aged_proportions
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dimensions of `obj_props_wgt_len_age_sex`
    expected_dimensions = tuple([4, 3])

    # ---- Expected dataframe output
    expected_output = pd.DataFrame(
        {
            "stratum_num": [0, 0, 1, 1],
            "sex": ["female", "male", "female", "male"],
            "proportion_weight_sex": [0.309550, 0.690450, 0.210486, 0.789514],
        }
    )

    # ----------------------------------
    # Run tests: `compute_index_aged_weight_proportions`
    # ----------------------------------
    # Process the specimen data
    # Check shape
    assert eval_proportions_unaged_weight_sex.shape == expected_dimensions

    # Check data value equality
    # ---- `stratum_num`
    assert np.all(eval_proportions_unaged_weight_sex.stratum_num == expected_output.stratum_num)


def test_calculate_unaged_biomass():

    # Mock data for `kriging_biomass_df`
    test_kriging_biomass_df = pd.DataFrame(
        {
            "centroid_latitude": [30.0, 40.0, 50.0],
            "centroid_longitude": [30.0, 31.0, 32.0],
            "fraction_cell_in_polygon": [1.0, 1.0, 1.0],
            "stratum_num": [0.0, 1.0, 2.0],
            "B_a_adult_mean": [1e3, 1e4, 1e5],
            "B_a_adult_prediction_variance": [1e1, 1e2, 1e3],
            "B_a_adult_sample_variance": [1e1, 1e2, 1e3],
            "cell_area_nmi2": [10.0, 10.0, 5.0],
            "B_adult_kriged": [1e4, 1e5, 5e5],
        }
    )

    # Mock data for `specimen_data`
    test_length_data = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1], 4),
            "sex": np.tile(["male", "female"], 4),
            "species_id": np.repeat([19350], 8),
            "length": [12, 12, 19, 19, 12, 12, 19, 19],
            "length_count": [5, 10, 15, 20, 20, 15, 10, 5],
        }
    )

    # Mock data for `regression_parameters`
    test_regression_parameters = pd.DataFrame(
        {
            "sex": ["female", "male", "all"],
            "rate": [1.0, 3.0, 2.0],
            "initial": [-2.0, -4.0, -3.0],
        }
    )

    # Mock data for `length_intervals`
    test_length_intervals = np.linspace(9, 21, 3)

    # Mock data for `length_weight_df`
    test_length_weight_df = pd.DataFrame(
        {
            "length_bin": pd.cut([12, 12, 12, 19, 19, 19], test_length_intervals),
            "sex": ["all", "male", "female", "all", "male", "female"],
            "weight_modeled": [1.0, 0.9, 1.1, 2.0, 1.9, 2.1],
        }
    )

    # Mock data for `aged_proportions`
    test_aged_proportions = pd.DataFrame(
        {
            "stratum_num": np.array([0.0, 1.0]).astype(np.int64),
            "weight": [20.0, 20.0],
            "weight_stratum_all": [100.0, 200.0],
            "proportion_aged_weight_all": [0.20, 0.10],
            "proportion_aged_weight_adult": [0.10, 0.05],
            "proportion_unaged_weight_all": [0.80, 0.90],
            "proportion_unaged_weight_adult": [0.40, 0.45],
        }
    )

    # Evaluate for later comparison
    eval_unaged_sex_biomass = calculate_unaged_biomass(
        test_kriging_biomass_df,
        test_length_data,
        test_length_intervals,
        test_length_weight_df,
        test_regression_parameters,
        test_aged_proportions,
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # `eval_aged_sex_biomass`
    # ---- Expected dimensions
    expected_dimensions = tuple([4, 5])
    # ---- Expected dataframe output
    expected_output = pd.DataFrame(
        {
            "species_id": np.repeat(19350, 4).astype(np.int64),
            "sex": np.repeat(["female", "male"], 2),
            "length_bin": pd.cut([12, 19, 12, 19], np.linspace(9, 21, 3)),
            "biomass_sexed_unaged_all": [10637.488113, 10782.651554, 39235.815054, 37344.045279],
            "biomass_sexed_unaged_adult": [5318.744056, 5391.325777, 19617.907527, 18672.022639],
        }
    )
    expected_output["length_bin"] = pd.IntervalIndex(expected_output["length_bin"])
    expected_output["length_bin"] = pd.Categorical(
        expected_output["length_bin"],
        categories=expected_output["length_bin"].unique(),
        ordered=True,
    )

    # ----------------------------------
    # Run tests: `distribute_aged_weight_proportions`
    # ----------------------------------
    # `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_unaged_sex_biomass.shape == expected_dimensions
    # ---- Dataframe equality
    assert np.allclose(
        eval_unaged_sex_biomass[["biomass_sexed_unaged_all", "biomass_sexed_unaged_adult"]],
        expected_output[["biomass_sexed_unaged_all", "biomass_sexed_unaged_adult"]],
    )
    # ---- Check expected summed values
    assert np.isclose(eval_unaged_sex_biomass["biomass_sexed_unaged_all"].sum(), 98e3)
    assert np.isclose(eval_unaged_sex_biomass["biomass_sexed_unaged_adult"].sum(), 49e3)


def test_apply_age_bins():

    # Mock data for `aged_sex_biomass`
    test_aged_sex_biomass = pd.DataFrame(
        {
            "species_id": np.repeat(19350, 8).astype(np.int64),
            "sex": np.repeat(["female", "male"], 4),
            "length_bin": pd.cut(np.tile([12, 12, 19, 19], 2), np.linspace(9, 21, 3)),
            "age_bin": pd.cut(np.tile([1, 2], 4), np.array([0.5, 1.5, 2.5])),
            "biomass_sexed_aged_all": [
                1000.0,
                2000.0,
                3000.0,
                5000.0,
                4000.0,
                3000.0,
                2000.0,
                200.0,
            ],
            "biomass_sexed_aged_adult": [0.0, 500.0, 0.0, 1250.0, 0.0, 750.0, 0.0, 50.0],
        }
    )
    test_aged_sex_biomass["length_bin"] = pd.IntervalIndex(test_aged_sex_biomass["length_bin"])
    test_aged_sex_biomass["length_bin"] = pd.Categorical(
        test_aged_sex_biomass["length_bin"],
        categories=test_aged_sex_biomass["length_bin"].unique(),
        ordered=True,
    )
    test_aged_sex_biomass["age_bin"] = pd.IntervalIndex(test_aged_sex_biomass["age_bin"])
    test_aged_sex_biomass["age_bin"] = pd.Categorical(
        test_aged_sex_biomass["age_bin"],
        categories=test_aged_sex_biomass["age_bin"].unique(),
        ordered=True,
    )

    # Mock data for `unaged_sex_biomass`
    test_unaged_sex_biomass = pd.DataFrame(
        {
            "species_id": np.repeat(19350, 4).astype(np.int64),
            "sex": np.repeat(["female", "male"], 2),
            "length_bin": pd.cut([12, 19, 12, 19], np.linspace(9, 21, 3)),
            "biomass_sexed_unaged_all": [5000.0, 2000.0, 4000.0, 0.0],
            "biomass_sexed_unaged_adult": [1250.0, 500.0, 1000.0, 0.0],
        }
    )
    test_unaged_sex_biomass["length_bin"] = pd.IntervalIndex(test_unaged_sex_biomass["length_bin"])
    test_unaged_sex_biomass["length_bin"] = pd.Categorical(
        test_unaged_sex_biomass["length_bin"],
        categories=test_unaged_sex_biomass["length_bin"].unique(),
        ordered=True,
    )

    # Evaluate for later comparison
    eval_redistrib_unaged_sex_biomass, eval_redistrib_unaged_biomass = apply_age_bins(
        test_aged_sex_biomass, test_unaged_sex_biomass
    )

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- `eval_redistrib_unaged_sex_biomass`
    # ---- Expected dimensions of `obj_props_wgt_len_age_sex`
    expected_dimensions_redistrib_unaged_sex_biomass = tuple([8, 6])
    # ---- Expected dataframe output
    expected_output_redistrib_unaged_sex_biomass = pd.DataFrame(
        {
            "species_id": np.repeat(19350, 8).astype(np.int64),
            "sex": np.repeat(["female", "male"], 4),
            "length_bin": pd.cut(np.tile([12, 12, 19, 19], 2), np.linspace(9, 21, 3)),
            "age_bin": pd.cut(np.tile([1, 2], 4), np.array([0.5, 1.5, 2.5])),
            "biomass_sexed_unaged_all": [
                1666.666668,
                3333.333333,
                750.000000,
                1250.000000,
                2285.714286,
                1714.285714,
                0.000000,
                0.000000,
            ],
            "biomass_sexed_unaged_adult": [0.0, 1250.0, 0.0, 500.0, 0.0, 1000.0, 0.0, 0.0],
        }
    )
    expected_output_redistrib_unaged_sex_biomass["length_bin"] = pd.IntervalIndex(
        expected_output_redistrib_unaged_sex_biomass["length_bin"]
    )
    expected_output_redistrib_unaged_sex_biomass["length_bin"] = pd.Categorical(
        expected_output_redistrib_unaged_sex_biomass["length_bin"],
        categories=expected_output_redistrib_unaged_sex_biomass["length_bin"].unique(),
        ordered=True,
    )
    expected_output_redistrib_unaged_sex_biomass["age_bin"] = pd.IntervalIndex(
        expected_output_redistrib_unaged_sex_biomass["age_bin"]
    )
    expected_output_redistrib_unaged_sex_biomass["age_bin"] = pd.Categorical(
        expected_output_redistrib_unaged_sex_biomass["age_bin"],
        categories=expected_output_redistrib_unaged_sex_biomass["age_bin"].unique(),
        ordered=True,
    )
    # ---- `eval_redistrib_unaged_biomass`
    # ---- Expected dimensions of `obj_props_wgt_len_age_sex`
    expected_dimensions_redistrib_unaged_biomass = tuple([4, 5])
    # ---- Expected dataframe output
    expected_output_redistrib_unaged_biomass = pd.DataFrame(
        {
            "length_bin": pd.cut(np.repeat([12, 19], 2), np.linspace(9, 21, 3)),
            "age_bin": pd.cut(np.tile([1, 2], 2), np.array([0.5, 1.5, 2.5])),
            "species_id": np.repeat(19350, 4).astype(np.int64),
            "biomass_unaged_all": [3952.4, 5047.6, 750.0, 1250.0],
            "biomass_unaged_adult": [0.0, 2250.0, 0.0, 500.0],
        }
    )
    expected_output_redistrib_unaged_biomass["length_bin"] = pd.IntervalIndex(
        expected_output_redistrib_unaged_biomass["length_bin"]
    )
    expected_output_redistrib_unaged_biomass["length_bin"] = pd.Categorical(
        expected_output_redistrib_unaged_biomass["length_bin"],
        categories=expected_output_redistrib_unaged_biomass["length_bin"].unique(),
        ordered=True,
    )
    expected_output_redistrib_unaged_biomass["age_bin"] = pd.IntervalIndex(
        expected_output_redistrib_unaged_biomass["age_bin"]
    )
    expected_output_redistrib_unaged_biomass["age_bin"] = pd.Categorical(
        expected_output_redistrib_unaged_biomass["age_bin"],
        categories=expected_output_redistrib_unaged_biomass["age_bin"].unique(),
        ordered=True,
    )
    # ----------------------------------
    # Run tests: `apply_age_bins`
    # ----------------------------------
    # Check shape
    assert (
        eval_redistrib_unaged_sex_biomass.shape == expected_dimensions_redistrib_unaged_sex_biomass
    )
    assert eval_redistrib_unaged_biomass.shape == expected_dimensions_redistrib_unaged_biomass
    # Check datatypes
    assert np.all(
        eval_redistrib_unaged_sex_biomass.dtypes
        == expected_output_redistrib_unaged_sex_biomass.dtypes
    )
    assert np.all(
        eval_redistrib_unaged_biomass.dtypes == expected_output_redistrib_unaged_biomass.dtypes
    )
    # ---- Dataframe equality
    assert np.allclose(
        eval_redistrib_unaged_sex_biomass[
            ["biomass_sexed_unaged_all", "biomass_sexed_unaged_adult"]
        ],
        expected_output_redistrib_unaged_sex_biomass[
            ["biomass_sexed_unaged_all", "biomass_sexed_unaged_adult"]
        ],
    )
    assert np.allclose(
        eval_redistrib_unaged_biomass[["biomass_unaged_all", "biomass_unaged_adult"]],
        expected_output_redistrib_unaged_biomass[["biomass_unaged_all", "biomass_unaged_adult"]],
    )
    # ---- Check expected summed values
    assert eval_redistrib_unaged_biomass["biomass_unaged_all"].sum() == 11e3
    assert eval_redistrib_unaged_biomass["biomass_unaged_adult"].sum() == 2.75e3
    # ---- Evaluate the sum and comparison between 'total' and 'sexed'
    assert np.isclose(
        eval_redistrib_unaged_sex_biomass["biomass_sexed_unaged_all"].sum(),
        eval_redistrib_unaged_biomass["biomass_unaged_all"].sum(),
    )
    assert np.isclose(
        eval_redistrib_unaged_sex_biomass["biomass_sexed_unaged_adult"].sum(),
        eval_redistrib_unaged_biomass["biomass_unaged_adult"].sum(),
    )
