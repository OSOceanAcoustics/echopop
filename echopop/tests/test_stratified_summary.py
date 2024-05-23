import numpy as np
import pandas as pd

from echopop.computation.statistics import confidence_interval, stratified_transect_statistic
from echopop.tests.conftest import assert_dictionary_equal


def test_stratified_transect_statistic():

    # Create mock data for `transect_summary`
    test_transect_summary = pd.DataFrame(
        {
            "transect_num": [1, 2, 3, 4],
            "minimum_longitude": [-5.0, -3.0, -1.0, 1.0],
            "maxmum_longitude": [-2.0, 5.0, 3.0, 7.0],
            "center_latitude": [10.0, 11.0, 12.5, 13.5],
            "transect_distance": [177.600950, 472.070493, 234.766275, 350.736855],
            "transect_spacing": [2.0, 2.0, 2.0, 2.0],
            "transect_area": [355.201900, 944.140986, 469.532550, 701.473710],
            "B_adult": [1e2, 1e3, 1e5, 1e4],
            "stratum_inpfc": [1, 1, 2, 2],
        },
    )

    # Create mock data for `strata_summary`
    test_strata_summary = pd.DataFrame(
        {
            "stratum_inpfc": [1, 2],
            "num_transects": [2, 2],
            "total_transect_area": [1299.342886, 1171.006260],
        },
    )

    # Evaluate for later comparison
    # ---- Replicates == 1
    # ---- Transect sample proportion == 100%
    test_transect_sample = 1.0
    test_transect_replicates = 1
    eval_single_stratified_results = stratified_transect_statistic(
        test_transect_summary,
        test_strata_summary,
        test_transect_sample,
        test_transect_replicates,
        parameter="B_adult",
    )
    # ---- Replicates == 10
    # ---- Transect sample proportion == 100%
    test_transect_sample = 1.0
    test_transect_replicates = 10
    eval_single_rep_stratified_results = stratified_transect_statistic(
        test_transect_summary,
        test_strata_summary,
        test_transect_sample,
        test_transect_replicates,
        parameter="B_adult",
    )

    # ---- Replicates == 1
    # ---- Transect sample proportion == 50%
    test_transect_sample = 0.5
    test_transect_replicates = 1
    np.random.seed(10)
    eval_single_sub_stratified_results = stratified_transect_statistic(
        test_transect_summary,
        test_strata_summary,
        test_transect_sample,
        test_transect_replicates,
        parameter="B_adult",
    )

    # ---- Replicates == 1
    # ---- Transect sample proportion == 50%
    test_transect_sample = 0.5
    test_transect_replicates = 10
    np.random.seed(1800)
    eval_single_sub_rep_stratified_results = stratified_transect_statistic(
        test_transect_summary,
        test_strata_summary,
        test_transect_sample,
        test_transect_replicates,
        parameter="B_adult",
    )

    # ++++ Bundle!
    eval_dictionary = {
        "single": eval_single_stratified_results,
        "single_rep": eval_single_rep_stratified_results,
        "single_sub": eval_single_sub_stratified_results,
        "single_rep_sub": eval_single_sub_rep_stratified_results,
    }

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dtypes
    expected_dtypes = {
        "single": {
            "biomass": {
                "mean": {
                    "estimate": np.floating,
                    "confidence_interval": np.array([54947653.27600001, 54947653.27600001]),
                },
                "variance": {
                    "estimate": 54846534.456292756,
                    "confidence_interval": np.array([54846534.45629276, 54846534.45629276]),
                },
                "CV": {
                    "estimate": 0.9981597245072626,
                    "confidence_interval": np.array([0.99815972, 0.99815972]),
                },
            },
        },
        "single_rep": {
            "biomass": {
                "mean": {
                    "estimate": 54947653.27600001,
                    "confidence_interval": np.array([54947653.27600001, 54947653.27600001]),
                },
                "variance": {
                    "estimate": 54846534.45629276,
                    "confidence_interval": np.array([54846534.45629275, 54846534.45629278]),
                },
                "CV": {
                    "estimate": 0.9981597245072626,
                    "confidence_interval": np.array([0.99815972, 0.99815972]),
                },
            },
        },
        "single_sub": {
            "biomass": {
                "mean": {
                    "estimate": 117230560.28860001,
                    "confidence_interval": np.array([1.1723056e08, 1.1723056e8]),
                },
                "variance": {
                    "estimate": 116601900.95605445,
                    "confidence_interval": np.array([1.16601901e8, 1.16601901e8]),
                },
                "CV": {
                    "estimate": 0.994637410833848,
                    "confidence_interval": np.array([0.99463741, 0.99463741]),
                },
            },
        },
        "single_rep_sub": {
            "biomass": {
                "mean": {
                    "estimate": 54463985.68756001,
                    "confidence_interval": np.array([-4.69233576e7, 1.55851329e8]),
                },
                "variance": {
                    "estimate": 53662832.43264915,
                    "confidence_interval": np.array([-4.70645276e7, 1.54390192e8]),
                },
                "CV": {
                    "estimate": 0.9710233886235905,
                    "confidence_interval": np.array([0.90408889, 1.03795788]),
                },
            },
        },
    }
    # ---- Expected output
    expected_output = {
        "single": {
            "biomass": {
                "mean": {
                    "estimate": 54947653.27600001,
                    "confidence_interval": np.array([54947653.27600001, 54947653.27600001]),
                },
                "variance": {
                    "estimate": 54846534.456292756,
                    "confidence_interval": np.array([54846534.45629276, 54846534.45629276]),
                },
                "CV": {
                    "estimate": 0.9981597245072626,
                    "confidence_interval": np.array([0.99815972, 0.99815972]),
                },
            },
        },
        "single_rep": {
            "biomass": {
                "mean": {
                    "estimate": 54947653.27600001,
                    "confidence_interval": np.array([54947653.27600001, 54947653.27600001]),
                },
                "variance": {
                    "estimate": 54846534.45629276,
                    "confidence_interval": np.array([54846534.45629275, 54846534.45629278]),
                },
                "CV": {
                    "estimate": 0.9981597245072626,
                    "confidence_interval": np.array([0.99815972, 0.99815972]),
                },
            },
        },
        "single_sub": {
            "biomass": {
                "mean": {
                    "estimate": 117230560.28860001,
                    "confidence_interval": np.array([1.1723056e08, 1.1723056e8]),
                },
                "variance": {
                    "estimate": 116601900.95605445,
                    "confidence_interval": np.array([1.16601901e8, 1.16601901e8]),
                },
                "CV": {
                    "estimate": 0.994637410833848,
                    "confidence_interval": np.array([0.99463741, 0.99463741]),
                },
            },
        },
        "single_rep_sub": {
            "biomass": {
                "mean": {
                    "estimate": 54463985.68756001,
                    "confidence_interval": np.array([-4.69233576e7, 1.55851329e8]),
                },
                "variance": {
                    "estimate": 53662832.43264915,
                    "confidence_interval": np.array([-4.70645276e7, 1.54390192e8]),
                },
                "CV": {
                    "estimate": 0.9710233886235905,
                    "confidence_interval": np.array([0.90408889, 1.03795788]),
                },
            },
        },
    }

    # ----------------------------------
    # Run tests: `stratified_transect_statistic`
    # ----------------------------------
    assert_dictionary_equal(eval_dictionary, expected_dtypes, expected_output)


def test_confidence_interval():

    # Mock values
    test_values = [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]

    # Evaluate for comparison later
    eval_ci_values = confidence_interval(test_values)

    # --------------------------------
    # Expected outcomes
    # --------------------------------
    # ---- Expected dimensions
    expected_dimensions = tuple(
        [
            2,
        ]
    )
    # --- Expected dtype
    # ---- Expected output
    expected_output = np.array([0.20104371, 5.35451185])

    # ----------------------------------
    # Run tests: `confidence_interval`
    # ----------------------------------
    # Check shape
    assert eval_ci_values.shape == expected_dimensions
    # Check dtype
    assert np.issubdtype(eval_ci_values.dtype, np.floating)
    # Check output
    assert np.allclose(eval_ci_values, expected_output)
