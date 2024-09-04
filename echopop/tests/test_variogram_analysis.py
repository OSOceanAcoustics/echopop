import re

import numpy as np
import pandas as pd
import pytest

import echopop.spatial.variogram as esv
from echopop.spatial.variogram import (
    VARIOGRAM_MODELS,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
    prepare_variogram_matrices,
    quantize_lags,
    semivariance,
    variogram_matrix_filter,
)
from echopop.tests.conftest import assert_dictionary_values_equal, load_json_data

___EXPECTED_OUTCOME_FILENAME__ = "variogram_analysis.json"


@pytest.mark.parametrize(
    "model, expected, exception",
    [
        (
            "exponential",
            ({"distance_lags", "sill", "nugget", "correlation_range"}, esv.exponential),
            None,
        ),
        (
            "gaussian",
            ({"distance_lags", "sill", "nugget", "correlation_range"}, esv.gaussian),
            None,
        ),
        ("jbessel", ({"distance_lags", "sill", "nugget", "hole_effect_range"}, esv.jbessel), None),
        ("kbessel", ({"distance_lags", "sill", "nugget", "hole_effect_range"}, esv.kbessel), None),
        ("linear", ({"distance_lags", "sill", "nugget"}, esv.linear), None),
        ("nugget", ({"distance_lags", "sill", "nugget"}, esv.nugget), None),
        ("sinc", ({"distance_lags", "sill", "nugget", "hole_effect_range"}, esv.sinc), None),
        (
            "spherical",
            ({"distance_lags", "sill", "nugget", "correlation_range"}, esv.spherical),
            None,
        ),
        (
            ("bessel", "exponential"),
            (
                {
                    "distance_lags",
                    "sill",
                    "nugget",
                    "correlation_range",
                    "hole_effect_range",
                    "decay_power",
                },
                esv.bessel_exponential,
            ),
            None,
        ),
        (
            ("bessel", "gaussian"),
            (
                {"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"},
                esv.bessel_gaussian,
            ),
            None,
        ),
        (
            ("cosine", "exponential"),
            (
                {
                    "distance_lags",
                    "sill",
                    "nugget",
                    "correlation_range",
                    "hole_effect_range",
                    "enhance_semivariance",
                },
                esv.cosine_exponential,
            ),
            None,
        ),
        (
            ("cosine", "gaussian"),
            (
                {"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"},
                esv.cosine_gaussian,
            ),
            None,
        ),
        (
            ("exponential", "linear"),
            (
                {
                    "distance_lags",
                    "sill",
                    "nugget",
                    "correlation_range",
                    "hole_effect_range",
                    "decay_power",
                },
                esv.exponential_linear,
            ),
            None,
        ),
        (
            ("gaussian", "linear"),
            (
                {"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"},
                esv.gaussian_linear,
            ),
            None,
        ),
        ("invalid", (None, None), "could not be matched to an existing variogram method."),
        (
            ("invalid", "tuple"),
            (None, None),
            "could not be matched to an existing variogram method.",
        ),
        (
            ["invalid", "list"],
            (None, None),
            "could not be matched to an existing variogram method.",
        ),
    ],
    ids=[
        "Exponential arguments",
        "Gaussian arguments",
        "J-Bessel arguments",
        "K-Bessel arguments",
        "Linear arguments",
        "Nugget arguments",
        "Sinc arguments",
        "Spherical arguments",
        "Bessel-Exponential arguments",
        "Bessel-Gaussian arguments",
        "Cosine-Exponential arguments",
        "Cosine-Gaussian arguments",
        "Exponential-Linear arguments",
        "Gaussian-Linear arguments",
        "Unknown single-family variogram model (str)",
        "Unknown composite variogram model (tuple)",
        "Invalid composite variogram model (list)",
    ],
)
def test_get_variogram_arguments(model, expected, exception):

    if exception is not None:
        with pytest.raises(LookupError, match=re.escape(exception)):
            assert esv.get_variogram_arguments(model)
    else:
        # Get arguments
        args, func = esv.get_variogram_arguments(model)
        # Assert argument equality
        assert set(args.keys()) == set(expected[0])
        # Assert function equality
        assert func["model_function"] == expected[1]


@pytest.mark.parametrize(
    "description",
    ["Test `initialize_optimization_config` translation"],
    ids=["Test `initialize_optimization_config` translation"],
)
def test_initialize_optimization_config(description):

    # -------------------------
    # Mock data [ TEST 1 ]: ASSESS CORRECT `lmfit` PARAMETER/ARGUMENT CONVERSION
    # ---- Create input dictionary
    ECHOPOP_OPTIMIZATION_INPUTS = {
        "max_fun_evaluations": 500,
        "cost_fun_tolerance": 1e-8,
        "solution_tolerance": 1e-8,
        "gradient_tolerance": 1e-8,
        "finite_step_size": 1e-8,
        "trust_region_solver": "exact",
        "x_scale": "jacobian",
        "jacobian_approx": "forward",
    }

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- Create output
    output = esv.initialize_optimization_config(ECHOPOP_OPTIMIZATION_INPUTS)
    # ---- Create translation dictionary
    TRANSLATION_DICT = {
        "max_fun_evaluations": {"name": "max_nfev", "value": 500},
        "cost_fun_tolerance": {"name": "ftol", "value": 1e-8},
        "solution_tolerance": {"name": "xtol", "value": 1e-8},
        "gradient_tolerance": {"name": "gtol", "value": 1e-8},
        "finite_step_size": {"name": "diff_step", "value": 1e-8},
        "trust_region_solver": {"name": "tr_solver", "value": "exact"},
        "x_scale": {"name": "x_scale", "value": "jac"},
        "jacobian_approx": {"name": "jac", "value": "2-point"},
    }
    # ---- Iterate through keys to assess equality across names [ EQUIVALENT KEYS ]
    assert all(
        [
            TRANSLATION_DICT[key]["name"] in output.keys()
            for key in ECHOPOP_OPTIMIZATION_INPUTS.keys()
        ]
    )
    # ---- Iterate through keys to assess equality across values [ EQUIVALENT VALUES ]
    assert all(
        [
            output[TRANSLATION_DICT[key]["name"]] == TRANSLATION_DICT[key]["value"]
            for key in ECHOPOP_OPTIMIZATION_INPUTS.keys()
        ]
    )


@pytest.fixture()
def transect_data():
    # Set a random seed for reproducibility
    np.random.seed(99)

    # Generate 20 random rows of data
    num_rows = 5
    data = {
        "transect_num": np.random.choice([1, 2], size=num_rows),
        "longitude": np.random.uniform(-121.5, -120.5, size=num_rows),
        "latitude": np.random.uniform(34.0, 35.0, size=num_rows),
        "stratum_num": np.random.choice([1, 2], size=num_rows),
        "transect_spacing": np.random.uniform(5.0, 15.0, size=num_rows),
        "biomass": np.random.uniform(0.0, 10.0, size=num_rows),
        "biomass_density": np.random.uniform(0.0, 2.0, size=num_rows),
        "x": np.random.uniform(-0.5, 0.5, size=num_rows),
        "y": np.random.uniform(-0.5, 0.5, size=num_rows),
    }

    return pd.DataFrame(data)


@pytest.mark.parametrize(
    "description",
    ["Test `prepare_variogram_matrices` azimuth and lag matrix outputs"],
    ids=["Test `prepare_variogram_matrices` azimuth and lag matrix outputs"],
)
def test_prepare_variogram_matrices(test_path, transect_data, description):

    # -------------------------
    # Load the expected outcomes
    expected_results = load_json_data(
        str(test_path["EXPECTED"] / ___EXPECTED_OUTCOME_FILENAME__),
        test_name="prepare_variogram_matrices",
    )

    # -------------------------
    # Mock additional parameters
    MOCK_LAG_RESOLUTION = 0.002

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- Run the function
    azimuth_matrix, lag_matrix = prepare_variogram_matrices(transect_data, MOCK_LAG_RESOLUTION)
    # ---- ASSERT azimuth_matrix equality
    assert np.array_equal(azimuth_matrix, expected_results["azimuth_matrix"], equal_nan=True)
    # ---- ASSERT azimuth_matrix equality
    assert np.array_equal(lag_matrix, expected_results["lag_matrix"])


@pytest.mark.parametrize(
    "data_matrix, mask_matrix, azimuth_matrix, azimuth_range, expected_output, exception",
    [
        # Test 2D matrix with valid mask and azimuth
        (
            np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
            0.5,
            np.array([1, 3, 5, 9]),
            None,
        ),
        # Test 1D data matrix broadcasting
        (
            np.array([1, 2, 3]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[0.1, -0.3, 0.7], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
            0.5,
            np.array([1, 2, 3]),
            None,
        ),
        # Test 2D matrix with NaN azimuths
        (
            np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[np.nan, -0.3, np.nan], [0.4, np.nan, 0.5], [-0.4, 0.3, np.nan]]),
            0.5,
            np.array([1, 3, 5, 9]),
            None,
        ),
        # Test 2D matrix with 0.0 azimuths
        (
            np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[0.0, -0.3, 0.0], [0.4, 0.0, 0.5], [-0.4, 0.3, 0.0]]),
            0.5,
            np.array([1, 3, 5, 9]),
            None,
        ),
        # Test dimension mismatch handling
        (
            np.array([[1, 2], [3, 4]]),
            np.array([[True, False], [False, True], [True, False]]),
            np.array([[0.1, -0.3], [0.4, -0.1], [0.0, 0.3]]),
            0.5,
            np.array([1]),
            # r"boolean index did not match indexed array"
            "boolean index did not match indexed array",
        ),
        # Test no matching values
        (
            np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
            np.array([[False, False, False], [False, False, False], [False, False, False]]),
            np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
            0.5,
            np.array([]),
            None,
        ),
    ],
    ids=[
        "2D matrix with valid mask and azimuth",
        "1D data matrix broadcasting",
        "2D matrix with NaN azimuths",
        "2D matrix with 0.0 azimuths",
        "Dimension mismatch handling",
        "No matching values",
    ],
)
def test_variogram_matrix_filter(
    data_matrix, mask_matrix, azimuth_matrix, azimuth_range, expected_output, exception
):
    if exception is not None:
        with pytest.raises(IndexError, match=re.escape(exception)):
            assert variogram_matrix_filter(data_matrix, mask_matrix, azimuth_matrix, azimuth_range)
    else:
        # Copmute the azimuth-filtered matrix
        output = variogram_matrix_filter(data_matrix, mask_matrix, azimuth_matrix, azimuth_range)
        # Assert argument equality
        assert np.array_equal(output, expected_output)


@pytest.fixture
def large_quantize_lags() -> tuple:

    # Set a random seed for reproducibility
    np.random.seed(99)

    # Set number of rows
    num_rows = int(1e3)

    # Simulate large number of rows
    (estimates, lag_matrix, mask_matrix, azimuth_matrix, azimuth_range, n_lags) = (
        np.random.rand(num_rows),
        np.random.randint(0, 10, size=(num_rows, num_rows)),
        np.random.choice([True, False], size=(num_rows, num_rows)),
        np.random.uniform(0, 360, size=(num_rows, num_rows)),
        180,
        3,
    )

    return estimates, lag_matrix, mask_matrix, azimuth_matrix, azimuth_range, n_lags


@pytest.fixture
def expected_quantize_lags(test_path):

    # -------------------------
    # Load the expected outcomes
    expected_results = load_json_data(
        str(test_path["EXPECTED"] / ___EXPECTED_OUTCOME_FILENAME__), test_name="quantize_lags"
    )

    return expected_results


@pytest.mark.parametrize(
    "estimates, lag_matrix, mask_matrix, azimuth_matrix, azimuth_range, n_lags, kwargs, json_key",
    [
        (
            np.array([1.0, 2.0, 3.0]),
            np.array([[1, 2, 3], [2, 3, 1], [3, 1, 2]]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
            0.5,
            3,
            {},
            "basic_inputs",
        ),
        (
            np.array([1.0, 2.0, 3.0]),
            np.array([[1, 2, 3], [2, 3, 1], [3, 1, 2]]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
            0.5,
            3,
            {"dummy": "variable"},
            "basic_inputs_plus_kwargs",
        ),
        (
            np.array([1.0, 1.0, 1.0, 1.0, 1.0]),
            np.array([0, 1, 2, 3, 4]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
            1.0,
            5,
            {},
            "invalid_1d_array",
        ),
        (
            np.array([1.0, 1.0, 1.0, 1.0, 1.0]),
            np.array([0, 1, 2, 3, 4]),
            np.array([True, True, True, True, True]),
            np.array([0.1, 0.2, 0.3, 0.4, 0.5]),
            1.0,
            5,
            {},
            "invalid_1d_arrays",
        ),
        (
            np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]),
            np.array([[1, 2, 3], [2, 3, 1], [3, 1, 2]]),
            np.array([[True, False, True], [False, True, False], [True, False, True]]),
            np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
            0.5,
            3,
            {},
            "invalid_estimates_input",
        ),
        (None, None, None, None, None, None, {}, "large_sample_size"),
    ],
    ids=[
        "Basic valid inputs",
        "Superfluous kwargs args",
        "Invalid 1D array",
        "Multiple invalid 1D arrays",
        "Invalid 2D array for estimates",
        "Large data input",
    ],
)
def test_quantize_lags(
    large_quantize_lags,
    expected_quantize_lags,
    estimates,
    lag_matrix,
    mask_matrix,
    azimuth_matrix,
    azimuth_range,
    n_lags,
    kwargs,
    json_key,
):

    # Get the key args
    expected_results = expected_quantize_lags[json_key]

    # Test for exceptions
    if expected_results["exceptions"] is not None:
        with pytest.raises(ValueError, match=re.escape(expected_results["exceptions"])):
            assert quantize_lags(
                estimates, lag_matrix, mask_matrix, azimuth_matrix, azimuth_range, n_lags, **kwargs
            )

    else:
        # Generate output from fixture if indicated
        if "fixture" in expected_results and expected_results["fixture"]:
            output = quantize_lags(*large_quantize_lags)
        else:
            # Copmute the azimuth-filtered matrix from the parameterized inputs
            output = quantize_lags(
                estimates, lag_matrix, mask_matrix, azimuth_matrix, azimuth_range, n_lags, **kwargs
            )

        # Convert `expected_results["expected"]` to a tuple
        expected_tuple = tuple(
            expected_results["expected"][key] for key in expected_results["expected"]
        )

        # Assert same shapes
        assert all([e.shape == o.shape for (e, o) in zip(expected_tuple, output)])

        # Assert value equality
        assert all([np.allclose(e, o) for (e, o) in zip(expected_tuple, output)])


# quantize_lags(estimates, lag_matrix, mask_matrix, azimuth_matrix, azimuth_range, n_lags)
# quantize_lags(*large_quantize_lags())


# new_entries = {
#     "basic_inputs": {
#         "expected": {
#             "lag_counts": np.array([1, 1]),
#             "lag_estimates": np.array([1.0, 3.0]),
#             "lag_estimates_squared": np.array([1.0, 9.0]),
#             "lag_deviations": np.array([0.0, 0.0])
#         },
#         "exceptions": None,
#         "id": "Basic valid inputs"
#     },
#     "basic_inputs_plus_kwargs": {
#         "expected": {
#             "lag_counts": np.array([1, 1]),
#             "lag_estimates": np.array([1.0, 3.0]),
#             "lag_estimates_squared": np.array([1.0, 9.0]),
#             "lag_deviations": np.array([0.0, 0.0])
#         },
#         "exceptions": None,
#         "id": "Superfluous kwargs args"
#     },
#     "invalid_1d_array": {
#         "expected": None,
#         "exceptions": "requires arrays to be 2D",
#         "id": "Invalid 1D array",
#     },
#     "invalid_1d_arrays": {
#         "expected": None,
#         "exceptions": "requires arrays to be 2D",
#         "id": "Multiple invalid 1D array",
#     },
#     "invalid_estimates_input": {
#         "expected": None,
#         "exceptions": "must be a 1D array",
#         "id": "Invalid 2D array for estimates",
#     },
#     "large_sample_size": {
#         "expected": {
#             "lag_counts": np.array([12635, 12519]),
#             "lag_estimates": np.array([6582.74123055, 6538.70487521]),
#             "lag_estimates_squared": np.array([4476.91982075, 4469.42487519]),
#             "lag_deviations": np.array([2088.67551011, 2026.94478262])
#         },
#         "exceptions": None,
#         "id": "Large data input",
#         "fixture": True,
#     },
#     # "is_numpy": ["lag_counts", "lag_estimates", "lag_estimates_squared", "lag_deviations"],
#     # "is_pandas": []
# }
# from pathlib import Path
# TEST_DATA_ROOT = Path("C:/Users/Brandyn/Documents/GitHub/echopop/echopop/test_data")
# test_path = {
#         "ROOT": TEST_DATA_ROOT,
#         "CONFIG": TEST_DATA_ROOT / "config_files",
#         "INPUT": TEST_DATA_ROOT / "input_files",
#         "EXPECTED": TEST_DATA_ROOT / "expected_outputs"
#     }

# filename = test_path["EXPECTED"] / ___EXPECTED_OUTCOME_FILENAME__
# mock_survey = Survey(
#         init_config_path=Path(test_path["CONFIG"] / "config_init.yml"),
#         survey_year_config_path=Path(test_path["CONFIG"] / "config_survey.yml"),
#     )

# mock_survey.load_acoustic_data()
# mock_survey.input
# # test_name = "quantize_lags"

# expected_counts = np.array([1, 1])
# expected_estimates = np.array([1.0, 3.0])
# expected_estimates_squared = np.array([1.0, 9.0])
# expected_deviations = np.array([0.0, 0.0])
# [

#     (
#         np.array([1.0, 2.0, 3.0]),
#         np.array([[1, 2, 3], [2, 3, 1], [3, 1, 2]]),
#         np.array([[True, False, True], [False, True, False], [True, False, True]]),
#         np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
#         0.5,
#         3
#     )
# ]


# @pytest.mark.parametrize(
#     "data_matrix, mask_matrix, azimuth_matrix, azimuth_range",
# )

# {
#     "variogram_matrix_filter": {
#         "INPUTS": {
#             "valid_mask_azimuth": {
#                 "data_matrix": np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
#                 "mask_matrix": np.array([[True, False, True], [False, True, False],
#                                          [True, False, True]]),
#                "azimuth_matrix": np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
#                 "azimuth_range": 0.5
#             },
#             "matrix_broadcasting": {
#                 "data_matrix": np.array([1, 2, 3]),
#                 "mask_matrix": np.array([[True, False, True], [False, True, False],
#                                          [True, False, True]]),
#                "azimuth_matrix": np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
#                 "azimuth_range": 0.5
#             },
#             "nan_azimuths": {
#                 "data_matrix": np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
#                 "mask_matrix": np.array([[True, False, True], [False, True, False],
#                                          [True, False, True]]),
#                 "azimuth_matrix": np.array([[np.nan, -0.3, np.nan], [0.4, np.nan, 0.5],
#                                             [-0.4, 0.3, np.nan]]),
#                 "azimuth_range": 0.5
#             },
#             "dimension_mismatch": {
#                 "data_matrix": np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
#                 "mask_matrix": np.array([[True, False, True], [False, True, False],
#                                          [True, False, True]]),
#                "azimuth_matrix": np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
#                 "azimuth_range": 0.5
#             },
#         },
#         "OUTPUTS": {
#             "valid_mask_azimuth": np.array([1, 9]),
#             "matrix_broadcasting": np.array([1, 3]),
#             "nan_azimuths": np.array([1, 9]),
#             "dimension_mismatch": np.array([1]),
#             "no_matching_values": np.array([]),
#         }
#     }
# }


# # Test 2D matrix with NaN azimuths

# # Test dimension mismatch handling
# (
#     np.array([[1, 2], [3, 4]]),
#     np.array([[True, False], [False, True], [True, False]]),
#     np.array([[0.1, -0.3], [0.4, -0.1], [0.0, 0.3]]),
#     0.5,

# ),
# # Test no matching values
# (
#     np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
#     np.array([[False, False, False], [False, False, False], [False, False, False]]),
#     np.array([[0.1, -0.3, 0.2], [0.4, -0.1, 0.5], [-0.4, 0.3, 0.0]]),
#     0.5,
# ),
# ]

# from echopop.survey import Survey
# import copy
# from pathlib import Path
# from typing import List, Literal, Optional, Union

# from echopop.analysis import (
#     acoustics_to_biology,
#     apportion_kriged_values,
#     krige,
#     process_transect_data,
#     stratified_summary,
#     variogram_analysis,
# )
# from echopop.utils.validate import VariogramBase, VariogramInitial, VariogramOptimize

# init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
# file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
# survey_data = Survey(init_config, file_config)
# survey_data.load_acoustic_data()
# survey_data.load_survey_data()
# survey_data.transect_analysis()
# survey_data.fit_variogram(optimization_parameters={"gradient_tolerance":
# 1e-8}, initialize_variogram=["sill", ])
# self = survey_data
# variable = "biomass"
# n_lags: int = 30
# azimuth_range: float = 360.0
# standardize_coordinates: bool = True
# force_lag_zero: bool = True
# model: Union[str, List[str]] = ["bessel", "exponential"]
# initialize_variogram: VariogramInitial = [
#     "nugget",
#     "sill",
#     "correlation_range",
#     "hole_effect_range",
#     "decay_power",
# ]
# verbose = True
# # Initialize Survey-class object
# self.analysis.update({"variogram": {}})

# # # Parameterize analysis settings that will be applied to the variogram fitting and analysis
# self.analysis["settings"].update(
#     {
#         "variogram": {
#             "azimuth_range": azimuth_range,
#             "fit_parameters": (
#                 initialize_variogram.keys()
#                 if isinstance(initialize_variogram, dict)
#                 else initialize_variogram
#             ),
#             "force_lag_zero": force_lag_zero,
#             "model": model,
#             "standardize_coordinates": standardize_coordinates,
#             "stratum_name": self.analysis["settings"]["transect"]["stratum_name"],
#             "variable": variable,
#             "verbose": verbose,
#         }
#     }
# )
# # Append `kriging_parameters` to the settings dictionary
# if standardize_coordinates:
#     self.analysis["settings"]["variogram"].update(
#         {"kriging_parameters": self.input["statistics"]["kriging"]["model_config"]}
#     )

# # Create a copy of the existing variogram settings
# default_variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
# # ---- Update model, n_lags
# default_variogram_parameters.update({"model": model, "n_lags": n_lags})
# optimization_parameters = VariogramOptimize.create()
# # # Create optimization settings dictionary
# # # ---- Add to settings
# # self.analysis["settings"]["variogram"].update({"optimization": optimization_parameters})
# variogram_parameters = {}
# default_variogram_parameters
# optimization_parameters
# initialize_variogram
# transect_dict = self.analysis["transect"]
# settings_dict = self.analysis["settings"]["variogram"]
# isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

# from echopop.utils.validate import VariogramEmpirical
# import echopop.spatial.variogram as esv
# from echopop.spatial.variogram import initialize_initial_optimization_values,
# initialize_optimization_config, initialize_variogram_parameters, empirical_variogram

# from echopop.spatial.transect import edit_transect_columns
# from echopop.spatial.projection import transform_geometry
# # Validate the relevant empirical variogram parameters
# empirical_variogram_params = VariogramEmpirical.create(**settings_dict)

# # # Initialize and validate the variogram model parameters
# valid_variogram_params = esv.initialize_variogram_parameters(
#     variogram_parameters, default_variogram_parameters
# )

# # # Initialize and validate the optimization parameters
# valid_optimization_params = esv.initialize_optimization_config(optimization_parameters)

# # # Initialize and validate the initial values/boundary inputs
# valid_initial_values = esv.initialize_initial_optimization_values(
#     initialize_variogram, valid_variogram_params
# )

# # # Prepare the transect data
# # ---- Create a copy of the transect dictionary
# transect_input = copy.deepcopy(transect_dict)
# # ---- Edit the transect data
# transect_data = edit_transect_columns(transect_input, settings_dict)

# # # Standardize the transect coordinates, if necessary
# if settings_dict["standardize_coordinates"]:
#     # ---- Transform geometry
#     transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)
#     # ---- Print message if verbose
#     if settings_dict["verbose"]:
#         # ---- Print alert
#         print(
#             "Longitude and latitude coordinates (WGS84) converted to standardized "
#             "coordinates (x and y)."
#         )
# else:
#     # ---- x
#     transect_data["x"] = "longitude"
#     # ---- y
#     transect_data["y"] = "latitude"

# # Compute the empirical variogram
# lags, gamma_h, lag_counts, _ = empirical_variogram(
#     transect_data, {**valid_variogram_params, **empirical_variogram_params}, settings_dict
# )

# #
# lag_counts = np.array([1, 20, 15, 35, 40, 80, 75, 100, 85, 120, 75])
# lags = np.linspace(0.0, 100.0, 11)
# gamma_h = np.array([0.000, 0.200, 0.300, 0.450, 0.775, 0.850, 0.952, 0.986, 0.976, 0.999, 0.984])
# initialize_variogram = ["nugget", "sill", "correlation_range", "hole_effect_range",
# "enhance_semivariance", "decay_power"]
# default_variogram_parameters = {
#     "nugget": 0.2, "sill": 0.25, "correlation_range": 20.0, "hole_effect_range": 1.00,
#     "decay_power": 2.00, "enhance_semivariance": True
# }
# variogram_parameters = {"model": "exponential",
#                         "lag_resolution": np.diff(lags).mean(),
#                         "n_lags": lags.size - 1,
#                         "range": lags.max()}
# valid_variogram_params = initialize_variogram_parameters(
#     variogram_parameters, default_variogram_parameters
# )
# optimization_parameters = {}
# valid_optimization_params = initialize_optimization_config(optimization_parameters)
# valid_initial_values = initialize_initial_optimization_values(
#     initialize_variogram, valid_variogram_params
# )
# optimization_settings = {
#     "parameters": valid_initial_values,
#     "config": valid_optimization_params,
# }
# variogram_parameters = valid_variogram_params

# from echopop.spatial.variogram import VARIOGRAM_MODELS
# template = {k: None for k in default_variogram_parameters}


@pytest.fixture
def optimize_variogram_data():

    # Define test cases
    data = {
        "normal_inputs": {
            "lag_counts": np.array([1, 20, 15, 35, 40, 80, 75, 100, 85, 120, 75]),
            "lags": np.linspace(0.0, 100.0, 11),
            "gamma_h": np.array(
                [0.000, 0.200, 0.300, 0.450, 0.775, 0.850, 0.952, 0.986, 0.976, 0.999, 0.984]
            ),
            "default_variogram_parameters": {
                "nugget": 0.2,
                "sill": 0.25,
                "correlation_range": 20.0,
                "hole_effect_range": 1.00,
                "decay_power": 2.00,
                "enhance_semivariance": True,
            },
            "initialize_variogram": [
                "nugget",
                "sill",
                "correlation_range",
                "hole_effect_range",
                "enhance_semivariance",
                "decay_power",
            ],
            "optimization_parameters": {},
            "exception": None,
        },
        "changed_default_parameters": {
            "lag_counts": np.array([1, 20, 15, 35, 40, 80, 75, 100, 85, 120, 75]),
            "lags": np.linspace(0.0, 100.0, 11),
            "gamma_h": np.array(
                [0.000, 0.200, 0.300, 0.450, 0.775, 0.850, 0.952, 0.986, 0.976, 0.999, 0.984]
            ),
            "default_variogram_parameters": {
                "nugget": 0.5,
                "sill": 0.75,
                "correlation_range": 10.0,
                "hole_effect_range": 5.00,
                "decay_power": 1.25,
                "enhance_semivariance": True,
            },
            "initialize_variogram": [
                "nugget",
                "sill",
                "correlation_range",
                "hole_effect_range",
                "enhance_semivariance",
                "decay_power",
            ],
            "optimization_parameters": {},
            "exception": None,
        },
        "adjusted_some_optimization_params": {
            "lag_counts": np.array([1, 20, 15, 35, 40, 80, 75, 100, 85, 120, 75]),
            "lags": np.linspace(0.0, 100.0, 11),
            "gamma_h": np.array(
                [0.000, 0.200, 0.300, 0.450, 0.775, 0.850, 0.952, 0.986, 0.976, 0.999, 0.984]
            ),
            "default_variogram_parameters": {
                "nugget": 0.5,
                "sill": 0.75,
                "correlation_range": 10.0,
                "hole_effect_range": 5.00,
                "decay_power": 1.25,
                "enhance_semivariance": True,
            },
            "initialize_variogram": [
                "nugget",
                "sill",
                "correlation_range",
                "hole_effect_range",
                "enhance_semivariance",
                "decay_power",
            ],
            "optimization_parameters": {
                "max_fun_evaluations": 100,
                "gradient_tolerance": 1e-10,
                "solution_tolerance": 1e-10,
                "cost_fun_tolerance": 1e-12,
            },
            "exception": None,
        },
        # "nugget_only": {
        #     "lag_counts": np.array([1, 20, 15, 35, 40, 80, 75, 100, 85, 120, 75]),
        #     "lags": np.linspace(0.0, 100.0, 11),
        #     "gamma_h": np.array([0.000, 0.200, 0.300, 0.450, 0.775, 0.850, 0.952, 0.986, 0.976,
        #                         0.999, 0.984]),
        #     "default_variogram_parameters": {
        #         "nugget": 0.5, "sill": 0.75, "correlation_range": 10.0, "hole_effect_range": 5.00,
        #         "decay_power": 1.25, "enhance_semivariance": True
        #     },
        #     "initialize_variogram": ["nugget"],
        #     "optimization_parameters": {"max_fun_evaluations": 100, "gradient_tolerance": 1e-10,
        #                                 "solution_tolerance": 1e-10, "cost_fun_tolerance": 1e-12},
        #     "exception": None
        # },
        "assign_init_bounds": {
            "lag_counts": np.array([1, 20, 15, 35, 40, 80, 75, 100, 85, 120, 75]),
            "lags": np.linspace(0.0, 100.0, 11),
            "gamma_h": np.array(
                [0.100, 0.200, 0.300, 0.450, 0.775, 0.850, 0.952, 0.986, 0.976, 0.999, 0.984]
            ),
            "default_variogram_parameters": {
                "nugget": 0.5,
                "sill": 0.75,
                "correlation_range": 5.0,
                "hole_effect_range": 5.00,
                "decay_power": 1.50,
                "enhance_semivariance": True,
            },
            "initialize_variogram": {
                "decay_power": {"min": 1.00, "max": 2.00},
                "sill": {"min": 0.50, "max": 1.00},
                "correlation_range": {"max": 5.0},
            },
            "optimization_parameters": {},
            "exception": None,
        },
        "assign_init_bounds_additional": {
            "lag_counts": np.array([1, 20, 15, 35, 40, 80, 75, 100, 85, 120, 75]),
            "lags": np.linspace(0.0, 100.0, 11),
            "gamma_h": np.array(
                [0.100, 0.200, 0.300, 0.450, 0.775, 0.850, 0.952, 0.986, 0.976, 0.999, 0.984]
            ),
            "default_variogram_parameters": {
                "nugget": 0.5,
                "sill": 1.00,
                "correlation_range": 5.0,
                "hole_effect_range": 5.00,
                "decay_power": 1.50,
                "enhance_semivariance": True,
            },
            "initialize_variogram": {
                "decay_power": {"min": 1.45, "max": 1.55},
                "sill": {"min": 0.95, "max": 1.05},
                "correlation_range": {"min": 5.0, "max": 10.0},
                "hole_effect_range": {"value": 1.50},
            },
            "optimization_parameters": {},
            "exception": None,
        },
    }

    yield data


@pytest.fixture
def expected_optimize_variogram(test_path):

    # -------------------------
    # Load the expected outcomes
    expected_results = load_json_data(
        str(test_path["EXPECTED"] / ___EXPECTED_OUTCOME_FILENAME__), test_name="optimize_variogram"
    )

    return expected_results


@pytest.mark.parametrize(
    "parameterization",
    [
        "normal_inputs",
        "changed_default_parameters",
        "adjusted_some_optimization_params",
        # "nugget_only",
        "assign_init_bounds",
        "assign_init_bounds_additional",
    ],
    ids=[
        "Standard inputs",
        "Change default variogram parameter values",
        "Adjust optimization parameters",
        # "Nugget-only fit",
        "Specify initial boundaries and values",
        "Modify additional initial values and boundaries",
    ],
)
def test_optimize_variogram(parameterization, expected_optimize_variogram, optimize_variogram_data):

    # Get the key args
    expected_results = expected_optimize_variogram[parameterization]

    # Get the input parameters
    test_data = optimize_variogram_data[parameterization]
    lag_counts = test_data["lag_counts"]
    lags = test_data["lags"]
    gamma_h = test_data["gamma_h"]
    initialize_variogram = test_data["initialize_variogram"]
    default_variogram_parameters = test_data["default_variogram_parameters"]
    optimization_parameters = test_data["optimization_parameters"]

    for mod in {**VARIOGRAM_MODELS["single"], **VARIOGRAM_MODELS["composite"]}.keys():
        if mod == ("cosine", "exponential") or mod == ("cosine", "gaussian"):
            pass
        else:
            # Amend model name if needed
            if isinstance(mod, tuple):
                mod = list(mod)
                model_name = "_".join(mod)
            else:
                model_name = mod

            # Get relevant arguments for model
            valid_args = expected_results[model_name]["initial_fit"][0]

            # Filter any defined initialized parameters
            if isinstance(initialize_variogram, list):
                initial_parameters = [
                    param for param in initialize_variogram if param in valid_args
                ]
            else:
                initial_parameters = {
                    k: v for k, v in initialize_variogram.items() if k in valid_args
                }

            # Filter the default parameters dictionary
            valid_defaults = {
                k: v for k, v in default_variogram_parameters.items() if k in valid_args
            }

            # Compute the variogram parameters
            variogram_parameters = {
                "model": mod,
                "lag_resolution": np.diff(lags).mean(),
                "n_lags": lags.size - 1,
                "range": lags.max(),
            }

            # Initialize the parameters
            valid_variogram_params = initialize_variogram_parameters(
                variogram_parameters, valid_defaults
            )

            # Initialize the optimization parameters
            valid_optimization_params = initialize_optimization_config(optimization_parameters)

            # Initialize the initial parameter values
            # ---- Reorder if needed
            if isinstance(initial_parameters, dict):
                initial_parameters = {k: initial_parameters[k] for k in sorted(initial_parameters)}
            else:
                initial_parameters = sorted(initial_parameters)
            # ---- And validate
            valid_initial_values = initialize_initial_optimization_values(
                initial_parameters, valid_variogram_params
            )

            # Bundle the optimization parameters and configuration
            optimization_settings = {
                "parameters": valid_initial_values,
                "config": valid_optimization_params,
            }

            # Compute the results
            best_fit_variogram, initial_fit, optimized_fit = optimize_variogram(
                lag_counts, lags, gamma_h, optimization_settings, **valid_variogram_params
            )
            # ---- Force order
            new_order = expected_results[model_name]["initial_fit"][0]
            # -------- Initialize new Parameters()
            old_order = [np.where(np.array([new_order]) == o)[1].tolist() for o in initial_fit[0]]
            # -------- Adjust the fit
            initial_fit = list(initial_fit)
            # ----
            if initial_fit[0] != new_order:
                initial_fit[0] = [
                    initial_fit[0][item] for item in np.concatenate(old_order).tolist()
                ]
                initial_fit[1] = [
                    initial_fit[1][item] for item in np.concatenate(old_order).tolist()
                ]
            optimized_fit = list(optimized_fit)
            # ----
            if optimized_fit[0] != new_order:
                optimized_fit[0] = [
                    optimized_fit[0][item] for item in np.concatenate(old_order).tolist()
                ]
                optimized_fit[1] = [
                    optimized_fit[1][item] for item in np.concatenate(old_order).tolist()
                ]
            # ---- Format as dictionary
            test_results = {
                "best_fit_variogram": best_fit_variogram,
                "initial_fit": initial_fit,
                "optimized_fit": optimized_fit,
            }

            # Get the associates results for the specific model
            assert_dictionary_values_equal(expected_results[model_name], test_results)


@pytest.fixture
def semivariance_data():

    # Set a random seed for reproducibility
    np.random.seed(100)

    # Generate mxn size arrays
    # ---- Rows (m)
    m = int(1e2)
    # ---- Columns (n)
    n = int(20)

    # Define test cases
    data = {
        "standard_inputs": (
            np.array([1.0] * m),
            np.array([2e09, 1e10, 5e09] + [1e09] * (n - 3)),
            np.array([2e09, 1e10, 5e09] + [1e09] * (n - 3)) ** 2,
            np.array([5000, 10000, 3000] + [10000] * (n - 3)),
            np.array([2e15, 4e16, 1e17] + [1e14] * (n - 3)),
            np.random.randint(0, 10, size=(m, n)),
        ),
        "uniform_inputs": (
            np.array([1.0] * m),
            np.array([1e09] * n),
            np.array([1e09] * n) ** 2,
            np.array([1000] * n),
            np.array([1e14] * n),
            np.ones((m, n)),
        ),
        "empty_inputs": (
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([[]]),
        ),
        "zero_counts": (
            np.array([5.0] * m),
            np.array([1e09] * n),
            np.array([1e09] * n),
            np.array([5000] * n),
            np.zeros(n),
            np.ones((m, n)),
        ),
        "nan_deviations": (
            np.array([5.0] * m),
            np.array([1e09] * n),
            np.array([1e09] * n),
            np.array([5000] * n),
            np.full(n, np.nan),
            np.ones((m, n)),
        ),
        "nan_deviations_estimates": (
            np.full(m, np.nan),
            np.array([1e09] * n),
            np.array([1e09] * n),
            np.array([5000] * n),
            np.full(n, np.nan),
            np.ones((m, n)),
        ),
        "all_nan": (
            np.full(m, np.nan),
            np.full(n, np.nan),
            np.full(n, np.nan),
            np.full(n, np.nan),
            np.full(n, np.nan),
            np.full((m, n), np.nan),
        ),
        "mismatched_shapes": (
            np.full(n, np.nan),
            np.full(n, np.nan),
            np.full(n, np.nan),
            np.full(n, np.nan),
            np.full(n, np.nan),
            np.full((m, n), np.nan),
        ),
    }

    yield data


@pytest.fixture
def expected_semivariance(test_path):

    # -------------------------
    # Load the expected outcomes
    expected_results = load_json_data(
        str(test_path["EXPECTED"] / ___EXPECTED_OUTCOME_FILENAME__), test_name="semivariance"
    )

    return expected_results


@pytest.mark.parametrize(
    "test_name",
    [
        "standard_inputs",
        "uniform_inputs",
        "empty_inputs",
        "zero_counts",
        "nan_deviations",
        "nan_deviations_estimates",
        "all_nan",
        "mismatched_shapes",
    ],
    ids=[
        "Standard (valid) inputs",
        "No variability in inputs",
        "Empty inputs",
        "All lag bins with zero-counts",
        "Lag deviations are all NaN",
        "Estimates and lag deviations are all NaN",
        "All inputs are NaN",
        "Mismatch in input array shapes",
    ],
)
def test_semivariance(test_name, semivariance_data, expected_semivariance):

    # Get the key args
    expected_results = expected_semivariance[test_name]

    # Get the input parameters
    estimates, lag_estimates, lag_estimates_squared, lag_counts, lag_deviations, head_index = (
        semivariance_data[test_name]
    )
    # ---- Package inputs together as a tuple
    inputs_tpl = (
        estimates,
        lag_estimates,
        lag_estimates_squared,
        lag_counts,
        lag_deviations,
        head_index,
    )

    # Test for exceptions
    if expected_results["exception"] is not None:
        with pytest.raises(ValueError, match=re.escape(expected_results["exception"])):
            assert semivariance(*inputs_tpl)

    else:
        # Compute the semivariance and mean lag covariance
        output_tpl = semivariance(*inputs_tpl)

        # Convert `expected_results["expected"]` to a tuple
        expected_tuple = tuple(
            expected_results[key] for key in expected_results if key != "exception"
        )

        # Assert same shapes (for just arrays)
        assert all(
            [
                e.shape == o.shape
                for (e, o) in zip(expected_tuple, output_tpl)
                if not isinstance(o, (int, float))
            ]
        )

        # Assert value equality
        assert [np.array_equal(e, o, equal_nan=True) for (e, o) in zip(expected_tuple, output_tpl)]
