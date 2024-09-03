import echopop.spatial.variogram as esv
import pytest
import re
import pandas as pd
import numpy as np

from echopop.spatial.variogram import prepare_variogram_matrices

from echopop.tests.conftest import load_json_data

___EXPECTED_OUTCOME_FILENAME__ = "variogram_analysis.json"

@pytest.mark.parametrize(
    "model, expected, exception",
    [
        (
            "exponential",
            ({"distance_lags", "sill", "nugget", "correlation_range"},
             esv.exponential), None
        ),
        (
             "gaussian",
             ({"distance_lags", "sill", "nugget", "correlation_range"},
             esv.gaussian), None
        ),
        (
             "jbessel",
             ({"distance_lags", "sill", "nugget", "hole_effect_range"},
             esv.jbessel), None
        ),
        (
             "kbessel",
             ({"distance_lags", "sill", "nugget", "hole_effect_range"}, 
             esv.kbessel), None
        ),
        (
             "linear",
             ({"distance_lags", "sill", "nugget"}, 
             esv.linear), None
        ),
        (
             "nugget",
             ({"distance_lags", "sill", "nugget"}, 
             esv.nugget), None
        ),
        (
             "sinc",
             ({"distance_lags", "sill", "nugget", "hole_effect_range"}, 
             esv.sinc), None
        ),
        (
             "spherical",
             ({"distance_lags", "sill", "nugget", "correlation_range"}, 
             esv.spherical), None
        ),
        (
             ("bessel", "exponential"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range", 
              "decay_power"}, esv.bessel_exponential), None
        ),
        (
             ("bessel", "gaussian"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"}, 
              esv.bessel_gaussian), None
        ),
        (
             ("cosine", "exponential"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range", 
              "enhance_semivariance"}, esv.cosine_exponential), None
        ),
        (
             ("cosine", "gaussian"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"}, 
              esv.cosine_gaussian), None
        ),
        (
             ("exponential", "linear"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range", 
              "decay_power"}, esv.exponential_linear), None
        ),
        (
             ("gaussian", "linear"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"}, 
              esv.gaussian_linear), None
        ),
        (
             "invalid",
             (None, None), "could not be matched to an existing variogram method."
        ),
        (
             ("invalid", "tuple"),
             (None, None), "could not be matched to an existing variogram method."
        ),
        (
             ["invalid", "list"],
             (None, None), "could not be matched to an existing variogram method."
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
        "Invalid composite variogram model (list)"
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
        "max_fun_evaluations": {
            "name": "max_nfev",
            "value": 500
        },
        "cost_fun_tolerance": {
            "name": "ftol",
            "value": 1e-8
        },
        "solution_tolerance": {
            "name": "xtol",
            "value": 1e-8
        },
        "gradient_tolerance": {
            "name": "gtol",
            "value": 1e-8
        },
        "finite_step_size": {
            "name": "diff_step",
            "value": 1e-8
        },
        "trust_region_solver": {
            "name": "tr_solver",
            "value": "exact"
        },
        "x_scale": {
            "name": "x_scale",
            "value": "jac"
        },
        "jacobian_approx": {
            "name": "jac",
            "value": "2-point"
        },
    }
    # ---- Iterate through keys to assess equality across names [ EQUIVALENT KEYS ]
    assert all([TRANSLATION_DICT[key]["name"] in output.keys()
                for key in ECHOPOP_OPTIMIZATION_INPUTS.keys()])
    # ---- Iterate through keys to assess equality across values [ EQUIVALENT VALUES ]
    assert all([output[TRANSLATION_DICT[key]["name"]] == TRANSLATION_DICT[key]["value"] 
                for key in ECHOPOP_OPTIMIZATION_INPUTS.keys()])

@pytest.fixture()
def transect_data():
    # Set a random seed for reproducibility
    np.random.seed(99)

    # Generate 20 random rows of data
    num_rows = 5
    data = {
        'transect_num': np.random.choice([1, 2], size=num_rows),
        'longitude': np.random.uniform(-121.5, -120.5, size=num_rows),
        'latitude': np.random.uniform(34.0, 35.0, size=num_rows),
        'stratum_num': np.random.choice([1, 2], size=num_rows),
        'transect_spacing': np.random.uniform(5.0, 15.0, size=num_rows),
        'biomass': np.random.uniform(0.0, 10.0, size=num_rows),
        'biomass_density': np.random.uniform(0.0, 2.0, size=num_rows),
        'x': np.random.uniform(-0.5, 0.5, size=num_rows),
        'y': np.random.uniform(-0.5, 0.5, size=num_rows)
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
    expected_results = load_json_data(str(test_path["EXPECTED"] / ___EXPECTED_OUTCOME_FILENAME__),
                                      test_name="prepare_variogram_matrices")
    
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

# # Parameterize analysis settings that will be applied to the variogram fitting and analysis
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
# # Create optimization settings dictionary
# # ---- Add to settings
# self.analysis["settings"]["variogram"].update({"optimization": optimization_parameters})
# variogram_parameters = {}
# default_variogram_parameters
# optimization_parameters
# initialize_variogram
# transect_dict = self.analysis["transect"]
# settings_dict = self.analysis["settings"]["variogram"]
# isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

# from echopop.utils.validate import VariogramEmpirical
# import echopop.spatial.variogram as esv
# from echopop.spatial.transect import edit_transect_columns
# from echopop.spatial.projection import transform_geometry
# # Validate the relevant empirical variogram parameters
# empirical_variogram_params = VariogramEmpirical.create(**settings_dict)

# # Initialize and validate the variogram model parameters
# valid_variogram_params = esv.initialize_variogram_parameters(
#     variogram_parameters, default_variogram_parameters
# )

# # Initialize and validate the optimization parameters
# valid_optimization_params = esv.initialize_optimization_config(optimization_parameters)

# # Initialize and validate the initial values/boundary inputs
# valid_initial_values = esv.initialize_initial_optimization_values(
#     initialize_variogram, valid_variogram_params
# )

# # Prepare the transect data
# # ---- Create a copy of the transect dictionary
# transect_input = copy.deepcopy(transect_dict)
# # ---- Edit the transect data
# transect_data = edit_transect_columns(transect_input, settings_dict)

# # Standardize the transect coordinates, if necessary
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
# # lags, gamma_h, lag_counts, _ = empirical_variogram(
# #     transect_data, {**valid_variogram_params, **empirical_variogram_params}, settings_dict
# # )
# import pandas as pd
# import numpy as np 
# from echopop.spatial.mesh import griddify_lag_distances
# from echopop.spatial.variogram import variogram_matrix_filter
# variogram_parameters = {**valid_variogram_params, **empirical_variogram_params}

# def prepare_variogram_matrices(transect_data: pd.DataFrame,
#                                lag_resolution: float,
#                                **kwargs):
    
#     # Calculate the lag distance matrix among transect data
#     transect_distance_matrix, transect_azimuth_matrix = griddify_lag_distances(
#         transect_data, transect_data, angles=True
#     )

#     # Compute the lag matrix
#     lag_matrix = np.round(transect_distance_matrix / lag_resolution).astype(int) + 1

#     return transect_azimuth_matrix, lag_matrix


# azimuth_matrix, lag_matrix = prepare_variogram_matrices(transect_data, **variogram_parameters)

# def compute_lag_deviations(
#     estimates: np.ndarray[float],
#     equivalent_lags: np.ndarray[int],
#     mask: np.ndarray[bool],
#     azimuth_matrix: np.ndarray[float],
#     n_lags: int,
#     azimuth_range: float
# ):
#     """
#     Compute the deviations within each lag
#     """

#     # Subset the lag array via a boolean bitmap
#     lag_bitmap = equivalent_lags < n_lags

#     # Create a dummy array that produces the row indices for the estimate matrix/array and then
#     # apply the triangle mask and azimuth filter
#     estimate_rows = variogram_matrix_filter(
#         np.arange(len(estimates))[:, np.newaxis],
#         mask,
#         azimuth_matrix,
#         azimuth_range,
#     )

#     # Calculate the deviations between indexed estimates and lag-specific ones
#     deviations = (estimates[estimate_rows][lag_bitmap] - estimates_filtered[lag_bitmap]) ** 2