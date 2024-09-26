# import numpy as np
# import pandas as pd

# from echopop.biology import (
#     age1_metric_proportions
# )

# def test_age1_metric_proportions():

#     # -------------------------
#     # Mock values
#     # ---- Input dictionary
#     mock_input_dict = {
#         "biology": {
#             "distributions": {
#                 "length_bins_df": pd.DataFrame(
#                     {
#                         "length_bins": np.array([12.5, 17.5, 22.5, 27.5]),
#                         "length_intervals": np.array(
#                             [
#                                 pd.Interval(left=10.0, right=15.0),
#                                 pd.Interval(left=15.0, right=20.0),
#                                 pd.Interval(left=20.0, right=25.0),
#                                 pd.Interval(left=25.0, right=30.0),
#                             ]
#                         ),
#                     }
#                 )
#             }
#         }
#     }
#     # ---- Proportions dictionary
#     mock_proportions_dict = {
#         "number": {
#             "aged_length_proportions_df": pd.DataFrame(
#                 {
#                     "stratum_num": np.repeat([1, 2], 24),
#                     "species_id": np.repeat([9933], 48),
#                     "length_bin": np.tile(
#                         [
#                             pd.Interval(left=10.0, right=15.0),
#                             pd.Interval(left=15.0, right=20.0),
#                             pd.Interval(left=20.0, right=25.0),
#                             pd.Interval(left=25.0, right=30.0),
#                         ],
#                         12,
#                     ),
#                     "age_bin": np.repeat(
#                         [
#                             pd.Interval(left=0.5, right=1.5),
#                             pd.Interval(left=1.5, right=2.5),
#                             pd.Interval(left=2.5, right=3.5),
#                         ],
#                         16,
#                     ),
#                     "sex": np.tile(
#                         [
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                         ],
#                         2,
#                     ),
#                     "proportion_number_aged": np.array(
#                         [
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.50,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                         ]
#                     ),
#                 }
#             ),
#             "unaged_length_proportions_df": pd.DataFrame(
#                 {
#                     "stratum_num": np.repeat([1, 2], 12),
#                     "species_id": np.repeat([9933], 24),
#                     "length_bin": np.tile(
#                         [
#                             pd.Interval(left=10.0, right=15.0),
#                             pd.Interval(left=15.0, right=20.0),
#                             pd.Interval(left=20.0, right=25.0),
#                             pd.Interval(left=25.0, right=30.0),
#                         ],
#                         6,
#                     ),
#                     "sex": np.tile(
#                         [
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                         ],
#                         2,
#                     ),
#                     "proportion_number_unaged": np.array(
#                         [
#                             0.50,
#                             0.50,
#                             0.25,
#                             0.25,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.50,
#                             0.50,
#                             0.50,
#                             0.50,
#                             0.50,
#                             0.50,
#                             0.50,
#                             0.50,
#                         ]
#                     ),
#                 }
#             ),
#         },
#         "weight": {
#             "aged_weight_proportions_df": pd.DataFrame(
#                 {
#                     "stratum_num": np.repeat([1, 2], 24),
#                     "species_id": np.repeat([9933], 48),
#                     "length_bin": np.tile(
#                         [
#                             pd.Interval(left=15.0, right=20.0),
#                             pd.Interval(left=20.0, right=25.0),
#                             pd.Interval(left=20.0, right=25.0),
#                             pd.Interval(left=25.0, right=30.0),
#                         ],
#                         12,
#                     ),
#                     "age_bin": np.repeat(
#                         [
#                             pd.Interval(left=0.5, right=1.5),
#                             pd.Interval(left=1.5, right=2.5),
#                             pd.Interval(left=2.5, right=3.5),
#                         ],
#                         16,
#                     ),
#                     "sex": np.tile(
#                         [
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "all",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "male",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                             "female",
#                         ],
#                         2,
#                     ),
#                     "weight_proportions": np.array(
#                         [
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.50,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.00,
#                             0.00,
#                             0.50,
#                             0.50,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                             0.25,
#                         ]
#                     ),
#                 }
#             )
#         },
#     }
#     # ---- Mock TS-length regression parameters
#     mock_TS_L_parameters = {
#         "number_code": 9933,
#         "TS_L_slope": 10.0,
#         "TS_L_intercept": -40.0,
#         "length_units": "cm",
#     }
#     # ---- Mock distributions dictionary
#     mock_distributions_dict = {
#         "length_bins_df": pd.DataFrame(
#             {
#                 "length_bins": np.array([12.5, 17.5, 22.5, 27.5]),
#                 "length_intervals": np.array(
#                     [
#                         pd.Interval(left=15.0, right=20.0),
#                         pd.Interval(left=20.0, right=25.0),
#                         pd.Interval(left=20.0, right=25.0),
#                         pd.Interval(left=25.0, right=30.0),
#                     ]
#                 ),
#             }
#         )
#     }
#     # ---- Mock settings dictionary
#     mock_settings_dict = {"transect": {"stratum_name": "stratum_num"}}

#     age1_metric_proportions(
#         mock_distributions_dict, mock_proportions_dict, mock_TS_L_parameters, mock_settings_dict
#     )

#     pd.DataFrame(
#         {
#             "stratum_num": np.repeat([1, 2], 12),
#             "species_id": np.repeat([9933], 24),
#             "length_bin": np.tile(
#                 [pd.Interval(left=15.0, right=20.0), pd.Interval(left=20.0, right=25.0)], 12
#             ),
#             "age_bin": np.repeat(
#                 [pd.Interval(left=0.5, right=1.5), pd.Interval(left=1.5, right=2.5)], 12
#             ),
#             "sex": np.tile(
#                 [
#                     "all",
#                     "all",
#                     "all",
#                     "all",
#                     "male",
#                     "male",
#                     "male",
#                     "male",
#                     "female",
#                     "female",
#                     "female",
#                     "female",
#                 ],
#                 2,
#             ),
#             "proportion_number_aged": np.array(
#                 [
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.50,
#                     0.00,
#                     0.00,
#                     0.50,
#                     0.50,
#                     0.00,
#                     0.00,
#                     0.50,
#                     0.50,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                     0.25,
#                 ]
#             ),
#         }
#     )
#     mock_length_array = np.array([1.0, 2.0, 3.0])
#     # ---- x values [ FLOAT input ]
#     mock_length_float = np.array([1.0])
#     # ---- Slope
#     mock_slope = 5.0
#     # ---- Intercept
#     mock_intercept = -2.0

#     # -------------------------
#     # Evaluate [ ARRAY ]
#     test_results_array = ts_length_regression(mock_length_array, mock_slope, mock_intercept)
#     # Evaluate [ FLOAT ]
#     test_results_float = ts_length_regression(mock_length_float, mock_slope, mock_intercept)

#     # -------------------------
#     # Expected outcomes
#     # Expected [ ARRAY ]
#     expected_array = np.array([-2.0, -0.49485002, 0.38560627])
#     # Expected [ FLOAT ]
#     expected_float = np.array([-2.0])

#     # -------------------------
#     # Run tests [ ARRAY ]
#     # ---- Type
#     assert type(test_results_array) == np.ndarray
#     # ---- Equality
#     assert np.allclose(test_results_array, expected_array)
#     # Run tests [ FLOAT ]
#     # ---- Type
#     assert type(test_results_float) == np.ndarray
#     # ---- Equality
#     assert np.allclose(test_results_float, expected_float)
