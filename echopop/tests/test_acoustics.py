import numpy as np
import pandas as pd

from echopop.acoustics import impute_missing_sigma_bs, to_dB, to_linear, ts_length_regression
from echopop.tests.conftest import assert_dataframe_equal


def test_ts_length_regression():

    # -------------------------
    # Mock values
    # ---- length values [ ARRAY input ]
    mock_length_array = np.array([1.0, 2.0, 3.0])
    # ---- x values [ FLOAT input ]
    mock_length_float = np.array([1.0])
    # ---- Slope
    mock_slope = 5.0
    # ---- Intercept
    mock_intercept = -2.0

    # -------------------------
    # Evaluate [ ARRAY ]
    test_results_array = ts_length_regression(mock_length_array, mock_slope, mock_intercept)
    # Evaluate [ FLOAT ]
    test_results_float = ts_length_regression(mock_length_float, mock_slope, mock_intercept)

    # -------------------------
    # Expected outcomes
    # Expected [ ARRAY ]
    expected_array = np.array([-2.0, -0.49485002, 0.38560627])
    # Expected [ FLOAT ]
    expected_float = np.array([-2.0])

    # -------------------------
    # Run tests [ ARRAY ]
    # ---- Type
    assert isinstance(test_results_array, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_array, expected_array)
    # Run tests [ FLOAT ]
    # ---- Type
    assert isinstance(test_results_float, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_float, expected_float)


def test_to_linear():

    # -------------------------
    # Mock values
    # ---- length values [ ARRAY input ]
    mock_db_array = np.array([-80.0, -60.0, -40.0])
    # ---- x values [ FLOAT input ]
    mock_db_float = np.array([-60.0])

    # -------------------------
    # Evaluate [ ARRAY ]
    test_results_array = to_linear(mock_db_array)
    # Evaluate [ FLOAT ]
    test_results_float = to_linear(mock_db_float)

    # -------------------------
    # Expected outcomes
    # Expected [ ARRAY ]
    expected_array = np.array([1e-8, 1e-6, 1e-4])
    # Expected [ FLOAT ]
    expected_float = np.array([1e-6])

    # -------------------------
    # Run tests [ ARRAY ]
    # ---- Type
    assert isinstance(test_results_array, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_array, expected_array)
    # Run tests [ FLOAT ]
    # ---- Type
    assert isinstance(test_results_float, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_float, expected_float)


def test_to_dB():

    # -------------------------
    # Mock values
    # ---- length values [ ARRAY input ]
    mock_linear_array = np.array([1e-8, 1e-6, 1e-4])
    # ---- x values [ FLOAT input ]
    mock_linear_float = np.array([1e-6])

    # -------------------------
    # Evaluate [ ARRAY ]
    test_results_array = to_dB(mock_linear_array)
    # Evaluate [ FLOAT ]
    test_results_float = to_dB(mock_linear_float)

    # -------------------------
    # Expected outcomes
    # Expected [ ARRAY ]
    expected_array = np.array([-80.0, -60.0, -40.0])
    # Expected [ FLOAT ]
    expected_float = np.array([-60.0])

    # -------------------------
    # Run tests [ ARRAY ]
    # ---- Type
    assert isinstance(test_results_array, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_array, expected_array)
    # Run tests [ FLOAT ]
    # ---- Type
    assert isinstance(test_results_float, np.ndarray)
    # ---- Equality
    assert np.allclose(test_results_float, expected_float)


def test_impute_missing_sigma_bs():

    # -------------------------
    # Mock values
    # ---- Stratified mean sigma_bs [ NO MISSING STRATA ]
    mock_sigma_bs_full = pd.DataFrame(
        {
            "stratum_num": [1, 2, 3, 4, 5],
            "species_id": np.repeat(94832, 5),
            "sigma_bs_mean": [1e-1, 1e-2, 1e-3, 1e-4, 1e-5],
        }
    )
    # ---- Stratified mean sigma_bs [ MISSING STRATA IN MIDDLE: STRATUM == 3 ]
    mock_sigma_bs_mid = pd.DataFrame(
        {
            "stratum_num": [1, 2, 4, 5],
            "species_id": np.repeat(94832, 4),
            "sigma_bs_mean": [1e-1, 1e-2, 1e-4, 1e-5],
        }
    )
    # ---- Stratified mean sigma_bs [ MISSING STRATA ON EDGES: STRATUM == [1, 5] ]
    mock_sigma_bs_edge = pd.DataFrame(
        {
            "stratum_num": [2, 3, 4],
            "species_id": np.repeat(94832, 3),
            "sigma_bs_mean": [1e-2, 1e-3, 1e-4],
        }
    )
    # ---- Strata definition dictionary
    mock_strata_dict = {"strata_df": pd.DataFrame({"stratum_num": [1, 2, 3, 4, 5]})}

    # -------------------------
    # Evaluate [ FULL ]
    test_results_full = impute_missing_sigma_bs(mock_strata_dict, mock_sigma_bs_full)
    # Evaluate [ MID ]
    test_results_mid = impute_missing_sigma_bs(mock_strata_dict, mock_sigma_bs_mid)
    # Evaluate [ EDGE ]
    test_results_edge = impute_missing_sigma_bs(mock_strata_dict, mock_sigma_bs_edge)

    # -------------------------
    # Expected outcomes
    # ---- Types [~ALL]
    expected_dtypes = {
        "stratum_num": np.integer,
        "species_id": np.integer,
        "sigma_bs_mean": np.floating,
    }
    # Expected [ FULL ]
    expected_full = mock_sigma_bs_full
    # Expected [ MID ]
    expected_mid = pd.DataFrame(
        {
            "stratum_num": [1, 2, 3, 4, 5],
            "species_id": np.repeat(94832, 5),
            "sigma_bs_mean": [1e-1, 1e-2, 5.05e-3, 1e-4, 1e-5],
        }
    )
    # Expected [ EDGE ]
    expected_edge = pd.DataFrame(
        {
            "stratum_num": [1, 2, 3, 4, 5],
            "species_id": np.repeat(94832, 5),
            "sigma_bs_mean": [1e-2, 1e-2, 1e-3, 1e-4, 1e-4],
        }
    )

    # -------------------------
    # Run tests [ FULL ]
    # ---- Type [ OUTPUT: FULL ]
    assert isinstance(test_results_full, pd.DataFrame)
    # ---- Shape, types, values [ OUTPUT: FULL ]
    assert_dataframe_equal(test_results_full, expected_dtypes, expected_full)
    # Run tests [ MID ]
    # ---- Type [ OUTPUT: MID ]
    assert isinstance(test_results_mid, pd.DataFrame)
    # ---- Shape, types, values [ OUTPUT: MID ]
    assert_dataframe_equal(test_results_mid, expected_dtypes, expected_mid)
    # Run tests [ EDGE ]
    # ---- Type [ OUTPUT: EDGE ]
    assert isinstance(test_results_edge, pd.DataFrame)
    # ---- Shape, types, values [ OUTPUT: EDGE ]
    assert_dataframe_equal(test_results_edge, expected_dtypes, expected_edge)


# def test_summarize_sigma_bs( ):

#     # -------------------------
#     # Mock values
#     # ---- Length data (unaged)
#     mock_length_data = pd.DataFrame( {
#             'stratum_num': np.repeat( [ 1 , 2 , 3 ] , 4 ) ,
#             'haul_num': np.repeat( [ 10 , 20 , 30 ] , 4 ) ,
#             'species_id': np.append( np.repeat( [ 19350 ] , 11 ) , 43130 ) ,
#             'length': np.linspace( 15 , 100 , 12 ) ,
#             'length_count': np.linspace( 10 , 105 , 12 ).astype( int ) ,
#             'group_sex': np.repeat( 'sexed' , 12 ) } )
#     # -------- > DROP SPECIES == 43130
#     mock_length_data_sub = mock_length_data[ mock_length_data[ 'species_id' ] != 43130 ]
#     # ---- Specimen data (aged)
#     mock_specimen_data = pd.DataFrame(
#         {
#             'stratum_num': np.repeat( [ 1 , 2 , 3 ] , 4 ) ,
#             'haul_num': np.tile( [ 5 , 10 , 15 , 20 , 25 , 30 ] , 2 ) ,
#             'species_id': np.append( np.repeat( [ 19350 ] , 11 ) , 43130 ) ,
#             'length': np.linspace( 10 , 100 , 12 ) ,
#             'group_sex': np.repeat( 'sexed' , 12 )
#         }
#     )
#     # -------- > DROP SPECIES == 43130
#     mock_specimen_data_sub = mock_specimen_data[ mock_specimen_data[ 'species_id' ] != 43130 ]
#     # -------- > DROP SPECIES == 43130
#     # ---- Strata definition dictionary
#     mock_strata_dict = {
#         'strata_df': pd.DataFrame( { 'stratum_num': [ 1 , 2 , 3 ] } )
#     }
#     # ---- Configuration dictionary
#     mock_configuration_dict = {
#         'TS_length_regression_parameters': {
#             'mega_guppies': {
#                 'number_code': 19350 ,
#                 'TS_L_slope': 10.0 ,
#                 'TS_L_intercept': -30.0 ,
#                 'length_units': 'cm'
#             } ,
#             'flying_unicorns': {
#                 'number_code': 43130 ,
#                 'TS_L_slope': 10.0 ,
#                 'TS_L_intercept': -40.0 ,
#                 'length_units': 'cm'
#             }
#         }
#     }
#     # ---- Settings dictionary
#     mock_settings_dict = {
#         'transect': {
#             'stratum_name': 'stratum_num'
#         }
#     }

#     # -------------------------
#     # Evaluate [ FULL ]
#     test_results_full = summarize_sigma_bs( mock_length_data ,
#                                             mock_specimen_data ,
#                                             mock_strata_dict ,
#                                             mock_configuration_dict ,
#                                             mock_settings_dict )
#     # Evaluate [ SUB ]
#     test_results_sub = summarize_sigma_bs( mock_length_data_sub ,
#                                            mock_specimen_data_sub ,
#                                            mock_strata_dict ,
#                                            mock_configuration_dict ,
#                                            mock_settings_dict )

#     # -------------------------
#     # Expected outcomes
#     # ---- Output type
#     expected_output_type = dict
#     # ---- Data types
#     expected_dtypes = {
#         'haul_mean_df': {
#             'species_id': np.integer ,
#             'haul_num': np.integer ,
#             'stratum_num': np.integer ,
#             'sigma_bs_mean': np.floating ,
#             'TS_mean': np.floating
#         } ,
#         'strata_mean_df': {
#             'stratum_num': np.integer ,
#             'species_id': np.integer ,
#             'sigma_bs_mean': np.floating ,
#             'TS_mean': np.floating
#         }
#     }
#     # Expected [ FULL ]
#     expected_full = {
#         'haul_mean_df': pd.DataFrame( {
#             'species_id': np.append( np.repeat( [ 19350 ] , 13 ) , 43130 ) ,
#             'haul_num': np.concatenate( [ np.repeat( [ 5 , 10 , 15 ] , 2 ) ,
#                                           np.repeat( [ 20 ] , 3 ) ,
#                                           np.repeat( [ 25 ] , 2 ) ,
#                                           np.repeat( [ 30 ] , 3 ) ] ) ,
#             'stratum_num': np.array( [ 1 , 2 , 1 , 2 , 1 , 3 , 1 ,
#                                        2 , 3 , 2 , 3 , 2 , 3 , 3 ] ) ,
#             'sigma_bs_mean': np.array( [ 0.01 , 0.059090909090909124 , 0.03006493506493507 ,
#                                          0.0672727272727273 , 0.026363636363636356 ,
#                                          0.07545454545454544 , 0.03454545454545455 ,
#                                          0.05895733652312599 , 0.08363636363636363 ,
#                                          0.042727272727272704 , 0.09181818181818176 ,
#                                          0.050909090909090904, 0.08504684247050656, 0.01 ] ) ,
#             'TS_mean': np.array( [-20.0, -12.284793285153693, -15.219397298185097,
#                                   -11.721609654272488, -15.78994687259269, -11.223145927821513,
#                                   -14.616090885414149, -12.294621445424438, -10.776048578126698,
#                                   -13.692948272225077, -10.370713113755826, -12.932046581520247,
#                                   -10.703418057797036, -20.0 ] ) } ) ,
#         'strata_mean_df': pd.DataFrame( {
#             'stratum_num': np.array( [ 1 , 2 , 3 , 3 ] ) ,
#             'species_id': np.append( np.repeat( 19350 , 3 ) , 43130 ) ,
#             'sigma_bs_mean': np.array( [0.025243506493506495, 0.05579146730462521,
#                                         0.08398898334489935, 0.01] ) ,
#             'TS_mean': np.array( [-15.978503188015692, -12.534322165799113,
#                                   -10.757776756796664, -20.0] )
#         } )
#     }
#     # Expected [ SUB ]
#     expected_sub = {
#         'haul_mean_df': pd.DataFrame( {
#             'species_id': np.repeat( [ 19350 ] , 13 ),
#             'haul_num': np.concatenate( [ np.repeat( [ 5 , 10 , 15 ] , 2 ) ,
#                                           np.repeat( [ 20 ] , 3 ) ,
#                                           np.repeat( [ 25 , 30 ] , 2 ) ] ) ,
#             'stratum_num': np.array( [ 1 , 2 , 1 , 2 , 1 , 3 , 1 ,
#                                        2 , 3 , 2 , 3 , 2 , 3 ] ) ,
#             'sigma_bs_mean': np.array( [ 0.01 , 0.059090909090909124 , 0.03006493506493507 ,
#                                          0.0672727272727273 , 0.026363636363636356 ,
#                                          0.07545454545454544 , 0.03454545454545455 ,
#                                          0.05895733652312599 , 0.08363636363636363 ,
#                                          0.042727272727272704 , 0.09181818181818176 ,
#                                          0.050909090909090904, 0.08504684247050656 ] ) ,
#             'TS_mean': np.array( [-20.0, -12.284793285153693, -15.219397298185097,
#                                   -11.721609654272488, -15.78994687259269, -11.223145927821513,
#                                   -14.616090885414149, -12.294621445424438, -10.776048578126698,
#                                   -13.692948272225077, -10.370713113755826, -12.932046581520247,
#                                   -10.703418057797036] ) } ) ,
#         'strata_mean_df': pd.DataFrame( {
#             'stratum_num': np.array( [ 1 , 2 , 3 ] ) ,
#             'species_id': np.repeat( 19350 , 3 ) ,
#             'sigma_bs_mean': np.array( [0.025243506493506495, 0.05579146730462521,
#                                         0.08398898334489935] ) ,
#             'TS_mean': np.array( [-15.978503188015692, -12.534322165799113,
#                                   -10.757776756796664] )
#         } )
#     }

#     # -------------------------
#     # Run tests [ FULL ]
#     # ---- Type [ OUTPUT: FULL ]
#     assert isinstance( test_results_full , expected_output_type )
#     # ---- Shape, types, values [ OUTPUT: FULL ]
#     assert_dataframe_equal( test_results_full , expected_dtypes , expected_full )
#     # Run tests [ SUB ]
#     # ---- Type [ OUTPUT: SUB ]
#     assert isinstance( test_results_sub , expected_output_type )
#     # ---- Shape, types, values [ OUTPUT: SUB ]
#     assert_dataframe_equal( test_results_sub , expected_dtypes , expected_sub )

# def test_nasc_to_biomass( ):

#     # -------------------------
#     # Mock values
#     # ---- Input dictionary
#     input_dict = {
#         'acoustics': {
#             'nasc_df': pd.DataFrame(
#                 { 'transect_num': np.array( [ 1 , 1 , 2 , 2 , 3 , 3 , 4 , 4 ] ) ,
#                   'vessel_log_start': np.linspace( 5.0 , 20.0 , 8 ) ,
#                   'vessel_log_end': np.linspace( 5.0 , 20.0 , 8 ) + 1.0 ,
#                   'latitude': np.linspace( 80.0 , 85.0 , 8 ) ,
#                   'longitude': np.linspace( -5.0 , 5.0 , 8 ) ,
#                   'stratum_num': np.array( [ 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 ] ) ,
#                   'transect_spacing': np.repeat( 5.0 , 8 ) ,
#                   'NASC_no_age1': np.array( [ 0.0 , 0.0 , 100.0 , 500.0 ,
#                                               1000.0 , 500.0 , 100.0 , 0.0 ] ) ,
#                   'haul_num': np.array( [ 5 , 10 , 15 , 15 , 20 , 25 , 30 , 30 ] ) ,
#                   'NASC_all_ages': np.array( [ 0.0 , 0.0 , 200.0 , 500.0 ,
#                                                2000.0 , 500.0 , 100.0 , 100.0 ] ) ,
#                   'stratum_inpfc': np.array( [ 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 ] ) } )
#         } ,
#         'biology': {
#             'distributions': {
#                 'length_bins_df': pd.DataFrame(
#                     { 'length_bins': np.array( [ 7.5 , 12.5 ] ) ,
#                       'length_intervals': np.array( [ pd.Interval( left = 5.0 , right = 10.0 ) ,
#                                                       pd.Interval( left = 10.0 , right = 15.0)])})
#             }
#         }
#     }
#     # ---- Analysis dictionary
#     analysis_dict = {
#         'acoustics': {
#             'sigma_bs': {
#                 'strata_mean_df': pd.DataFrame(
#                     { 'stratum_num': np.array( [ 1 , 2 ] ) ,
#                       'species_id': np.repeat( [ 9933 ] , 2 ) ,
#                       'sigma_bs_mean': np.array( [ 1e-5 , 1e-6 ] ) ,
#                       'TS_mean': np.array( [ -50.0 , -60.0 ] ) } )
#             }
#         } ,
#         'biology': {
#             'distributions': {

#             } ,
#             'proportions': {
#                 'number': {
#                     'aged_length_proportions_df': pd.DataFrame(
#                         { 'stratum_num': np.repeat( [ 1 , 2 ] , 12 ) ,
#                           'species_id': np.repeat( [ 9933 ] , 24 ) ,
#                           'length_bin': np.tile( [ pd.Interval( left = 5.0 , right = 10.0 ) ,
#                                                    pd.Interval(left =10.0 ,right = 15.0 ) ] ,12) ,
#                           'age_bin': np.tile( [ pd.Interval( left = 0.5 , right = 1.5 ) ,
#                                                 pd.Interval( left = 1.5 , right = 2.5 ) ] , 12 ) ,
#                           'sex': np.tile( [ 'all' , 'all' , 'all' , 'all' ,
#                                             'male' , 'male' , 'male' , 'male' ,
#                                             'female' , 'female' , 'female' , 'female' ] , 2 ) ,
#                           'proportion_number_aged': np.array( [ 0.00 , 0.00 , 0.50 , 0.50 ,
#                                                                 0.00 , 0.00 , 0.50 , 0.50 ,
#                                                                 0.00 , 0.00 , 0.50 , 0.50 ,
#                                                                 0.25 , 0.25 , 0.25 , 0.25 ,
#                                                                 0.25 , 0.25 , 0.25 , 0.25 ,
#                                                                 0.25 , 0.25 , 0.25 , 0.25 ] ) } ),
#                     'unaged_length_proportions_df': pd.DataFrame(
#                         { 'stratum_num': np.repeat( [ 1 , 2 ] , 6 ) ,
#                           'species_id': np.repeat( [ 9933 ] , 12 ) ,
#                           'length_bin': np.tile( [ pd.Interval( left = 5.0 , right = 10.0 ) ,
#                                                    pd.Interval( left =10.0,right=15.0 )  ] , 6 ) ,
#                           'sex': np.tile( [ 'all' ,'all','male','male','female','female'] , 2 ) ,
#                           'propotion_number_aged': np.array( [ 0.00 , 0.00 , 0.25 , 0.25  ,
#                                                                0.00 , 0.00 , 0.50 , 0.50 ,
#                                                                0.00 , 0.00 , 0.50 , 0.50 ,
#                                                                0.25 , 0.25 , 0.25 , 0.25 ,
#                                                                0.25 , 0.25 , 0.25 , 0.25 ,
#                                                                0.25 , 0.25 , 0.25 , 0.25 ] ) } )
#                 }
#             }
#             'weight': {
#                 'weight_stratum_df': pd.DataFrame(
#                     { 'stratum_num': np.tile( [ 1 , 2 , 3 ] , 3 ) ,
#                       'sex': np.repeat( [ 'all' , 'male' , 'female' ] , 3 ) ,
#                       'average_weight': np.array( [ 1.0 , 2.0 , 3.0 ,
#                                                     0.5 , 1.5 , 2.5 ,
#                                                     1.5 , 2.5 , 3.5 ] ) }
#                 )
#             }
#         }
#     }
#     # ---- Configuration dictionary
#     mock_configuration_dict = {
#         'TS_length_regression_parameters': {
#             'pacific_hake': {
#                 'number_code': 9933 ,
#                 'TS_L_slope': 10.0 ,
#                 'TS_L_intercept': -40.0 ,
#                 'length_units': 'cm'
#             }
#         }
#     }
#     # ---- Settings dictionary
#     settings_dict = {
#         'transect': {
#             'stratum_name': 'stratum_num' ,
#             'exclude_age1': True
#         }
#     }
#     # ---- Length data (unaged)
#     mock_length_data = pd.DataFrame( {
#             'stratum_num': np.repeat( [ 1 , 2 , 3 ] , 4 ) ,
#             'haul_num': np.repeat( [ 10 , 20 , 30 ] , 4 ) ,
#             'species_id': np.append( np.repeat( [ 19350 ] , 11 ) , 43130 ) ,
#             'length': np.linspace( 15 , 100 , 12 ) ,
#             'length_count': np.linspace( 10 , 105 , 12 ).astype( int ) ,
#             'group_sex': np.repeat( 'sexed' , 12 ) } )
#     # -------- > DROP SPECIES == 43130
#     mock_length_data_sub = mock_length_data[ mock_length_data[ 'species_id' ] != 43130 ]
#     # ---- Specimen data (aged)
#     mock_specimen_data = pd.DataFrame(
#         {
#             'stratum_num': np.repeat( [ 1 , 2 , 3 ] , 4 ) ,
#             'haul_num': np.tile( [ 5 , 10 , 15 , 20 , 25 , 30 ] , 2 ) ,
#             'species_id': np.append( np.repeat( [ 19350 ] , 11 ) , 43130 ) ,
#             'length': np.linspace( 10 , 100 , 12 ) ,
#             'group_sex': np.repeat( 'sexed' , 12 )
#         }
#     )
#     # -------- > DROP SPECIES == 43130
#     mock_specimen_data_sub = mock_specimen_data[ mock_specimen_data[ 'species_id' ] != 43130 ]
#     # -------- > DROP SPECIES == 43130
#     # ---- Strata definition dictionary
#     mock_strata_dict = {
#         'strata_df': pd.DataFrame( { 'stratum_num': [ 1 , 2 , 3 ] } )
#     }
#     # ---- Configuration dictionary
#     mock_configuration_dict = {
#         'TS_length_regression_parameters': {
#             'mega_guppies': {
#                 'number_code': 19350 ,
#                 'TS_L_slope': 10.0 ,
#                 'TS_L_intercept': -30.0 ,
#                 'length_units': 'cm'
#             } ,
#             'flying_unicorns': {
#                 'number_code': 43130 ,
#                 'TS_L_slope': 10.0 ,
#                 'TS_L_intercept': -40.0 ,
#                 'length_units': 'cm'
#             }
#         }
#     }
#     # ---- Settings dictionary
#     mock_settings_dict = {
#         'transect': {
#             'stratum_name': 'stratum_num'
#         }
#     }
