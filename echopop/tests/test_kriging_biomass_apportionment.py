import pandas as pd
import numpy as np
from echopop.computation.biology import sum_strata_weight , compute_index_aged_weight_proportions

def test_sum_strata_weight( mock_survey ):

    #### Pull in mock Survey object
    objS = mock_survey

    ### Re-parameterize `specimen_df` with dummy data 
    objS.biology[ 'specimen_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 , 2 , 4 , 5 ] , 4 ) ,
            'haul_num': np.repeat( [ 1 , 2 , 3 , 4 , 5 ] , 4 ) ,
            'species_id': np.repeat( [ 19350 ] , 20 ) ,
            'length': np.linspace( 10 , 100 , 20 ) ,
            'weight': np.linspace( 1 , 5 , 20 ) ,     
        }
    )

    ### Re-parameterize `length_df` with dummy data
    objS.biology[ 'length_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 , 2 , 4 , 5 ] , 4 ) ,
            'haul_num': np.repeat( [ 1 , 2 , 3 , 4 , 5 ] , 4 ) ,
            'species_id': np.repeat( [ 19350 ] , 20 ) ,
            'length': np.linspace( 10 , 100 , 20 ) ,
            'length_count': np.linspace( 10 , 100 , 20 ) ,     
        }
    )

    ### Re-parameterize `catch_df` with dummy data
    objS.biology[ 'catch_df' ] = pd.DataFrame(
        {
            'stratum_num': [ 0 , 1 , 2 , 4 , 5 ] ,
            'haul_num': [ 1 , 2 , 3 , 4 , 5 ] ,
            'haul_weight': [ 51.4 , 0.6 , 81.7 , 16.2 , 12.9 ] ,
        }
    )

    ### Evaluate object for later comparison 
    object_weight_strata_aged_unaged , object_weight_strata = sum_strata_weight( objS.biology[ 'catch_df' ] ,
                                                                                objS.biology[ 'specimen_df' ] )

    #----------------------------------
    ### Run tests: `sum_strata_weight`
    #----------------------------------
    ### Evaluate shape
    # ---- `object_weight_strata`
    assert object_weight_strata.shape == tuple( [ 5 , 2 ] )
    # ---- `object_weight_strata_aged_unaged`
    assert object_weight_strata_aged_unaged.shape == tuple( [ 10 , 3 ] )

    ### Evaluate value equality
    # ---- `object_weight_strata`
    check_values = np.array( [ 56.663158 , 9.231579 , 93.700000 , 31.568421 , 31.636842 ] )
    assert np.allclose( check_values , object_weight_strata.weight_stratum_all )

    # ---- `object_weight_strata_aged_unaged`
    check_values = np.array( [ 51.400000 , 0.600000 , 81.700000 , 16.200000 , 12.900000 ,
                            5.263158 , 8.631579 , 12.000000 , 15.368421 , 18.736842  ] )
    assert np.allclose( check_values , object_weight_strata_aged_unaged.stratum_weight )

def test_compute_index_aged_weight_proportions( mock_survey ):
    #### Pull in mock Survey object
    objS = mock_survey

    ### Re-parameterize `specimen_df` with dummy data 
    objS.biology[ 'specimen_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'haul_num': np.tile( [ 1 , 2 ] , 4 ) ,
            'species_id': np.repeat( [ 19350 ] , 8 ) ,
            'length': [ 10 , 20 , 19 , 11 , 12 , 18 , 17 , 13 ] ,
            'weight': np.linspace( 1 , 5 , 8 ) ,
            'age': [ 1 , 2 , 1 , 2 , 1 , 2 , 1 , 2 ]    
        }
    )

    ### Length interval
    objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] = np.linspace( 9 , 21 , 3 )
    
    ### Age interval
    objS.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ] = np.array( [ 0.5 , 1.5 , 2.5 ] )

    ### Evaluate object for later comparison 
    obj_props_wgt_len_age_sex = compute_index_aged_weight_proportions( objS.biology[ 'specimen_df' ] ,
                                                                       objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] ,
                                                                       objS.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ] )
    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions of `obj_props_wgt_len_age_sex`   
    expected_dimensions = tuple( [ 16 , 11 ] )
    
    # ---- Expected dataframe output
    expected_output = pd.DataFrame( {
        'stratum_num': np.repeat( [ 0 , 1 ] , 8 ).astype( np.int64 ) ,
        'species_id': np.repeat( 19350 , 16 ).astype( np.int64 ) ,
        'sex': np.tile( [ 'female' , 'male' ] , 8 ).astype( object ) ,
        'length_bin': pd.IntervalIndex.from_arrays( np.tile( [9.0 , 15.0 ] , 8 ) , 
                                                    np.tile( [ 15.0 , 21.0 ] , 8 ) , 
                                                    closed = 'right' ) ,
        'age_bin': pd.IntervalIndex.from_arrays( np.tile( [ 0.5, 1.5 ] , 8 ) , 
                                                 np.tile( [ 1.5 , 2.5 ] , 8 ) , 
                                                 closed = 'right' ) ,
        'weight_all': [ 0.0, 2.714286, 0.0, 1.571429, 1.0, 0.0, 2.142857, 0.0, 0.0, 
                        5.0, 0.0, 3.857143, 3.285714, 0.0, 4.428571, 0.0 ] ,
        'weight_adult': [ 0.0, 2.714286, 0.0, 1.571429, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 
                          0.0, 3.857143, 0.0, 0.0, 0.0, 0.0 ] ,
        'total_weight_sex_all': [ 4.285714, 4.285714, 4.285714, 4.285714, 3.142857, 
                                  3.142857, 3.142857, 3.142857, 8.857143, 8.857143, 
                                  8.857143, 8.857143, 7.714286, 7.714286, 7.714286, 
                                  7.714286 ] ,
        'total_weight_sex_adult': [ 4.285714, 4.285714, 4.285714, 4.285714, 0.0, 0.0, 
                                   0.0, 0.0, 8.857143, 8.857143, 8.857143, 8.857143, 
                                   0.0, 0.0, 0.0, 0.0 ] ,
        'proportion_weight_sex_all': [ 0.0, 0.633333, 0.0, 0.366667, 0.318182, 0.0, 
                                      0.681818, 0.0, 0.0, 0.564516, 0.0, 0.435484, 
                                      0.425926, 0.0, 0.574074, 0.0 ] ,
        'proportion_weight_sex_adult': [ 0.0, 0.633333, 0.0, 0.366667, 0.0, 0.0, 0.0, 
                                        0.0, 0.0, 0.564516, 0.0, 0.435484, 0.0, 0.0, 
                                        0.0, 0.0 ]
    } )
    expected_output[ 'length_bin' ] = pd.IntervalIndex( expected_output[ 'length_bin' ] )
    expected_output[ 'length_bin' ] = pd.Categorical( expected_output[ 'length_bin' ] , 
                                                      categories =  expected_output[ 'length_bin' ].unique( ) , 
                                                      ordered = True )
    expected_output[ 'age_bin' ] = pd.IntervalIndex( expected_output[ 'age_bin' ] )
    expected_output[ 'age_bin' ] = pd.Categorical( expected_output[ 'age_bin' ] , 
                                                 categories = expected_output[ 'age_bin' ].unique( ) , 
                                                 ordered=True)

    #----------------------------------
    ### Run tests: `compute_index_aged_weight_proportions`
    #----------------------------------
    ### Process the specimen data 
    ### Check shape 
    assert obj_props_wgt_len_age_sex.shape == expected_dimensions

    ### Check datatypes
    assert np.all( obj_props_wgt_len_age_sex.dtypes == expected_output.dtypes )
    
    ### Check data value equality
    assert np.all( obj_props_wgt_len_age_sex.dtypes == expected_output.dtypes )

def test_compute_summed_aged_proportions( mock_survey ):
    #### Pull in mock Survey object
    objS = mock_survey

    ### Mock data for `proportions_weight_length_age_sex`
    test_proportions_weight_length_age_sex = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0.0 , 1.0 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'weight_all': [ 1.0 , 3.0 , 2.0 , 4.0 , 0.0 , 6.0 , 3.0 , 3.0 ] ,
            'weight_adult': [ 0.0 , 1.5, 2.0 , 4.0 , 0.0 , 4.0 , 3.0 , 3.0 ]
        }
    )

    ### Mock data for `weight_strata`
    test_weight_strata = pd.DataFrame(
        {
            'stratum_num': [ 0.0 , 1.0 ] ,
            'weight_stratum_all': [ 10.0 , 100.0 ]
        }
    )
    
    ### Evaluate for later comparison 
    eval_aged_sex_proportions , eval_aged_proportions = compute_summed_aged_proportions( test_proportions_weight_length_age_sex ,
                                                                                         test_weight_strata )
    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    ### `eval_aged_sex_proportions`
    # ---- Expected dimensions
    expected_dimensions_aged_sex_proportions = tuple( [ 4 , 4 ] )
    # ---- Expected dataframe output
    expected_output_aged_sex_proportions = pd.DataFrame( {
        'stratum_num': np.repeat( [ 0.0 , 1.0 ] , 2 ) ,
        'sex': np.tile( [ 'female' , 'male' ] , 2 ) ,
        'proportion_weight_all': [ 0.70 , 0.30 , 0.09 , 0.03 ] ,
        'proportion_weight_adult': [ 0.55 , 0.20 , 0.07 , 0.03 ]
    } )

    ### `eval_aged_proportions`
    # ---- Expected dimensions
    expected_dimensions_aged_proportions = tuple( [ 2 , 3 ] )
    # ---- Expected dataframe output
    expected_output_aged_proportions = pd.DataFrame( {
        'stratum_num': [ 0.0 , 1.0 ] ,
        'proportion_weight_all': [ 1.00 , 0.12 ] ,
        'proportion_weight_adult': [ 0.75 , 0.10 ]
    } )

    #----------------------------------
    ### Run tests: `compute_index_aged_weight_proportions`
    #----------------------------------
    ### `eval_aged_sex_proportions`
    # ---- Shape
    assert eval_aged_sex_proportions.shape == expected_dimensions_aged_sex_proportions
    # ---- Datatypes
    assert np.all( eval_aged_sex_proportions.dtypes == expected_output_aged_sex_proportions.dtypes )
    # ---- Dataframe equality
    assert np.all( eval_aged_sex_proportions == expected_output_aged_sex_proportions )

    ### `eval_aged_proportions`
    # ---- Shape
    assert eval_aged_proportions.shape == expected_dimensions_aged_proportions
    # ---- Datatypes
    assert np.all( eval_aged_proportions.dtypes == expected_output_aged_proportions.dtypes )
    # ---- Dataframe equality
    assert np.all( eval_aged_proportions == expected_output_aged_proportions )