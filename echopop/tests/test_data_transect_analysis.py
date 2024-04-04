import pandas as pd
import numpy as np
from echopop.survey import Survey

def test_fit_binned_length_weight_relationship( mock_survey ):
    
    #### Pull in mock Survey object
    objS = mock_survey
    
    ### Initialize objS for `length_weight`
    objS.statistics[ 'length_weight' ] = { }
    
    ### Re-parameterize `specimen_df` with dummy data 
    objS.biology[ 'specimen_df' ] = pd.DataFrame(
        {
            'stratum_num': [ 0 , 0 , 1 , 1 , 2 , 2 , 3 , 3 ] ,
            'haul_num': [ 1 , 1 , 2 , 2 , 3 , 3 , 4 , 4 ] ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'group': np.repeat( 'sexed' , 8 ) ,
            'species_id': np.repeat( [ 8675309 ] , 8 ) ,
            'length': [ 2.0 , 3.0 , 4.0 , 5.0 , 6.0 , 7.0 , 8.0 , 9.0  ] ,
            'weight': [ 4.0 , 9.0 , 16.0 , 25.0 , 36.0 , 49.0 , 64.0 , 81.0 ] ,
        }
    )
    
    ### Re-parameterize `length_bins` with dummy data
    objS.biology[ 'distributions' ][ 'length' ][ 'length_bins_arr' ] = (
        [ 2.0 , 5.0 , 8.0 , 11.0 ]
    )
    
    ### Re-parameterize `length_interval` with dummy data
    objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] = (
        [ 0.5 , 3.5 , 6.5 , 9.5 , 12.5 ]
    )
    
    ### Evaluate object for later comparison 
    objS.fit_binned_length_weight_relationship( species_id = 8675309 )
    
    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- `objS.statistics[ 'length_weight' ][ 'regression_parameters' ]`
    # ---- Expected dimensions
    expected_dimensions_regression_parameters = tuple( [ 3 , 3 ] )
    # ---- Expected output
    expected_output_regression_parameters = pd.DataFrame(
        {
            'sex': [ 'all' , 'female' , 'male' ] ,
            'rate': [ 2.0 , 2.0 , 2.0 ] ,
            'initial': [ 4.7e-16 , -2.2e-16 , 1.1e-15 ]
        }
    )
    # ---- `objS.statistics[ 'length_weight' ][ 'length_weight_df' ]`
    # ---- Expected dimensions
    expected_dimensions_length_weight_df = tuple( [ 12 , 10 ] )
    # ---- Expected output
    expected_output_length_weight_df = pd.DataFrame(
        {
            'length_bin': pd.cut( np.repeat( [ 1 , 4 , 7 , 10 ] , 3 ) ,
                                  np.array( [ 0.5 , 3.5 , 6.5 , 9.5 , 12.5 ] ) ) ,
            'sex': np.tile( [ 'all' , 'female' , 'male' ] , 4 ) ,
            'mean_length': [ 2.5 , 3.0 , 2.0 , 5.0 , 5.0 , 5.0 ,
                             8.0 , 8.0 , 8.0 , 0.0 , 0.0 , 0.0 ] ,
            'n_length': [ 2 , 1 , 1 , 3 , 1 , 2 ,
                          3 , 2 , 1 , 0 , 0 , 0 ] ,
            'mean_weight': [ 6.50 , 9.00 , 4.00 , 25.67 , 25.00 , 26.00 ,
                             64.67 , 65.00 , 64.00 , 0.00 , 0.00 , 0.00 ] ,
            'n_weight': [ 2 , 1 , 1 , 3 , 1 , 2 ,
                          3 , 2 , 1 , 0 , 0 , 0 ] ,
            'rate': np.repeat( 2.0 , 12 ) ,
            'initial': np.tile( [ 4.7e-16 , -2.2e-16 , 1.1e-15 ] , 4 ) ,
            'weight_fitted': [ 4.0 , 4.0 , 4.0 , 25.0 , 25.0 , 25.0 ,
                               64.0 , 64.0 , 64.0 , 121.0 , 121.0 , 121.0 ] ,
            'weight_modeled': [ 4.0 , 4.0 , 4.0 , 25.0 , 25.0 , 25.0 ,
                                64.0 , 64.0 , 64.0 , 121.0 , 121.0 , 121.0 ]
        }
    )
    expected_output_length_weight_df[ 'length_bin' ] = pd.IntervalIndex( expected_output_length_weight_df[ 'length_bin' ] )
    expected_output_length_weight_df[ 'length_bin' ] = pd.Categorical( expected_output_length_weight_df[ 'length_bin' ] , 
                                                                       categories = expected_output_length_weight_df[ 'length_bin' ].unique( ) , 
                                                                       ordered = True )
    #----------------------------------
    ### Run tests: `fit_binned_length_weight_relationship`
    #----------------------------------
    eval_regression_parameters = objS.statistics[ 'length_weight' ][ 'regression_parameters' ]
    eval_length_weight_df = objS.statistics[ 'length_weight' ][ 'length_weight_df' ]
    ### Check shape 
    assert eval_regression_parameters.shape == expected_dimensions_regression_parameters
    assert eval_length_weight_df.shape == expected_dimensions_length_weight_df
    ### Check datatypes
    assert np.all( eval_regression_parameters.dtypes == expected_output_regression_parameters.dtypes )
    assert np.all( eval_length_weight_df.dtypes == expected_output_length_weight_df.dtypes )
    ### Dataframe equality
    assert np.allclose( eval_regression_parameters[ [ 'rate' , 'initial' ] ] , 
                        expected_output_regression_parameters[ [ 'rate' , 'initial' ] ] ,
                        rtol = 1e-1 )
    # ---- Non-float/high-precision
    assert eval_length_weight_df[ [ 'length_bin' , 'sex' , 'mean_length' , 'n_length' , 'n_weight' ] ].equals(
        expected_output_length_weight_df[ [ 'length_bin' , 'sex' , 'mean_length' , 'n_length' , 'n_weight' ] ]
    )
    # ---- Float/high-precision
    assert np.allclose( eval_length_weight_df[ [ 'mean_weight' , 'rate' , 'initial' ,'weight_fitted' , 'weight_modeled' ] ] ,
                        expected_output_length_weight_df[ [ 'mean_weight' , 'rate' , 'initial' ,'weight_fitted' , 'weight_modeled' ] ] ,
                        rtol = 1e-1 )
    
def test_strata_sex_weight_proportions( mock_survey ):

    #### Pull in mock Survey object
    objS = mock_survey
    
    ### Initialize objS for `weight`
    objS.biology[ 'weight' ] = { }
    
    ### Initialize objS for `length_weight`
    objS.statistics[ 'length_weight' ] = { }
    
    ### Re-parameterize `specimen_df` with dummy data 
    objS.biology[ 'specimen_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'group': np.repeat( 'sexed' , 8 ) ,
            'haul_num': np.tile( [ 1 , 2 ] , 4 ) ,
            'species_id': np.repeat( [ 8675309 ] , 8 ) ,
            'length': [ 12.0 , 12.0 , 19.0 , 19.0 , 12.0 , 12.0 , 19.0 , 19.0 ] ,
            'weight': [ 2.0 , 3.0 , 3.0 , 2.0 , 2.0 , 3.0 , 2.0 , 3.0 ] ,
            'age': [ 1 , 1 , 2 , 2 , 1 , 1 , 2 , 2 ]    
        }
    )    
    
    ### Re-parameterize `length_df` with dummy data 
    objS.biology[ 'length_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'group': np.repeat( 'sexed' , 8 ) ,
            'species_id': np.repeat( [ 8675309 ] , 8 ) ,
            'length': [ 12 , 12 , 19 , 19 , 12 , 12 , 19 , 19 ] ,
            'length_count': [ 5 , 10 , 15 , 20 , 20 , 15 , 10 , 5 ]
        }
    )
    
    ### Re-parameterize `fitted_weight` with dummy data 
    objS.statistics[ 'length_weight' ][ 'length_weight_df' ] = pd.DataFrame(
        {
            'length_bin': pd.cut( np.repeat( [ 12 , 18 ] , 3 ) ,
                                  np.linspace( 9 , 21 , 3 ) ) ,
            'sex': np.repeat( [ 'all' , 'female' , 'male' ] , 2 ) ,
            'n_length': [ 4 , 2 , 2 , 4 , 2 , 2 ] ,
            'mean_weight': [ 2.5 , 3.5 , 1.5 , 7.5 , 6.5 , 8.5 ] ,
            'n_weight': [ 4 , 2 , 2 , 4 , 2 , 2 ] ,
            'rate': [ 2.63 , 1.36 , 3.90 , 2.63 , 1.36 , 3.90 ] ,
            'initial': [ -2.49 , -0.93 , -4.06 , -2.49 , -0.93 , -4.06 ] ,
            'weight_fitted': [ 2.21 , 3.46 , 1.41 , 6.43 , 6.02 , 6.87 ] ,
            'weight_modeled': [ 2.21 , 3.46 , 1.41 , 6.43 , 6.02 , 6.87 ]
        }
    )


    ### Re-parameterize `length_df` with dummy data 
    objS.biology[ 'length_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'group': np.repeat( 'sexed' , 8 ) ,
            'species_id': np.repeat( [ 8675309 ] , 8 ) ,
            'length': [ 12 , 12 , 19 , 19 , 12 , 12 , 19 , 19 ] ,
            'length_count': [ 5 , 10 , 15 , 20 , 20 , 15 , 10 , 5 ]
        }
    )
    
    ### Re-parameterize `fitted_weight` with dummy data 
    objS.statistics[ 'length_weight' ][ 'length_weight_df' ] = pd.DataFrame(
        {
            'length_bin': pd.cut( np.repeat( [ 12 , 18 ] , 3 ) ,
                                  np.linspace( 9 , 21 , 3 ) ) ,
            'sex': np.repeat( [ 'all' , 'female' , 'male' ] , 2 ) ,
            'n_length': [ 4 , 2 , 2 , 4 , 2 , 2 ] ,
            'mean_weight': [ 2.5 , 3.5 , 1.5 , 7.5 , 6.5 , 8.5 ] ,
            'n_weight': [ 4 , 2 , 2 , 4 , 2 , 2 ] ,
            'rate': [ 2.63 , 1.36 , 3.90 , 2.63 , 1.36 , 3.90 ] ,
            'initial': [ -2.49 , -0.93 , -4.06 , -2.49 , -0.93 , -4.06 ] ,
            'weight_fitted': [ 2.21 , 3.46 , 1.41 , 6.43 , 6.02 , 6.87 ] ,
            'weight_modeled': [ 2.21 , 3.46 , 1.41 , 6.43 , 6.02 , 6.87 ]
        }
    )

    ### Re-parameterize `length_bins` with dummy data
    objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] = np.linspace( 9 , 21 , 3 )

    ### Evaluate object for later comparison 
    objS.strata_sex_weight_proportions( species_id = 8675309 )
    
    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions = tuple( [ 2 , 8 ] )
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            'stratum_num': np.array( [ 0 , 1 ] ).astype( np.int32 ) ,
            'proportion_female': [ 0.59 , 0.41 ] ,
            'proportion_male': [ 0.41 , 0.59 ] ,
            'proportion_station_1': [ 0.93 , 0.93 ] ,
            'proportion_station_2': [ 0.07 , 0.07 ] ,
            'average_weight_female': [ 4.72 , 2.71 ] ,
            'average_weight_male': [ 6.64 , 6.30 ] ,
            'average_weight_total': [ 3.07 , 2.60 ]
        }
    )
    #----------------------------------
    ### Run tests: `strata_sex_weight_proportions`
    #----------------------------------
    eval_weight_strata_df = objS.biology[ 'weight' ][ 'weight_strata_df' ]
    ### Check shape 
    assert eval_weight_strata_df.shape == expected_dimensions
    ### Check datatypes
    assert np.all( eval_weight_strata_df.dtypes == expected_output.dtypes )
    ### Dataframe equality
    assert np.allclose( eval_weight_strata_df , expected_output , rtol = 1e-1 )
    
def test_strata_age_binned_weight_proportions( mock_survey ):
    
    #### Pull in mock Survey object
    objS = mock_survey
    
    ### Initialize objS for `weight`
    objS.biology[ 'weight' ] = { }
    
    ### Re-parameterize `specimen_df` with dummy data 
    objS.biology[ 'specimen_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'group': np.repeat( 'sexed' , 8 ) ,
            'haul_num': np.tile( [ 1 , 2 ] , 4 ) ,
            'species_id': np.repeat( [ 8675309 ] , 8 ) ,
            'length': [ 12.0 , 12.0 , 19.0 , 19.0 , 12.0 , 12.0 , 19.0 , 19.0 ] ,
            'weight': [ 2.0 , 3.0 , 3.0 , 2.0 , 2.0 , 3.0 , 2.0 , 3.0 ] ,
            'age': [ 1 , 1 , 2 , 2 , 1 , 1 , 2 , 2 ]    
        }
    )

    ### Re-parameterize `length_bins` with dummy data
    objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] = np.linspace( 9 , 21 , 3 )

    ### Evaluate object for later comparison 
    objS.strata_age_binned_weight_proportions( species_id = 8675309 )
    
    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions = { 'age_proportions': tuple( [ 4 , 4 ] ) ,
                            'age_weight_proportions': tuple( [ 4 , 4 ] ) ,
                            'sex_age_weight_proportions': tuple( [ 12 , 5 ] ) ,
                            'length_sex_age_weight_proportions': tuple( [ 24 , 9 ] ) }
    
    # ---- Expected output
    expected_output = {
        'age_proportions': pd.DataFrame( {
           'stratum_num': np.repeat( [ 0 , 1 ] , 2 ).astype( np.int64 ) ,
           'age': np.tile( [ 1 , 2 ] , 2 ).astype( np.int64 ) ,
           'count_age_proportion_all': np.repeat( 0.5 , 4 ) ,
           'count_age_proportion_adult': [ 0.0 , 1.0 , 0.0 , 1.0 ]
        } ) ,
        'age_weight_proportions': pd.DataFrame( {
           'stratum_num': np.repeat( [ 0 , 1 ] , 2 ).astype( np.int64 ) ,
           'age': np.tile( [ 1 , 2 ] , 2 ).astype( np.int64 ) ,
           'weight_age_proportion_all': [ 0.50 , 0.50 , 0.50 , 0.50 ] ,
           'weight_age_proportion_adult': [ 0.0 , 1.0 , 0.0 , 1.0 ]
        } ) ,
        'sex_age_weight_proportions': pd.DataFrame( {
            'stratum_num': np.repeat( [ 0 , 1 ] , 6 ).astype( np.int64 ) ,
            'age': np.tile( [ 1 , 1 , 1 , 2 , 2 , 2 ] , 2 ).astype( np.int64 ) ,
            'sex': np.tile( [ 'all' , 'female' , 'male' ] , 4 ) ,
            'weight_sex_proportion_all': [ 0.5 , 0.6 , 0.4 , 0.5 , 0.4 , 0.6 ,
                                           0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ] ,
            'weight_sex_proportion_adult': np.tile( [ 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 1.0 ] , 2 ) 
        } ) ,
        'length_sex_age_weight_proportions': pd.DataFrame( {
            'stratum_num': np.repeat( [ 0 , 1 ] , 12 ).astype( np.int64 ) ,
            'age': np.tile( [ 1 , 1 , 1 , 1 , 1 , 1 , 
                              2 , 2 , 2 , 2 , 2 , 2 ] , 2 ).astype( np.int64 ) ,
            'length_bin': pd.cut( np.tile( [ 12.0 , 12.0 , 12.0 , 18.0 , 18.0 , 18.0 ] , 4 ) ,
                                  np.linspace( 9 , 21 , 3 ) ) ,
            'sex': np.tile( [ 'all' , 'female' , 'male' ] , 8 ) ,
            'count': [ 5.0 , 3.0 , 2.0 , 0.0 , 0.0 , 0.0 ,
                       0.0 , 0.0 , 0.0 , 5.0 , 2.0 , 3.0 ,
                       5.0 , 3.0 , 2.0 , 0.0 , 0.0 , 0.0 ,
                       0.0 , 0.0 , 0.0 , 5.0 , 3.0 , 2.0 ] ,
            'weight_total_all': [ 10.0 , 5.0 , 5.0 , 10.0 , 5.0 , 5.0 ,
                                  10.0 , 5.0 , 5.0 , 10.0 , 5.0 , 5.0 ,
                                  10.0 , 6.0 , 4.0 , 10.0 , 6.0 , 4.0 ,
                                  10.0 , 6.0 , 4.0 , 10.0 , 6.0 , 4.0 ] ,
            'weight_total_adult': [ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
                                    5.0 , 2.0 , 3.0 , 5.0 , 2.0 , 3.0 ,
                                    0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
                                    5.0 , 3.0 , 2.0 , 5.0 , 3.0 , 2.0 ] ,
            'weight_length_sex_proportion_all': [ 0.5 , 0.6 , 0.4 , 0.0 , 0.0 , 0.0 ,
                                                  0.0 , 0.0 , 0.0 , 0.5 , 0.4 , 0.6 ,
                                                  0.5 , 0.5 , 0.5 , 0.0 , 0.0 , 0.0 ,
                                                  0.0 , 0.0 , 0.0 , 0.5 , 0.5 , 0.5 ] ,
            'weight_length_sex_proportion_adult': np.tile( [ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
                                                             0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 1.0 ] , 2 )            
        } )
    }

    #----------------------------------
    ### Run tests: `strata_age_binned_weight_proportions`
    #----------------------------------
    eval_age_proportions_df = objS.biology[ 'weight' ][ 'proportions' ][ 'age_proportions_df' ]
    eval_age_weight_proportions_df = objS.biology[ 'weight' ][ 'proportions' ][ 'age_weight_proportions_df' ]
    eval_sex_age_weight_proportions_df = objS.biology[ 'weight' ][ 'proportions' ][ 'sex_age_weight_proportions_df' ]
    eval_length_sex_age_weight_proportions_df = objS.biology[ 'weight' ][ 'proportions' ][ 'length_sex_age_weight_proportions_df' ]
    ### Check shape 
    assert eval_age_proportions_df.shape == expected_dimensions[ 'age_proportions' ]
    assert eval_age_weight_proportions_df.shape == expected_dimensions[ 'age_weight_proportions' ]
    assert eval_sex_age_weight_proportions_df.shape == expected_dimensions[ 'sex_age_weight_proportions' ]
    assert eval_length_sex_age_weight_proportions_df.shape == expected_dimensions[ 'length_sex_age_weight_proportions' ]
    ### Check datatypes
    assert np.all( eval_age_proportions_df.dtypes == expected_output[ 'age_proportions' ].dtypes )
    assert np.all( eval_age_weight_proportions_df.dtypes == expected_output[ 'age_weight_proportions' ].dtypes )
    assert np.all( eval_sex_age_weight_proportions_df.dtypes == expected_output[ 'sex_age_weight_proportions' ].dtypes )
    assert np.all( eval_length_sex_age_weight_proportions_df.dtypes == expected_output[ 'length_sex_age_weight_proportions' ].dtypes )
    ### Dataframe equality
    assert eval_age_proportions_df.equals( expected_output[ 'age_proportions' ] )    
    assert eval_age_weight_proportions_df.equals( expected_output[ 'age_weight_proportions' ] )    
    assert eval_sex_age_weight_proportions_df.equals( expected_output[ 'sex_age_weight_proportions' ] )
    assert eval_length_sex_age_weight_proportions_df.equals( expected_output[ 'length_sex_age_weight_proportions' ] )