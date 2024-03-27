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

def test_compute_summed_aged_proportions( mock_survey ):
    #### Pull in mock Survey object
    objS = mock_survey

    ### Re-parameterize `specimen_df` with dummy data 
    objS.biology[ 'specimen_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 , 2 , 4 , 5 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 10 ) ,
            'haul_num': np.repeat( [ 1 , 2 , 3 , 4 , 5 ] , 4 ) ,
            'species_id': np.repeat( [ 19350 ] , 20 ) ,
            'length': np.linspace( 10 , 20 , 20 ) ,
            'weight': np.linspace( 1 , 5 , 20 ) ,
            'age': [ 1 , 1 , 1 , 2 , 2 , 2 , 2 , 3 , 3 , 3 ,
                     4 , 4 , 5 , 5 , 5 , 6 , 6 , 6 , 7 , 8 ]     
        }
    )

    ### Length interval
    objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] = np.linspace( 9 , 21 , 5 )
    
    ### Age interval
    objS.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ] = (
        np.array( [ 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5 ] )
    ) 

    ### Evaluate object for later comparison 
    obj_props_wgt_len_age_sex = compute_index_aged_weight_proportions( objS.biology[ 'specimen_df' ] ,
                                                                       objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] ,
                                                                       objS.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ] )
    
    #----------------------------------
    ### Run tests: `compute_index_aged_weight_proportions`
    #----------------------------------
    ### Process the specimen data 
    specimen_binned = (
         objS.biology[ 'specimen_df' ]
        # ---- Drop unaged fish 
        # !!! TODO: pending what FEAT says, weights associated with 
        # missing ages should be added into the 'unaged' category. 
        # This would further mean apportioning these into `weight_strata_aged_uanged`
        .dropna( how = 'any' )
        # ---- Remove unsexed fish
        .loc[ lambda df: df.sex != 'unsexed' ]
        # ---- Bin length
        .bin_variable( objS.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ] , 'length' )
        # ---- Age bin
        .bin_variable( objS.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ] , 'age' )
    )

    ### Sum weights within each length and age bin for each sex within each stratum
    specimen_binned_weight = (
        specimen_binned
        # ---- Group weight summations across distributions of length and age 
        # ---- for each sex within each stratum
        .groupby( [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] )
        # ---- Sum the weights 
        .apply( lambda df: pd.Series( { 'weight_all': df.weight.sum( ) ,
                                        'weight_adult': df.loc[ df.age > 1 ].weight.sum( ) } ) ) 
        # ---- Fill empty/non-existent values with 0's
        .fillna( 0 )
        .reset_index( )
    )

    ### Calculate the relative weight proportions of each length-age bin for each sex within each stratum
    proportions_weight_length_age_sex = (
        specimen_binned_weight
        # ---- Calculate total sex-specific weights for each stratum
        .assign( total_weight_sex_all = lambda df: df.groupby( [ 'stratum_num' , 'sex' ] )[ 'weight_all' ].transform( sum ) ,
                 total_weight_sex_adult = lambda df: df.groupby( [ 'stratum_num' , 'sex' ] )[ 'weight_adult' ].transform( sum ) )
        # ---- Calculate the weight proportions within each sex: from Matlab --> Len_age_key_wgt_*n
        .assign( proportion_weight_sex_all = lambda df: df.weight_all / df.total_weight_sex_all ,
                 proportion_weight_sex_adult = lambda df: df.weight_adult / df.total_weight_sex_adult )
    )

    ### Fill empty/non-existent values with 0's
    proportions_weight_length_age_sex[ 'proportion_weight_sex_all' ] = (
        proportions_weight_length_age_sex[ 'proportion_weight_sex_all' ].fillna( 0 )
    )

    proportions_weight_length_age_sex[ 'proportion_weight_sex_adult' ] = (
        proportions_weight_length_age_sex[ 'proportion_weight_sex_adult' ].fillna( 0 )
    )

    ### Check shape 
    test_dimensions = (
        len( np.unique( specimen_binned.age_bin ) ) *
        len( np.unique( specimen_binned.length_bin ) ) *
        len( np.unique( specimen_binned.sex ) ) *
        len( np.unique( specimen_binned.stratum_num ) )
    )

    assert obj_props_wgt_len_age_sex.shape == tuple( [ test_dimensions , 10 ] )

    ### Check values
    assert obj_props_wgt_len_age_sex.equals( proportions_weight_length_age_sex )