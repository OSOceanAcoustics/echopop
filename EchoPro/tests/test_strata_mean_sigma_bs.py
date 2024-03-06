import pandas as pd
import numpy as np
from EchoPro.computation.acoustics import to_linear , ts_length_regression

def test_strata_mean_sigma_bs( mock_survey ):
    #### Pull in mock Survey object
    objS = mock_survey

    ### Re-parameterize `specimen_df` with dummy data 
    objS.biology[ 'specimen_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 , 2 , 4 , 5 ] , 4 ) ,
            'haul_num': np.repeat( [ 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 ] , 2 ) ,
            'species_id': np.append( np.repeat( [ 19350 ] , 19 ) , 43130 ) ,
            'length': np.linspace( 10 , 100 , 20 ) ,
            'weight': np.linspace( 1 , 5 , 20 ) ,     
        }
    )

    ### Re-parameterize `length_df` with dummy data
    objS.biology[ 'length_df' ] = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 0 , 1 , 2 , 4 , 5 ] , 4 ) ,
            'haul_num': np.repeat( [ 1 , 2 , 3 , 4 , 5 ] , 4 ) ,
            'species_id': np.append( np.repeat( [ 19350 ] , 19 ) , 43130 ) ,
            'length': np.linspace( 10 , 100 , 20 ) ,
            'length_count': np.linspace( 10 , 100 , 20 ) ,     
        }
    )

    ### Re-parameterize `strata_df` with dummy data 
    objS.spatial[ 'strata_df' ] = pd.DataFrame(
        {
            'stratum_num': [ 0 , 1 , 2 , 3 , 4 , 5 , 6 ]
        }
    )

    ### Dummy parameters
    objS.config[ 'TS_length_regression_parameters' ][ 'pacific_hake' ] = { 'species_code': 19350 ,
                                                                           'TS_L_slope': 20.0 ,
                                                                           'TS_L_intercept': -84.1 ,
                                                                           'length_units': 'cm' }
    
    ### Define dummy `species_id` code
    species_id = 19350

    #----------------------------------
    ### Run tests: `strata_mean_sigma_bs`
    #----------------------------------
    ### Evaluate whether non-specified `species_id` code are removed 
    # `specimen_df`: 20 rows -> 19 rows
    specimen_df_copy = objS.biology[ 'specimen_df' ].copy( )
    specimen_df_copy = specimen_df_copy[ specimen_df_copy.species_id == species_id ]
    assert specimen_df_copy.shape[ 0 ] == 19

    # `length_df`: 20 rows -> 19 rows
    length_df_copy = objS.biology[ 'length_df' ].copy( )
    length_df_copy = length_df_copy[ length_df_copy.species_id == species_id ]
    assert length_df_copy.shape[ 0 ] == 19

    #----------------------------------
    ### Next step: Concatenate the two dataframes
    # Re-bin `specimen_df_copy`
    spec_df_reframed = (
        specimen_df_copy
        .groupby( [ 'haul_num' , 'stratum_num' , 'species_id', 'length' ] )
        .apply( lambda x: len( x[ 'length' ] ) )
        .reset_index( name = 'length_count' )
    )

    # Concatenate
    all_length_df = pd.concat( [ spec_df_reframed , length_df_copy ] , 
                               join = 'inner' )
    assert all_length_df.shape[ 0 ] == 38

    #----------------------------------
    ### TS-length parameterization & modeling
    # Fit parameters
    ts_length_parameters = objS.config[ 'TS_length_regression_parameters' ][ 'pacific_hake' ]
    slope = ts_length_parameters[ 'TS_L_slope' ]
    intercept = ts_length_parameters[ 'TS_L_intercept' ]
    assert slope == 20.0
    assert intercept == -84.1

    # Predict `TS`
    all_length_df[ 'TS' ] = ts_length_regression( all_length_df[ 'length' ] , slope , intercept )
    assert np.isclose( np.median( all_length_df[ 'TS' ] ) , -49.675 )

    # Linearize to `sigma_bs`
    all_length_df[ 'sigma_bs' ] = to_linear( all_length_df[ 'TS' ] )
    np.isclose( all_length_df[ 'sigma_bs' ].mean( ) , 1.3395e-5 )

    #----------------------------------
    ### Calculate mean `sigma_bs` per `haul_num` and then `stratum_num`
    # `haul_num`
    mean_haul_sigma_bs = (
        all_length_df
        .groupby( [ 'haul_num' , 'stratum_num' , 'species_id' ] )[ [ 'sigma_bs' , 'length_count' ] ]
        .apply(lambda x: np.average( x[ 'sigma_bs' ] , weights= x[ 'length_count' ] ) )
        .to_frame( 'sigma_bs_mean' )
        .reset_index( )
    )
    assert mean_haul_sigma_bs.shape[ 0 ] == 14

    # `stratum_num`
    mean_strata_sigma_bs = (
        mean_haul_sigma_bs
        .groupby( [ 'stratum_num' , 'species_id' ] )[ 'sigma_bs_mean' ]
        .mean( )
        .reset_index( )
    )
    assert mean_strata_sigma_bs.shape[ 0 ] == 5
    assert np.allclose( mean_strata_sigma_bs.sigma_bs_mean , 
                        np.array( [ 1.659e-6 , 5.238e-6 , 1.195e-5 , 2.145e-5 , 3.254e-5 ] ) )
    assert any( ~ mean_strata_sigma_bs.stratum_num.isin( [ 3 , 6 ] ) )

    #----------------------------------
    ### Add to object as dictionary
    objS.acoustics[ 'sigma_bs' ] = { 
        'length_binned': all_length_df ,
        'haul_mean': mean_haul_sigma_bs ,
        'strata_mean': mean_strata_sigma_bs ,
    }
    assert 'sigma_bs' in objS.acoustics.keys( )

    #----------------------------------
    ### `impute_missing_sigma_bs`
    # Run function -- missing strata (3, 6)
    objS.impute_missing_sigma_bs( species_id )

    # Pull results
    sigma_bs_imputed_means = objS.acoustics[ 'sigma_bs' ][ 'strata_mean' ].reset_index( drop = True )

    assert any( sigma_bs_imputed_means.stratum_num.isin( [ 3 , 6 ] ) )
    assert np.allclose( sigma_bs_imputed_means.loc[ lambda x: x.stratum_num.isin( [ 3 , 6 ] ) ][ 'sigma_bs_mean' ] ,
                        np.array( [ 1.670e-5 , 3.254e-5 ] ) )