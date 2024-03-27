import pandas as pd
import numpy as np
from echopop.computation.biology import sum_strata_weight

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