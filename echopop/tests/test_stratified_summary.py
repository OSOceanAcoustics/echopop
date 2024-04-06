import numpy as np
import pandas as pd
from echopop.tests.utility_testing_functions import dictionary_shape_equal
from echopop.computation.statistics import stratified_transect_statistic

def test_stratified_transect_statistic( ):

    ### Create mock data for `transect_summary`
    test_transect_summary = pd.DataFrame(
        {
            'transect_num': [ 1 , 2 , 3 , 4 ] ,
            'minimum_longitude': [ -5.0 , -3.0 , -1.0 , 1.0 ] ,
            'maxmum_longitude': [ -2.0 , 5.0 , 3.0 , 7.0 ] ,
            'center_latitude': [ 10.0 , 11.0 , 12.5 , 13.5 ] ,
            'transect_distance': [ 177.600950 , 472.070493 , 234.766275 , 350.736855 ] ,
            'transect_spacing': [ 2.0 , 2.0 , 2.0 , 2.0 ] ,
            'transect_area': [ 355.201900 , 944.140986 , 469.532550 , 701.473710 ] ,
            'B_adult': [ 1e2 , 1e3 , 1e5 , 1e4 ] ,
            'stratum_inpfc': [ 1 , 1 , 2 , 2 ]
        }
    )

    ### Create mock data for `strata_summary`
    test_strata_summary = pd.DataFrame(
        {
            'stratum_inpfc': [ 1 , 2 ] ,
            'num_transects': [ 2 , 2 ] ,
            'total_transect_area': [ 1299.342886 , 1171.006260 ] ,
        }
    )

    ### Evaluate for later comparison
    # ---- Replicates == 1
    # ---- Transect sample proportion == 100%
    test_transect_sample = 1.0
    test_transect_replicates = 1
    eval_single_stratified_results = stratified_transect_statistic( test_transect_summary ,
                                                                    test_strata_summary ,
                                                                    test_transect_sample ,
                                                                    test_transect_replicates ,
                                                                    parameter = 'B_adult' )
    # ---- Replicates == 10
    # ---- Transect sample proportion == 100%
    test_transect_sample = 1.0
    test_transect_replicates = 10
    eval_single_rep_stratified_results = stratified_transect_statistic( test_transect_summary ,
                                                                        test_strata_summary ,
                                                                        test_transect_sample ,
                                                                        test_transect_replicates ,
                                                                        parameter = 'B_adult' )
    
    # ---- Replicates == 1 
    # ---- Transect sample proportion == 50%
    test_transect_sample = 0.5
    test_transect_replicates = 1
    np.random.seed( 10 )
    eval_single_sub_stratified_results = stratified_transect_statistic( test_transect_summary ,
                                                                        test_strata_summary ,
                                                                        test_transect_sample ,
                                                                        test_transect_replicates ,
                                                                        parameter = 'B_adult' )
    
    # ---- Replicates == 1 
    # ---- Transect sample proportion == 50%
    test_transect_sample = 0.5
    test_transect_replicates = 10
    np.random.seed( 1800 )
    eval_single_sub_rep_stratified_results = stratified_transect_statistic( test_transect_summary ,
                                                                            test_strata_summary ,
                                                                            test_transect_sample ,
                                                                            test_transect_replicates ,
                                                                            parameter = 'B_adult' )
    
    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    expected_output = {
        'biomass': {
            'mean': {
                'estimate': 1 ,
                'confidence_interval': np.array( [ 1 , 1 ] ) ,
            } ,
            'variance': {
                'estimate': 1 ,
                'confidence_interval': np.array( [ 1 , 1 ] ) ,
            } ,
            'CV': {
                'estimate': 1 ,
                'confidence_interval': np.array( [ 1 , 1 ] ) ,
            } ,            
        }
    }

    #----------------------------------
    ### Run tests: `stratified_transect_statistic`
    #----------------------------------
    ### Dictionary structure
    # !!! TODO: based on the original data structure -- will need to be updated once the core data structure is also updated 
    # ---- Check attributes
    assert set( eval_single_stratified_results[ 'biomass' ].keys( ) ) == ( set( [ 'mean' , 'variance' , 'CV' ] ) )
    assert set( eval_single_rep_stratified_results[ 'biomass' ].keys( ) ) == ( set( [ 'mean' , 'variance' , 'CV' ] ) )
    assert set( eval_single_sub_stratified_results[ 'biomass' ].keys( ) ) == ( set( [ 'mean' , 'variance' , 'CV' ] ) )
    assert set( eval_single_sub_rep_stratified_results[ 'biomass' ].keys( ) ) == ( set( [ 'mean' , 'variance' , 'CV' ] ) )
    # ---- Check sub-directory keys and structure
    assert dictionary_shape_equal( eval_single_stratified_results , expected_output )
    assert dictionary_shape_equal( eval_single_rep_stratified_results , expected_output )
    assert dictionary_shape_equal( eval_single_sub_stratified_results , expected_output )
    assert dictionary_shape_equal( eval_single_sub_rep_stratified_results , expected_output )
    ### Data outputs
    # ++++ mean
    # ---- > estimate
    assert np.isclose( eval_single_stratified_results[ 'biomass' ][ 'mean' ][ 'estimate' ] , 
                       54947653.0 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_rep_stratified_results[ 'biomass' ][ 'mean' ][ 'estimate' ] , 
                       54947653.0 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_sub_stratified_results[ 'biomass' ][ 'mean' ][ 'estimate' ] , 
                       117230560.0 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_sub_rep_stratified_results[ 'biomass' ][ 'mean' ][ 'estimate' ] , 
                       54463985.0 ,
                       rtol = 1e-2 )
    # ---- > confidence interval
    assert np.allclose( eval_single_stratified_results[ 'biomass' ][ 'mean' ][ 'confidence_interval' ] ,
                        np.array( [ 54947653.28 , 54947653.28 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_rep_stratified_results[ 'biomass' ][ 'mean' ][ 'confidence_interval' ] ,
                        np.array( [ 54947653.28 , 54947653.28 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_sub_stratified_results[ 'biomass' ][ 'mean' ][ 'confidence_interval' ] ,
                        np.array( [ 1.17e8 , 1.172e8 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_sub_rep_stratified_results[ 'biomass' ][ 'mean' ][ 'confidence_interval' ] ,
                        np.array( [ -4.69e7 , 1.57e8 ] ) ,
                        rtol = 1e-2 )
    # ++++ variance
    assert np.isclose( eval_single_stratified_results[ 'biomass' ][ 'variance' ][ 'estimate' ] , 
                       54846534.0 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_rep_stratified_results[ 'biomass' ][ 'variance' ][ 'estimate' ] , 
                       54846534.0 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_sub_stratified_results[ 'biomass' ][ 'variance' ][ 'estimate' ] , 
                       116601900.0 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_sub_rep_stratified_results[ 'biomass' ][ 'variance' ][ 'estimate' ] , 
                       53662832.0 ,
                       rtol = 1e-2 )
    # ---- > confidence interval
    assert np.allclose( eval_single_stratified_results[ 'biomass' ][ 'variance' ][ 'confidence_interval' ] ,
                        np.array( [ 54846534.0 , 54846534.0 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_rep_stratified_results[ 'biomass' ][ 'variance' ][ 'confidence_interval' ] ,
                        np.array( [ 54846534.0 , 54846534.0 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_sub_stratified_results[ 'biomass' ][ 'variance' ][ 'confidence_interval' ] ,
                        np.array( [ 1.17e8 , 1.17e8] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_sub_rep_stratified_results[ 'biomass' ][ 'variance' ][ 'confidence_interval' ] ,
                        np.array( [ -4.71e7 , 1.53e8 ] ) ,
                        rtol = 1e-2 )
    # ++++ CV
    assert np.isclose( eval_single_stratified_results[ 'biomass' ][ 'CV' ][ 'estimate' ] , 
                       0.998 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_rep_stratified_results[ 'biomass' ][ 'CV' ][ 'estimate' ] , 
                       0.998 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_sub_stratified_results[ 'biomass' ][ 'CV' ][ 'estimate' ] , 
                       0.995 ,
                       rtol = 1e-2 )
    assert np.isclose( eval_single_sub_rep_stratified_results[ 'biomass' ][ 'CV' ][ 'estimate' ] , 
                       0.971 ,
                       rtol = 1e-2 )
    # ---- > confidence interval
    assert np.allclose( eval_single_stratified_results[ 'biomass' ][ 'CV' ][ 'confidence_interval' ] ,
                        np.array( [ 0.998 , 0.998 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_rep_stratified_results[ 'biomass' ][ 'CV' ][ 'confidence_interval' ] ,
                        np.array( [ 0.998 , 0.998 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_sub_stratified_results[ 'biomass' ][ 'CV' ][ 'confidence_interval' ] ,
                        np.array( [ 0.995 , 0.995 ] ) ,
                        rtol = 1e-2 )
    assert np.allclose( eval_single_sub_rep_stratified_results[ 'biomass' ][ 'CV' ][ 'confidence_interval' ] ,
                        np.array( [ 0.904 , 1.038 ] ) ,
                        rtol = 1e-2 )
    

    
