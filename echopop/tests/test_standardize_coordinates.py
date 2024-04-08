import pandas as pd
import numpy as np
from echopop.computation.spatial import transform_geometry

def test_transform_geometry( ):

    ### Mock data for `dataframe`
    test_dataframe = pd.DataFrame(
        {
            'transect_num': [ 1 , 2 , 3 , 4 , 5 ] ,
            'latitude': [ -5.0 , -2.5 , 0 , 2.5 , 5.0 ] ,
            'longitude': [ -180.0 , -120.0 , -60.0 , 0.0 , 60.0 ] ,
            'stratum_num': [ 1 , 1 , 1 , 2 , 2 ] ,
            'metric_value': [ 0.0 , 1.0 , 2.0 , 3.0 , 4.0 ] ,
        } ,
    )

    ### Mock data for `reference_grid`
    test_reference_grid = pd.DataFrame(
        {
            'latitude': [ -3.0 , -1.0 , 1.0 , 3.0 ] ,
            'longitude': [ -20.0 , -10.0 , 10.0 , 20.0 ] ,
        } ,
    )

    ### Mock data for `kriging_grid_parameters`
    test_kriging_grid_parameters = {
        'longitude_reference': 10.0 , 
        'longitude_offset': 10.0 ,
        'latitude_offset': 1.0 ,
    }

    ### Mock data for `projection`
    test_projection = 'epsg:4326'

    ### Evaluate when `range_output == True` and no `d_longitude`/`d_latitude` input
    # ---- Expect: tuple output
    output1 , output2 , output3 = transform_geometry( test_dataframe ,
                                                      test_reference_grid ,
                                                      test_kriging_grid_parameters ,
                                                      test_projection ,
                                                      range_output = True )
    # ---- Check output type
    assert isinstance( output1 , pd.DataFrame )
    assert isinstance( output2 , float )
    assert isinstance( output3 , float )
    # ---- Check output shape
    assert output1.shape == tuple( [ 5 , 10 ] )
    # ---- Check data value equality
    non_na_values = ~np.isnan( output1.longitude_transformed )
    assert np.all( output1.longitude_transformed[ non_na_values ] == np.array( [ -92.5 , -50.0 , -7.5 ] ) )
    assert np.isnan( output1.longitude_transformed[ 0 ] ) & np.isnan( output1.longitude_transformed[ 4 ] )
    assert np.allclose( output1.x_transformed.values[ non_na_values ] , 
                        np.array( [ -1.205 , -0.706 , -0.206 ] ) ,
                        rtol = 1e-1 )
    assert np.allclose( output1.y_transformed.values[ non_na_values ] , 
                        np.array( [ -1.1 , -0.6 , -0.1 ] ) ,
                        rtol = 1e-1 )
    assert output2 == 85.0
    assert output3 == 10.0

    ### Evaluate when `range_output == False` and no `d_longitude`/`d_latitude` input
    # ---- Expect: single dataframe output
    output1 = transform_geometry( test_dataframe ,
                                  test_reference_grid ,
                                  test_kriging_grid_parameters ,
                                  test_projection ,
                                  range_output = False )
    # ---- Check output type
    assert isinstance( output1 , pd.DataFrame )
    # ---- Check output shape
    assert output1.shape == tuple( [ 5 , 10 ] )
    # ---- Check data value equality
    non_na_values = ~np.isnan( output1.longitude_transformed )
    assert np.all( output1.longitude_transformed[ non_na_values ] == np.array( [ -92.5 , -50.0 , -7.5 ] ) )
    assert np.isnan( output1.longitude_transformed[ 0 ] ) & np.isnan( output1.longitude_transformed[ 4 ] )
    assert np.allclose( output1.x_transformed.values[ non_na_values ] , 
                        np.array( [ -1.205 , -0.706 , -0.206 ] ) ,
                        rtol = 1e-1 )
    assert np.allclose( output1.y_transformed.values[ non_na_values ] , 
                        np.array( [ -1.1 , -0.6 , -0.1 ] ) ,
                        rtol = 1e-1 )    
    
    ### Evaluate when `range_output == False` and `d_longitude`/`d_latitude` input
    # ---- Expect: tuple
    output1 = transform_geometry( test_dataframe ,
                                  test_reference_grid ,
                                  test_kriging_grid_parameters ,
                                  test_projection ,
                                  d_longitude = 60.0 ,
                                  d_latitude = 5.0 ,
                                  range_output = False )
    # ---- Check output type
    assert isinstance( output1 , pd.DataFrame )
    # ---- Check output shape
    assert output1.shape == tuple( [ 5 , 10 ] )
    # ---- Check data value equality
    non_na_values = ~np.isnan( output1.longitude_transformed )
    assert np.all( output1.longitude_transformed[ non_na_values ] == np.array( [ -92.5 , -50.0 , -7.5 ] ) )
    assert np.isnan( output1.longitude_transformed[ 0 ] ) & np.isnan( output1.longitude_transformed[ 4 ] )
    assert np.allclose( output1.x_transformed.values[ non_na_values ] , 
                        np.array( [ -1.707 , -1.000 , -0.291 ] ) ,
                        rtol = 1e-1 )
    assert np.allclose( output1.y_transformed.values[ non_na_values ] , 
                        np.array( [ -1.558 , -0.850 , -0.142 ] ) ,
                        rtol = 1e-1 )
    
    ### Evaluate when `range_output == True` and `d_longitude`/`d_latitude` input
    # ---- Expect: single dataframe output
    output11 , output22 , output33 = transform_geometry( test_dataframe ,
                                                         test_reference_grid ,
                                                         test_kriging_grid_parameters ,
                                                         test_projection ,
                                                         d_longitude = 60.0 ,
                                                         d_latitude = 5.0 ,
                                                         range_output = True )
    # ---- Check output type
    assert isinstance( output11 , pd.DataFrame )
    assert isinstance( output22 , float )
    assert isinstance( output33 , float )
    # ---- Check output shape
    assert output11.shape == tuple( [ 5 , 10 ] )
    # ---- Check data value equality
    non_na_values = ~np.isnan( output11.longitude_transformed )
    assert np.all( output11.longitude_transformed[ non_na_values ] == np.array( [ -92.5 , -50.0 , -7.5 ] ) )
    assert np.isnan( output11.longitude_transformed[ 0 ] ) & np.isnan( output11.longitude_transformed[ 4 ] )
    assert np.allclose( output11.x_transformed.values[ non_na_values ] , 
                        np.array( [ -1.707 , -1.000 , -0.291 ] ) ,
                        rtol = 1e-1 )
    assert np.allclose( output11.y_transformed.values[ non_na_values ] , 
                        np.array( [ -1.558 , -0.850 , -0.142 ] ) ,
                        rtol = 1e-1 )
    assert output22 == 60.0
    assert output33 == 5.0