import pytest
from typing import Union
import numpy as np
import pandas as pd
from pathlib import Path
from echopop import Survey
from _pytest.assertion.util import assertrepr_compare

### Set up path to the `test_data` folder
HERE = Path(__file__).parent.absolute()
TEST_DATA_ROOT = HERE.parent / "test_data"

### Fixtures
# ---- Test root/config/input file paths
@pytest.fixture( scope = "session" )
def test_path( ) :

    return {
        "ROOT" : TEST_DATA_ROOT ,
        "CONFIG" : TEST_DATA_ROOT / "config_files" ,
        "INPUT" : TEST_DATA_ROOT / "input_files" , 
    }

# ---- Mock `Survey` class object 
@pytest.fixture( scope = "session") 
def mock_survey( test_path ) -> Survey :

    return Survey(
        init_config_path =Path( test_path[ "CONFIG" ] / "config_init.yml" ) ,
        survey_year_config_path =Path( test_path[ "CONFIG" ] / "config_survey.yml" ) ,
    )

### Hook functions
def pytest_assertrepr_compare( config , op , left , right ):
    """
    Hook function that always shows the full `diff` on assertion
    failures by increasing the verbosity (`config.option.verbose`)
    """

    ### Adjust configuration `diff` verbosity
    config.option.verbose = 2

    return assertrepr_compare( config , op , left , right) 

### Utility functions
# ---- Dictionary shape comparison utility function
def dictionary_shape( dictionary: dict ) :
    """
    A utility test function that extracts the shape of a nested dictionary
    """

    if isinstance( dictionary , dict ) :
        return( { i: dictionary_shape( dictionary[ i ] ) for i in dictionary } )
    else: 
        return None

# ---- Test for comparing Dictionary shapes/dimensions
def dictionary_shape_equal( dictionary1: dict ,
                            dictionary2: dict ):
    """
    Tests equality between the shapes of two nested dictionaries
    """
    result = dictionary_shape( dictionary1 ) == dictionary_shape( dictionary2 )
    
    if result :
        return result 
    else:
        if set( dictionary_shape( dictionary1 ) ) <= set( dictionary_shape( dictionary2 ) ) :
            tracked_true = [ ]
            
            for j in dictionary2.keys( ) :
                test = set( dictionary1[ j ].keys( ) ) <= ( dictionary2[ j ].keys( ) )
                tracked_true.append( test )
                
            if np.all( tracked_true ) :
                return True 
            else : 
                return result 
        else :  
            return result

# ---- Test for dataframe shape equality
def dataframe_shape_equal( input: Union[ pd.DataFrame , dict ] ,
                           reference: Union[ tuple , dict ] ):

    ### DataFrame
    if ( isinstance( input , pd.DataFrame ) ) & ( isinstance( reference , tuple ) ) :
        assert input.shape == reference
    
    ### Dictionary
    elif ( isinstance( input , dict ) ) & ( isinstance( reference , dict ) ):
        assert extract_dataframe_shape( input ) == extract_dataframe_shape( reference )


# ---- Test for comparing Dictionary equality (including nested DataFrames)
def dictionary_equality( dictionary1: dict , 
                         dictionary2: dict ):
    
    ### Iterate through nested DataFrames within each dictionary
    for key , expected_df in dictionary2.items( ) :

        if isinstance( dictionary1[ key ] , pd.DataFrame ) :
            dataframe_equality( dictionary1[ key ] , expected_df )
        
        else :
            for sub_key , _ in dictionary2[ key ].items( ):
                dataframe_equality( dictionary1[ key ][ sub_key ] , 
                                    expected_df[ sub_key ] )

# ---- Extract dataframe shape
def extract_dataframe_shape( input: Union[ pd.DataFrame , dict ] ):

    ### DataFrame
    if isinstance( input , pd.DataFrame ) :

        return input.shape
    
    ### Dictionary
    elif isinstance( input , dict ) :
        dataframe_shapes = { }

        for key , value in input.items( ):
            if isinstance( value , pd.DataFrame ) :
                dataframe_shapes[ key ] = value.shape
            elif isinstance( value , dict ) :
                dataframe_shapes[ key ] = extract_dataframe_shape( value )

        return dataframe_shapes

# ---- Extract dataframe dtypes
def dataframe_dtypes_equal( dataframe: pd.DataFrame ,
                            reference_dictionary: dict ):
    
    ### Separate evaluation for categorical-type
    # ---- Parse expected categorical variables
    categorical_columns = [ k for k , v in reference_dictionary.items( ) if isinstance( v , pd.CategoricalDtype ) ]

    # ---- Assert that all categorical columns in the reference dictionary match the categorical
    # ----- columns in the tested dataframe
    assert np.all( dataframe.select_dtypes( include = [ 'category' ] ).columns.isin( categorical_columns ) )

    # ---- Remove categorical columns from the dataframe
    dataframe = dataframe.copy( ).drop( categorical_columns , axis = 1 )

    ### Loop through columns to assert that dtypes from the tested dataframe
    ### match those expected in a reference dictionary
    for column , dtype in dataframe.dtypes.items( ):
        assert np.issubdtype( dtype , reference_dictionary.get( column , object ) ) , \
            f"Data type mismatch for column '{ column }'" 


# ---- Test for evaluating differences in DataFrame `dtypes`
def dataframe_dtypes_equality( input: Union[ pd.DataFrame , dict ] ,
                               reference: dict ):
    
    ### DataFrame
    if isinstance( input , pd.DataFrame ) :
        dataframe_dtypes_equal( input , reference )

    ### Dictionary
    elif isinstance( input , dict ) :
        for category , data in reference.items( ) :

            # ---- Single Dictionary layer
            if isinstance( input[ category ] , pd.DataFrame ):
                dataframe_dtypes_equal( input[ category ] , 
                                        reference[ category ] )
            
            # ---- Nested Dictionary layers
            else:
                for df_name , _ in data.items( ):
                    dataframe_dtypes_equal( input[ category ][ df_name ] , reference[ category ][ df_name ] )

# ---- Test for evaluating equality between two dataframes
def dataframe_equality( dataframe1: pd.DataFrame ,
                        dataframe2: pd.DataFrame ):

    ### Evaluate equality between numerical values
    assert np.allclose( dataframe1.select_dtypes( include = [ 'number' ] ) ,
                        dataframe2.select_dtypes( include = [ 'number' ] ) )
    
    ### Evaluate equality between non-numerical values
    assert np.all( dataframe1.select_dtypes( exclude = [ 'number' ] ) == dataframe2.select_dtypes( exclude = [ 'number' ] ) )
