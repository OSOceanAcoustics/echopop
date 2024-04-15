import pytest
import numpy as np
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

# ---- Dictionary shape comparison utility function
@pytest.fixture( scope = "session") 
def dictionary_shape( dictionary: dict ) :
    """
    A utility test function that extracts the shape of a nested dictionary
    """

    if isinstance( dictionary , dict ) :
        return( { i: dictionary_shape( dictionary[ i ] ) for i in dictionary } )
    else: 
        return None

# ---- Test for comparing Dictionary shapes/dimensions
@pytest.fixture( scope = "session") 
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

### Hook functions
def pytest_assertrepr_compare( config , op , left , right ):
    """
    Hook function that always shows the full `diff` on assertion
    failures by increasing the verbosity (`config.option.verbose`)
    """

    ### Adjust configuration `diff` verbosity
    config.option.verbose = 2

    return assertrepr_compare( config , op , left , right) 