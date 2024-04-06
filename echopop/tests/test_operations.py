import numpy as np
import pandas as pd
from echopop.computation.operations import bin_variable , bin_stats , count_variable , meld , stretch , group_merge

def test_bin_variable( ):

    ### Mock dataframe
    test_dataframe = pd.DataFrame( 
        {
            'animal': [ 'pretty pink pony' , 'big blue bass' , 'silly silver silkworm' ] ,
            'length': [ 2.0 , 4.0 , 8.0 ] ,
        } ,
    )

    ### Mock bin_values
    test_bin_values = np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] )

    ### Evaluate for later comparison
    # ---- Monkey patch method (TEMPORARY)
    eval_dataframe_monkey = test_dataframe.bin_variable( test_bin_values , 'length' )
    # ---- Normal function
    eval_dataframe_function = bin_variable( test_dataframe , test_bin_values , 'length' )

    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions = tuple( [ 3 , 3 ] )
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            'animal': [ 'pretty pink pony' , 'big blue bass' , 'silly silver silkworm' ] ,
            'length': [ 2.0 , 4.0 , 8.0 ] ,
            'length_bin': pd.cut( [ 2.0 , 4.0 , 8.0 ] ,
                                  np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
        } ,
    )

    #----------------------------------
    ### Run tests: `bin_variable`
    #----------------------------------
    ### Check shape
    assert eval_dataframe_monkey.shape == expected_dimensions
    assert eval_dataframe_function.shape == expected_dimensions
    ### Check output
    assert eval_dataframe_monkey.equals( expected_output )
    assert eval_dataframe_function.equals( expected_output )