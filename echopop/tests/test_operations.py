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

def test_bin_stats( ):

    ### Mock dataframe
    test_dataframe = pd.DataFrame( 
        {
            'animal': [ 'pretty pink pony' , 'big blue bass' , 'silly silver silkworm' ,
                        'gnarly green grouse' , 'roudy red rabbit' , 'magenta mad manatee' ] ,
            'length': [ 2.0 , 4.0 , 8.0 , 3.0 , 6.0 , 7.0 ] ,
            'weight': [ 100.0 , 200.0 , 300.0 , 300.0 , 200.0 , 100.0 ] ,
            'location': [ 'timbuktu' , 'timbuktu' , 'timbuktu' ,
                          'lost city of z' , 'lost city of z' , 'lost city of z' ] ,
        } ,
    )

    ### Mock bin_values
    test_bin_values = np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] )

    ### Evaluate for later comparison
    # ++++ No contrast | length + weight
    # ---- Monkey patch method (TEMPORARY)
    eval_dataframe_monkey_lwnc = test_dataframe.bin_stats( 'length' , test_bin_values )
    # ---- Normal function
    eval_dataframe_function_lwnc = bin_stats( test_dataframe , 'length' , test_bin_values )
    # ++++ No contrast | length
    eval_dataframe_monkey_lnc = test_dataframe.bin_stats( 'length' , test_bin_values , variables = 'length' )
    # ---- Normal function
    eval_dataframe_function_lnc = bin_stats( test_dataframe , 'length' , test_bin_values , variables = 'length' )
    # ++++ No contrast | length ~ function: just mean
    eval_dataframe_monkey_lncm = test_dataframe.bin_stats( 'length' , test_bin_values , variables = 'length' , functions = [ 'mean' ] )
    # ---- Normal function
    eval_dataframe_function_lncm = bin_stats( test_dataframe , 'length' , test_bin_values , variables = 'length' , functions = [ 'mean' ] )
    # ++++ No contrast | length ~ function: just mean
    eval_dataframe_monkey_lwc = test_dataframe.bin_stats( 'length' , test_bin_values , contrasts = [ 'location' ] , variables = 'length' )
    # ---- Normal function
    eval_dataframe_function_lwc = bin_stats( test_dataframe , 'length' , test_bin_values , contrasts = [ 'location' ] , variables = 'length' )

    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions_lwnc = tuple( [ 4 , 5 ] )
    expected_dimensions_lnc = tuple( [ 4 , 3 ] )
    expected_dimensions_lncm = tuple( [ 4 , 2 ] )
    expected_dimensions_lwc = tuple( [ 8 , 4 ] )
    # ---- Expected output
    expected_output_lwnc = pd.DataFrame(
        {   
            'length_bin': pd.cut( [ 2.0 , 4.0 , 6.0 , 8.0 ] ,
                                  np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
            'mean_length': [ 2.5 , 4.0 , 6.5 , 8.0 ] ,
            'n_length': [ 2 , 1 , 2 , 1 ] ,
            'mean_weight': [ 200.0 , 200.0 , 150.0 , 300.0 ] ,
            'n_weight': [ 2 , 1 , 2 , 1 ] ,
        } ,
    )
    expected_output_lnc = pd.DataFrame(
        {   
            'length_bin': pd.cut( [ 2.0 , 4.0 , 6.0 , 8.0 ] ,
                                  np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
            'mean_length': [ 2.5 , 4.0 , 6.5 , 8.0 ] ,
            'n_length': [ 2 , 1 , 2 , 1 ] ,
        } ,
    )
    expected_output_lncm = pd.DataFrame(
        {   
            'length_bin': pd.cut( [ 2.0 , 4.0 , 6.0 , 8.0 ] ,
                                  np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
            'mean_length': [ 2.5 , 4.0 , 6.5 , 8.0 ]
        } ,
    )
    expected_output_lwc = pd.DataFrame(
        {   
            'length_bin': pd.cut( np.repeat( [ 2.0 , 4.0 , 6.0 , 8.0 ] , 2 ) ,
                                  np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
            'location': np.tile( [ 'lost city of z' , 'timbuktu' ] , 4 ) ,
            'mean_length': [ 3.0 , 2.0 , 0.0 , 4.0 , 6.5 , 0.0 , 0.0 , 8.0 ] ,
            'n_length': [ 1 , 1 , 0 , 1 , 2 , 0 , 0 , 1 ] ,
        } ,
    )

    #----------------------------------
    ### Run tests: `bin_stats`
    #----------------------------------
    ### Check shape
    assert eval_dataframe_monkey_lwnc.shape == expected_dimensions_lwnc
    assert eval_dataframe_function_lwnc.shape == expected_dimensions_lwnc
    assert eval_dataframe_monkey_lnc.shape == expected_dimensions_lnc
    assert eval_dataframe_function_lnc.shape == expected_dimensions_lnc
    assert eval_dataframe_monkey_lncm.shape == expected_dimensions_lncm
    assert eval_dataframe_function_lncm.shape == expected_dimensions_lncm
    assert eval_dataframe_monkey_lwc.shape == expected_dimensions_lwc
    assert eval_dataframe_function_lwc.shape == expected_dimensions_lwc
    ### Check output
    assert eval_dataframe_monkey_lwnc.equals( expected_output_lwnc )
    assert eval_dataframe_function_lwnc.equals( expected_output_lwnc )
    assert eval_dataframe_monkey_lnc.equals( expected_output_lnc )
    assert eval_dataframe_function_lnc.equals( expected_output_lnc )
    assert eval_dataframe_monkey_lncm.equals( expected_output_lncm )
    assert eval_dataframe_function_lncm.equals( expected_output_lncm )
    assert eval_dataframe_monkey_lwc.equals( expected_output_lwc )
    assert eval_dataframe_function_lwc.equals( expected_output_lwc )

def test_count_variable( ):

    ### Mock dataframe
    test_dataframe = pd.DataFrame( 
        {
            'animal': [ 'pretty pink pony' , 'big blue bass' , 'silly silver silkworm' ,
                        'gnarly green grouse' , 'roudy red rabbit' , 'magenta mad manatee' ,
                        'pretty pink pony' , 'big blue bass' , 'silly silver silkworm' ,
                        'gnarly green grouse' , 'roudy red rabbit' , 'magenta mad manatee' ] ,
            'length': [ 2.0 , 4.0 , 8.0 , 3.0 , 6.0 , 7.0 ,
                        2.0 , 4.0 , 8.0 , 3.0 , 6.0 , 7.0 ] ,            
            'location': [ 'timbuktu' , 'timbuktu' , 'timbuktu' ,
                         'timbuktu' , 'timbuktu' , 'timbuktu' ,
                         'lost city of z' , 'lost city of z' , 'lost city of z' ,
                         'lost city of z' , 'lost city of z' , 'lost city of z' ] ,
            'length_count': [ 10 , 20 , 30 , 40 , 50 , 60 ,
                              60 , 50 , 40 , 30 , 20 , 10 ] ,
        } ,
    )

    ### Evaluate for later comparison
    # ---- Monkey patch method (TEMPORARY)
    eval_dataframe_monkey = test_dataframe.count_variable( [ 'location' , 'animal' ] , 'length_count' , 'sum' )
    # ---- Normal function
    eval_dataframe_function = count_variable( test_dataframe , [ 'location' , 'animal' ] , 'length_count' , 'sum' )

    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions = tuple( [ 12 , 3 ] )
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            'location': np.repeat( [ 'lost city of z' , 'timbuktu'  ] , 6 ) ,
            'animal': np.tile( [ 'big blue bass' , 'gnarly green grouse' , 'magenta mad manatee' ,
                                 'pretty pink pony' , 'roudy red rabbit' , 'silly silver silkworm' ] , 2 ) ,
            'count': [ 50 , 30 , 10 , 60 , 20 , 40 ,
                       20 , 40 , 60 , 10 , 50 , 30 ] ,
        } ,
    )

    #----------------------------------
    ### Run tests: `count_variable`
    #----------------------------------
    ### Check shape
    assert eval_dataframe_monkey.shape == expected_dimensions
    assert eval_dataframe_function.shape == expected_dimensions
    ### Check dataframe equality
    assert eval_dataframe_monkey.equals( expected_output )
    assert eval_dataframe_function.equals( expected_output )