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

def test_meld( ):

    ### Mock specimen dataframe
    test_specimen_dataframe = pd.DataFrame(
        {
            'stratum_num': np.repeat( 1 , 12 ) ,
            'species_id': np.tile( [ 'big blue bass' , 'pretty pink pony' , 'silly silver silkworm' ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 6 ) ,
            'group': np.repeat( 'sexed' , 12 ) ,
            'station': np.repeat( 'clouds' , 12 ) ,
            'length': [ 5.0 , 4.0 , 6.0 , 5.0 , 5.0 , 4.0 ,
                        5.0 , 5.0 , 7.0 , 4.0 , 5.0 , 6.0 ] ,
            'length_bin': pd.cut( [ 5.0 , 4.0 , 6.0 , 5.0 , 5.0 , 4.0 ,
                        5.0 , 5.0 , 7.0 , 4.0 , 5.0 , 6.0 ] ,
                                  np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
        } ,
    )

    ### Mock length dataframe
    test_length_dataframe = pd.DataFrame(
        {
            'stratum_num': np.repeat( 1 , 6 ) ,
            'species_id': np.tile( [ 'big blue bass' , 'pretty pink pony' , 'silly silver silkworm' ] , 2 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 3 ) ,
            'group': np.repeat( 'sexed' , 6 ) ,
            'station': np.repeat( 'waves' , 6 ) ,
            'length': [ 2.0 , 4.0 , 3.0 , 2.0 , 4.0 , 3.0 ] ,
            'length_bin': pd.cut(  [ 2.0 , 4.0 , 3.0 , 2.0 , 4.0 , 3.0 ] ,
                                  np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
            'length_count': [ 10 , 20 , 30 , 30 , 20 , 10 ] ,
        } ,
    ) 

    ### Evaluate for later comparison
    # ---- Monkey patch method (TEMPORARY)
    eval_dataframe_monkey = test_specimen_dataframe.meld( test_length_dataframe )
    # ---- Normal function
    eval_dataframe_function = meld( test_specimen_dataframe , test_length_dataframe )

    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions = tuple( [ 16 , 8 ] )
    # ---- Expected output 
    expected_output = pd.DataFrame(
        {
            'stratum_num': np.repeat( 1 , 16 ) ,
            'species_id': np.concatenate( [ np.repeat( 'big blue bass' , 3 ) ,
                                            np.repeat( 'pretty pink pony' , 3 ) ,
                                            np.repeat( 'silly silver silkworm' , 4 ) ,
                                            np.tile( [ 'big blue bass' , 'pretty pink pony' , 
                                                       'silly silver silkworm' ] , 2 ) ] ) ,
            'sex': [ 'female' , 'female' , 'male' , 'female' , 'female' , 'male' , 'female' , 'female' ,
                     'male' , 'male' , 'male' , 'female' , 'male' , 'female' , 'male' , 'female' ] ,
            'group': np.repeat( 'sexed' , 16 ) ,
            'station': np.concatenate( [ np.repeat( 'clouds' , 10 ) ,
                                         np.repeat( 'waves' , 6 ) ] ) ,
            'length': [ 4.0 , 5.0 , 5.0 , 4.0 , 5.0 , 5.0 , 4.0 , 6.0 , 
                        6.0 , 7.0 , 2.0 , 4.0 , 3.0 , 2.0 , 4.0 , 3.0 ] ,
            'length_bin': pd.cut( [ 4.0 , 5.0 , 5.0 , 4.0 , 5.0 , 5.0 , 4.0 , 6.0 , 
                                    6.0 , 7.0 , 2.0 , 4.0 , 3.0 , 2.0 , 4.0 , 3.0 ] ,
                                 np.array( [ 1.0 , 3.0 , 5.0 , 7.0 , 9.0 ] ) ) ,
            'length_count': [ 1 , 1 , 2 , 1 , 1 , 2 , 1 , 1 , 1 , 1 ,
                              10 , 20 , 30 , 30 , 20 , 10 ] ,
        } ,
    )

    #----------------------------------
    ### Run tests: `count_variable`
    #----------------------------------
    ### Check shape
    assert eval_dataframe_monkey.shape == expected_dimensions
    assert eval_dataframe_function.shape == expected_dimensions
    ### Check output
    assert np.all( eval_dataframe_monkey == expected_output )
    assert np.all( eval_dataframe_function == expected_output )

def test_stretch( ):

    ### Create mock dataframe
    test_dataframe = pd.DataFrame(
        {
            'stratum_num': [ 1 , 1 , 2 , 2 ] ,
            'transect_num': [ 1 , 2 , 3 , 4 ] ,
            'latitude': [ 0.0 , 1.0 , 3.0 , 4.0 ] ,
            'longitude': [ -1.0 , 0.0 , 1.0 , 2.0 ] ,
            'load_a_male': [ 5.0 , 4.0 , 2.0 , 1.0 ] ,
            'load_a_female': [ 10.0 , 3.0 , 5.0 , 6.0 ] ,
        } ,
    )

    ### Eval for later comparison 
    # ---- Monkey patch method (TEMPORARY)
    eval_dataframe_monkey = test_dataframe.stretch( variable = 'load_a' )
    # ---- Normal function
    eval_dataframe_function = stretch( test_dataframe , variable = 'load_a' )

    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions = tuple( [ 8 , 6 ] )
    # ---- Expected output
    expected_output = pd.DataFrame(
        {
            'transect_num': np.repeat( [ 1 , 2 , 3 , 4 ] , 2 ) ,
            'latitude': np.repeat( [ 0.0 , 1.0 , 3.0 , 4.0 ] , 2 ) ,
            'longitude': np.repeat( [ -1.0 , 0.0 , 1.0 , 2.0 ] , 2 ) ,
            'stratum_num': np.repeat( [ 1 , 2 ] , 4 ) ,
            'sex': np.tile( [ 'male' , 'female' ] , 4 ) ,
            'load_a': [ 5.0 , 10.0 , 4.0 , 3.0 ,
                        2.0 , 5.0 , 1.0 , 6.0 ] ,
        } ,
    )

    #----------------------------------
    ### Run tests: `count_variable`
    #----------------------------------
    ### Check shape
    assert eval_dataframe_monkey.shape == expected_dimensions
    assert eval_dataframe_function.shape == expected_dimensions
    ### Check output
    assert np.all( eval_dataframe_monkey == expected_output )
    assert np.all( eval_dataframe_function == expected_output )

def test_group_merge( ):

    ### Create mock dataframe 1
    test_dataframe_a = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 1 , 2 ] , 6 ) ,
            'animal': np.tile( [ 'big blue bass' , 'gnarly green grouse' , 'magenta mad manatee' ,
                                 'pretty pink pony' , 'roudy red rabbit' , 'silly silver silkworm' ] , 2 ) ,
            'insert_metric_here': [ 1.00 , 1.00 , 1.00 , 0.75 , 0.75 , 0.75 ,
                                    0.50 , 0.50 , 0.50 , 0.75 , 0.75 , 1.00 ] ,

        } ,
    )

    ### Create mock dataframe 2
    test_dataframe_b = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 1 , 2 ] , 6 ) ,
            'animal': np.tile( [ 'big blue bass' , 'gnarly green grouse' , 'magenta mad manatee' ,
                                 'pretty pink pony' , 'roudy red rabbit' , 'silly silver silkworm' ] , 2 ) ,
            'group': np.repeat( [ 'sleepy' , 'alert' ] , 6 ) ,
            'new_metric_here': [ 0.1 , 0.1 , 0.2 , 0.2 , 0.3 , 0.3 ,
                                 0.5 , 0.2 , 0.2 , 0.4 , 0.4 , 0.5 ] ,

        } ,
    )

    ### Create mock dataframe 3
    test_dataframe_c = pd.DataFrame(
        {   
            'stratum_num': np.repeat( [ 1 , 2 ] , 2 ) ,
            'group': np.tile( [ 'sleepy' , 'alert' ] , 2 ) ,
            'categorical_metric': np.tile( [ 'zippity' , 'doo' ] , 2 ) ,
        } ,
    )

    ### Evaluate for later comparison
    # ++++ Drop NA
    # ---- Monkey patch method (TEMPORARY)
    eval_dataframe_monkey_dropna = test_dataframe_a.group_merge( [ test_dataframe_b , test_dataframe_c ] , 
                                                                 inner_on = 'group' , 
                                                                 outer_on = [ 'stratum_num' ] ,
                                                                 drop_na = True )
    # ---- Normal function
    eval_dataframe_function_dropna = group_merge( test_dataframe_a ,
                                                  [ test_dataframe_b , test_dataframe_c ] , 
                                                  inner_on = 'group' , 
                                                  outer_on = [ 'stratum_num' ] ,
                                                  drop_na = True )
    # ++++ Don't drop NA
    # ---- Monkey patch method (TEMPORARY)
    eval_dataframe_monkey_keepna = test_dataframe_a.group_merge( [ test_dataframe_b , test_dataframe_c ] , 
                                                                 inner_on = 'group' , 
                                                                 outer_on = [ 'stratum_num' ] ,
                                                                 drop_na = False )
    # ---- Normal function
    eval_dataframe_function_keepna = group_merge( test_dataframe_a ,
                                                  [ test_dataframe_b , test_dataframe_c ] , 
                                                  inner_on = 'group' , 
                                                  outer_on = [ 'stratum_num' ] ,
                                                  drop_na = False )

    ###--------------------------------
    ### Expected outcomes
    ###--------------------------------
    # ---- Expected dimensions
    expected_dimensions_dropna = tuple( [ 12 , 6 ] )
    expected_dimensions_keepna = tuple( [ 14 , 6 ] )
    # ---- Expected output
    expected_output_dropna = pd.DataFrame(
        {
            'stratum_num': np.repeat( [ 1 , 2 ] , 6 ) ,
            'animal': np.tile( [ 'big blue bass' , 'gnarly green grouse' , 'magenta mad manatee' ,
                                 'pretty pink pony' , 'roudy red rabbit' , 'silly silver silkworm' ] , 2 ) ,
            'insert_metric_here': [ 1.00 , 1.00 , 1.00 , 0.75 , 0.75 , 0.75 ,
                                    0.50 , 0.50 , 0.50 , 0.75 , 0.75 , 1.00 ] ,
            'group': np.repeat( [ 'sleepy' , 'alert' ] , 6 ) ,
            'new_metric_here': [ 0.1 , 0.1 , 0.2 , 0.2 , 0.3 , 0.3 ,
                                 0.5 , 0.2 , 0.2 , 0.4 , 0.4 , 0.5 ] ,
            'categorical_metric': np.repeat( [ 'zippity' , 'doo' ] , 6 ) ,
        } ,
    )
    expected_output_keepna = pd.DataFrame(
        {
            'stratum_num': np.concatenate( [ np.repeat( [ 1 , 2 ] , 6 ) , [ 1 , 2 ] ] ) ,
            'animal': np.concatenate( [ np.tile( [ 'big blue bass' , 'gnarly green grouse' , 'magenta mad manatee' ,
                                                  'pretty pink pony' , 'roudy red rabbit' , 'silly silver silkworm' ] , 2 ) ,
                                      np.repeat( np.nan , 2 ).astype( 'object' ) ] ) ,
            'insert_metric_here': [ 1.00 , 1.00 , 1.00 , 0.75 , 0.75 , 0.75 ,
                                    0.50 , 0.50 , 0.50 , 0.75 , 0.75 , 1.00 ,
                                    np.nan , np.nan ] ,
            'group': np.concatenate( [ np.repeat( [ 'sleepy' , 'alert' ] , 6 ) ,
                                     [ 'alert' , 'sleepy' ] ] ) ,
            'new_metric_here': [ 0.1 , 0.1 , 0.2 , 0.2 , 0.3 , 0.3 ,
                                 0.5 , 0.2 , 0.2 , 0.4 , 0.4 , 0.5 , 
                                 np.nan , np.nan ] ,
            'categorical_metric': np.concatenate( [ np.repeat( [ 'zippity' , 'doo' ] , 6 ) ,
                                                  [ 'doo' , 'zippity' ] ] ) ,
        } ,
    )

    #----------------------------------
    ### Run tests: `count_variable`
    #----------------------------------
    ### Check shape
    # ++++ NaN removed
    assert eval_dataframe_monkey_dropna.shape == expected_dimensions_dropna
    assert eval_dataframe_function_dropna.shape == expected_dimensions_dropna
    # ++++ NaN kept
    assert eval_dataframe_monkey_keepna.shape == expected_dimensions_keepna
    assert eval_dataframe_function_keepna.shape == expected_dimensions_keepna
    ### Check output
    # ++++ NaN removed
    assert np.all( eval_dataframe_monkey_dropna == expected_output_dropna )
    assert np.all( eval_dataframe_function_dropna == expected_output_dropna )
    # ++++ NaN kept
    assert eval_dataframe_function_keepna.equals( expected_output_keepna )
    assert eval_dataframe_function_keepna.equals( expected_output_keepna )