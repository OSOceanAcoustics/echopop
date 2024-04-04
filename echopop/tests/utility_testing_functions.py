import numpy as np

def dictionary_shape( dictionary: dict ):
    """
    A utility test function that extracts the shape of a nested dictionary
    """
    if isinstance( dictionary , dict ):
        return( { i: dictionary_shape( dictionary[ i ] ) for i in dictionary } )
    else: 
        return None

def dictionary_shape_equal( dictionary1: dict ,
                            dictionary2: dict ):
    """
    Tests equality between the shapes of two nested dictionaries
    """
    result = dictionary_shape( dictionary1 ) == dictionary_shape( dictionary2 )
    
    if result:
        return result 
    else:
        if set( dictionary_shape( dictionary1 ) ) <= set( dictionary_shape( dictionary2 ) ):
            tracked_true = [ ]
            
            for j in dictionary2.keys( ):
                test = set( dictionary1[ j ].keys( ) ) <= ( dictionary2[ j ].keys( ) )
                tracked_true.append( test )
                
            if np.all( tracked_true ):
                return True 
            else: 
                return result 
        else:  
            return result