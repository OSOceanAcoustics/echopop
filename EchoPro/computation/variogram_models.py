
import numpy as np
from scipy import special
from typing import Union

### Dictionary containing available variogram models for user input
VARIOGRAM_MODELS = {
    'single': {
        'bessel': bessel ,
        'exponential': exponential ,
        'gaussian': gaussian ,
        'linear': linear ,
        'spherical': spherical ,
    } ,
    'composite': {
        ( 'bessel' , 'exponential' ): bessel_exponential ,
        ( 'bessel' , 'gaussian' ): bessel_gaussian ,
        ( 'cosine' , 'exponential' ): cosine_exponential ,
        ( 'cosine' , 'gaussian' ): cosine_gaussian ,
        ( 'exponential' , 'linear' ): exponential_linear ,
        ( 'gaussian' , 'linear' ): gaussian_linear ,
    }
}

def variogram( model: Union[ str , list ] = [ 'bessel' , 'exponential' ] , 
               **kwargs ):

    ### Convert to lowercase to match reference model dictionary
    ### And then convert to lowercase
    if isinstance( model , str ):
        model_input = model.lower( )
    elif isinstance( model , list ) & len( model ) == 1:
        model_input = ''.join( model ).lower( )
    else:
        model_input = [ name.lower( ) for name in model ]

         ### Alphabetic sort 
        model_input.sort( )   

    ### Parse user input from reference model dictionary
    # Check against VARIOGRAM_MODELS options to ensure model exists
    if ( len( model_input ) > 1 ) & ( tuple( model_input ) in VARIOGRAM_MODELS[ 'composite' ] ):

        ### Parse model function
        model_function = VARIOGRAM_MODELS[ 'composite' ][ tuple( model_input ) ]

    elif ( len( [ model_input ] ) == 1 ) & ( model_input in VARIOGRAM_MODELS[ 'single' ] ):

        ### Parse model function
        model_function = VARIOGRAM_MODELS[ 'single' ][ model_input ]
    
    ### Pass the additional user arguments (kwargs) to the child function
    return model_function( **kwargs )

#################################################################
# Single-family models
#################################################################
def bessel( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay = 1.0 - special.j0( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * decay + kwargs[ 'nugget' ]

def exponential( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay= 1.0 - np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) )
    ###
    return partial_sill * decay + kwargs[ 'nugget' ]

def gaussian( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay = 1.0 - np.exp( - ( kwargs[ 'distance_lags' ] ** 2 / kwargs[ 'correlation_range' ] ** 2.0 ) )

    ###
    return partial_sill * decay + kwargs[ 'nugget' ]

def linear( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    return partial_sill * kwargs[ 'distance_lags' ] + kwargs[ 'nugget' ]

def spherical( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay = (
        ( 3.0 * kwargs[ 'distance_lags' ] ) / ( 2.0 * kwargs[ 'correlation_range' ] ) -
        kwargs[ 'distance_lags' ] ** 3.0 / ( 2.0 * kwargs[ 'correlation_range' ] ** 3.0 )        
    )

    ###
    return np.where( kwargs[ 'distance_lags' ] <= kwargs[ 'correlation_range' ] ,
                     partial_sill * decay + kwargs[ 'nugget' ] ,
                     partial_sill )

#################################################################
# Composite-family models
#################################################################
def bessel_gaussian( **kwargs ):
    
    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay = 1.0 - np.exp( - ( ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** 2 ) )

    ### 
    hole_effect = special.j0( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * ( decay * hole_effect ) + kwargs[ 'nugget' ]

def bessel_exponential( decay_power: np.float64 = 1.0 , **kwargs ):
    
    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay = 1.0 - np.exp( - ( ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** decay_power ) )

    ### 
    hole_effect = special.j0( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * ( decay * hole_effect ) + kwargs[ 'nugget' ]

def cosine_exponential( enhance_semivariance = False , **kwargs ):   

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_modifier = - 1.0 if enhance_semivariance else 1.0
    decay = (
        decay_modifier *        
        np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) )
    )

    ### 
    hole_effect = np.cos( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * ( 1.0 - decay * hole_effect ) + kwargs[ 'nugget' ]

def cosine_gaussian(  **kwargs ):   

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay = (
        np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** 2 )
    )

    ### 
    hole_effect = np.cos( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * ( decay * hole_effect ) + kwargs[ 'nugget' ]

def exponential_linear( decay_power: np.float64 = 1.0 , **kwargs ):
    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_function = 1.0 - np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** decay_power )
 
    ### 
    hole_effect = 1.0 - kwargs[ 'hole_effect_range' ] * kwargs[ 'distance_lags' ] ** decay_power

    ###
    return partial_sill * ( decay_function * hole_effect ) + kwargs[ 'nugget' ]

def gaussian_linear( decay_power: np.float64 = 1.0 , **kwargs ):
    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay = 1.0 - np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** 2 )
 
    ### 
    hole_effect = 1.0 - kwargs[ 'hole_effect_range' ] * kwargs[ 'distance_lags' ] ** decay_power

    ###
    return partial_sill * ( decay * hole_effect ) + kwargs[ 'nugget' ]