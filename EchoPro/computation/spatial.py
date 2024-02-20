import pandas as pd
import numpy as np
import geopandas as gpd
import geopy.distance
from typing import Optional , Tuple
from scipy import interpolate , special
                                
def correct_transect_intervals( dataframe: pd.DataFrame ,
                                threshold: np.float64 = 0.05 ):
    """
    Calculate along-transect intervals and impute erroneous values

    Parameters
    ----------
    dataframe: pd.DataFrame
        DataFrame

    Notes
    -----
    This function calculates the along-track transect interval length and areas.
    It then 'searches' for possible erroneous values at the end of each line 
    and replaces/imputes with alternative lengths/areas.
    """
    return (
            dataframe
                        # Calculate along-transect interval distances
                .pipe( lambda df: df.assign( interval = df[ 'vessel_log_start' ].diff( periods = -1 ).abs( ) ) 
                                    .replace( np.nan , df[ 'vessel_log_end' ].iloc[ -1 ] - df[ 'vessel_log_start' ].iloc[ -1 ] ) )
                # Replace likely erroneous interval lengths associated with the edges of each transect
                .pipe
                ( lambda df: df.assign( median_interval = np.median( df[ 'interval' ] ) )
                                    .assign( interval = lambda x: np.where( np.abs( x[ 'interval' ] - x[ 'median_interval' ] > threshold ) ,
                                                                            x.vessel_log_end - x.vessel_log_start ,
                                                                            x.interval ) ) )
                # Calculate interval area
                .pipe( lambda df: df.assign( interval_area = df[ 'interval' ] * df[ 'transect_spacing' ] ) )                            
                # Keep dataframe tidy by only retaining the necessary columns/variables
                .loc[ : , [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' , 'haul_num' , 
                            'interval' , 'interval_area' , 'NASC_all_ages' , 'NASC_no_age1' ] ]
            )


def calculate_start_end_coordinates( group ,
                                     contrast ):
    """
    Calculates start and end latitude/longitude   

    Parameters
    ----------
    group: pd.DataFrameGroupBy
        Grouped DataFrame
    contrast: List
        Target contrast for grouping

    Notes
    -----
    This function calculates the bounding rectangle surrounding the latitude/longitude values
    for each grouped value (e.g. transect)
    """ 

    return (
        group
        .groupby( contrast )
        .apply( lambda x: pd.Series( { 'minimum_longitude': x['longitude'].min() , 
                                       'maximum_longitude': x['longitude'].max() ,
                                       'center_latitude': x[ 'latitude' ].mean() } ) )
        .reset_index( )
    )

def calculate_transect_distance( dataframe ,
                                 contrast = 'transect_num' ):
    """
    Calculates spatial features of each transect    

    Parameters
    ----------
    dataframe: pd.DataFrame
        DataFrame
    contrast: List
        Target contrast for grouping

    Notes
    -----
    This function calculates the bounding rectangle surrounding the latitude/longitude values
    for each transect and stratum, the average spacing between transects, approximate areas 
    relative to each transect, and the distance for each transect
    """ 
    
    ### Calculate mean transect spacinge
    transect_spacing = dataframe.groupby( contrast )[ 'transect_spacing' ].mean().reset_index()
    
    ### Calculate minimum/maximum longitude, mean latitude, spacing, and area
    return (
            dataframe
            .pipe( lambda df: calculate_start_end_coordinates( df , [ contrast ] ) )
            .assign(transect_distance=lambda x: x.apply( lambda row: geopy.distance.distance( 
                    ( row[ 'center_latitude' ] , row[ 'minimum_longitude' ] ) , 
                    ( row[ 'center_latitude' ] , row[ 'maximum_longitude' ] ) ).nm , 
                    axis=1 ) )
            .merge( transect_spacing , on = [ contrast ] )
            .assign( transect_area = lambda x: x.transect_distance * x.transect_spacing )
    )

def georeference( dataframe: pd.DataFrame ,
                  projection: str = 'epsg:4326' ):
    """
    Converts a DataFrame to a GeoDataFrame

    Parameters
    ----------
    dataframe: pd.DataFrame
        DataFrame
    projection: str
        Projection (EPSG) string

    Notes
    -----
    This function converts a DataFrame into a GeoDataFrame by searching the column
    names for 'longitude' and 'latitude'. It then converts these two columns into a
    geometry Point string that can be used for geospatial operations
    """

    ### Parse longitude and latitude column names
    lat_col = [ col for col in dataframe.columns if 'latitude' in col.lower( ) ][ 0 ]
    lon_col = [ col for col in dataframe.columns if 'longitude' in col.lower( ) ][ 0 ]
    
    ### Create GeoDataFrame
    gdf = gpd.GeoDataFrame( dataframe.copy( ) , 
                            geometry = gpd.points_from_xy( dataframe.copy( )[ lon_col ] , 
                                                           dataframe.copy( )[ lat_col ] ) ,
                            crs = projection )
    
    ### Carriage return
    return gdf

def transform_geometry( geodataframe: gpd.GeoDataFrame ,
                        geodataframe_reference: gpd.GeoDataFrame ,
                        longitude_reference: np.float64 ,
                        longitude_offset: np.float64 ,
                        latitude_offset: np.float64 ,
                        d_longitude: Optional[ np.float64 ] = None ,
                        d_latitude: Optional[ np.float64 ] = None , ) -> Tuple[ pd.DataFrame , 
                                                                                np.float64 , 
                                                                                np.float64 ]:
    
    ### Interpolate coordinates
    # This transformation will be applied the appropriate datasets
    coordinate_interp = interpolate.interp1d(
        geodataframe_reference.geometry.y ,
        geodataframe_reference.geometry.x ,
        kind = 'linear' ,
        bounds_error = False
    )

    ### Apply transformation to pre-existing coordinates
    geodataframe_transformed = (
        geodataframe
        .copy( )
        .assign( longitude_transformed = lambda df: df.geometry.x - coordinate_interp( df.geometry.y ) + longitude_reference )
        .pipe( lambda df: df.assign( geometry = gpd.points_from_xy( df.longitude_transformed , df.geometry.y ) ) )
    )

    ### Calculate geospatial distances in longitudinal and latitudinal axes
    # Based on 'd_longitude' and 'd_latitude' inputs
    # If 'd_longitude' and 'd_latitude' is not defined
    if ( d_longitude is None ) & ( d_latitude is None ):
        d_longitude = geodataframe_transformed.geometry.x.max( ) - geodataframe_transformed.geometry.x.min( )
        d_latitude = geodataframe_transformed.geometry.y.max( ) - geodataframe_transformed.geometry.y.min( )

    ### Transform coordinates from longitude/latitude to distance
    # Generate copies
    x_copy = geodataframe_transformed.geometry.x 
    y_copy = geodataframe_transformed.geometry.y 

    # Apply transformation
    x = np.cos( np.pi / 180.0 * y_copy ) * ( x_copy - longitude_offset ) / d_longitude
    y = ( y_copy - latitude_offset ) / d_latitude

    ### Assign columns to original geodataframe
    geodataframe_transformed[ 'x_transformed' ] = x
    geodataframe_transformed[ 'y_transformed' ] = y 

    ### Carriage return
    return geodataframe_transformed

#################################
### variogram models
## From the Matlab code 
## 01 = spherical
## 02 = exponential
## 03 = gaussian
## 04 = linear
## 05 = sinusoidal
## 06 = combined exponential-cosine
## --- The cosine term reduces the overall semivariogram values via exponential decay
## 07 = combined exponential-cosine
## --- The cosine term increases the overall semivariogram values
## 08 = combined exponential-squared cosine
## --- The cosine term increases the overall semivariogram values
## --- The exponential decay term was squared introduces curvature that can 
## --- yield a more complex decay pattern
## 09 = Bessel
## 10 = combined exponential-Bessel 
## 11 = combined exponential-squared Bessel
## --- Squared exponential decay 
## 12 = combined exponential-squared linear 
## --- Squared exponential decay 
## 13 = combined exponential-power Bessel
## --- Exponential decay is exponentiated by power p
# y(h) = o^2(1 - e^(-ah))Jv(Bh)
# y(h) = C(1 - e^(-ah))Jv(Bh) 
# y(h) = (sill - nugget)(1 - e^(-ah))Jv(Bh)
# y(h) = (sill - nugget)(1 - e^(-(1/ls)^exp_pow h))Jv(Bh)  
# y(h) = (sill - nugget)(1 - e^(-(1/ls)^exp_pow lag_vec))Jv(B lag_vec) 
# y(h) = (sill - nugget)(1 - e^(-(1/ls)^exp_pow lag_vec))Jv(ls_hole_eff lag_vec) 
# y(h) = (sill - nugget)(1 - e^(-(1/ls)^exp_pow lag_vec))Jv(ls_hole_eff lag_vec) + nugget -- accounts for discontinuity at zero lag distance (ie nugget effect)
### --- 
# ls = correlation range (where a = range and therefore: a = 1 / ls )
# h = lag_vec = lag distances
# ls_hole_eff = hole effect correlation range
# decay_power = decay function exponent


def composite_family( family: list = [ 'bessel' , 'exponential' ] , **kwargs ):
    ### Allowed family combinations
    composite_options = {
        ( 'bessel' , 'exponential' ): bessel_exponential ,
        ( 'cosine', 'exponential' ): cosine_exponential ,
        ( 'exponential' , 'linear' ): exponential_linear ,
    }

    ###
    families = [ model.lower( ) for model in family ]

    ### Sort -- alphabetical
    families.sort( )

    ###
    if tuple( families ) in composite_options:
        model_name = composite_options[ tuple( families ) ]

        return model_name( **kwargs )

def bessel_exponential( decay_power: np.float64 = 1.0 , **kwargs ):
    
    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_function = 1.0 - np.exp( - ( ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** decay_power ) )

    ### 
    hole_effect = special.j0( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * ( decay_function * hole_effect ) + kwargs[ 'nugget' ]

def cosine_exponential( decay_power: np.float64 = 1.0 , enhance_semivariance = False , **kwargs ):   

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_modifier = 1.0 if enhance_semivariance else -1.0
    decay_function = (
        1.0 + decay_modifier *
        np.exp( - ( ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** decay_power ) )
    )

    ### 
    hole_effect = np.cos( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * ( decay_function * hole_effect ) + kwargs[ 'nugget' ]

def exponential_linear( decay_power: np.float64 = 1.0 , **kwargs ):
    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_function = 1.0 - np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** decay_power )
 
    ### 
    hole_effect = 1.0 - kwargs[ 'hole_effect_range' ] * kwargs[ 'distance_lags' ] ** decay_power

    ###
    return partial_sill * ( decay_function * hole_effect ) + kwargs[ 'nugget' ]

def exponential( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_function = 1.0 - np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) )

    ###
    return partial_sill * decay_function + kwargs[ 'nugget' ]

def gaussian( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_function = 1.0 - np.exp( - ( kwargs[ 'distance_lags' ] / kwargs[ 'correlation_range' ] ) ** 2.0 )

    ###
    return partial_sill * decay_function + kwargs[ 'nugget' ]

def linear( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    return partial_sill * kwargs[ 'distance_lags' ] + kwargs[ 'nugget' ]

def sinusoid( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_function = 1.0 - np.sin( kwargs[ 'hole_effect_range' ] * kwargs[ 'distance_lags' ] ) / kwargs[ 'distance_lags' ]

    ###
    return partial_sill * decay_function + kwargs[ 'nugget' ]

def bessel( **kwargs ):

    ###
    partial_sill = kwargs[ 'sill' ] - kwargs[ 'nugget' ]

    ###
    decay_function = 1.0 - special.j0( kwargs[ 'hole_effect_range' ] *  kwargs[ 'distance_lags' ] )

    ###
    return partial_sill * decay_function + kwargs[ 'nugget' ]