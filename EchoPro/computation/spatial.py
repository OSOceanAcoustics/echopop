import pandas as pd
import numpy as np
import geopandas as gpd
import geopy.distance
from typing import Union
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

def transform_geometry( dataframe: pd.DataFrame ,
                        dataframe_reference: pd.DataFrame ,
                        longitude_reference: np.float64 ,
                        longitude_offset: np.float64 ,
                        latitude_offset: np.float64 ,
                        projection: str = 'epsg:4326' ,
                        range_output: bool = True ,
                        d_longitude: Optional[ np.float64 ] = None ,
                        d_latitude: Optional[ np.float64 ] = None , ) -> Union[ Tuple[ pd.DataFrame , 
                                                                                       np.float64 ,
                                                                                       np.float64 ] ,
                                                                                pd.DataFrame ]:
    """
    Transforms the geometry of a GeoDataFrame to reference coordinates

    Parameters
    ----------
    dataframe: pd.DataFrame
        DataFrame
    geodataframe_reference: gpd.GeoDataFrame
        Reference GeoDataFrame
    longitude_reference: np.float64
        Reference longitude used to 'center' all coordinates
    longitude_offset: np.float64
        Longitude value used to shift longitude coordinates
    latitude_offset: np.float64
        Latitude value used to shift latitude coordinates
    projection: str
        Projection (EPSG) string
    d_longitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates
    d_latitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates   
    """

    ### Parse longitude and latitude column names
    ref_lat_col = [ col for col in dataframe_reference.columns if 'latitude' in col.lower( ) ][ 0 ]
    ref_lon_col = [ col for col in dataframe_reference.columns if 'longitude' in col.lower( ) ][ 0 ]
    
    ### Interpolate coordinates
    # This transformation will be applied the appropriate datasets
    coordinate_interp = interpolate.interp1d(
        dataframe_reference[ ref_lat_col ] ,
        dataframe_reference[ ref_lon_col ] ,
        kind = 'linear' ,
        bounds_error = False
    )
    
    ### Parse longitude and latitude column names
    lat_col = [ col for col in dataframe.columns if 'latitude' in col.lower( ) ][ 0 ]
    lon_col = [ col for col in dataframe.columns if 'longitude' in col.lower( ) ][ 0 ]
    
    ### Reduce coordinates to prevent issues with unwieldly dataframes
    dataframe_geometry = (
        # dataframe.loc[ : , [ lon_col , lat_col ] ]
        dataframe[ [ lon_col , lat_col ] ]
        .drop_duplicates( )
        .assign( geometry = lambda x: gpd.points_from_xy( x[ lon_col ] , x[ lat_col ] ) )
    )
    
    ### Convert to GeoDataFrame
    geodataframe_geometry = gpd.GeoDataFrame( dataframe_geometry ,
                                              geometry = 'geometry' ,
                                              crs = projection )

    ### Apply transformation to pre-existing coordinates
    geodataframe_transformed = (
        geodataframe_geometry
        .assign( longitude_transformed = lambda df: ( 
                    df.geometry.x - coordinate_interp( df.geometry.y ) + longitude_reference ) 
                )
        .pipe( lambda df: df.assign( geometry = gpd.points_from_xy( df.longitude_transformed , df.geometry.y ) ) )
    )

    ### Calculate geospatial distances in longitudinal and latitudinal axes
    # Based on 'd_longitude' and 'd_latitude' inputs
    # If 'd_longitude' and 'd_latitude' is not defined
    if ( d_longitude is None ) & ( d_latitude is None ):
        d_longitude = geodataframe_transformed.geometry.x.max( ) - geodataframe_transformed.geometry.x.min( )
        d_latitude = geodataframe_transformed.geometry.y.max( ) - geodataframe_transformed.geometry.y.min( )

    ### Apply transformation to coordinates and assign to reduced geodataframe
    geodataframe_result = (
        geodataframe_transformed
        .assign( x_transformed = lambda dfs: ( np.cos( np.pi / 180.0 * dfs.geometry.y ) *
                                               ( dfs.geometry.x - longitude_offset ) /
                                               d_longitude ) ,
                 y_transformed = lambda dfs: ( dfs.geometry.y - latitude_offset ) /
                                               d_latitude ) )  

    ### Merge back with original input data in case of reduction
    if range_output:
        return geodataframe_result.merge( dataframe , on = [ lon_col , lat_col ] ) , d_longitude , d_latitude
    else:
        return geodataframe_result.merge( dataframe , on = [ lon_col , lat_col ] )
     
from scipy import spatial

def transform_geometry1( dataframe: gpd.GeoDataFrame ,
                         geodataframe_reference: gpd.GeoDataFrame ,
                         longitude_reference: np.float64 ,
                         longitude_offset: np.float64 ,
                         latitude_offset: np.float64 ,
                         projection: str = 'epsg:4326' ,
                         d_longitude: Optional[ np.float64 ] = None ,
                         d_latitude: Optional[ np.float64 ] = None ) -> Tuple[ pd.DataFrame,
                                                                               np.float64,
                                                                          np.float64]:
    """
    Transforms the geometry of a GeoDataFrame to reference coordinates

    Parameters
    ----------
    dataframe: gpd.GeoDataFrame
        DataFrame
    geodataframe_reference: gpd.GeoDataFrame
        Reference GeoDataFrame
    longitude_reference: np.float64
        Reference longitude used to 'center' all coordinates
    longitude_offset: np.float64
        Longitude value used to shift longitude coordinates
    latitude_offset: np.float64
        Latitude value used to shift latitude coordinates
    projection: str
        Projection (EPSG) string
    d_longitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates
    d_latitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates
    """

    # Build a spatial index for the reference GeoDataFrame
    reference_tree = spatial.cKDTree(geodataframe_reference.geometry.apply(lambda x: (x.x, x.y)).tolist())

    # Parse longitude and latitude column names
    lat_col = next(col for col in dataframe.columns if 'latitude' in col.lower())
    lon_col = next(col for col in dataframe.columns if 'longitude' in col.lower())

    # Reduce coordinates to prevent issues with unwieldly dataframes
    unique_coordinates = dataframe[[lon_col, lat_col]].drop_duplicates()

    # Create geometry column using Numpy
    unique_coordinates['geometry'] = gpd.points_from_xy(unique_coordinates[lon_col], unique_coordinates[lat_col])

    # Convert to GeoDataFrame
    geodataframe_geometry = gpd.GeoDataFrame(unique_coordinates,
                                              geometry='geometry',
                                              crs=projection)

    # Query the spatial index to interpolate coordinates
    _, indices = reference_tree.query(geodataframe_geometry.geometry.apply(lambda x: (x.x, x.y)).tolist())

    interpolated_lon = geodataframe_reference.geometry.x.iloc[indices].values
    interpolated_lat = geodataframe_reference.geometry.y.iloc[indices].values

    # Apply transformation to pre-existing coordinates
    geodataframe_transformed = geodataframe_geometry.assign(
        longitude_transformed=lambda df: df.geometry.x - interpolated_lon + longitude_reference
    )

    # Transform coordinates from longitude/latitude to distance using Numpy operations
    if (d_longitude is None) & (d_latitude is None):
        d_longitude = geodataframe_transformed.geometry.x.max() - geodataframe_transformed.geometry.x.min()
        d_latitude = geodataframe_transformed.geometry.y.max() - geodataframe_transformed.geometry.y.min()

    geodataframe_transformed['x_transformed'] = (
        np.cos(np.radians(geodataframe_transformed.geometry.y)) *
        (geodataframe_transformed.geometry.x - longitude_offset) / d_longitude
    )

    geodataframe_transformed['y_transformed'] = (
        (geodataframe_transformed.geometry.y - latitude_offset) / d_latitude
    )

    # Merge back with original input data in case of reduction
    return geodataframe_transformed.merge(dataframe, on=[lon_col, lat_col])


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


def composite_family( models: list = [ 'bessel' , 'exponential' ] , **kwargs ):
    ### Allowed family combinations
    composite_options = {
        ( 'bessel' , 'exponential' ): bessel_exponential ,
        ( 'cosine', 'exponential' ): cosine_exponential ,
        ( 'exponential' , 'linear' ): exponential_linear ,
    }

    ###
    families = [ model.lower( ) for model in models ]

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

def sine( **kwargs ):

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





def composite_family( models: list = [ 'bessel' , 'exponential' ] , **kwargs ):
    ### Allowed family combinations
    composite_options = {
        ( 'bessel' , 'exponential' ): bessel_exponential ,
        ( 'cosine', 'exponential' ): cosine_exponential ,
        ( 'exponential' , 'linear' ): exponential_linear ,
    }

    ###
    families = [ model.lower( ) for model in models ]

    ### Sort -- alphabetical
    families.sort( )

    ###
    if tuple( families ) in composite_options:
        model_name = composite_options[ tuple( families ) ]

        return model_name( **kwargs )