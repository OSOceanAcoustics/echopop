import geopandas as gpd
import pandas as pd
import numpy as np

from typing import Union

from ..spatial.projection import wgs84_to_utm
from ..spatial.transect import transect_extent

def crop_mesh( transect_data: pd.DataFrame ,
               mesh_data: pd.DataFrame ,               
               settings_dict: dict ):
    
    # Extract the analysis settings
    # ---- Number of nearest transects
    num_nearest_transects = settings_dict[ 'num_nearest_transect' ]
    # ---- Grid buffer distance (nmi)
    mesh_buffer_distance = settings_dict[ 'mesh_buffer_distance' ]

    # Rename the mesh coordinate names, if necessary
    # ---- Longitude
    mesh_longitude = (
        [ col for col in mesh_data.columns if 'lon' in col.lower( ) ][ 0 ]
    )  
    # ---- Latitude
    mesh_latitude = (
        [ col for col in mesh_data.columns if 'lat' in col.lower( ) ][ 0 ]
    )
    # ---- Rename the dataframe
    mesh = mesh_data.copy( ).rename( columns = { f"{mesh_longitude}": 'longitude' ,
                                                 f"{mesh_latitude}": 'latitude' } )
    
    # Convert the mesh dataframes to a geodataframe
    mesh_gdf = gpd.GeoDataFrame( mesh ,
                                 geometry = gpd.points_from_xy( mesh[ 'longitude' ] ,
                                                                mesh[ 'latitude' ] ) ,
                                 crs = settings_dict[ 'projection' ] )
    
    # Convert the mesh projection to UTM (m) -- In place
    wgs84_to_utm( mesh_gdf )
    
    # Determine the survey extent by generating the border polygon    
    survey_polygon = transect_extent( transect_data , 
                                      settings_dict[ 'projection' ] , 
                                      num_nearest_transects )
    
    # Find the mesh coordinates that fall within the buffered polygon
    # ---- Convert `grid_buffer` (nmi) to m and add buffer to polygon
    survey_polygon_buffered = survey_polygon.buffer( mesh_buffer_distance * 1852 ) 
    # ---- Inclusion/union filter mask
    within_polygon_mask = mesh_gdf.geometry.within( survey_polygon_buffered )
    # ---- Apply mask to the mesh grid
    mesh_gdf_masked = mesh_gdf[ within_polygon_mask ]

    # Return the masked mesh dataframe
    return mesh_gdf_masked.drop( columns = 'geometry' )

def griddify_lag_distances( coordinates_1: Union[pd.DataFrame , np.ndarray] ,
                            coordinates_2: Union[pd.DataFrame , np.ndarray] ):
    """
    Calculate point-to-point distances between two gridded dataframes

    Parameters
    ----------
    mesh_grid: pd.DataFrame
        Background dataframe mesh that represents the "complete" field 
        of values
    spatial_data: pd.DataFrame
        Georeferenced dataframe
        
    Notes
    ----------
    This is used to effectively create a matrix comprising gridded 
    distance values (i.e. 'griddify').
    """

    # If the inputs are dataframes
    if isinstance( coordinates_1 , pd.DataFrame ) and isinstance( coordinates_2 , pd.DataFrame ):
        # ---- Differences across x-coordinates
        x_distance = np.subtract.outer( coordinates_1[ 'x' ].to_numpy( ) , 
                                        coordinates_2[ 'x' ].to_numpy( ) )
        # ---- Differences across y-coordinates
        y_distance = np.subtract.outer( coordinates_1[ 'y' ].to_numpy( ) , 
                                        coordinates_2[ 'y' ].to_numpy( ) )
    # If the inputs are arrays
    elif isinstance( coordinates_1 , np.ndarray ) and isinstance( coordinates_2 , np.ndarray ):
        # ---- Differences across x-coordinates
        x_distance = np.subtract.outer( coordinates_1 , coordinates_1 )
        # ---- Differences across y-coordinates
        y_distance = np.subtract.outer( coordinates_2 , coordinates_2 )    

    # Return Euclidean distances
    return np.sqrt( x_distance * x_distance + y_distance * y_distance )

def stratify_mesh( input_dict: dict ,
                   kriged_mesh: pd.DataFrame ,
                   settings_dict: dict ) -> pd.DataFrame:
        
    # Extract the geographic-delimited strata
    if settings_dict[ 'stratum' ].lower( ) == 'ks': 
        geo_strata = input_dict[ 'spatial' ][ 'geo_strata_df' ]
    elif settings_dict[ 'stratum' ].lower( ) == 'inpfc': 
        geo_strata = input_dict[ 'spatial' ][ 'inpfc_strata_df' ]  

    # Define the latitude bin array
    latitude_bins = np.concatenate( [ [ -90.0 ] , geo_strata[ 'northlimit_latitude' ] , [ 90.0 ] ] )

    # Append the stratum variable to the kriging mesh
    # ---- Extract the stratum column name
    stratum_col = settings_dict[ 'stratum_name' ]
    # ---- Append a stratum column to the kriging mesh
    kriged_mesh[ f"{stratum_col}"] = pd.cut( kriged_mesh[ 'latitude' ] ,
                                             latitude_bins ,
                                             labels = list( geo_strata[ f"{stratum_col}" ] ) + [ 1 ] ,
                                             ordered = False )
    
    # Return output
    return kriged_mesh
