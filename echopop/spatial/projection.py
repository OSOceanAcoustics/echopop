import numpy as np
import pandas as pd
import geopandas as gpd

from scipy import interpolate
from shapely.geometry import Point , Polygon
from shapely.ops import unary_union
from typing import List, Union, Optional , Tuple

def utm_string_generator( longitude: float , latitude: float ) :
    """
    Converts projection string from longitude/latitude (WGS84) to equivalent UTM
    
    Parameters
    ----------
    longitude: float
        Longitude coordinate
    latitude: float
        Latitude coordinate
    """

    # Calculate UTM band value
    utm_value = str( ( np.floor( ( longitude + 180 ) / 6 ) % 60 + 1 ).astype( int ) )
    
    # Construct string to create equivalent EPSG code
    if len( utm_value ) == 1 :
        utm_value = '0' + utm_value
        
    if latitude >= 0.0 :
        epsg = '326' + utm_value
    else :
        epsg = '327' + utm_value
    
    return epsg    

def wgs84_to_utm( geodataframe: gpd.GeoDataFrame ) :
    """
    Changes the coordinate reference system (CRS) from WGS84 to UTM
    
    Parameters
    ----------
    geodataframe: float
        Longitude coordinate
    latitude: float
        Latitude coordinate
    """    
    
    # Detect the correct longitude and latitude coordinates
    # ---- Latitude
    lat_col = [ col for col in geodataframe.columns if 'lat' in col.lower( ) ][ 0 ]
    # ---- Longitude
    lon_col = [ col for col in geodataframe.columns if 'long' in col.lower( ) ][ 0 ]    

    # Generate the equivalent UTM EPSG string 
    utm_code = utm_string_generator( np.median( geodataframe[ lon_col ] ) , 
                                     np.median( geodataframe[ lat_col ] ) )
    
    # Apply the CRS change
    geodataframe.to_crs( f"epsg:{utm_code}" , inplace = True )

def mesh_crop( transect_data: pd.DataFrame ,
               mesh_grid: pd.DataFrame ,
               projection: str = 'epsg:4326' ,
               nearest_neighbors: int = 4 ,
               grid_buffer: float = 1.25 , ) :
    
    # Calculate the polygon extents of the entire survey
    survey_polygon = transect_extent( transect_data , projection , nearest_neighbors )

    # Create distance buffer that limits extrapolation
    # ---- Convert `grid_buffer` (nmi) to m
    survey_polygon_buffered = survey_polygon.buffer( grid_buffer * 1852 )

    # Detect the correct longitude and latitude coordinates of the mesh
    # ---- Latitude
    lat_col = [ col for col in mesh_grid.columns if 'lat' in col.lower( ) ][ 0 ]
    # ---- Longitude
    lon_col = [ col for col in mesh_grid.columns if 'long' in col.lower( ) ][ 0 ]  

    # Convert to GeoDataFrame
    mesh_grid_gdf = gpd.GeoDataFrame(
        mesh_grid ,
        geometry = gpd.points_from_xy( mesh_grid[ lon_col ] ,
                                       mesh_grid[ lat_col ] ) ,
        crs = projection ,       
    )

    # Convert from WGS84 to UTM
    wgs84_to_utm( mesh_grid_gdf )
    
    # Find the mesh coordinates contained within the buffered survey polygon
    within_polygon_index = mesh_grid_gdf.geometry.within( survey_polygon_buffered )
    
    # Filter the kriging mesh grid contained within the buffer survey polygon
    # ---- Apply within-polygon index to mesh grid
    mesh_grid_gdf_within = mesh_grid_gdf[ within_polygon_index ]
    # ---- Select necessary columns
    return mesh_grid_gdf_within[ [ 'fraction_cell_in_polygon' , 'centroid_longitude' , 'centroid_latitude' ] ]

def transform_geometry( dataframe: pd.DataFrame ,
                        reference_grid: pd.DataFrame ,
                        settings_dict: dict ,
                        delta_longitude: Optional[ np.float64 ] = None ,
                        delta_latitude: Optional[ np.float64 ] = None ) -> Tuple[ pd.DataFrame , float , float ] :
    """
    Transforms the geometry of a GeoDataFrame to reference coordinates

    Parameters
    ----------
    dataframe: pd.DataFrame
        DataFrame
    reference_grid: gpd.GeoDataFrame
        Reference GeoDataFrame
    settings_dict: dict
        Dictionary containing parameters defining longitudinal and latitudinal 
        reference and offset values
    delta_longitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates
    delta_latitude: np.float64
        Total longitudinal distance (degrees) used for standardizing coordinates   
    """

    # Parse the longitude and latitude columns for the reference grid
    # ---- Longitude
    lon_col = [ col for col in reference_grid.columns if 'lon' in col.lower( ) ][ 0 ]
    # ---- Latitude
    lat_col = [ col for col in reference_grid.columns if 'lat' in col.lower( ) ][ 0 ]
    # ---- Rename columns, if necessary
    reference = reference_grid.copy( ).rename( columns = { f"{lon_col }": 'longitude' ,
                                                           f"{lat_col}": 'latitude' } )
    
    # Extract the stratum and variable names
    # ---- Stratum
    stratum_col = settings_dict[ 'stratum_name' ]
    # ---- Biological variable
    variable_col = settings_dict[ 'variable' ]
    
    # Create interpolation function from reference grid coordinates (to interpolate longitude)
    reference_interp = interpolate.interp1d(
        reference[ 'latitude' ] ,
        reference[ 'longitude' ] ,
        kind = 'linear' ,
        bounds_error = False
    )

    # Apply offset to the longitudae and latitude coordinates
    # ---- Create copy
    dataframe_copy = dataframe.copy( )
    # ---- Longitude
    dataframe_copy[ 'longitude_transformed' ] = (
        dataframe_copy[ 'longitude' ] - reference_interp( dataframe_copy[ 'latitude' ] )
        + settings_dict[ 'kriging_parameters' ][ 'longitude_reference' ]
    )

    # Calculate the geospatial distances along the longitudinal and latitudinal axes
    if delta_longitude is None and delta_latitude is None:
        # ---- Longitude
        delta_longitude = (
            dataframe_copy[ 'longitude_transformed' ].max( ) 
            - dataframe_copy[ 'longitude_transformed' ].min( )
        )
        # ---- Latitude
        delta_latitude = dataframe_copy[ 'latitude' ].max( ) - dataframe_copy[ 'latitude' ].min( )

    # Standardize the x- and y-coordinates
    # ---- longitude --> x
    dataframe_copy[ 'x' ] = (
        np.cos( np.pi / 180.0 * dataframe_copy[ 'latitude' ] )
        * ( dataframe_copy[ 'longitude_transformed' ] 
           - settings_dict[ 'kriging_parameters' ][ 'longitude_offset' ] ) / delta_longitude
    )
    # ---- latitude --> y
    dataframe_copy[ 'y' ] = (
        ( dataframe_copy[ 'latitude' ] 
        - settings_dict[ 'kriging_parameters' ][ 'latitude_offset' ] )
        / delta_latitude
    )

    # Return the output tuple 
    return (
        dataframe_copy.filter( regex = "^(?!.*(_transformed))") ,
        delta_longitude ,
        delta_latitude
    )
 