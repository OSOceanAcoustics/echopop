import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.geometry import Point , Polygon
from shapely.ops import unary_union

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

def group_centroid( spatial_grouped: gpd.GeoSeries ) :
    """
    Calculate the centroid of a given spatial group
    
    Parameters
    ----------
    spatial_grouped: gpd.GeoSeries
        A GeoSeries comprising coordinates (i.e. points)
    """    

    # Compute the union of all coordinates within `spatial_grouped`
    centroid_point = spatial_grouped.unary_union.centroid
    
    # Return output
    return Point( centroid_point )

def transect_extent( spatial_data: pd.DataFrame ,
                     projection: str ,
                     nearest_neighbors: int , ) :

    # Detect the correct longitude and latitude coordinates
    # ---- Latitude
    lat_col = [ col for col in spatial_data.columns if 'lat' in col.lower( ) ][ 0 ]
    # ---- Longitude
    lon_col = [ col for col in spatial_data.columns if 'long' in col.lower( ) ][ 0 ]    

    # Drop duplicate coordinates
    # ---- Remove extraneous columns
    spatial_data_reduced = spatial_data.copy( )[ [ 'transect_num' , lon_col , lat_col ] ]
    # ---- Drop duplicates
    spatial_data_reduced.drop_duplicates( [ lon_col , lat_col ] , inplace = True )

    # Convert to GeoDataFrame
    spatial_data_gdf = gpd.GeoDataFrame(
        spatial_data_reduced ,
        geometry = gpd.points_from_xy( spatial_data_reduced[ lon_col ] ,
                                       spatial_data_reduced[ lat_col ] ) ,
        crs = projection ,
    )

    # Convert from WGS84 to UTM
    wgs84_to_utm( spatial_data_gdf )

    # Calculate the centroid of each transect line
    transect_centroid = spatial_data_gdf.groupby( 'transect_num' )[ 'geometry' ].apply( group_centroid )

    # Generate grouped polygons around each transect line
    # ---- Initialize polygon list
    transect_polygons = [ ]
    # ---- Iterate through each transect
    for transect in transect_centroid.index :
        # ---- Extract coordinates of the transect
        coord_centroid = transect_centroid[ transect ]
        # ---- Extract all remaining centroids
        other_centroids = transect_centroid[ transect_centroid.index != transect ].to_frame( )
        # ---- Calculate the distance between centroids
        other_centroids[ 'distance_centroid' ] = (
            other_centroids.geometry.apply( lambda g: coord_centroid.distance( g ) )
        )
        # ---- Find the 'n' nearest transect centroids
        nearest_centroids = other_centroids.distance_centroid.nsmallest( nearest_neighbors ) 
        # ---- Filter the transect centroids
        nearest_transects = other_centroids[ other_centroids.distance_centroid.isin( nearest_centroids ) ]
        # ---- Parse the coordinates of the relevant transect numbers
        unique_transects = np.append( nearest_transects.index , transect )
        # ---- Get the full coordinates of the relevant transects
        transect_coords = spatial_data_gdf[ spatial_data_gdf.transect_num.isin( unique_transects ) ]
        # ---- Generate the local polygon
        polygon = Polygon( list( transect_coords.geometry ) )
        # ---- Append the convex hull of the transect polygon to `transect_polygons`
        transect_polygons.append( polygon.convex_hull )
    
    # Merge the polygons via the union of each set
    return unary_union( transect_polygons )

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
