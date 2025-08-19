from typing import Tuple

import geopandas as gpd
import pandas as pd

from .projection import wgs84_to_utm
from .spatial import transect_extent

def hull_crop(
    transect_df: pd.DataFrame,
    mesh_df: pd.DataFrame,
    num_nearest_transects: int = 3,
    mesh_buffer_distance: float = 2.5,
    projection: str = "epsg:4326",
    coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
):
    """
    Crop the kriging mesh using convex hull polygons generated from survey transects.

    This function creates a survey boundary by generating convex hulls around each transect
    and its nearest neighbors, then filters the mesh to include only cells within the
    buffered survey area. This approach provides a more flexible alternative to
    region-based cropping.

    Parameters
    ----------
    transect_df : pd.DataFrame
        Georeferenced survey transect data used for defining the spatial extent for the kriging
        mesh grid. Must contain columns: 'longitude', 'latitude', 'transect_num'.
    mesh_df : pd.DataFrame
        Complete kriging mesh DataFrame that is subsequently cropped. Must contain columns:
        'longitude', 'latitude'.
    num_nearest_transects : int, default=3
        The number of nearest-neighbor transects used for defining the local extent around each
        transect. These convex hulls are then combined to generate the full survey extent hull.
        Higher values create more inclusive boundaries.
    mesh_buffer_distance : float, default=2.5
        Buffer distance in nautical miles applied to the survey polygon before filtering
        mesh cells. This ensures adequate coverage around the survey boundary.
    projection : str, default='epsg:4326'
        EPSG projection code for the input coordinate system. Default is WGS84.
    coordinate_names : Tuple[str, str], default=("longitude", "latitude")
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).

    Returns
    -------
    pd.DataFrame
        Cropped mesh DataFrame containing only cells within the buffered survey extent.
        The 'geometry' column is removed from the output.

    Examples
    --------
    >>> cropped_mesh = hull_crop(
    ...     transect_df, mesh_df,
    ...     num_nearest_transects=5,
    ...     mesh_buffer_distance=3.0
    ... )
    >>> print(f"Original mesh size: {len(mesh_df)}")
    >>> print(f"Cropped mesh size: {len(cropped_mesh)}")

    Notes
    -----
    The function performs the following steps:
    1. Converts the mesh DataFrame to a GeoDataFrame with point geometries
    2. Transforms coordinates from WGS84 to UTM for accurate distance calculations
    3. Generates survey extent polygon using transect_extent() function
    4. Applies buffer distance (converted from nautical miles to meters)
    5. Filters mesh cells to those within the buffered polygon
    6. Returns the filtered mesh without geometry column

    The UTM transformation ensures accurate distance calculations for the convex hull
    generation and buffering operations. The buffer distance helps ensure adequate
    mesh coverage around the survey boundary.

    This method is particularly useful for irregularly shaped survey areas where
    region-based cropping may be too restrictive or complex.
    """

    # Get coordinate names
    x_coord, y_coord = coordinate_names

    # Convert mesh DataFrame into a GeoDataframe
    mesh_gdf = gpd.GeoDataFrame(
        mesh_df,
        geometry=gpd.points_from_xy(mesh_df[x_coord], mesh_df[y_coord]),
        crs=projection,
    )

    # Convert the mesh projection to UTM (m)
    wgs84_to_utm(mesh_gdf)

    # Determine the survey extent by generating the border polygon
    survey_polygon = transect_extent(transect_df, projection, num_nearest_transects)

    # Find the mesh coordinates that fall within the buffered polygon
    # ---- Convert `grid_buffer` (nmi) to m and add buffer to polygon
    survey_polygon_buffered = survey_polygon.buffer(mesh_buffer_distance * 1852)
    # ---- Inclusion/union filter mask
    within_polygon_mask = mesh_gdf.geometry.within(survey_polygon_buffered)
    # ---- Apply mask to the mesh grid
    mesh_gdf_masked = mesh_gdf[within_polygon_mask]

    # Return the masked DataFrame
    return mesh_gdf_masked.drop(columns="geometry")
