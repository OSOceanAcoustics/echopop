import warnings
from typing import Optional, Tuple, Union

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from ..validators.spatial import ValidateHullCropArgs
from .projection import wgs84_to_utm

# Set warnings filter
warnings.simplefilter("always")


def transform_coordinates(
    data: pd.DataFrame,
    x_offset: float = 0.0,
    y_offset: float = 0.0,
    coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
    reference: Optional[pd.DataFrame] = None,
    delta_x: Optional[float] = None,
    delta_y: Optional[float] = None,
) -> Tuple[pd.DataFrame, Union[float, None], Union[float, None]]:
    """
    Transform the x- and y-coordinates of a georeferenced dataset.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame with coordinates
    x_offset : float, default=0.
        Offset to apply to the x-coordinates that corresponds to `coordinate_names[0]`
    y_offset : float, default=0.
        Offset to apply to the y-coordinates that corresponds to `coordinate_names[0]`
    coordinate_names : Tuple[str, str], default=("longitude", "latitude")
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    reference : pd.DataFrame, optional
        Reference DataFrame with x and y coordinates for interpolation that is
        used as an additional offset to the x-axis.
    delta_x : float, optional
        Total x-axis distance used for standardizing coordinates. Will use the full range of the
        original x-axis if not provided.
    delta_y : float, optional
        Total y-axis distance used for standardizing coordinates. Will use the full range of the
        original y-axis if not provided.

    Returns
    -------
    pd.DataFrame
        DataFrame with the new transformed coordinates 'x' and 'y'.
    float or None
        Distance of the pre-transformed x-axis coordinates that can be used to transform other
        georeferenced datasets (assuming shared projections).
    float or None
        Distance of the pre-transformed y-axis coordinates that can be used to transform other
        georeferenced datasets (assuming shared projections).
    """

    # Get the coordinate names
    x_coord, y_coord = coordinate_names

    # Create interpolation function from reference grid coordinates (to interpolate longitude)
    if reference is not None:
        reference_interp = interpolate.interp1d(
            reference[y_coord], reference[x_coord], kind="linear", bounds_error=False
        )
        reference_offset = reference_interp(data[y_coord])
    else:
        reference_offset = 0.0

    # Transform longitude
    transformed_x = data[x_coord] - reference_offset + x_offset

    # Calculate the geospatial distances along the x- and y-axes [if missing]
    # ---- Longitude
    if delta_x is None:
        delta_x = transformed_x.max() - transformed_x.min()
    # ---- Latitude
    if delta_y is None:
        delta_y = data[y_coord].max() - data[y_coord].min()

    # Transform the x- and y-coordinates
    # ---- x
    data["x"] = np.cos(np.pi / 180.0 * data[y_coord]) * (transformed_x - x_offset) / delta_x
    # ---- y
    data["y"] = (data[y_coord] - y_offset) / delta_y

    # Return the output tuple
    return (data, delta_x, delta_y)


def transect_coordinate_centroid(spatial_grouped: gpd.GeoSeries):
    """
    Calculate the centroid of a given spatial group.

    This function computes the geometric centroid of a collection of spatial points,
    which is useful for determining the center point of transect lines or other
    spatial groupings.

    Parameters
    ----------
    spatial_grouped: gpd.GeoSeries
        A GeoSeries comprising coordinates (i.e. points). Each element should be
        a Point geometry representing spatial locations.

    Returns
    -------
    Point
        A shapely Point object representing the centroid of all input coordinates.

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point
    >>> points = gpd.GeoSeries([Point(0, 0), Point(1, 1), Point(2, 0)])
    >>> centroid = transect_coordinate_centroid(points)
    >>> print(f"Centroid: ({centroid.x:.1f}, {centroid.y:.1f})")
    Centroid: (1.0, 0.3)

    Notes
    -----
    The function uses the union_all() method to combine all geometries before
    calculating the centroid, which ensures proper handling of the spatial
    reference system.
    """

    # Compute the union of all coordinates within `spatial_grouped`
    centroid_point = spatial_grouped.union_all().centroid

    # Return output
    return Point(centroid_point)


def transect_extent(transects: pd.DataFrame, projection: str, num_nearest_transects: int, **kwargs):
    """
    Compute the spatial extent of survey transects using convex hull generation.

    This function creates a polygon representing the spatial extent of survey transects
    by generating convex hulls around each transect and its nearest neighbors, then
    unioning all hulls to create the overall survey boundary.

    Parameters
    ----------
    transects : pd.DataFrame
        Dataframe containing survey transect data with columns:
        - 'longitude': Longitude coordinates
        - 'latitude': Latitude coordinates
        - 'transect_num': Transect identifier numbers
    projection : str
        EPSG projection code string (e.g., 'epsg:4326' for WGS84)
    num_nearest_transects : int
        Number of nearest neighbor transects to include when generating
        the convex hull around each transect

    Returns
    -------
    shapely.geometry.base.BaseGeometry
        A shapely geometry object representing the union of all transect convex hulls,
        defining the overall spatial extent of the survey area.

    Examples
    --------
    >>> import pandas as pd
    >>> transect_data = pd.DataFrame({
    ...     'longitude': [-125.0, -125.1, -125.2],
    ...     'latitude': [48.0, 48.1, 48.2],
    ...     'transect_num': [1, 2, 3]
    ... })
    >>> extent = transect_extent(transect_data, 'epsg:4326', 2)
    >>> print(f"Extent type: {type(extent)}")
    Extent type: <class 'shapely.geometry.polygon.Polygon'>

    Notes
    -----
    The function performs the following steps:
    1. Converts the DataFrame to a GeoDataFrame with point geometries
    2. Transforms coordinates from WGS84 to UTM for accurate distance calculations
    3. Calculates centroids for each transect
    4. For each transect, finds the nearest neighbor transects
    5. Generates convex hulls around each transect and its neighbors
    6. Returns the union of all convex hulls

    The resulting polygon can be used for spatial filtering, mesh cropping,
    or defining survey boundaries for analysis.
    """

    # Copy
    transect_df = transects.copy()

    # Convert to GeoDataFrame
    transect_gdf = gpd.GeoDataFrame(
        transects,
        geometry=gpd.points_from_xy(transects["longitude"], transect_df["latitude"]),
        crs=projection,
    )

    # Convert from WGS84 to UTM
    wgs84_to_utm(transect_gdf)

    # Calculate the centroid of each transect line
    transect_centroid = transect_gdf.groupby("transect_num")["geometry"].apply(
        transect_coordinate_centroid
    )

    # Generate grouped polygons around each transect line
    # ---- Initialize polygon list
    transect_polygons = []
    # ---- Iterate through each transect
    for transect in transect_centroid.index:
        # ---- Extract coordinates of the transect
        coord_centroid = transect_centroid[transect]
        # ---- Extract all remaining centroids
        other_centroids = transect_centroid[transect_centroid.index != transect].to_frame()
        # ---- Handle case where there's only one transect (no other centroids)
        if len(other_centroids) == 0:
            # -------- Just use the current transect to create a polygon
            unique_transects = np.array([transect])
            transect_coords = transect_gdf[transect_gdf.transect_num.isin(unique_transects)]
            polygon = Polygon(list(transect_coords.geometry))
            transect_polygons.append(polygon.convex_hull)
            continue
        # ---- Calculate the distance between centroids
        other_centroids["distance_centroid"] = other_centroids.geometry.apply(
            lambda g: coord_centroid.distance(g)
        )
        # ---- Find the 'n' nearest transect centroids
        nearest_centroids = other_centroids.distance_centroid.nsmallest(num_nearest_transects)
        # ---- Filter the transect centroids
        nearest_transects = other_centroids[
            other_centroids.distance_centroid.isin(nearest_centroids)
        ]
        # ---- Parse the coordinates of the relevant transect numbers
        unique_transects = np.append(nearest_transects.index, transect)
        # ---- Get the full coordinates of the relevant transects
        transect_coords = transect_gdf[transect_gdf.transect_num.isin(unique_transects)]
        # ---- Generate the local polygon
        polygon = Polygon(list(transect_coords.geometry))
        # ---- Append the convex hull of the transect polygon to `transect_polygons`
        transect_polygons.append(polygon.convex_hull)

    # Merge the polygons via the union of each set
    return unary_union(transect_polygons)


def hull_crop(
    transects: pd.DataFrame,
    mesh: pd.DataFrame,
    num_nearest_transects: int = 3,
    mesh_buffer_distance: float = 2.5,
    projection: str = "epsg:4326",
    coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
) -> pd.DataFrame:
    """
    Crop the kriging mesh using convex hull polygons generated from survey transects.

    This function creates a survey boundary by generating convex hulls around each transect
    and its nearest neighbors, then filters the mesh to include only cells within the
    buffered survey area. This approach provides a more flexible alternative to
    region-based cropping.

    Parameters
    ----------
    transects : pd.DataFrame
        Georeferenced survey transect data used for defining the spatial extent for the kriging
        mesh grid. Must contain columns: 'longitude', 'latitude', 'transect_num'.
    mesh : pd.DataFrame
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
    ...     transects, mesh,
    ...     num_nearest_transects=5,
    ...     mesh_buffer_distance=3.0
    ... )
    >>> print(f"Original mesh size: {len(mesh)}")
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

    # Validate parameters
    try:
        valid_params = ValidateHullCropArgs.create(
            **dict(
                transects=transects,
                mesh=mesh,
                num_nearest_transects=num_nearest_transects,
                mesh_buffer_distance=mesh_buffer_distance,
                projection=projection,
                coordinate_names=coordinate_names,
            )
        )
    except Exception as e:
        raise e from None

    # Get coordinate names
    x_coord, y_coord = valid_params["coordinate_names"]

    # Convert mesh DataFrame into a GeoDataframe
    mesh_gdf = gpd.GeoDataFrame(
        valid_params["mesh"],
        geometry=gpd.points_from_xy(valid_params["mesh"][x_coord], valid_params["mesh"][y_coord]),
        crs=valid_params["projection"],
    )

    # Convert the mesh projection to UTM (m)
    wgs84_to_utm(mesh_gdf)

    # Determine the survey extent by generating the border polygon
    survey_polygon = transect_extent(**valid_params)

    # Find the mesh coordinates that fall within the buffered polygon
    # ---- Convert `grid_buffer` (nmi) to m and add buffer to polygon
    survey_polygon_buffered = survey_polygon.buffer(valid_params["mesh_buffer_distance"] * 1852)
    # ---- Inclusion/union filter mask
    within_polygon_mask = mesh_gdf.geometry.within(survey_polygon_buffered)
    # ---- Apply mask to the mesh grid
    mesh_gdf_masked = mesh_gdf[within_polygon_mask]

    # Return the masked DataFrame
    return mesh_gdf_masked.drop(columns="geometry")
