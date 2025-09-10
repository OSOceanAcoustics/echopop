from typing import Tuple

import geopandas as gpd
import numpy as np
import pandas as pd


def utm_string_generator(longitude: float, latitude: float):
    """
    Generate UTM EPSG projection string from longitude/latitude coordinates.

    This function converts WGS84 coordinates to the appropriate UTM zone EPSG code
    string, automatically determining the correct UTM zone based on longitude and
    the hemisphere (north/south) based on latitude.

    Parameters
    ----------
    longitude : float
        Longitude coordinate in decimal degrees (WGS84)
    latitude : float
        Latitude coordinate in decimal degrees (WGS84)

    Returns
    -------
    str
        EPSG code string for the appropriate UTM zone (e.g., '32610' for UTM zone 10N)

    Examples
    --------
    >>> utm_code = utm_string_generator(-125.0, 48.0)
    >>> print(f"UTM EPSG code: {utm_code}")
    UTM EPSG code: 32610

    >>> utm_code = utm_string_generator(-125.0, -48.0)
    >>> print(f"UTM EPSG code: {utm_code}")
    UTM EPSG code: 32710

    Notes
    -----
    UTM zones are numbered from 1 to 60, with each zone spanning 6 degrees of longitude.
    The zone number is calculated as: floor((longitude + 180) / 6) + 1

    EPSG codes follow the pattern:
    - 326XX for northern hemisphere (where XX is the zero-padded zone number)
    - 327XX for southern hemisphere (where XX is the zero-padded zone number)
    """

    # Calculate UTM band value
    utm_value = str((np.floor((longitude + 180) / 6) % 60 + 1).astype(int))

    # Construct string to create equivalent EPSG code
    if len(utm_value) == 1:
        utm_value = "0" + utm_value

    if latitude >= 0.0:
        epsg = "326" + utm_value
    else:
        epsg = "327" + utm_value

    return epsg


def wgs84_to_utm(geodataframe: gpd.GeoDataFrame):
    """
    Transform a GeoDataFrame from WGS84 to the appropriate UTM coordinate system.

    This function automatically determines the correct UTM zone based on the median
    longitude and latitude of the GeoDataFrame and transforms the coordinate reference
    system (CRS) in place. The transformation improves accuracy for distance and area
    calculations.

    Parameters
    ----------
    geodataframe : gpd.GeoDataFrame
        GeoDataFrame containing spatial data with WGS84 coordinates. Must contain
        columns with 'lat' and 'long' in their names (case-insensitive).

    Returns
    -------
    None
        The function modifies the GeoDataFrame in place by changing its CRS.

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point
    >>> df = gpd.GeoDataFrame({
    ...     'longitude': [-125.0, -125.1],
    ...     'latitude': [48.0, 48.1],
    ...     'geometry': [Point(-125.0, 48.0), Point(-125.1, 48.1)]
    ... }, crs='epsg:4326')
    >>> print(f"Original CRS: {df.crs}")
    Original CRS: EPSG:4326
    >>> wgs84_to_utm(df)
    >>> print(f"Transformed CRS: {df.crs}")
    Transformed CRS: EPSG:32610

    Notes
    -----
    The function performs the following steps:
    1. Automatically detects longitude and latitude columns by name
    2. Calculates the median coordinates to determine the appropriate UTM zone
    3. Generates the UTM EPSG code using utm_string_generator()
    4. Transforms the GeoDataFrame to the new CRS in place

    The transformation uses the median coordinates to ensure the UTM zone is
    appropriate for the entire dataset, which is particularly important for
    datasets spanning multiple UTM zones.
    """

    # Detect the correct longitude and latitude coordinates
    # ---- Latitude
    lat_col = [col for col in geodataframe.columns if "lat" in col.lower()][0]
    # ---- Longitude
    lon_col = [col for col in geodataframe.columns if "long" in col.lower()][0]

    # Generate the equivalent UTM EPSG string
    utm_code = utm_string_generator(
        np.median(geodataframe[lon_col]), np.median(geodataframe[lat_col])
    )

    # Apply the CRS change
    geodataframe.to_crs(f"epsg:{utm_code}", inplace=True)


def reproject_dataset(
    data_df: pd.DataFrame,
    crs_out: str,
    coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
    projection: str = "epsg:4326",
) -> pd.DataFrame:
    """
    Transform coordinates using a new projection via GeoPandas.

    This function converts coordinates from one coordinate reference system (CRS) to another
    using GeoPandas for accurate cartographic projections. It creates a temporary GeoDataFrame
    to perform the transformation and returns the results as a regular DataFrame with new
    'x' and 'y' columns containing the projected coordinates.

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame containing coordinate data to be reprojected.
    crs_out : str
        Target Coordinate Reference System (CRS) string (e.g., 'epsg:32610' for UTM Zone 10N).
    coordinate_names : Tuple[str, str], default=("longitude", "latitude")
        Names of the coordinate columns in the input DataFrame. Expected format: (x_col, y_col).
    projection : str, default='epsg:4326'
        Input Coordinate Reference System (CRS) string representing the original projection
        of the coordinate data (default is WGS84 geographic coordinates).

    Returns
    -------
    pd.DataFrame
        DataFrame with the original data plus new 'x' and 'y' columns containing the
        reprojected coordinates in the target CRS.

    Examples
    --------
    >>> # Project from WGS84 to UTM Zone 10N
    >>> df_projected = reproject_dataset(
    ...     data_df=survey_data,
    ...     crs_out='epsg:32610',
    ...     coordinate_names=('longitude', 'latitude')
    ... )
    >>> print(df_projected[['x', 'y']].head())

    >>> # Project from UTM back to WGS84
    >>> df_geo = reproject_dataset(
    ...     data_df=utm_data,
    ...     crs_out='epsg:4326',
    ...     coordinate_names=('x', 'y'),
    ...     projection='epsg:32610'
    ... )

    Notes
    -----
    This function uses GeoPandas for accurate coordinate transformations, which is more
    reliable than simple mathematical conversions for cartographic projections.
    The geometry column created during processing is automatically dropped from the
    output DataFrame.
    """

    # Get the coordinate names
    x_coord, y_coord = coordinate_names

    # Convert DataFrame into GeoDataFrame
    gdf = gpd.GeoDataFrame(
        data_df, geometry=gpd.points_from_xy(data_df[x_coord], data_df[y_coord]), crs=projection
    )

    # Project to new CRS
    gdf_proj = gdf.to_crs(crs_out)

    # Add projected x/y columns
    df_out = gdf_proj.copy()
    df_out["x"] = gdf_proj.geometry.x
    df_out["y"] = gdf_proj.geometry.y

    # Return the reprojected data
    return df_out.drop(columns="geometry")
