import geopandas as gpd
import numpy as np


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
