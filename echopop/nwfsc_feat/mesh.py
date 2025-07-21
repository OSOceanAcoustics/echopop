from typing import Callable, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate

from .projection import wgs84_to_utm
from .spatial import transect_extent


# Helper functions for setting E-W/N-S assignments
def region_13_conditions(x, boundary_column: str, position: str):
    """
    Map mesh regions 1 and 3 boundaries based on encoded transect boundary values.

    This helper function interprets the encoded boundary values for regions 1 and 3,
    which use latitude-parallel transects. The encoding uses decimal fractions to
    indicate spatial positions (e.g., .1 for west, .4 for east, .6 for south, .9 for north).

    Parameters
    ----------
    x : pd.DataFrame
        DataFrame subset containing transect data for a specific transect number
    boundary_column : str
        Name of the column containing the encoded boundary values
    position : str
        Position identifier ('east', 'west', 'north', 'south') used for output column naming

    Returns
    -------
    pd.Series
        Series containing longitude and latitude coordinates for the specified position,
        with column names formatted as 'longitude_{position}' and 'latitude_{position}'

    Notes
    -----
    The encoding scheme:
    - x.1: Western boundary (minimum longitude)
    - x.4: Eastern boundary (maximum longitude)
    - x.6: Southern boundary (minimum latitude)
    - x.9: Northern boundary (maximum latitude)
    """
    # Calculate the floor
    x_floor = np.round((x[boundary_column].iloc[0] % 1) * 10)

    # Assign output based on conditions
    if x_floor == 1.0:
        return pd.Series(
            {
                f"longitude_{position}": x["longitude"].min(),
                f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmin()],
            }
        )
    elif x_floor == 4.0:
        return pd.Series(
            {
                f"longitude_{position}": x["longitude"].max(),
                f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmax()],
            }
        )
    elif x_floor == 6.0:
        return pd.Series(
            {
                f"longitude_{position}": x["longitude"].iloc[
                    x["latitude"].argmin(), f"latitude_{position}" : x["latitude"].min()
                ]
            }
        )
    elif x_floor == 9.0:
        return pd.Series(
            {
                f"longitude_{position}": x["longitude"].iloc[
                    x["latitude"].argmax(), f"latitude_{position}" : x["latitude"].max()
                ]
            }
        )


def region_2_conditions(x, boundary_column: str, position: str):
    """
    Map mesh region 2 boundaries based on encoded transect boundary values.

    This helper function interprets the encoded boundary values for region 2,
    which uses longitude-parallel transects. The encoding uses decimal fractions to
    indicate spatial positions.

    Parameters
    ----------
    x : pd.DataFrame
        DataFrame subset containing transect data for a specific transect number
    boundary_column : str
        Name of the column containing the encoded boundary values
    position : str
        Position identifier ('north', 'south') used for output column naming

    Returns
    -------
    pd.Series
        Series containing longitude and latitude coordinates for the specified position,
        with column names formatted as 'longitude_{position}' and 'latitude_{position}'

    Notes
    -----
    The encoding scheme for region 2:
    - x.1: Western boundary (minimum longitude)
    - x.4: Eastern boundary (maximum longitude)
    - x.6: Southern boundary (minimum latitude)
    - x.9: Northern boundary (maximum latitude)

    Region 2 transects run parallel to longitudes (north-south orientation).
    """
    # Calculate the floor
    x_floor = np.round((x[boundary_column].iloc[0] % 1) * 10)

    # Assign output based on conditions
    if x_floor == 1.0:
        return pd.Series(
            {
                f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmin()],
                f"longitude_{position}": x["longitude"].min(),
            }
        )
    elif x_floor == 4.0:
        return pd.Series(
            {
                f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmax()],
                f"longitude_{position}": x["longitude"].max(),
            }
        )
    elif x_floor == 6.0:
        return pd.Series(
            {
                f"latitude_{position}": x["latitude"].min(),
                f"longitude_{position}": x["longitude"].iloc[x["latitude"].argmin()],
            }
        )
    elif x_floor == 9.0:
        return pd.Series(
            {
                f"latitude_{position}": x["latitude"].max(),
                f"longitude_{position}": x["longitude"].iloc[x["latitude"].argmax()],
            }
        )


def transect_ends_crop(
    transect_df: pd.DataFrame,
    mesh_df: pd.DataFrame,
    latitude_resolution: float,
    transect_mesh_region_function: Callable,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Crop the kriging mesh by interpolating the eastern and western extents of survey transects
    partitioned into discrete regions via user-defined sorting functions.

    This function processes survey transect data to define spatial boundaries for mesh cropping.
    It divides transects into three regions, interpolates their boundaries, and identifies
    mesh cells that fall within the survey extent.

    Parameters
    ----------
    transect_df : pd.DataFrame
        Georeferenced survey transect data used for defining the spatial extent for the kriging
        mesh grid. Must contain columns: 'transect_num', 'longitude', 'latitude'.
    mesh_df : pd.DataFrame
        Complete kriging mesh DataFrame that is subsequently cropped. Must contain columns:
        'longitude', 'latitude'.
    latitude_resolution : float
        The latitudinal resolution (in degrees) used for the interpolation. This determines
        the spacing between interpolation points and affects the precision of boundary detection.
    transect_mesh_region_function : Callable
        A sorting function that maps specific transect numbers to their respective discretized
        mesh regions. The outputs of this function are expected to be:
        - 'transect_start' (np.number): the first transect number of a particular region,
        - 'transect_end' (np.number): the last transect number of a particular region,
        - 'transect_lower_bound' (List[np.number]): a list of encoded transect numbers indicating
        whether the transect is north, south, west, or east,
        - 'transect_upper_bound' (List[np.number]): a list of encoded transect numbers indicating
        whether the transect is north, south, west, or east.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        A tuple comprising:
        - Cropped kriging mesh DataFrame containing only cells within the survey extent
        - Annotated transect data with each transect number's respective mesh region assignment

    Examples
    --------
    >>> from echopop.nwfsc_feat.FEAT import transect_mesh_region_2019
    >>> cropped_mesh, annotated_transects = transect_ends_crop(
    ...     transect_df, mesh_df, 0.05, transect_mesh_region_2019
    ... )
    >>> print(f"Original mesh size: {len(mesh_df)}")
    >>> print(f"Cropped mesh size: {len(cropped_mesh)}")
    """

    # Create dictionary from user-defined transect-mesh region assignment
    mesh_region_dict = {
        region: {
            name: value
            for name, value in zip(
                ["start", "end", "upper", "lower"], transect_mesh_region_function(region)
            )
        }
        for region in [1, 2, 3]
    }

    # Create an equivalent DataFrame for merging
    mesh_region_df = pd.concat(
        [
            pd.DataFrame(
                {
                    "transect_num": np.arange(values["start"], values["end"] + 1),
                    "mesh_region": region,
                    "transect_lower_bound": values["lower"],
                    "transect_upper_bound": values["upper"],
                }
            )
            for region, values in mesh_region_dict.items()
        ]
    )
    # ---- Merge
    transect_df = transect_df.merge(mesh_region_df, on="transect_num", how="inner")
    # ---- Set indices
    transect_df.set_index(["mesh_region"], inplace=True)
    mesh_region_df.set_index(["mesh_region"], inplace=True)

    # Generate arrays that will be used for interpolation for each region
    # ---- Region 1
    region_1_latitude = np.arange(
        transect_df.loc[1, "latitude"].min(),
        transect_df.loc[1, "latitude"].max(),
        latitude_resolution,
    )
    # ---- Region 2
    # -------- Compute longitude resolution
    longitude_resolution_region_2 = latitude_resolution * np.cos(
        np.radians(transect_df.loc[2, "latitude"].mean())
    )
    # -------- Complete the array
    region_2_longitude = np.arange(
        transect_df.loc[2, "longitude"].min(),
        transect_df.loc[2, "longitude"].max(),
        longitude_resolution_region_2,
    )
    # ---- Region 3
    # -------- Get start and end transect numbers
    t3_start = mesh_region_dict[3]["start"]
    t3_end = mesh_region_dict[3]["end"]
    # -------- Complete the array
    region_3_latitude = np.arange(
        transect_df.loc[3]
        .loc[lambda x: x["transect_num"].isin([t3_start, t3_end]), "latitude"]
        .min(),
        transect_df.loc[3]
        .loc[lambda x: x["transect_num"].isin([t3_start, t3_end]), "latitude"]
        .max(),
        latitude_resolution,
    )

    # Process region 1
    # ---- Get eastern and western extents
    region_1_extents = pd.concat(
        [
            transect_df.loc[1]
            .groupby(["transect_num"])
            .apply(region_13_conditions, "transect_lower_bound", "east", include_groups=False),
            transect_df.loc[1]
            .groupby(["transect_num"])
            .apply(region_13_conditions, "transect_upper_bound", "west", include_groups=False),
        ],
        axis=1,
    )
    # ---- Generate interpolated coordinates
    region_1_interp = pd.DataFrame(
        {
            "longitude_west": interpolate.interp1d(
                region_1_extents["latitude_west"],
                region_1_extents["longitude_west"],
                kind="linear",
                bounds_error=False,
            )(region_1_latitude),
            "longitude_east": interpolate.interp1d(
                region_1_extents["latitude_east"],
                region_1_extents["longitude_east"],
                kind="linear",
                bounds_error=False,
            )(region_1_latitude),
        }
    )
    # ---- Calculate the longitudinal interval
    delta_longitude_region_1 = latitude_resolution * np.cos(np.radians(region_1_latitude))
    # ---- Initialize Region 1 index
    region_1_index = []
    # ---- Iterate through
    for i in range(len(region_1_latitude)):
        # -------- Find the mesh indices that are within the survey extent
        idx = np.where(
            (
                mesh_df["longitude"]
                >= region_1_interp.loc[i, "longitude_west"] - delta_longitude_region_1[i]
            )
            & (
                mesh_df["longitude"]
                <= region_1_interp.loc[i, "longitude_east"] + delta_longitude_region_1[i]
            )
            & (mesh_df["latitude"] >= region_1_latitude[i] - latitude_resolution)
            & (mesh_df["latitude"] < region_1_latitude[i] + latitude_resolution)
        )
        # -------- Append the indices
        region_1_index.append(idx[0])
    # ---- Assign the unique indices
    region_1_index_unique = np.unique(np.concatenate(region_1_index))

    # Process region 2
    # ---- Get northern and southern extents
    region_2_extents = pd.concat(
        [
            transect_df.loc[2]
            .groupby(["transect_num"])
            .apply(region_2_conditions, "transect_upper_bound", "south", include_groups=False),
            transect_df.loc[2]
            .groupby(["transect_num"])
            .apply(region_2_conditions, "transect_lower_bound", "north", include_groups=False),
        ],
        axis=1,
    )
    # ---- Generate interpolated coordinates
    region_2_interp = pd.DataFrame(
        {
            "latitude_south": interpolate.interp1d(
                region_2_extents["longitude_south"],
                region_2_extents["latitude_south"],
                kind="linear",
                bounds_error=False,
            )(region_2_longitude),
            "latitude_north": interpolate.interp1d(
                region_2_extents["longitude_north"],
                region_2_extents["latitude_north"],
                kind="linear",
                bounds_error=False,
            )(region_2_longitude),
        }
    )
    # ---- Calculate the longitudinal interval
    delta_longitude_region_2 = latitude_resolution * np.cos(
        np.radians(transect_df.loc[2, "latitude"].mean())
    )
    # ---- Initialize Region 2 index
    region_2_index = []
    # ---- Iterate through
    for i in range(len(region_2_longitude)):
        # -------- Find the mesh indices that are within the survey extent: southern limit
        if ~np.isnan(region_2_interp.loc[i, "latitude_south"]) or ~np.isnan(
            region_2_interp.loc[i, "latitude_north"]
        ):
            # -------- Compute the indices for the northern- and southernmost coordinates
            # -------- North
            lon_n_min = np.argmin(
                np.abs(region_2_longitude[i] - region_2_extents["longitude_north"])
            )
            # -------- South
            lon_s_min = np.argmin(
                np.abs(region_2_longitude[i] - region_2_extents["longitude_south"])
            )
            # -------- Slope
            slope = (
                region_2_extents["latitude_north"].iloc[lon_n_min]
                - region_2_extents["latitude_south"].iloc[lon_s_min]
            ) / (
                region_2_extents["longitude_north"].iloc[lon_n_min]
                - region_2_extents["longitude_south"].iloc[lon_s_min]
            )
            # -------- Set a new border threshold
            latitude_slope_i = (
                slope
                * (region_2_longitude[i] - region_2_extents["longitude_south"].iloc[lon_s_min])
                + region_2_extents["latitude_south"].iloc[lon_s_min]
            )
            if np.isnan(region_2_interp.loc[i, "latitude_south"]):
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_df["longitude"] >= region_2_longitude[i] - delta_longitude_region_2)
                    & (mesh_df["longitude"] <= region_2_longitude[i] + delta_longitude_region_2)
                    & (mesh_df["latitude"] >= latitude_slope_i - latitude_resolution)
                    & (
                        mesh_df["latitude"]
                        < region_2_interp.loc[i, "latitude_north"] + latitude_resolution
                    )
                )
            elif np.isnan(region_2_interp.loc[i, "latitude_north"]):
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_df["longitude"] >= region_2_longitude[i] - delta_longitude_region_2)
                    & (mesh_df["longitude"] <= region_2_longitude[i] + delta_longitude_region_2)
                    & (
                        mesh_df["latitude"]
                        >= region_2_interp.loc[i, "latitude_south"] - latitude_resolution
                    )
                    & (mesh_df["latitude"] < latitude_slope_i + latitude_resolution)
                )
            else:
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_df["longitude"] >= region_2_longitude[i] - delta_longitude_region_2)
                    & (mesh_df["longitude"] <= region_2_longitude[i] + delta_longitude_region_2)
                    & (
                        mesh_df["latitude"]
                        >= region_2_interp.loc[i, "latitude_south"] - latitude_resolution
                    )
                    & (
                        mesh_df["latitude"]
                        < region_2_interp.loc[i, "latitude_north"] + latitude_resolution
                    )
                )
            # -------- Append the indices
            region_2_index.append(idx[0])
    # ---- Assign the unique indices
    region_2_index_unique = np.unique(np.concatenate(region_2_index))

    # Process region 3
    # ---- Get eastern and western extents
    region_3_extents = pd.concat(
        [
            transect_df.loc[3]
            .groupby(["transect_num"])
            .apply(region_13_conditions, "transect_upper_bound", "west", include_groups=False),
            transect_df.loc[3]
            .groupby(["transect_num"])
            .apply(region_13_conditions, "transect_lower_bound", "east", include_groups=False),
        ],
        axis=1,
    )
    # ---- Generate interpolated coordinates
    region_3_interp = pd.DataFrame(
        {
            "longitude_west": interpolate.interp1d(
                region_3_extents["latitude_west"],
                region_3_extents["longitude_west"],
                kind="linear",
                bounds_error=False,
            )(region_3_latitude),
            "longitude_east": interpolate.interp1d(
                region_3_extents["latitude_east"],
                region_3_extents["longitude_east"],
                kind="linear",
                bounds_error=False,
            )(region_3_latitude),
        }
    )
    # ---- Calculate the longitudinal interval
    delta_longitude_region_3 = latitude_resolution * np.cos(np.radians(region_3_latitude))
    # ---- Initialize Region 1 index
    region_3_index = []
    # -------- Compute the indices for the eastern- and western-most coordinates
    # -------- West
    lat_w_max = np.argmax(region_3_extents["latitude_west"])
    # -------- East
    lat_e_max = np.argmax(region_3_extents["latitude_east"])
    # -------- Slope
    slope = (
        region_3_extents.iloc[lat_w_max]["longitude_west"]
        - region_3_extents.iloc[lat_e_max]["longitude_east"]
    ) / (
        region_3_extents.iloc[lat_w_max]["latitude_west"]
        - region_3_extents.iloc[lat_e_max]["latitude_east"]
    )
    # ---- Iterate through
    for i in range(len(region_3_latitude)):
        if ~np.isnan(region_3_interp.loc[i, "longitude_west"]) or ~np.isnan(
            region_3_interp.loc[i, "longitude_east"]
        ):
            # -------- Set a new border threshold
            longitude_slope_i = (
                slope * (region_3_latitude[i] - region_3_extents.iloc[lat_e_max]["latitude_east"])
                + region_3_extents.iloc[lat_w_max]["longitude_east"]
            )
            if np.isnan(region_3_interp.loc[i, "longitude_west"]):
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_df["longitude"] >= longitude_slope_i - delta_longitude_region_3[i])
                    & (
                        mesh_df["longitude"]
                        <= region_3_interp.loc[i, "longitude_east"] + delta_longitude_region_3[i]
                    )
                    & (mesh_df["latitude"] >= region_3_latitude[i] - latitude_resolution)
                    & (mesh_df["latitude"] < region_3_latitude[i] + latitude_resolution)
                )
            elif np.isnan(region_3_interp.loc[i, "longitude_east"]):
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (
                        mesh_df["longitude"]
                        >= region_3_interp.loc[i, "longitude_west"] - delta_longitude_region_3[i]
                    )
                    & (mesh_df["longitude"] <= longitude_slope_i + delta_longitude_region_3[i])
                    & (mesh_df["latitude"] >= region_3_latitude[i] - latitude_resolution)
                    & (mesh_df["latitude"] < region_3_latitude[i] + latitude_resolution)
                )
            else:
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (
                        mesh_df["longitude"]
                        >= region_3_interp.loc[i, "longitude_west"] - delta_longitude_region_3[i]
                    )
                    & (
                        mesh_df["longitude"]
                        <= region_3_interp.loc[i, "longitude_east"] + delta_longitude_region_3[i]
                    )
                    & (mesh_df["latitude"] >= region_3_latitude[i] - latitude_resolution)
                    & (mesh_df["latitude"] < region_3_latitude[i] + latitude_resolution)
                )
            # -------- Append the indices
            region_3_index.append(idx[0])
    # ---- Assign the unique indices
    region_3_index_unique = np.unique(np.concatenate(region_3_index))

    # Concatenate the indices
    mesh_indices = np.concatenate(
        [region_1_index_unique, region_2_index_unique, region_3_index_unique]
    )

    # Return the DataFrames
    return mesh_df.loc[mesh_indices], transect_df.reset_index()


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