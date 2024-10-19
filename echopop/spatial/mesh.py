from typing import Union

import geopandas as gpd
import geopy.distance
import numpy as np
import pandas as pd
from scipy import interpolate

from ..spatial.projection import wgs84_to_utm
from ..spatial.transect import transect_bearing, transect_extent


def crop_mesh(transect_data: pd.DataFrame, mesh_data: pd.DataFrame, cropping_parameters: dict):
    """
    Crop survey kriging mesh.

    Parameters
    ----------
    transect_data: pd.DataFrame
        Georeferenced transect data.
    mesh_data: pd.DataFrame
        Kriging mesh.
    cropping_parameters: dict
        Dictionary containing relevant algorithm variables and arguments.
    """
    # Rename the mesh coordinate names, if necessary
    # ---- Longitude
    mesh_longitude = [col for col in mesh_data.columns if "lon" in col.lower()][0]
    # ---- Latitude
    mesh_latitude = [col for col in mesh_data.columns if "lat" in col.lower()][0]
    # ---- Rename the dataframe
    mesh = mesh_data.copy().rename(
        columns={f"{mesh_longitude}": "longitude", f"{mesh_latitude}": "latitude"}
    )

    # Select and return cropped mesh depending on the cropping method
    # ---- Interpolation
    if cropping_parameters["crop_method"] == "transect_ends":
        return transect_ends_crop_method(transect_data.copy(), mesh, cropping_parameters)
    # ---- Convex hull
    elif cropping_parameters["crop_method"] == "convex_hull":
        return hull_crop_method(transect_data.copy(), mesh, cropping_parameters)


def hull_crop_method(transect_data: pd.DataFrame, mesh_data: pd.DataFrame, settings_dict: dict):
    """
    Crop the kriging mesh via convex hull polygons.

    Parameters
    ----------
    transect_data: pd.DataFrame
        Georeferenced transect data.
    mesh_data: pd.DataFrame
        Kriging mesh.
    settings_dict: dict
        Dictionary containing relevant algorithm variables and arguments.
    """
    # Extract the analysis settings
    # ---- Number of nearest transects
    num_nearest_transects = settings_dict["num_nearest_transect"]
    # ---- Grid buffer distance (nmi)
    mesh_buffer_distance = settings_dict["mesh_buffer_distance"]

    # Convert the mesh dataframes to a geodataframe
    mesh_gdf = gpd.GeoDataFrame(
        mesh_data,
        geometry=gpd.points_from_xy(mesh_data["longitude"], mesh_data["latitude"]),
        crs=settings_dict["projection"],
    )

    # Convert the mesh projection to UTM (m) -- In place
    wgs84_to_utm(mesh_gdf)

    # Determine the survey extent by generating the border polygon
    survey_polygon = transect_extent(
        transect_data, settings_dict["projection"], num_nearest_transects
    )

    # Find the mesh coordinates that fall within the buffered polygon
    # ---- Convert `grid_buffer` (nmi) to m and add buffer to polygon
    survey_polygon_buffered = survey_polygon.buffer(mesh_buffer_distance * 1852)
    # ---- Inclusion/union filter mask
    within_polygon_mask = mesh_gdf.geometry.within(survey_polygon_buffered)
    # ---- Apply mask to the mesh grid
    mesh_gdf_masked = mesh_gdf[within_polygon_mask]

    # Return the masked mesh dataframe
    return mesh_gdf_masked.drop(columns="geometry")


def transect_ends_crop_method(
    transect_data: pd.DataFrame, mesh_data: pd.DataFrame, cropping_parameters: dict
):
    """
    Crop the kriging mesh by interpolating the eastern and western boundaries of the survey
    partitioned into three geographical regions.

    Parameters
    ----------
    transect_data: pd.DataFrame
        Georeferenced transect data.
    mesh_data: pd.DataFrame
        Kriging mesh.
    cropping_parameters: dict
        Dictionary containing relevant algorithm variables and arguments.
    """

    # Extract the analysis settings
    # ---- Number of nearest transects
    latitude_resolution = cropping_parameters["latitude_resolution"]
    # ---- Grid buffer distance (nmi)
    bearing_tolerance = cropping_parameters["bearing_tolerance"]

    # Convert latitude resolution to degrees latitude
    latitude_resolution_deg = latitude_resolution / 60.0

    # Calculate the transect bearings
    transect_headings = transect_bearing(transect_data)

    # Find the transects that face north-to-south (Region 2)
    transect_headings_ns = transect_headings[
        (transect_headings["heading"] < bearing_tolerance)
        | (360.0 - transect_headings["heading"] < bearing_tolerance)
    ]
    # ---- Sub-sample the transect coordinates belonging to Region 2
    transect_data["mesh_region"] = np.where(
        transect_data["transect_num"] < transect_headings_ns["transect_num"].min(),
        1,
        np.where(transect_data["transect_num"] > transect_headings_ns["transect_num"].max(), 3, 2),
    )

    # Compute the transect extents across each region
    # ---- Mean latitude
    transect_data["latitude_mean"] = transect_data.groupby(["transect_num", "mesh_region"])[
        "latitude"
    ].transform("mean")
    # ---- Northernmost extent
    transect_data["latitude_north"] = transect_data.groupby(["transect_num", "mesh_region"])[
        "latitude"
    ].transform("max")
    # ---- Southernmost extent
    transect_data["latitude_south"] = transect_data.groupby(["transect_num", "mesh_region"])[
        "latitude"
    ].transform("min")
    # ---- Eastern extent
    transect_data["longitude_east"] = transect_data.groupby(["transect_num", "mesh_region"])[
        "longitude"
    ].transform("max")
    # ---- Westernmost extent
    transect_data["longitude_west"] = transect_data.groupby(["transect_num", "mesh_region"])[
        "longitude"
    ].transform("min")
    # ---- Index by region
    transect_data.set_index("mesh_region", inplace=True)

    # Generate arrays that will be used for interpolation for each region
    # ---- Region 1
    region_1_latitude = np.arange(
        transect_data.loc[1, "latitude"].min(),
        transect_data.loc[1, "latitude"].max(),
        latitude_resolution_deg,
    )
    # ---- Region 2
    # -------- Compute the requisite longitudinal resolution
    longitude_resolution_deg = latitude_resolution_deg * np.cos(
        np.radians(
            transect_data.loc[
                transect_data["transect_num"] == transect_headings_ns["transect_num"].max(),
                "latitude_mean",
            ].mean()
        )
    )
    # -------- Compute the array
    region_2_longitude = np.arange(
        transect_data.loc[2, "longitude_west"].min(),
        transect_data.loc[2, "longitude_east"].max(),
        longitude_resolution_deg,
    )
    # ---- Region 3
    region_3_latitude = np.arange(
        transect_data.loc[3, "latitude_south"].min(),
        transect_data.loc[3, "latitude_north"].max(),
        latitude_resolution_deg,
    )

    # Generate the new paired interpolated coordinates
    # ---- Region 1
    region_1_extents = interpolate_survey_extent(
        region_1_latitude, transect_data.loc[1], "latitude", "longitude"
    )
    # ---- Region 2
    region_2_extents = interpolate_survey_extent(
        region_2_longitude, transect_data.loc[2], "longitude", "latitude"
    )
    # ---- Region 3
    region_3_extents = interpolate_survey_extent(
        region_3_latitude, transect_data.loc[3], "latitude", "longitude"
    )

    # Iterate through each region to crop the mesh
    # ---- Region 1
    region_1_index = []
    # -------- Compute the change in longitude (degrees)
    delta_longitude = latitude_resolution_deg * np.cos(np.radians(region_1_latitude))
    # -------- Iterate through
    for i in range(len(delta_longitude)):
        # -------- Find the mesh indices that are within the survey extent
        idx = np.where(
            (mesh_data["longitude"] >= region_1_extents[0][i] - delta_longitude[i])
            & (mesh_data["longitude"] <= region_1_extents[1][i] + delta_longitude[i])
            & (mesh_data["latitude"] >= region_1_latitude[i] - latitude_resolution_deg)
            & (mesh_data["latitude"] < region_1_latitude[i] + latitude_resolution_deg)
        )
        # -------- Append the indices
        region_1_index.append(idx[0])
    # ---- Region 2
    region_2_index = []
    # -------- Extract the northern and southern components separately
    # -------- North
    transect_north = transect_data.loc[
        transect_data["latitude"] == transect_data["latitude_north"]
    ].loc[2]
    # -------- South
    transect_south = transect_data.loc[
        transect_data["latitude"] == transect_data["latitude_south"]
    ].loc[2]
    # -------- Iterate through
    for i in range(len(region_2_extents[0])):
        # -------- Find the mesh indices that are within the survey extent: southern limit
        if np.isnan(region_2_extents[0][i]) | np.isnan(region_2_extents[1][i]):
            # -------- Compute the indices for the northern- and southernmost coordinates
            # -------- North
            lon_n_min = np.argmin(np.abs(region_2_longitude[i] - transect_north["longitude"]))
            # -------- South
            lon_s_min = np.argmin(np.abs(region_2_longitude[i] - transect_south["longitude"]))
            # -------- Slope
            slope = (
                transect_north["latitude"].iloc[lon_n_min]
                - transect_south["latitude"].iloc[lon_s_min]
            ) / (
                transect_north["longitude"].iloc[lon_n_min]
                - transect_south["longitude"].iloc[lon_s_min]
            )
            # -------- Set a new border threshold
            latitude_slope_i = (
                slope * (region_2_longitude[i] - transect_south["longitude"].iloc[lon_s_min])
                + transect_south["latitude"].iloc[lon_s_min]
            )
            if np.isnan(region_2_extents[0][i]):
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_data["longitude"] >= region_2_longitude[i] - longitude_resolution_deg)
                    & (mesh_data["longitude"] <= region_2_longitude[i] + longitude_resolution_deg)
                    & (mesh_data["latitude"] >= latitude_slope_i - latitude_resolution_deg)
                    & (mesh_data["latitude"] < region_2_extents[1][i] + latitude_resolution_deg)
                )
            else:
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_data["longitude"] >= region_2_longitude[i] - longitude_resolution_deg)
                    & (mesh_data["longitude"] <= region_2_longitude[i] + longitude_resolution_deg)
                    & (mesh_data["latitude"] >= region_2_extents[0][i] - latitude_resolution_deg)
                    & (mesh_data["latitude"] < latitude_slope_i + latitude_resolution_deg)
                )
        else:
            # -------- Find the mesh indices that are within the survey extent
            idx = np.where(
                (mesh_data["longitude"] >= region_2_longitude[i] - longitude_resolution_deg)
                & (mesh_data["longitude"] <= region_2_longitude[i] + longitude_resolution_deg)
                & (mesh_data["latitude"] >= region_2_extents[0][i] - latitude_resolution_deg)
                & (mesh_data["latitude"] < region_2_extents[1][i] + latitude_resolution_deg)
            )
        # -------- Append the indices
        region_2_index.append(idx[0])
    # ---- Region 3
    region_3_index = []
    # -------- Compute the change in longitude (degrees)
    delta_longitude = latitude_resolution_deg * np.cos(np.radians(region_3_latitude))
    # -------- Extract the northern and southern components separately
    # -------- West
    transect_west = transect_data.loc[
        transect_data["longitude"] == transect_data["longitude_west"]
    ].loc[3]
    # -------- East
    transect_east = transect_data.loc[
        transect_data["longitude"] == transect_data["longitude_east"]
    ].loc[3]
    # -------- Iterate through
    for i in range(len(region_3_extents[0])):
        # -------- Find the mesh indices that are within the survey extent: southern limit
        if np.isnan(region_3_extents[0][i]) | np.isnan(region_3_extents[1][i]):
            # -------- Compute the indices for the northern- and southernmost coordinates
            # -------- North
            # lat_w_max = np.argmax(transect_west["latitude"])
            lat_w_max = np.argmin(transect_west["latitude"])
            # -------- South
            # lat_e_max = np.argmax(transect_east["latitude"])
            lat_e_max = np.argmin(transect_east["latitude"])
            # -------- Slope
            # slope = (
            #     transect_west["longitude"].iloc[lat_w_max]
            #     - transect_east["longitude"].iloc[lat_e_max]
            # ) / (transect_west["latitude"].min() - transect_east["latitude"].min())
            slope = (transect_west["latitude"].min() - transect_east["latitude"].min()) / (
                transect_west["longitude"].iloc[lat_w_max]
                - transect_east["longitude"].iloc[lat_e_max]
            )
            # -------- Set a new border threshold
            longitude_slope_i = (
                # slope * (region_3_latitude[i] - transect_east["latitude"].max())
                slope * (region_3_latitude[i] - transect_east["latitude"].min())
                + transect_east["longitude"].iloc[lat_e_max]
            )
            if np.isnan(region_3_extents[0][i]):
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_data["longitude"] >= longitude_slope_i - delta_longitude[i])
                    & (mesh_data["longitude"] <= region_3_extents[1][i] + delta_longitude[i])
                    & (mesh_data["latitude"] >= region_3_latitude[i] - latitude_resolution_deg)
                    & (mesh_data["latitude"] < region_3_latitude[i] + latitude_resolution_deg)
                )
            else:
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_data["longitude"] >= region_3_extents[0][i] - delta_longitude[i])
                    & (mesh_data["longitude"] <= longitude_slope_i + delta_longitude[i])
                    & (mesh_data["latitude"] >= region_3_latitude[i] - latitude_resolution_deg)
                    & (mesh_data["latitude"] < region_3_latitude[i] + latitude_resolution_deg)
                )
        else:
            # -------- Find the mesh indices that are within the survey extent
            idx = np.where(
                (mesh_data["longitude"] >= region_3_extents[0][i] - delta_longitude[i])
                & (mesh_data["longitude"] <= region_3_extents[1][i] + delta_longitude[i])
                & (mesh_data["latitude"] >= region_3_latitude[i] - latitude_resolution_deg)
                & (mesh_data["latitude"] < region_3_latitude[i] + latitude_resolution_deg)
            )
        # -------- Append the indices
        region_3_index.append(idx[0])
    # ---- Concatenate the region indices
    interpolated_indices = np.unique(
        np.concatenate(
            [
                np.concatenate(region_1_index),
                np.concatenate(region_2_index),
                np.concatenate(region_3_index),
            ]
        )
    )

    # Crop the mesh data and return the output
    return mesh_data.loc[interpolated_indices]


def griddify_lag_distances(
    coordinates_1: Union[pd.DataFrame, np.ndarray],
    coordinates_2: Union[pd.DataFrame, np.ndarray],
    angles: bool = False,
):
    """
    Calculate point-to-point distances between two gridded dataframes

    Parameters
    ----------
    coordinates_1: pd.DataFrame
        Background dataframe mesh that represents the "complete" field
        of values
    coordinates_2: pd.DataFrame
        Georeferenced dataframe
    angles: bool
        A boolean flag determining whether the point-to-point azimuth angles are also computed
        alongside the distances.

    Notes
    ----------
    This is used to effectively create a matrix comprising gridded
    distance values (i.e. 'griddify').
    """

    # If the inputs are dataframes
    if isinstance(coordinates_1, pd.DataFrame) and isinstance(coordinates_2, pd.DataFrame):
        # ---- Differences across x-coordinates
        x_distance = np.subtract.outer(coordinates_1["x"].to_numpy(), coordinates_2["x"].to_numpy())
        # ---- Differences across y-coordinates
        y_distance = np.subtract.outer(coordinates_1["y"].to_numpy(), coordinates_2["y"].to_numpy())
    # If the inputs are arrays
    elif isinstance(coordinates_1, np.ndarray) and isinstance(coordinates_2, np.ndarray):
        # ---- Differences across x-coordinates
        x_distance = np.subtract.outer(coordinates_1, coordinates_1)
        # ---- Differences across y-coordinates
        y_distance = np.subtract.outer(coordinates_2, coordinates_2)

    # Return Euclidean distances (and angles if appropriate)
    if angles:
        # ---- Copy x-array
        x_angles = x_distance.copy()
        # ---- Replace the self-points with NaN
        np.fill_diagonal(x_angles, np.nan)
        # ---- Copy y-array
        y_angles = y_distance.copy()
        # ---- Replace the self-points with NaN
        np.fill_diagonal(y_angles, np.nan)
        # ---- Calculate the azimuth angle grid
        angularity = np.arctan(y_angles / x_angles) * 180.0 / np.pi + 180 % 180
        # ---- Return output
        return np.sqrt(x_distance * x_distance + y_distance * y_distance), angularity
    else:
        # Return Euclidean distances
        return np.sqrt(x_distance * x_distance + y_distance * y_distance)


def stratify_mesh(input_dict: dict, kriged_mesh: pd.DataFrame, settings_dict: dict) -> pd.DataFrame:
    """
    Partition the kriging mesh into separate strata.

    Parameters
    ----------
    input_dict: dict
        Dictionary comprising data inputs.
    kriged_mesh: pd.DataFrame
        Kriging mesh.
    settings_dict: dict
    """

    # Extract the geographic-delimited strata
    if settings_dict["stratum"].lower() == "ks":
        geo_strata = input_dict["spatial"]["geo_strata_df"]
    elif settings_dict["stratum"].lower() == "inpfc":
        geo_strata = input_dict["spatial"]["inpfc_strata_df"]

    # Define the latitude bin array
    latitude_bins = np.concatenate([[-90.0], geo_strata["northlimit_latitude"], [90.0]])

    # Append the stratum variable to the kriging mesh
    # ---- Extract the stratum column name
    stratum_col = settings_dict["stratum_name"]
    # ---- Append a stratum column to the kriging mesh
    kriged_mesh[f"{stratum_col}"] = pd.cut(
        kriged_mesh["latitude"],
        latitude_bins,
        labels=list(geo_strata[f"{stratum_col}"]) + [1],
        ordered=False,
    )

    # Return output
    return kriged_mesh


def mesh_to_transects(
    kriging_dict: dict, spatial_dict: dict, settings_dict: dict
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Synthesize virtual transects from the kriging mesh.

    Parameters
    ----------
    kriging_dict: dict
        Dictionary comprising kriged (interpolated) data.
    spatial_dict: dict
        Dictionary containing spatial coordinate and stratum definitions.
    settings_dict: dict
    """
    # Extract mesh data
    mesh_data = kriging_dict["mesh_results_df"].copy()

    # Extract settings values
    # ---- Number of virtual transects per degree latitude
    transects_per_latitude = settings_dict["mesh_transects_per_latitude"]
    # ---- Stratum column name
    stratum_col = settings_dict["stratum_name"]
    # ---- Biological variable name
    variable_col = settings_dict["variable"]

    # Extract the appropriate stratum definitions/delimiters
    if settings_dict["stratum"].lower() == "inpfc":
        # ---- INPFC
        strata = spatial_dict["inpfc_strata_df"]
    elif settings_dict["stratum"].tolower() == "ks":
        # ---- KS
        strata = spatial_dict["geo_strata_df"]
    # ---- Create latitude bins
    latitude_bins = latitude_bins = np.concatenate([[-90.0], strata["northlimit_latitude"], [90.0]])

    # Partition the mesh into virtual/sytnthetic transects
    # ---- Get the unique latitudes of the virtual transects
    mesh_data.loc[:, "latitude"] = (
        np.round(mesh_data.loc[:, "latitude"] * transects_per_latitude + 0.5)
        / transects_per_latitude
    )
    # ---- Cut the virtual latitudes into their unique strata
    mesh_data.loc[:, f"{stratum_col}"] = pd.cut(
        mesh_data["latitude"],
        latitude_bins,
        labels=list(strata[f"{stratum_col}"]) + [1],
        ordered=False,
    )
    # ---- Get the unique latitude values
    unique_latitude_transect_key = pd.DataFrame(
        {
            "latitude": np.unique(mesh_data["latitude"]),
            "transect_num": np.arange(0, len(np.unique(mesh_data["latitude"])), 1),
        }
    ).set_index("latitude")
    # ---- Temporarily set `mesh_data` index
    mesh_data.set_index("latitude", inplace=True)
    # ---- Append the transect numkbers
    mesh_data["transect_num"] = unique_latitude_transect_key

    # Create equivalent transect dataframe needed for the stratified summary analysis
    virtual_transect_data = mesh_data.reset_index()[
        ["transect_num", "longitude", "latitude", stratum_col, "area", variable_col]
    ]
    # ---- Create density column
    virtual_transect_data[f"{variable_col}_density"] = (
        virtual_transect_data[variable_col] / virtual_transect_data["area"]
    )

    # Calculate the total virtual transect distance and area
    # ---- Initialize the dataframe
    virtual_transect_summary = (
        virtual_transect_data.drop_duplicates([stratum_col, "transect_num", "latitude"])[
            [stratum_col, "transect_num", "latitude"]
        ].set_index("transect_num")
    ).sort_index()
    # ---- Calculate the minimum longitude
    virtual_transect_summary["longitude_min"] = virtual_transect_data.groupby(["transect_num"])[
        "longitude"
    ].min()
    # ---- Calculate the mean longitude
    virtual_transect_summary["longitude_mean"] = virtual_transect_data.groupby(["transect_num"])[
        "longitude"
    ].mean()
    # ---- Calculate the minimum longitude
    virtual_transect_summary["longitude_min"] = virtual_transect_data.groupby(["transect_num"])[
        "longitude"
    ].min()
    # ---- Calculate the maximum longitude
    virtual_transect_summary["longitude_max"] = virtual_transect_data.groupby(["transect_num"])[
        "longitude"
    ].max()
    # ---- Calculate virtual transect distances
    virtual_transect_summary["transect_distance"] = virtual_transect_summary.apply(
        lambda row: geopy.distance.distance(
            (row["latitude"], row["longitude_min"]), (row["latitude"], row["longitude_max"])
        ).nm,
        axis=1,
    )
    # ---- Calculate the difference in latitude across the virtual transects (edge cases are handled
    # ---- using an additional step)
    virtual_transect_summary["d_latitude"] = np.where(
        ~virtual_transect_summary.index.isin([0, virtual_transect_summary.index[-1]]),
        np.concatenate([[np.nan], np.diff(virtual_transect_summary["latitude"])]),
        np.nanmean(np.diff(virtual_transect_summary["latitude"])),
    )
    # ---- Calculate the mean latitudinal distances to represent the mean transect spacing
    virtual_transect_summary["mean_spacing"] = virtual_transect_summary.apply(
        lambda row: geopy.distance.distance(
            (row["latitude"], row["longitude_mean"]),
            (row["latitude"] + row["d_latitude"], row["longitude_mean"]),
        ).nm,
        axis=1,
    )
    # ---- Calculate area
    virtual_transect_summary["transect_area"] = np.where(
        ~virtual_transect_summary.index.isin([0, virtual_transect_summary.index[-1]]),
        virtual_transect_summary["transect_distance"]
        * virtual_transect_summary["mean_spacing"]
        / 2.0,
        virtual_transect_summary["transect_distance"] * virtual_transect_summary["mean_spacing"],
    )
    # ---- Reset index
    virtual_transect_summary.reset_index(inplace=True)

    # Calculate the virtual strata summary
    # ---- Initialize the dataframe
    virtual_strata_summary = pd.DataFrame(
        {f"{stratum_col}": np.unique(virtual_transect_summary[stratum_col])},
    )
    # ---- Set index
    virtual_strata_summary.set_index(stratum_col, inplace=True)
    # ---- Calculate transect counts
    virtual_strata_summary["transect_count"] = virtual_transect_summary.groupby(
        [stratum_col], observed=False
    ).size()
    # ---- Calculate the total area
    virtual_strata_summary["transect_area_total"] = virtual_transect_summary.groupby(
        [stratum_col], observed=False
    )["transect_area"].sum()

    # Return the dataframes
    return virtual_transect_data, virtual_transect_summary, virtual_strata_summary.reset_index()


def interpolate_survey_extent(
    new_coords: np.ndarray, coordinate_data: pd.DataFrame, coordinates_x: str, coordinates_y: str
) -> tuple[np.ndarray, np.ndarray]:
    """
    Interpolate the eastern and western survey extent boundaries.

    Parameters
    ----------
    new_coords: np.ndarray
        New coordinates for interpolation.
    coordinate_data: pd.DataFrame
        Georeferenced points from the original dataset.
    coordinates_x: str
        'longitude' or 'latitude'
    coordinates_y: str
        'longitude' or 'latitude'
    """

    # Remove case-dependency
    coordinates_x = coordinates_x.lower()
    coordinates_y = coordinates_y.lower()

    # Error check
    if coordinates_x == coordinates_y:
        raise ValueError("Name for `coordinates_x` cannot be the same as `coordinates_y.")

    # Generate string that will be appended to the input strings
    if coordinates_y in ["longitude"]:
        add_string = ["_west", "_east"]
    else:
        add_string = ["_south", "_north"]

    # Generate the column strings
    # ---- South or West
    lower_col = coordinates_y + add_string[0]
    # ---- North or East
    upper_col = coordinates_y + add_string[1]

    # Reduce the dataframe
    # ---- South or West coordinates
    lower_coords = coordinate_data.loc[coordinate_data[coordinates_y] == coordinate_data[lower_col]]
    # ---- North or East coordinates
    upper_coords = coordinate_data.loc[coordinate_data[coordinates_y] == coordinate_data[upper_col]]

    # 1D interpolators
    # ---- South/West
    interpolator_lower = interpolate.interp1d(
        lower_coords[coordinates_x], lower_coords[coordinates_y], kind="linear", bounds_error=False
    )
    # ---- North/East
    interpolator_upper = interpolate.interp1d(
        upper_coords[coordinates_x], upper_coords[coordinates_y], kind="linear", bounds_error=False
    )

    # Apply the interpolators to the new coordinates and return the outputs
    return interpolator_lower(new_coords), interpolator_upper(new_coords)
