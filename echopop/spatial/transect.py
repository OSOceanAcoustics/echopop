import geopandas as gpd
import geopy.distance
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from ..spatial.projection import wgs84_to_utm


def correct_transect_intervals(transect_data: pd.DataFrame, interval_threshold: float = 0.05):
    """
    Calculate along-transect intervals and impute erroneous values

    Parameters
    ----------
    transect_data: pd.DataFrame
        Dataframe containing transect data
    interval_threshold: float
        Along-transect interval threshold

    Notes
    -----
    This function calculates the along-track transect interval length and areas.
    It then 'searches' for possible erroneous values at the end of each line
    and replaces/imputes with alternative lengths/areas.
    """

    # Create dataframe copy
    transect_data_copy = transect_data.copy()

    # Calculate the along-transect interval distance
    # ---- Entire dataframe
    transect_data_copy["interval"] = transect_data_copy["vessel_log_start"].diff(periods=-1).abs()
    # ---- Replace the trailing NaN
    transect_data_copy["interval"] = transect_data_copy["interval"].replace(
        np.nan,
        transect_data_copy["vessel_log_end"].iloc[-1]
        - transect_data_copy["vessel_log_start"].iloc[-1],
    )

    # Replace (likely) erroneous interval lengths associated at the ends of each transect
    # ---- Calculate median interval
    median_interval = np.median(transect_data_copy["interval"])
    # ---- Find indices where interval deviations from the median exceed the difference threshold
    transect_data_copy["interval"] = np.where(
        np.abs(transect_data_copy["interval"] - median_interval) > interval_threshold,
        transect_data_copy["vessel_log_end"] - transect_data_copy["vessel_log_start"],
        transect_data_copy["interval"],
    )

    # Calculate the interval area
    transect_data_copy["interval_area"] = (
        transect_data_copy["interval"] * transect_data_copy["transect_spacing"]
    )

    # Filter out unnecessary columns and return output
    return transect_data_copy.filter(
        regex="^(?=transect_num|latitude|longitude|stratum_|haul_|interval_area|NASC_).*"
    )


def save_transect_coordinates(transect_data: pd.DataFrame):

    # Extract transect numbers, coordinates, and strata
    return transect_data[
        [
            "transect_num",
            "stratum_num",
            "stratum_inpfc",
            "longitude",
            "latitude",
            "transect_spacing",
        ]
    ]


def edit_transect_columns(transect_dict: dict, settings_dict: dict):

    # Define the current stratum definition used for the transect data
    new_stratum = settings_dict["stratum_name"]

    # Prepare the transect data to be updated
    # ---- Call transect data
    transect_data = transect_dict["acoustics"]["adult_transect_df"].copy()
    # -------- Set index
    transect_data.set_index(["transect_num", "longitude", "latitude"], inplace=True)
    # ---- Call full coordinate information
    transect_info = transect_dict["coordinates"].copy()
    # -------- Set index
    transect_info.set_index(["transect_num", "longitude", "latitude"], inplace=True)

    # Adjust stratum values
    # ---- Change strata values
    transect_info = transect_info[[new_stratum, "transect_spacing"]]
    # ---- Select analysis variable
    transect_info[f"{settings_dict['variable']}"] = transect_data[f"{settings_dict['variable']}"]
    # ---- Additional density values, if relevant
    if settings_dict["variable"] == "biomass":
        transect_info["biomass_density"] = transect_data["biomass_density"]
    elif settings_dict["variable"] == "abundance":
        transect_info["number_density"] = transect_data["number_density"]

    # Return the output
    return transect_info.reset_index()


def transect_spatial_features(transect_data: pd.DataFrame):
    """
    Calculates spatial features of each transect

    Parameters
    ----------
    transect_data: pd.DataFrame
        Dataframe comprising georeferenced transect data

    Notes
    -----
    This function calculates the bounding rectangle surrounding the latitude/longitude values
    for each transect and stratum, the average spacing between transects, approximate areas
    relative to each transect, and the distance for each transect
    """

    # Get the name of the stratum column
    stratum_col = [col for col in transect_data.columns if "stratum" in col.lower()][0]

    # Calculate the minimum and maximum longitude, mean latitude, and area of each transect
    # ---- Initialize the dataframe
    transect_spatial = transect_data.drop_duplicates([stratum_col, "transect_num"]).set_index(
        "transect_num"
    )
    # ---- Mean latitude
    transect_spatial["latitude_mean"] = transect_data.groupby(["transect_num"])["latitude"].mean()
    # ---- Min longitude
    transect_spatial["longitude_min"] = transect_data.groupby(["transect_num"])["longitude"].min()
    # ---- Max longitude
    transect_spatial["longitude_max"] = transect_data.groupby(["transect_num"])["longitude"].max()
    # ---- Calculate mean distance (nmi)
    transect_spatial["transect_distance"] = transect_spatial.apply(
        lambda row: geopy.distance.distance(
            (row["latitude_mean"], row["longitude_min"]),
            (row["latitude_mean"], row["longitude_max"]),
        ).nm,
        axis=1,
    )
    # ---- Calculate mean spacing
    transect_spatial["mean_spacing"] = transect_data.groupby(["transect_num"])[
        "transect_spacing"
    ].mean()
    # ---- Calculate transect area
    transect_spatial["transect_area"] = (
        transect_spatial["transect_distance"] * transect_spatial["mean_spacing"]
    )

    # Return output
    return transect_spatial.reset_index()


def summarize_transect_strata(transect_summary: pd.DataFrame):
    """
    Calculate the total number of transects and area coverage within each stratum.
    """

    # Get the name of the stratum column
    stratum_col = [col for col in transect_summary.columns if "stratum" in col.lower()][0]

    # Summarize transect spatial information within each stratum
    # ---- Initialize the dataframe
    strata_summary = (
        pd.DataFrame({f"{stratum_col}": np.unique(transect_summary[f"{stratum_col}"])})
    ).set_index(stratum_col)
    # ---- Calculate the number of transects within each stratum
    strata_summary["transect_count"] = transect_summary.groupby([stratum_col])[
        "transect_num"
    ].size()
    # ---- Sum the total transect area within each stratum
    strata_summary["transect_area_total"] = transect_summary.groupby([stratum_col])[
        "transect_area"
    ].sum()

    # Return output
    return strata_summary.reset_index()


def transect_array(index_matrix: np.ndarray, data_series: pd.Series):
    """
    Helper function for indexing transect arrays.
    """
    return [data_series[idx] for idx in index_matrix]


def transect_coordinate_centroid(spatial_grouped: gpd.GeoSeries):
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
    return Point(centroid_point)


def transect_extent(transect_data: pd.DataFrame, projection: str, num_nearest_transects: int):
    """
    Compute the extent of each transect line.
    """

    # Detect the correct longitude and latitude coordinates
    # ---- Latitude
    lat_col = [col for col in transect_data.columns if "lat" in col.lower()][0]
    # ---- Longitude
    lon_col = [col for col in transect_data.columns if "lon" in col.lower()][0]
    # ---- Rename the dataframe
    transect = transect_data.copy().rename(
        columns={f"{lon_col}": "longitude", f"{lat_col}": "latitude"}
    )

    # Convert to GeoDataFrame
    transect_gdf = gpd.GeoDataFrame(
        transect,
        geometry=gpd.points_from_xy(transect["longitude"], transect["latitude"]),
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


def define_western_extent(transect_data: pd.DataFrame, latitude_threshold: float = 51.0):
    """
    Find the western extent of the survey transect data required for extrapolating
    semivariogram ranges during ordinary kriging

    Parameters
    ----------
    transect_data: pd.DataFrame
        A DataFrame that includes the georeferenced coordinates of survey transect intervals
    latitude_threshold: float
        A threshold that is applied to the georeferenced coordinates that further constrains
        any extrapolation that occurs during the kriging analysis.
    """

    # Apply latitude threshold
    transect_thresholded = transect_data[transect_data.latitude < latitude_threshold]

    # Reduce the dataframe to only the necessary columns
    transect_thresholded = transect_thresholded[["transect_num", "x", "y"]]

    # Parse the western-most coordinates of each transect
    western_extent_idx = transect_thresholded.groupby(["transect_num"])[["x"]].idxmin()

    # Use indices from `western_extent_idx` to extract the actual dataframe rows
    transect_western_extent = transect_thresholded.loc[western_extent_idx["x"]]

    # Return output
    return transect_western_extent


def transect_bearing(transect_data: pd.DataFrame):
    """
    Compute the average transect bearing.
    """

    # Calculate the average heading per transect line
    # ---- Initialize the dataframe
    transect_direction = pd.DataFrame({"transect_num": np.unique(transect_data["transect_num"])})
    # ---- Set index
    transect_direction.set_index("transect_num", inplace=True)
    # ---- Calculate the difference between the minimum and maximum longitudes
    transect_direction["longitude_min"] = np.radians(
        transect_data.groupby(["transect_num"])["longitude"].min()
    )
    transect_direction["longitude_max"] = np.radians(
        transect_data.groupby(["transect_num"])["longitude"].max()
    )
    # ---- Calculate the difference between the minimum and maximum latitudes
    transect_direction["latitude_min"] = np.radians(
        transect_data.groupby(["transect_num"])["latitude"].min()
    )
    transect_direction["latitude_max"] = np.radians(
        transect_data.groupby(["transect_num"])["latitude"].max()
    )
    # ---- Compute the change in longitude
    delta_longitude = transect_direction["longitude_max"] - transect_direction["longitude_min"]
    # ---- Compute the change in y-direction
    delta_y = np.sin(delta_longitude) * np.cos(transect_direction["latitude_max"])
    # ---- Compute the change in x-direction
    delta_x = np.cos(transect_direction["latitude_min"]) * np.sin(
        transect_direction["latitude_max"]
    ) - np.sin(transect_direction["latitude_min"]) * np.cos(
        transect_direction["latitude_max"]
    ) * np.cos(
        delta_longitude
    )
    # ---- Compute the bearings
    transect_direction["heading"] = np.degrees(np.arctan2(delta_y, delta_x)) + 360 % 360

    # Convert other columns back to degrees
    transect_direction[["longitude_min", "longitude_max", "latitude_min", "latitude_max"]] = (
        transect_direction[
            ["longitude_min", "longitude_max", "latitude_min", "latitude_max"]
        ].apply(lambda x: np.degrees(x))
    )

    # Return the output
    return transect_direction.reset_index()
