import geopandas as gpd
import pandas as pd
import numpy as np
from geopy.distance import distance
from ..spatial.projection import utm_string_generator
import shapely.geometry

def apply_inpfc_definitions(acoustic_data: dict, biology_data: dict, spatial_config: dict):

    # Extract the INPFC definitions
    inpfc_definitions = spatial_config["inpfc"]

    # Create latitude bins
    latitude_bins = np.concatenate([[-90.0], inpfc_definitions["latitude_max"], [90.0]])
    # ---- Append 1 more stratum layer
    bin_names = np.concatenate([inpfc_definitions["stratum_names"],
                                [np.max(inpfc_definitions["stratum_names"]) + 1]])
    
    # Create spatial key
    spatial_config["spatial_key"] = pd.DataFrame({
        "latitude_limit": inpfc_definitions["latitude_max"],
    })
    # ---- Cut
    spatial_config["spatial_key"]["stratum"] = (
        pd.cut(inpfc_definitions["latitude_max"],
               latitude_bins,
               right = True,
               labels = bin_names)
    )

    # Get the `prc_nasc_df` values, if they exist, and apply stratification information
    if not acoustic_data["prc_nasc_df"].empty:
        # ---- Bin the latitude data
        acoustic_data["prc_nasc_df"]["stratum"] = pd.cut(
            acoustic_data["prc_nasc_df"]["latitude"],
            latitude_bins,
            right = True,
            labels = bin_names,
        )

    # Get the `trawl_info_df` values, if they exist, and apply stratification information
    if not biology_data["trawl_info_df"].empty:
        # ---- Bin the latitude data
        biology_data["trawl_info_df"]["stratum"] = pd.cut(
            biology_data["trawl_info_df"]["latitude"],
            latitude_bins,
            right = True,
            labels = bin_names,
        )

def define_boundary_box(boundary_dict: dict, projection: str):
    
    # Get x-coordinates
    if "longitude" in boundary_dict.keys():
        x = np.array(boundary_dict["longitude"])
    else:
        x = np.array(boundary_dict["northings"])

    # Get y-coordinates
    if "latitude" in boundary_dict.keys():
        y = np.array(boundary_dict["latitude"])
    else:
        y = np.array(boundary_dict["eastings"])

    # Create a boundary DataFrame
    bound_df = pd.DataFrame({
        "x": np.array([x.min(), x.max(), x.max(), x.min(), x.min()]),
        "y":np.array([y.min(), y.max(), y.max(), y.min(), y.min()]),
    })

    # Convert to a GeoDataFrame and return the GeoDataFrame
    return gpd.GeoDataFrame(
        data=bound_df,
        geometry=gpd.points_from_xy(bound_df["x"], bound_df["y"]),
        crs=projection,
    )

def apply_griddify_definitions(acoustic_data: dict, biology_data: dict, spatial_config: dict):

    # Extract the griddification definitions
    griddify_definitions = spatial_config["griddify"]

    # Get the projection definition
    projection = spatial_config["projection"]

    # Compute the boundary box GeoDataFrame
    boundary_box = define_boundary_box(griddify_definitions["bounds"], projection)

    # Convert the coordinates, if needed
    if not set(["northings", "eastings"]).intersection(set(griddify_definitions["bounds"].keys())):
        # ---- Compute the equivalent UTM string
        utm_num = int(utm_string_generator(np.median(boundary_box.loc[0:3, "x"]), 
                                           np.median(boundary_box.loc[0:3, "y"])))
        # ---- Compute the boundary box GeoDataFrame with the new projection
        boundary_box = boundary_box.to_crs(utm_num)
        # ---- Create a new projection for later
        projection_new = f"epsg:{utm_num}"
    else:
        projection_new = projection

    # Define the step sizes
    # ---- Define x step size
    x_step = distance(nautical=griddify_definitions["grid_resolution"]["x_distance"]).meters
    # ---- Define y step size
    y_step = distance(nautical=griddify_definitions["grid_resolution"]["y_distance"]).meters

    # Get the boundary tuple
    xmin, ymin, xmax, ymax = boundary_box.total_bounds

    # Generate the cells
    grid_cells = []
    # ---- Iterate through
    for y0 in np.arange(ymin, ymax+y_step, y_step):
        for x0 in np.arange(xmin, xmax+x_step, x_step):
            x1 = x0-x_step
            y1 = y0+y_step
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))

    # Convert to a GeoDataFrame
    cells_gdf = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=projection_new)

    # Get the centroids
    cells_gdf["cell_centroid"] = cells_gdf["geometry"].centroid

    # Get the `prc_nasc_df` values, if they exist, and apply stratification information
    if not acoustic_data["prc_nasc_df"].empty:

        #
        prc_nasc_df = acoustic_data["prc_nasc_df"]

        # to GDF
        prc_nasc_gdf = gpd.GeoDataFrame(
            data=prc_nasc_df,
            geometry=gpd.points_from_xy(prc_nasc_df["longitude"], prc_nasc_df["latitude"]),
            crs=projection,
        )
        # to UTM
        prc_nasc_new = prc_nasc_gdf.to_crs(projection_new)

        prc_nasc_new["x"] = prc_nasc_new["geometry"].x
        prc_nasc_new["y"] = prc_nasc_new["geometry"].y

        # ---- Bin the latitude data
        prc_nasc_new["stratum_x"] = pd.cut(
            prc_nasc_new["x"],
            np.arange(xmin, xmax+x_step, x_step),
            right = True,
            labels = range(len(np.arange(xmin, xmax+x_step, x_step)) - 1),
        ).astype(int) + 1

        prc_nasc_new["stratum_y"] = pd.cut(
            prc_nasc_new["y"],
            np.arange(ymin, ymax+y_step, y_step),
            right = True,
            labels = range(len(np.arange(ymin, ymax+y_step, y_step)) - 1),
        ).astype(int) + 1

        #
        acoustic_data["prc_nasc_df"]["stratum"] = (
            prc_nasc_new["stratum_x"].astype(str) + "-" + prc_nasc_new["stratum_y"].astype(str)
        )

    if not biology_data["trawl_info_df"].empty:

        #
        trawl_info_df = biology_data["trawl_info_df"]

        # to GDF
        trawl_info_gdf = gpd.GeoDataFrame(
            data=trawl_info_df,
            geometry=gpd.points_from_xy(trawl_info_df["longitude"], trawl_info_df["latitude"]),
            crs=projection,
        )
        # to UTM
        trawl_info_new = trawl_info_gdf.to_crs(projection_new)

        trawl_info_new["x"] = trawl_info_new["geometry"].x
        trawl_info_new["y"] = trawl_info_new["geometry"].y

        # ---- Bin the latitude data
        trawl_info_new["stratum_x"] = pd.cut(
            trawl_info_new["x"],
            np.arange(xmin, xmax+x_step, x_step),
            right = True,
            labels = range(len(np.arange(xmin, xmax+x_step, x_step)) - 1),
        ).astype(int) + 1

        trawl_info_new["stratum_y"] = pd.cut(
            trawl_info_new["y"],
            np.arange(ymin, ymax+y_step, y_step),
            right = True,
            labels = range(len(np.arange(ymin, ymax+y_step, y_step)) - 1),
        ).astype(int) + 1

        #
        biology_data["trawl_info_df"]["stratum"] = (
            trawl_info_new["stratum_x"].astype(str) + "-" + trawl_info_new["stratum_y"].astype(str)
        )
