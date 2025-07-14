import abc
import numpy as np
import pandas as pd
from typing import Union, Dict, List, Optional, Any
from functools import reduce
import pytest

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs
from echopop.nwfsc_feat import utils
from typing import Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate

#######

########
from functools import partial
from echopop.spatial.variogram import variogram
variable: str = "biomass_density"
anisotropy: float = dict_kriging_params["anisotropy"]
search_radius: float = dict_kriging_params["search_radius"]
# OR: 
correlation_range: float = dict_best_fit_variogram_params["correlation_range"]
k_min: int = int(dict_kriging_params["kmin"])
k_max: int = int(dict_kriging_params["kmax"])
kriging_mesh: pd.DataFrame = df_mesh.copy()
transect_df: pd.DataFrame = df_nasc_all_ages.copy()
model = ["bessel", "exponential"]
coordinates = ("x", "y")
coordinate_names = coordinates
variogram_parameters = dict_best_fit_variogram_params
# KrigingParameterInputs.create(**{"correlation_range": correlation_range})
# KrigingAnalysis.create(**{})

##########
transect_df = df_nasc_all_ages.copy()
##########
tx = transect_df["transect_num"]
uniq_tx = tx.unique()
nt = len(uniq_tx)

from typing import List, Tuple
###################
def transect_mesh_region_2019(
    region
) -> Tuple[np.number, np.number, List[np.number], List[np.number]]:
    
    # Initialize variables
    transect_start = None
    transect_end = None
    transect_lower_bound = [] # W/S
    transect_upper_bound = [] # E/N

    # Region 1: parallel transects to latitudes from south of SCB to west of Gwaii Haida
    if region == 1:
        # ---- Southern-most transect
        transect_start = 1
        # ---- Northern-most transect
        transect_end = 119
        # ---- Western boundary
        transect_lower_bound = [i + 0.1 for i in range(transect_start, transect_end + 1)]
        # ---- Eastern boundary
        transect_upper_bound = [i + 0.4 for i in range(transect_start, transect_end + 1)] 
    # Region 2: transects parallel to longitudes north of Gwaii Haida
    elif region == 2:
        # ---- Western-most transect
        transect_start = 121
        # ---- Eastern-most transect
        transect_end = 127
        # ---- Southern boundary
        transect_lower_bound = [i + 0.6 for i in range(transect_start, transect_end + 1)]
        # ---- Northern boundary
        transect_upper_bound = [i + 0.9 for i in range(transect_start, transect_end + 1)] 
    # Region 3: parallel transects to latitudes west of Gwaii Haida
    else:
        # ---- Southern-most transect
        transect_start = 129
        # ---- Northern-most transect
        transect_end = 145
        # ---- Western boundary
        transect_lower_bound = [i + 0.1 for i in range(transect_start, transect_end + 1)]
        # ---- Eastern boundary
        transect_upper_bound = [i + 0.4 for i in range(transect_start, transect_end + 1)] 

    return transect_start, transect_end, transect_lower_bound, transect_upper_bound

###################
latitude_resolution = 1.25 / 60 # degrees
transect_mesh_region_function = transect_mesh_region_2019
###################

# Assign the mesh region definitions to the transect data
# ---- Initialize
transect_df["mesh_region"] = None
# ---- Iterate through 

    
mesh_region_info = [transect_mesh_region_2019(region) for region in [1, 2, 3]]

mesh_region_df = pd.concat([
    pd.DataFrame({
        "transect_num": np.arange(start, end + 1),
        "mesh_region": region,
        "transect_lower_bound": lower,
        "transect_upper_bound": upper
    })
    for region, (start, end, lower, upper) in zip([1, 2, 3], mesh_region_info)
])

mesh_region_dict = {
    region: {
        name: value for name, value in zip(["start", "end", "upper", "lower"], 
                                           transect_mesh_region_2019(region))
    }
    for region in [1, 2, 3]
}

mesh_region_df = pd.concat([
    pd.DataFrame({
        "transect_num": np.arange(start, end + 1),
        "mesh_region": region,
        "transect_lower_bound": lower,
        "transect_upper_bound": upper
    })
    for region, (start, end, lower, upper) in zip([1, 2, 3], mesh_region_info)
])

transect_df = transect_df.merge(mesh_region_df, 
                                on="transect_num", 
                                how="inner")
mesh_region_df.set_index(["mesh_region"], inplace=True)
# # Compute the transect extents across each region
# # ---- Mean latitude
# transect_df["latitude_mean"] = transect_df.groupby(["transect_num", "mesh_region"])[
#     "latitude"
# ].transform("mean")
# # ---- Northernmost extent
# transect_df["latitude_north"] = transect_df.groupby(["transect_num", "mesh_region"])[
#     "latitude"
# ].transform("max")
# # ---- Southernmost extent
# transect_df["latitude_south"] = transect_df.groupby(["transect_num", "mesh_region"])[
#     "latitude"
# ].transform("min")
# # ---- Eastern extent
# transect_df["longitude_east"] = transect_df.groupby(["transect_num", "mesh_region"])[
#     "longitude"
# ].transform("max")
# # ---- Westernmost extent
# transect_df["longitude_west"] = transect_df.groupby(["transect_num", "mesh_region"])[
#     "longitude"
# ].transform("min")
# # ---- Index by region
transect_df.set_index("mesh_region", inplace=True)



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
    transect_df.loc[3].loc[lambda x: x["transect_num"].isin([t3_start, t3_end]), "latitude"].min(),
    transect_df.loc[3].loc[lambda x: x["transect_num"].isin([t3_start, t3_end]), "latitude"].max(),
    latitude_resolution,
)


# ---- Get the minimum and maximum coordinates for each region
# ----- Region 1 conditions
def region_1_conditions(x, boundary_column: str, position: str):
    # Calculate the floor 
    x_floor = np.round((x[boundary_column].iloc[0] % 1) * 10)
          
    # Assign output based on conditions
    if   x_floor == 1.: return pd.Series(
        {f"longitude_{position}": x["longitude"].min(), 
         f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmin()]}
    )
    elif x_floor == 4.: return pd.Series(
        {f"longitude_{position}": x["longitude"].max(), 
         f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmax()]}
    )
    elif x_floor == 6.: return pd.Series(
        {f"longitude_{position}": x["longitude"].iloc[x["latitude"].argmin(), 
         f"latitude_{position}": x["latitude"].min()]}
    )
    elif x_floor == 9.: return pd.Series(
        {f"longitude_{position}": x["longitude"].iloc[x["latitude"].argmax(), 
         f"latitude_{position}": x["latitude"].max()]}
    )

region_1_west = transect_df.loc[1].groupby(["transect_num"]).apply(region_1_conditions, 
                                                                   "transect_lower_bound", 
                                                                   "west",
                                                                   include_groups=False)
region_1_east = transect_df.loc[1].groupby(["transect_num"]).apply(region_1_conditions, 
                                                                   "transect_upper_bound", 
                                                                   "east",
                                                                   include_groups=False)

# Generate the new paired interpolated coordinates
# ---- Region 1 west
interpolator_region1_west = interpolate.interp1d(    
    region_1_west["latitude_west"], 
    region_1_west["longitude_west"], 
    kind="linear", bounds_error=False
)(region_1_latitude)
# ---- Region 1 east
interpolator_region1_east = interpolate.interp1d(     
    region_1_east["latitude_east"], 
    region_1_east["longitude_east"],
    kind="linear", bounds_error=False
)(region_1_latitude)

region_1_index = []

delta_longitude_region_1 = latitude_resolution * np.cos(np.radians(region_1_latitude))
delta_longitude_region_2 = latitude_resolution * np.cos(np.radians(transect_df.loc[2, "latitude"].mean()))

for i in range(len(interpolator_region1_west)):
    # -------- Find the mesh indices that are within the survey extent
    idx = np.where(
        (mesh_df["longitude"] >= interpolator_region1_west[i] - delta_longitude_region_1[i])
        & (mesh_df["longitude"] <= interpolator_region1_east[i] + delta_longitude_region_1[i])
        & (mesh_df["latitude"] >= region_1_latitude[i] - latitude_resolution)
        & (mesh_df["latitude"] < region_1_latitude[i] + latitude_resolution)
    )
    # -------- Append the indices
    region_1_index.append(idx[0])

region_1_index_unique = np.unique(np.concatenate(region_1_index))

def region_2_conditions(x, boundary_column: str, position: str):
    # Calculate the floor 
    x_floor = np.round((x[boundary_column].iloc[0] % 1) * 10)
          
    # Assign output based on conditions
    if x_floor == 1.:
        return pd.Series({
            f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmin()],
            f"longitude_{position}": x["longitude"].min()
        })
    elif x_floor == 4.:
        return pd.Series({
            f"latitude_{position}": x["latitude"].iloc[x["longitude"].argmax()],
            f"longitude_{position}": x["longitude"].max()
        })
    elif x_floor == 6.:
        return pd.Series({
            f"latitude_{position}": x["latitude"].min(),
            f"longitude_{position}": x["longitude"].iloc[x["latitude"].argmin()]
        })
    elif x_floor == 9.:
        return pd.Series({
            f"latitude_{position}": x["latitude"].max(),
            f"longitude_{position}": x["longitude"].iloc[x["latitude"].argmax()]
        })


region_2_south = transect_df.loc[2].groupby(["transect_num"]).apply(region_2_conditions, 
                                                                    "transect_lower_bound", 
                                                                    "south",
                                                                    include_groups=False)
region_2_north = transect_df.loc[2].groupby(["transect_num"]).apply(region_2_conditions, 
                                                                    "transect_upper_bound", 
                                                                    "north",
                                                                    include_groups=False)

# Generate the new paired interpolated coordinates
# ---- Region 1 west
interpolator_region2_south = interpolate.interp1d(  
    region_2_south["longitude_south"], 
    region_2_south["latitude_south"], 
    kind="linear", bounds_error=False
)(region_2_longitude)
# ---- Region 1 east
interpolator_region2_north = interpolate.interp1d(         
    region_2_north["longitude_north"],
    region_2_north["latitude_north"], 
    kind="linear", bounds_error=False
)(region_2_longitude)

region_2_index = []
delta_longitude = latitude_resolution * np.cos(np.radians(region_1_latitude))

for i in range(len(region_2_longitude)):
    # -------- Find the mesh indices that are within the survey extent: southern limit
    if ~np.isnan(interpolator_region2_south[i]) or ~np.isnan(interpolator_region2_north[i]):
        # -------- Compute the indices for the northern- and southernmost coordinates
        # -------- North
        lon_n_min = np.argmin(np.abs(region_2_longitude[i] - region_2_north["longitude_north"]))
        # -------- South
        lon_s_min = np.argmin(np.abs(region_2_longitude[i] - region_2_south["longitude_south"]))  
        # -------- Slope
        slope = (
            region_2_north["latitude_north"].iloc[lon_n_min]
            - region_2_south["latitude_south"].iloc[lon_s_min]
        ) / (
            region_2_north["longitude_north"].iloc[lon_n_min]
            - region_2_south["longitude_south"].iloc[lon_s_min]
        )
        # -------- Set a new border threshold
        latitude_slope_i = (
            slope * (region_2_longitude[i] - region_2_south["longitude_south"].iloc[lon_s_min])
            + region_2_south["latitude_south"].iloc[lon_s_min]
        )
        if np.isnan(interpolator_region2_south[i]):
            # -------- Find the mesh indices that are within the survey extent
            idx = np.where(
                    (mesh_df["longitude"] >= region_2_longitude[i] - delta_longitude_region_2)
                    & (mesh_df["longitude"] <= region_2_longitude[i] + delta_longitude_region_2)
                    & (mesh_df["latitude"] >= latitude_slope_i - latitude_resolution)
                    & (mesh_df["latitude"] < interpolator_region2_north[i] + latitude_resolution)
                )
        elif np.isnan(interpolator_region2_north[i]):
                # -------- Find the mesh indices that are within the survey extent
                idx = np.where(
                    (mesh_df["longitude"] >= region_2_longitude[i] - delta_longitude_region_2)
                    & (mesh_df["longitude"] <= region_2_longitude[i] + delta_longitude_region_2)
                    & (mesh_df["latitude"] >= interpolator_region2_south[i] - latitude_resolution)
                    & (mesh_df["latitude"] < latitude_slope_i + latitude_resolution)
                )
        else:
            # -------- Find the mesh indices that are within the survey extent
            idx = np.where(
                (mesh_df["longitude"] >= region_2_longitude[i] - delta_longitude_region_2)
                & (mesh_df["longitude"] <= region_2_longitude[i] + delta_longitude_region_2)
                & (mesh_df["latitude"] >= interpolator_region2_south[i] - latitude_resolution)
                & (mesh_df["latitude"] < interpolator_region2_north[i] + latitude_resolution)
            ) 
        # -------- Append the indices
        region_2_index.append(idx[0])        
        
    # -------- Find the mesh indices that are within the survey extent
    idx = np.where(
        (mesh_df["longitude"] >= interpolator_region1_west[i] - delta_longitude[i])
        & (mesh_df["longitude"] <= interpolator_region1_east[i] + delta_longitude[i])
        & (mesh_df["latitude"] >= region_1_latitude[i] - latitude_resolution)
        & (mesh_df["latitude"] < region_1_latitude[i] + latitude_resolution)
    )
    # -------- Append the indices
    region_1_index.append(idx[0])

region_2_index_unique = np.unique(np.concatenate(region_2_index))

region_3_west = transect_df.loc[3].groupby(["transect_num"]).apply(region_1_conditions, 
                                                                   "transect_lower_bound", 
                                                                   "west",
                                                                   include_groups=False)
region_3_east = transect_df.loc[3].groupby(["transect_num"]).apply(region_1_conditions, 
                                                                   "transect_upper_bound", 
                                                                   "east",
                                                                   include_groups=False)

# Generate the new paired interpolated coordinates
# ---- Region 3 west
interpolator_region3_west = interpolate.interp1d(    
    region_3_west["latitude_west"], 
    region_3_west["longitude_west"], 
    kind="linear", bounds_error=False
)(region_3_latitude)
# ---- Region 3 east
interpolator_region1_east = interpolate.interp1d(     
    region_3_east["latitude_east"], 
    region_3_east["longitude_east"],
    kind="linear", bounds_error=False
)(region_3_latitude)

region_1_index = []

delta_longitude_region_1 = latitude_resolution * np.cos(np.radians(region_1_latitude))
delta_longitude_region_2 = latitude_resolution * np.cos(np.radians(transect_df.loc[2, "latitude"].mean()))

for i in range(len(interpolator_region1_west)):
    # -------- Find the mesh indices that are within the survey extent
    idx = np.where(
        (mesh_df["longitude"] >= interpolator_region1_west[i] - delta_longitude_region_1[i])
        & (mesh_df["longitude"] <= interpolator_region1_east[i] + delta_longitude_region_1[i])
        & (mesh_df["latitude"] >= region_1_latitude[i] - latitude_resolution)
        & (mesh_df["latitude"] < region_1_latitude[i] + latitude_resolution)
    )
    # -------- Append the indices
    region_1_index.append(idx[0])

region_1_index_unique = np.unique(np.concatenate(region_1_index))

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





# ---- Region 2
# -------- Compute the requisite longitudinal resolution
longitude_resolution_deg = latitude_resolution * np.cos(
    np.radians(
        transect_df.loc[
            transect_df["transect_num"] == transect_df.loc["transect_num"].max(),
            "latitude_mean",
        ].mean()
    )
)
# Create interpolation grids
region_1_lat = np.arange(transect_df.loc[1, "latitude"].min(), transect_df.loc[1, "latitude"].max(), latitude_resolution_deg)
region_3_lat = np.arange(transect_df.loc[3, "latitude_south"].min(), transect_df.loc[3, "latitude_north"].max(), latitude_resolution_deg)

lon_res_deg = latitude_resolution_deg * np.cos(np.radians(
    transect_df.loc[transect_df["transect_num"] == transect_df.loc[2, "transect_num"].max(), "latitude_mean"].mean()
))
region_2_lon = np.arange(transect_df.loc[2, "longitude_west"].min(), transect_df.loc[2, "longitude_east"].max(), lon_res_deg)

# Interpolate extents
region_1_extents = interpolate_survey_extent(region_1_lat, transect_df.loc[1], "latitude", "longitude")
region_2_extents = interpolate_survey_extent(region_2_lon, transect_df.loc[2], "longitude", "latitude")
region_3_extents = interpolate_survey_extent(region_3_lat, transect_df.loc[3], "latitude", "longitude")

def crop_region(lat_or_lon, extents, axis):
    idxs = []
    res_deg = latitude_resolution_deg if axis == "lat" else lon_res_deg
    delta = latitude_resolution_deg * np.cos(np.radians(lat_or_lon)) if axis == "lat" else res_deg
    for i in range(len(lat_or_lon)):
        lo, hi = extents[0][i], extents[1][i]
        if np.isnan(lo) or np.isnan(hi):
            idxs.append(np.array([], dtype=int))
            continue
        if axis == "lat":
            idx = np.where(
                (mesh_df["longitude"] >= lo - delta[i])
                & (mesh_df["longitude"] <= hi + delta[i])
                & (mesh_df["latitude"] >= lat_or_lon[i] - latitude_resolution_deg)
                & (mesh_df["latitude"] < lat_or_lon[i] + latitude_resolution_deg)
            )
        else:
            idx = np.where(
                (mesh_df["longitude"] >= lat_or_lon[i] - lon_res_deg)
                & (mesh_df["longitude"] <= lat_or_lon[i] + lon_res_deg)
                & (mesh_df["latitude"] >= lo - latitude_resolution_deg)
                & (mesh_df["latitude"] < hi + latitude_resolution_deg)
            )
        idxs.append(idx[0])
    return idxs

region_1_idx = crop_region(region_1_lat, region_1_extents, "lat")
region_2_idx = crop_region(region_2_lon, region_2_extents, "lon")
region_3_idx = crop_region(region_3_lat, region_3_extents, "lat")

all_indices = np.unique(np.concatenate([
    np.concatenate(region_1_idx),
    np.concatenate(region_2_idx),
    np.concatenate(region_3_idx)
]))

transect_mesh_regions = transect_df.reset_index().loc[:, "mesh_region":"latitude"]
return mesh_df.loc[all_indices], transect_mesh_regions



# Assign mesh_region using fixed region definitions
transect_df["mesh_region"] = None
for region in [1, 2, 3]:
    transect_start, transect_end, _, _ = transect_mesh_region_2019(region)
    mask = transect_data["transect_num"].between(transect_start, transect_end)
    transect_data.loc[mask, "mesh_region"] = region
transect_bounds = {
    region: transect_mesh_region_2019(region)[:2] for region in [1, 2, 3]
}
conditions = [
    transect_data["transect_num"].between(start, end)
    for region, (start, end) in transect_bounds.items()
]
choices = list(transect_bounds.keys())
transect_data["mesh_region"] = np.select(conditions, choices, default=np.nan)

transect_data = transect_data.dropna(subset=["mesh_region"]).copy()
transect_data["mesh_region"] = transect_data["mesh_region"].astype(int)

transect_data = transect_data.dropna(subset=["mesh_region"]).copy()
transect_data["mesh_region"] = transect_data["mesh_region"].astype(int)

# Compute spatial extent summaries
transect_data["latitude_mean"] = transect_data.groupby(["transect_num", "mesh_region"])["latitude"].transform("mean")
transect_data["latitude_north"] = transect_data.groupby(["transect_num", "mesh_region"])["latitude"].transform("max")
transect_data["latitude_south"] = transect_data.groupby(["transect_num", "mesh_region"])["latitude"].transform("min")
transect_data["longitude_east"] = transect_data.groupby(["transect_num", "mesh_region"])["longitude"].transform("max")
transect_data["longitude_west"] = transect_data.groupby(["transect_num", "mesh_region"])["longitude"].transform("min")
transect_data.set_index("mesh_region", inplace=True)

region_1_lat = np.arange(transect_data.loc[1, "latitude"].min(), transect_data.loc[1, "latitude"].max(), latitude_resolution_deg)
region_3_lat = np.arange(transect_data.loc[3, "latitude_south"].min(), transect_data.loc[3, "latitude_north"].max(), latitude_resolution_deg)

lon_res_deg = latitude_resolution_deg * np.cos(np.radians(
    transect_data.loc[transect_data["transect_num"] == transect_data.loc[2, "transect_num"].max(), "latitude_mean"].mean()
))
region_2_lon = np.arange(transect_data.loc[2, "longitude_west"].min(), transect_data.loc[2, "longitude_east"].max(), lon_res_deg)

region_1_extents = interpolate_survey_extent(region_1_lat, transect_data.loc[1], "latitude", "longitude")
region_2_extents = interpolate_survey_extent(region_2_lon, transect_data.loc[2], "longitude", "latitude")
region_3_extents = interpolate_survey_extent(region_3_lat, transect_data.loc[3], "latitude", "longitude")

def crop_region(lat_or_lon, extents, axis, data_slice):
    idxs = []
    res_deg = latitude_resolution_deg if axis == "lat" else lon_res_deg
    delta = latitude_resolution_deg * np.cos(np.radians(lat_or_lon)) if axis == "lat" else res_deg
    for i in range(len(lat_or_lon)):
        lo, hi = extents[0][i], extents[1][i]
        if np.isnan(lo) or np.isnan(hi):
            idxs.append(np.array([], dtype=int))
            continue
        if axis == "lat":
            idx = np.where(
                (mesh_data["longitude"] >= lo - delta[i])
                & (mesh_data["longitude"] <= hi + delta[i])
                & (mesh_data["latitude"] >= lat_or_lon[i] - latitude_resolution_deg)
                & (mesh_data["latitude"] < lat_or_lon[i] + latitude_resolution_deg)
            )
        else:
            idx = np.where(
                (mesh_data["longitude"] >= lat_or_lon[i] - lon_res_deg)
                & (mesh_data["longitude"] <= lat_or_lon[i] + lon_res_deg)
                & (mesh_data["latitude"] >= lo - latitude_resolution_deg)
                    & (mesh_data["latitude"] < hi + latitude_resolution_deg)
                )
            idxs.append(idx[0])
        return idxs

    region_1_idx = crop_region(region_1_lat, region_1_extents, "lat", transect_data.loc[1])
    region_2_idx = crop_region(region_2_lon, region_2_extents, "lon", transect_data.loc[2])
    region_3_idx = crop_region(region_3_lat, region_3_extents, "lat", transect_data.loc[3])

    all_indices = np.unique(np.concatenate([
        np.concatenate(region_1_idx),
        np.concatenate(region_2_idx),
        np.concatenate(region_3_idx)
    ]))

    transect_mesh_regions = transect_data.reset_index().loc[:, "mesh_region":"latitude"]
    # return mesh_data.loc[all_indices], transect_mesh_regions



####################################################################################################
from echopop.survey import Survey
from echopop.utils.validate_dict import KrigingParameterInputs, KrigingAnalysis
from echopop.spatial.transect import correct_transect_intervals
from echopop.acoustics import aggregate_sigma_bs, nasc_to_biomass
# from echopop.biology import (
#     # age1_metric_proportions,
#     # distribute_length_age,
#     # filter_species,
#     # fit_length_weight_relationship,
#     # fit_length_weights,
#     # impute_kriged_values,
#     # # number_proportions,
#     # # partition_transect_age,
#     # quantize_number_counts,
#     # quantize_weights,
#     # reallocate_kriged_age1,
#     # weight_proportions,
# )
from echopop.spatial.krige import kriging
from echopop.spatial.mesh import crop_mesh, mesh_to_transects, stratify_mesh
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import (
    edit_transect_columns,
    save_transect_coordinates,
    summarize_transect_strata,
    transect_spatial_features,
)
from echopop.analysis import (
    acoustics_to_biology,
    apportion_kriged_values,
    krige,
    process_transect_data,
    stratified_summary,
    variogram_analysis,
)
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import edit_transect_columns
from echopop.utils import load as el, load_nasc as eln, message as em
from echopop.utils.load import dataset_integrity
from echopop.spatial.variogram import (
    empirical_variogram,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
)
from echopop.spatial.mesh import griddify_lag_distances
from echopop.spatial.transect import define_western_extent
from echopop.spatial.krige import kriging
from echopop.statistics import stratified_transect_statistic
from echopop.utils.validate_dict import (
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
)
from echopop.spatial.krige import griddify_lag_distances, define_western_extent, adaptive_search_radius, count_within_radius, kriging_interpolation, kriging_lambda, kriging_matrix
from echopop.spatial.krige import kriging_lambda, kriging_matrix
from echopop.nwfsc_feat.spatial import kriging_lambda
survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data(ingest_exports="echoview")
survey.load_survey_data()
survey.transect_analysis(exclude_age1=False)
survey.fit_variogram()
survey.kriging_analysis(extrapolate=True)

self = survey
input_dict, analysis_dict, settings_dict = self.input, self.analysis, self.analysis["settings"]["kriging"]
transect_data, mesh_data, settings_dict = analysis_dict["kriging"]["transect_df"], analysis_dict["kriging"]["mesh_df"], settings_dict

survey.analysis["settings"]["variogram"][""]
self = survey
input_dict, analysis_dict, configuration_dict, settings_dict = self.input, self.analysis["transect"], self.config, self.analysis["settings"]
# Extract the necessary correct strata mean sigma_bs
sigma_bs_strata = analysis_dict["acoustics"]["sigma_bs"]["strata_mean_df"]

# Pull out the length-weight conversion for each stratum
length_weight_strata = analysis_dict["biology"]["weight"]["weight_stratum_df"]

# Get the name of the stratum column
stratum_col = settings_dict["transect"]["stratum_name"]

# Get group-specific columns
age_group_cols = settings_dict["transect"]["age_group_columns"]

# Extract the correct strata dataframe
# ---- Define `strata_df` if KS
if settings_dict["transect"]["stratum"] == "ks":
    strata_df = input_dict["spatial"]["strata_df"].copy()
# Define `inpfc_strata_df` if INPFC
elif settings_dict["transect"]["stratum"] == "inpfc":
    strata_df = input_dict["spatial"]["inpfc_strata_df"].copy()

# Get group-specific column names and create conversion key
name_conversion_key = {age_group_cols["haul_id"]: "haul_num", age_group_cols["nasc_id"]: "nasc"}
# ---- Update if the stratum is not equal to INPFC
if settings_dict["transect"]["stratum"] != "inpfc":
    name_conversion_key.update({age_group_cols["stratum_id"]: stratum_col})

# Rename columns
# ---- Extract NASC data
nasc_data = input_dict["acoustics"]["nasc_df"].copy()
# ---- Change names
nasc_data.rename(columns=name_conversion_key, inplace=True)

# Correct the acoustic survey transect intervals
nasc_interval_df = correct_transect_intervals(nasc_data)

distributions_dict, proportions_dict, TS_L_parameters, settings_dict = (
    input_dict["biology"]["distributions"],
    analysis_dict["biology"]["proportions"],
    configuration_dict["TS_length_regression_parameters"]["pacific_hake"],
    settings_dict
)