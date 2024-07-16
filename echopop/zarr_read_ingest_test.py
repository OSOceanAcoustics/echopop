import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union, Tuple
from pathlib import Path
import copy
import yaml
import glob
from datetime import datetime
import geopandas as gpd

####################################################################################################
# * Functionality for a) loading YAML configuration file, b) search defined directory for 
# * input files, c) ingest *.zarr/*.csv
# TODO: Incorporate complete YAML file validator
# TODO: Documentation
def live_configuration(live_init_config_path: Union[str, Path], 
                       live_file_config_path: Union[str, Path]):
    
    # Validate file existence
    # ---- str-to-Path conversion, if necessary
    live_init_config_path = Path(live_init_config_path)
    live_file_config_path = Path(live_file_config_path)
    # ---- Create list of both config paths
    config_files = [live_init_config_path, live_file_config_path]
    # ---- List of file existence checks
    config_existence = [live_init_config_path.exists(), live_file_config_path.exists()]
    # ---- Error evaluation and print message (if applicable)
    if not all(config_existence):
        missing_config = [
            files for files, exists in zip(config_files, config_existence) if not exists
        ]
        raise FileNotFoundError(f"The following configuration files do not exist: {missing_config}")

    # Read the YAML configuration/recipe file to parameterize the `LiveSurvey` class
    # ---- Initialization settings
    init_config = yaml.safe_load(Path(live_init_config_path).read_text())
    # ---- Filepath/directory settings
    file_config = yaml.safe_load(Path(live_file_config_path).read_text())
    
    # Check for intersecting/duplicative configuration keys
    # ---- Compare sets of keys from each dictionary
    config_intersect = set(init_config.keys()).intersection(set(file_config.keys()))
    # ---- Raise error if needed
    if config_intersect:
        raise ValueError(
            f"The initialization and file configuration files comprise the following intersecting "
            f"keys: {' ,'.join(config_intersect)}. Key names must be unique for each configuration "
            f"file."
        )
    
    # Combine both into a dictionary output that can be added to the `LiveSurvey` class object
    return {**init_config, **file_config}
####################################################################################################  
# TEST: YAML FILE CONFIGURATION
# ---- Define filepaths
live_init_config_path = "C:/Users/15052/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
live_file_config_path = "C:/Users/15052/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"
# ---- Run function: `live_configuration`
file_configuration = live_configuration(live_init_config_path, live_file_config_path)
file_configuration
####################################################################################################
# * Accessory function for tuning the acoustic transmit frequency units/scaling
# TODO: Documentation
def configure_transmit_frequency(frequency_values: pd.Series,
                                 transmit_settings: dict, 
                                 current_units: str):
    
    # Extract transmit frequency units defined in configuration file
    configuration_units = transmit_settings["units"]
    
    # Transform the units, if necessary
    # ---- Hz to kHz
    if current_units == "Hz" and configuration_units == "kHz":
        return frequency_values * 1e-3
    # ---- kHz to Hz
    elif current_units == "kHz" and configuration_units == "Hz":
        return frequency_values * 1e3
    # ---- No change
    else:
        return frequency_values
####################################################################################################
# * Define `LIVE_INPUT_FILE_CONFIG_MAP` configuration mapping (this will be in an equivalent 
# * `core.py`)
# TODO: Update structure with additional information (as needed)
# TODO: Documentation
LIVE_INPUT_FILE_CONFIG_MAP = {
    "acoustics": {
        "xarray_coordinates": {
            "distance": float,
            "depth": float,
        },
        "xarray_variables": {
            "NASC": float,
            "frequency_nominal": float, 
            "latitude": float,
            "longitude": float,
            "ping_time": "datetime64[ns]",
        }
    }
}
####################################################################################################
# * Functionality for reading in processed acoustic data
# TODO: Expand data validator and limit cases to '*.zarr' (for now)
# TODO: Refactor "extra" components such as the validation steps, xarray-to-dataframe piping, etc.
# TODO: Documentation
def load_acoustic_data(file_configuration: dict) -> Tuple[pd.DataFrame, xr.Dataset]:
    # Get acoustic directory and initialization settings
    # ---- Files
    acoustic_file_settings = file_configuration["input_directories"]["acoustic"]
    # ---- General settings
    acoustic_analysis_settings = file_configuration["acoustics"]
    
    # Create full filepath
    acoustic_directory_path = (
        Path(file_configuration["data_root_dir"]) / acoustic_file_settings["directory"]
    )
    
    # Validate filepath, columns, datatypes
    # ---- Directory check
    directory_existence = acoustic_directory_path.exists()
    # ---- Error evaluation (if applicable)
    if not directory_existence:
        raise FileNotFoundError(
            f"The acoustic data directory [{acoustic_directory_path}] does not exist."
        )
    # ---- Get the defined file extension
    file_extension = acoustic_file_settings["extension"]
    # ---- In the case of a *.zarr file
    if file_extension == "zarr":
        # ---- Create Path.glob generator object
        file_path_obj = acoustic_directory_path.glob(f"*{'.'+file_extension}")
        # ---- Find all zarr files
        zarr_files = list(file_path_obj)
        # ---- Ensure files exist or raise error otherwise
        if len(zarr_files) < 1:
            raise FileNotFoundError(
                f"No `*.zarr` files found in [{acoustic_directory_path}]!"
            )
        # ---- Get defined columns and datatypes from `LIVE_INPUT_FILE_CONFIG_MAP`
        acoustics_config_map = LIVE_INPUT_FILE_CONFIG_MAP["acoustics"]
        # ---- Create list of coordinate data variables
        specified_vars = list(acoustics_config_map["xarray_variables"].keys())
        # ---- Create set of coordinate variables
        specified_coords = list(acoustics_config_map["xarray_coordinates"].keys())      
        # ---- Concatenate into a full configuration map
        full_config_map = {**acoustics_config_map["xarray_coordinates"],
                           **acoustics_config_map["xarray_variables"]}          
        # ! [REQUIRES DASK] ---- Read in all listed files 
        # TODO: The sliding/overlapping windows makes this annoying -- in theory, only a single new zarr file will be ingested
        # TODO: So this needs to be replaced w/ `open_dataset` instead
        zarr_data_ds = xr.open_mfdataset(zarr_files, 
                                         engine="zarr",
                                         chunks="auto",
                                         data_vars=specified_vars,
                                         coords=specified_coords)
        # ---- Extract coordinate metadata
        coordinate_metadata = zarr_data_ds[["longitude", "latitude"]]
        # ---- Convert to a DataFrame
        zarr_data_df = zarr_data_ds.to_dataframe().reset_index()
        # ---- Check for any missing columns
        missing_columns = (
            [key for key in full_config_map.keys() if key not in zarr_data_df.columns]
        )
        # ---- Raise Error, if needed
        if missing_columns: 
            raise ValueError(
                f"The following columns are missing from at least one *.{file_extension} file in "
                f"[{acoustic_directory_path}]: {', '.join(missing_columns)}!"    
            )
        # ---- Select defined columns
        zarr_data_df_filtered = zarr_data_df[full_config_map.keys()]
        # ---- Validate data types
        zarr_data_df_filtered = (
            zarr_data_df_filtered
            .apply(lambda col: col.astype(full_config_map[col.name]) 
                   if col.name in full_config_map else col)
        )
        
    # Extract defined acoustic frequency
    # ---- From the configuration 
    transmit_settings = acoustic_analysis_settings["transmit"]
    # ---- Transform `frequency_nominal`, if necessary
    zarr_data_df_filtered["frequency_nominal"] = (
        configure_transmit_frequency(zarr_data_df_filtered["frequency_nominal"],
                                     transmit_settings,
                                     zarr_data_ds["frequency_nominal"].units)
    )
    # ---- Filter out any unused frequency coordinates
    zarr_data_df_output = (
        zarr_data_df_filtered
        [zarr_data_df_filtered["frequency_nominal"] == transmit_settings["frequency"]]
    )
    
    # Remaining adjustments to the acoustic data prior to being passed to the `LiveSurvey` object
    # ---- Replace NASC `NaN` values with `0.0`
    zarr_data_df_output.loc[:, "NASC"] = zarr_data_df_output.loc[:, "NASC"].fillna(0.0)
    # ---- Drop frequency column and return the output
    return zarr_data_df_output.drop(columns = ["frequency_nominal"]), coordinate_metadata
####################################################################################################  
# TEST: ACOUSTIC ZARR FILE INGESTION CONFIGURATION
# NOTE: 
# ---- Run function: `load_validated_acoustic_data` using previously defined `file_configuration`
acoustic_data, coordinate_metadata = load_acoustic_data(file_configuration)
acoustic_data
coordinate_metadata
####################################################################################################
def load_spatial_data(file_configuration: dict,
                      acoustic_data: pd.DataFrame,
                      coordinate_metadata: xr.Dataset):
    
    # Extract spatial strata *only* if spatial information from the configuration settings
    # ---- Extract the projection
    projection = file_configuration["geospatial"]["projection"]
    # ---- Extract the biology-acoustics linking method options
    acoustics_biology_link = file_configuration["geospatial"]["link_biology_acoustics"]

    # Convert the DataFrame to a GeoDataFrame
    acoustic_data_gdf = gpd.GeoDataFrame(
        data=acoustic_data,
        geometry=gpd.points_from_xy(acoustic_data["longitude"], acoustic_data["latitude"]),
        crs=projection
    )

    # Validate the spatial biology-acoustics linking method
    # ---- Get the biology-acoustics linking method
    link_method = next(key for key, value in acoustics_biology_link.items() if value)
    # ---- Flag Error if unexpected method
    if link_method not in ["global", "closest_haul", "INPFC", "weighted_haul"]:
        raise ValueError(
            f"Unexpected biology-acoustic linking parameter ([{link_method}]). Valid options "
            f"include: 'global', 'closest_haul', 'weighted_haul', and 'INPFC'."
        )
    
    # Create INPFC stratum dataframe
    # ---- Extract 
        
    # Validate projection information
    # ---- Create a dummy GeoDataFrame to extract CRS information
    # geo_crs = gpd.GeoDataFrame(geometry=[], crs=projection)
    # ---- Extract coordinate limits from the acoustic data
    # lat_min = coordinate_metadata.attrs['geospatial_lat_min']
    # lat_max = coordinate_metadata.attrs['geospatial_lat_max']
    # lon_min = coordinate_metadata.attrs['geospatial_lon_min']
    # lon_max = coordinate_metadata.attrs['geospatial_lon_max']
    # # ---- Create boundary box string
    # boundary_box_str = (
    #     f"POLYGON(({lon_min} {lat_min}, {lon_max} {lat_min}, {lon_max} {lat_max}, "
    #     f"{lon_min} {lat_max}, {lon_min} {lat_min}))"
    # )
    
    # data_gdf = gpd.GeoDataFrame(acoustic_data, geometry=gpd.points_from_xy(acoustic_data["longitude"], acoustic_data["latitude"]),crs=f"epsg:{utm_string_generator(lon_min, lat_min)}")
    # gpd.GeoDataFrame(acoustic_data, geometry=gpd.points_from_xy(acoustic_data["longitude"], acoustic_data["latitude"]),crs=f"epsg:4326").to_crs("epsg:32610")
    
    # from pyproj import CRS
    # from pyproj.aoi import AreaOfInterest
    # from pyproj.database import query_utm_crs_info
    
    # utm_crs_list = query_utm_crs_info(
    #     datum_name="WGS 84",
    #     area_of_interest=AreaOfInterest(
    #         west_lon_degree=lon_min,
    #         south_lat_degree=lat_min,
    #         east_lon_degree=-lon_max,
    #         north_lat_degree=lat_max,
    #     ),
    # )
    # CRS.from_epsg(utm_crs_list[0].code).to_epsg("+proj=latlon")
    
####################################################################################################
def live_data(file_configuration: dict): 
    
    # Extract the file directories (or from the configuration) containing acoustic, biological, and 
    # spatial definitions/data/parameters
    # ---- Acoustic data
    acoustic_data = load_validated_acoustic_data(file_configuration)
    # ---- Biological data 
    # ---- Spatial data
    


####################################################################################################
# * Define `LIVE_DATA_STRUCTURE` configuration mapping (this will be in an equivalent `core.py`)
# TODO: Update structure with additional information (as needed)
# TODO: Documentation
LIVE_DATA_STRUCTURE = {
    "meta": {
        "provenance": dict(),
        "date": list(),
    },
    "input": {
        "acoustics": {
            "nasc_df": pd.DataFrame(),
        },
        "biology": {
            "catch_df": pd.DataFrame(),
            "distributions": {
                "length_bins_df": pd.DataFrame(),
            },
            "length_df": pd.DataFrame(),
            "specimen_df": pd.DataFrame(),
        },
    },
    "results": {
        "acoustics": dict(),
        "biology": dict(),
        "stratified": dict(),        
    },
}
####################################################################################################
# * Define `LiveSurvey` class structure
# TODO: Incorporate validators
# TODO: Scope out full structure including accessors, attributes, and methods
# TODO: Configure input arguments (for initialization)
# TODO: Documentation
class LiveSurvey:
    """
    A real-time processing version of the `echopop` base `Survey` class that ingests biological, 
    acoustic, and event meta data to provide population estimates when generated.
    """

    def __init__(
        self,
        live_init_config_path: Union[str, Path], 
        live_file_config_path: Union[str, Path],
    ):
        # Initialize `meta` attribute
        self.meta = copy.deepcopy(LIVE_DATA_STRUCTURE["meta"])

        # Loading the configuration settings and definitions that are used for defining the 
        # configuration settings
        self.config = live_configuration(live_file_config_path, live_file_config_path)

        # Loading the datasets defined in the configuration files
        self.input = el.load_survey_data(self.config)

        # Initialize the `results` data attribute
        self.results = copy.deepcopy(LIVE_DATA_STRUCTURE["results"])

current_units = zarr_data_ds["frequency_nominal"].units
acoustic_analysis_settings["transmit"]
file_configuration

specimen_df = pd.DataFrame(
    {
        "haul_num": np.repeat([1,2,3], 4),
        "station": "specimen",
        "sex": np.tile(["male", "female"], 6),
        "length": np.array([11, 11, 11, 18, 21, 23, 13, 11, 19, 25, 18, 9]), 
        "weight": np.array([11, 14, 16, 18, 21, 23, 13, 11, 19, 25, 18, 9]) / 3.5,
    },
)

length_df = pd.DataFrame(
    {
        "haul_num": np.repeat([1,2,3], 4),
        "station": "length",
        "sex": np.tile(["male", "female"], 6),
        "length": np.array([16, 15, 19, 14, 9, 10, 18, 15, 16, 22, 17, 11]), 
        "length_count": np.array([103, 123, 257, 106, 52, 329, 131, 72, 101, 212, 93, 81]),
    },
)

catch_df = pd.DataFrame(
    {
        "haul_num": np.array([1, 2, 3]),
        "weight": np.array([503.12, 684.32, 978.54])
    }
)

TS_SLOPE = 20.0
TS_INTERCEPT = -68.0

####
# CONCATENATE FILE SOURCES
specimen_reframed = specimen_df.groupby(["haul_num", "station", "sex", "length"])["length"].value_counts().to_frame("length_count").reset_index()
specimen_reframed
# MELD
all_lengths = pd.concat([length_df, specimen_reframed])
# COMBINE 
comb_lengths = all_lengths.groupby(["haul_num", "sex", "length"])["length_count"].sum().to_frame("length_count").reset_index()


# CONVERT TO TS
comb_lengths["ts"] = TS_SLOPE * np.log10(comb_lengths["length"]) + TS_INTERCEPT
# TO SIGMA_BS
comb_lengths["sigma_bs"] = 10 ** (comb_lengths["ts"] / 10)
# WEIGHTED MEAN SIGMA_BS
sigma_mean = np.average(comb_lengths["sigma_bs"], weights=comb_lengths["length_count"])

### 
# INTEGRATE NASC
path2file = "C:/Users/15052/Downloads/win_1720457505_1720460000_NASC.zarr"

Path(path2file).exists()
xds = xr.open_dataset(path2file, engine="zarr")
xds
xdf = xds.to_dataframe().reset_index()
xdf["NASC"] = xdf["NASC"].fillna(0.0)
# convert frequency
xdf["frequency_nominal"] = (xdf["frequency_nominal"] * 1e-3).astype(int)
# filter
xdf_38 = xdf[xdf["frequency_nominal"] == nasc_frequency]

xdf_38.plot.scatter(x="distance", y="depth", c="NASC")
plt.show()

xdf_int = xdf_38.groupby(["distance", "longitude", "latitude"])["NASC"].sum().reset_index()

plt.scatter(xdf_int["longitude"], xdf_int["latitude"], c=xdf_int["NASC"])
plt.plot(xdf_int["longitude"], xdf_int["latitude"])
plt.show()

# CONVERT TO NUMBER DENSITY
xdf_int["number_density"] = xdf_int["NASC"] / (4.0 * np.pi * sigma_mean)


###################
from geopy.distance import distance
from shapely.geometry import Polygon, Point, box
import geopandas as gpd
from shapely.ops import unary_union
import pyproj


grid_settings = file_configuration["geospatial"]["griddify"]
grid = []
lat_step = distance(nautical=grid_settings["grid_resolution"]["x"]).meters
lon_step = distance(nautical=grid_settings["grid_resolution"]["y"]).meters
lat_min = grid_settings["bounds"]["latitude"][0]
lat_max = grid_settings["bounds"]["latitude"][1]
lon_min = grid_settings["bounds"]["longitude"][0]
lon_max = grid_settings["bounds"]["longitude"][1]

utm_str = utm_string_generator((lon_max + lon_min)/2, (lat_max + lat_min)/2)
utm_proj = pyproj.Proj(f"epsg:{utm_str}")
x_min, y_min = utm_proj(lon_min, lat_min)
x_max, y_max = utm_proj(lon_max, lat_max)

num_lon_steps = int((x_max - x_min) / lon_step)
num_lat_steps = int((y_max - y_min) / lat_step)

lon1 = np.linspace(x_min, x_max - lon_step, num_lon_steps)
lat1 = np.linspace(y_min, y_max - lat_step, num_lat_steps)
lon2 = lon1 + lon_step
lat2 = lat1 + lat_step

# Convert UTM coordinates back to degrees
lon_min_grid, lat_min_grid = np.meshgrid(lon1, lat1)
lon_max_grid, lat_max_grid = np.meshgrid(lon2, lat2)

# Convert UTM coordinates back to degrees with adjusted resolution
lon1_deg, lat1_deg = utm_proj(lon_min_grid.ravel(), lat_min_grid.ravel(), inverse=True)
lon2_deg, lat2_deg = utm_proj(lon_max_grid.ravel(), lat_max_grid.ravel(), inverse=True)

polygons = [box(lon1, lat1, lon2, lat2) for lon1, lat1, lon2, lat2 in zip(lon1_deg, lat1_deg, lon2_deg, lat2_deg)]
grid_gdf = gpd.GeoDataFrame({'geometry': polygons}, crs="epsg:4326")

world = gpd.read_file("C:/Users/15052/Documents/GitHub/echopop_data/live_2019_files/coastline/ne_110m_land/ne_110m_land.shp")
bbox = box(lon_min - 0.25, lat_min - 0.25, lon_max + 0.25, lat_max + 0.25)
shapefile = world
clipped_shapefile = gpd.clip(shapefile, bbox).to_crs(utm_proj.srs)
clipped_shapefile.to_crs(utm_proj.srs)
# clipped_geometry = bbox.intersection(world.union_all())
# clipped_gdf = gpd.GeoDataFrame(geometry=[clipped_geometry], crs=world.crs)

from shapely.geometry import MultiPolygon
# Create an empty list to store clipped geometries
# clipped_geometries = []

# # Iterate over each grid polygon
# for index, row in grid_gdf.iterrows():
#     # Intersect grid polygon with land shape
#     intersection = row['geometry'].intersection(clipped_shapefile.unary_union)

#     # If intersection is a MultiPolygon, get the difference with the land shape
#     if isinstance(intersection, MultiPolygon):
#         clipped = row['geometry'].difference(clipped_shapefile.unary_union)
#         if clipped.is_empty:
#             continue
#         clipped_geometries.append(clipped)
#     else:
#         # If intersection is a single Polygon, directly add to clipped geometries
#         clipped_geometries.append(intersection)

# clipped_grid = gpd.GeoDataFrame(geometry=clipped_geometries, crs=grid_gdf.crs)

clipped_geometries = grid_gdf['geometry'].to_crs(utm_proj.srs).difference(clipped_shapefile.geometry.union_all())
clipped_gdf = gpd.GeoDataFrame(geometry=clipped_geometries)
clipped_gdf.to_crs(epsg=32610)

invalid_geometries = clipped_gdf[~clipped_gdf.is_valid]
clipped_gdf = clipped_gdf.buffer(0.001)
clipped_gdf['area_sqm'] = clipped_gdf.area / 46300.00000000001**2

clipped_gdf.area

fig, ax = plt.subplots(figsize=(10, 8))
clipped_gdf.plot(ax=ax, facecolor="none", edgecolor="black")
clipped_shapefile.plot(ax=ax, edgecolor='black', linewidth=0.5)
plt.tight_layout()
plt.show()


bbox.crs = {"init": "epsg:4326"}
intersection = gpd.overlay(bbox, world, how='intersection')

world_cut = gpd.sjoin(world, gpd.GeoDataFrame(geometry=[bbox]), how='inner', op='intersects')

world_cut = world[world.geometry.intersects(bbox)]
world_cut.to_crs("epsg:4326")

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 10))
grid_gdf.plot(ax=ax, facecolor="none", edgecolor="black")
world_cut.plot(ax=ax, linewidth=2, color='blue')
plt.show()

for cell in grid_gdf:

    x, y = cell.exterior.xy  # Extract x and y coordinates of the cell
    ax.fill(x, y, facecolor='none', edgecolor='black')  # Plot the cell as a polygon patch
# Plot coastline
# world.plot(ax=ax, linewidth=2, color='blue')
plt.show()


bbox = (lat_min, lon_min, lat_max, lon_max)
G = ox.graph_from_bbox(bbox[2], bbox[3], bbox[0], bbox[1], network_type='none', simplify=False)
G = ox.geometries_from_bbox(north=bbox[2], south=bbox[0], east=bbox[3], west=bbox[1], tags={'natural': ['coastline']})



latitudes = range(int(lat_min), int(lat_max) + 1, int(lat_step))
longitudes = range(int(lon_min), int(lon_max) + 1, int(lon_step))

# Initialize `meta` attribute
meta = copy.deepcopy(LIVE_DATA_STRUCTURE["meta"])

# Loading the configuration settings and definitions that are used to
# initialize the Survey class object
config = yaml.safe_load(Path(initialization_config).read_text())

nasc_frequency = config["acoustics"]["nasc_frequency"]