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
import os
import re
import contextlib
from sqlalchemy import create_engine, text, Engine, inspect
from echopop.live.live_core import LIVE_DATA_STRUCTURE, LIVE_FILE_FORMAT_MAP, LIVE_INPUT_FILE_CONFIG_MAP
from echopop.live import live_data_processing as eldp

####################################################################################################  
# TEST: YAML FILE CONFIGURATION
# ---- Define filepaths
live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"
# ---- Run function: `live_configuration`
file_configuration = live_configuration(live_init_config_path, live_file_config_path)
file_configuration.update({"database": {"acoustics": None, "biology": None}})
####################################################################################################
# * Accessory function for tuning the acoustic transmit frequency units/scaling
def format_vlaue(x):
    pass

def format_value(x):
    if isinstance(x, str):
        return "'{}'".format(x.replace("'", "''"))
    elif isinstance(x, pd.Timestamp):
        return "'{}'".format(x)
    elif x is None:
        return 'NULL'
    else:
        return str(x)

data_str = ", ".join(
    "({})".format(", ".join(format_value(x) for x in row))
    for row in data_tuple
)



####################################################################################################
# * Functionality for reading in processed acoustic data
# TODO: Expand data validator and limit cases to '*.zarr' (for now)
# TODO: Refactor "extra" components such as the validation steps, xarray-to-dataframe piping, etc.
# TODO: Documentation
file_settings = file_configuration["input_directories"]["acoustics"]
root_directory = file_configuration["data_root_dir"]


####################################################################################################  
# TEST: ACOUSTIC ZARR FILE INGESTION CONFIGURATION
# NOTE: 
# ---- Run function: `load_validated_acoustic_data` using previously defined `file_configuration`
acoustic_data = load_acoustic_data(file_configuration)
acoustic_data
file_configuration["database"]

def estimate_echometrics(acoustic_data_df: pd.DataFrame):

    # Create copy
    acoustic_df = acoustic_data_df.copy().reset_index(drop=True)

    # Pre-compute the change in depth
    acoustic_df["dz"] = acoustic_df["depth"].diff()

    # Initialize echometrics dictionary
    echometrics = {}

    # Compute the metrics center-of-mass
    if acoustic_df["NASC"].sum() == 0.0:
        echometrics.update({
            "n_layers": 0,
            "mean_Sv": -999,
            "max_Sv": -999,
            "nasc_db": np.nan,
            "center_of_mass": np.nan,
            "dispersion": np.nan,
            "evenness": np.nan,
            "aggregation": np.nan,    
            "occupied_area": 0.0,        
        })
    else:
        
        # Compute the number of layers
        echometrics.update({
            "n_layers": acoustic_df["depth"][acoustic_df["NASC"] > 0.0].size
        })

        # Compute ABC
        # ---- Convert NASC to ABC
        acoustic_df["ABC"] = acoustic_df["NASC"] / (4 * np.pi * 1852 ** 2)
        # ---- Estimate mean Sv
        echometrics.update({
            "mean_Sv": 10.0 * np.log10(acoustic_df["ABC"].sum() / acoustic_df["depth"].max())
        })
        # --- Estimate max Sv (i.e. )
        echometrics.update({
            "max_Sv": 10 * np.log10(acoustic_df["ABC"].max() 
                                    / acoustic_df.loc[np.argmax(acoustic_df["ABC"]), "dz"])
        })

        # Compute (acoustic) abundance
        echometrics.update({
            "nasc_db": 10 * np.log10(acoustic_df["ABC"].sum())
        })

        # Compute center of mass
        echometrics.update({
            "center_of_mass": (
                (acoustic_df["depth"] * acoustic_df["NASC"]).sum()
                / (acoustic_df["NASC"]).sum()
            )
        })

        # Compute the dispersion
        echometrics.update({
            "dispersion": (
                ((acoustic_df["depth"] - echometrics["center_of_mass"]) ** 2 
                * acoustic_df["NASC"]).sum() / (acoustic_df["NASC"]).sum()                
            )
        })

        # Compute the evenness
        echometrics.update({
            "evenness": (acoustic_df["NASC"] **2).sum() / ((acoustic_df["NASC"]).sum()) ** 2
        })

        # Compute the index of aggregation
        echometrics.update({
            "aggregation": 1 / echometrics["evenness"]
        })

        # Get the occupied area
        echometrics.update({
            "occupied_area": (
                acoustic_df["dz"][acoustic_df["ABC"] > 0.0].sum() / acoustic_df["depth"].max()
            )
        })

    # Return the dictionary
    return echometrics

def integrate_nasc(acoustic_data_df: pd.DataFrame, echometrics: bool = True):

    # Vertically integrate PRC NASC
    nasc_dict = {"nasc": acoustic_data_df["NASC"].sum()}
    
    # Horizontally concatenate `echometrics`, if `True`
    if echometrics:
        # ---- Compute values
        # NOTE: This uses NASC instead of linear `sv`
        echometrics_dict = estimate_echometrics(acoustic_data_df)
        # ---- Merge
        nasc_dict.update(echometrics_dict)

    # Convert `nasc_dict` to a DataFrame and return the output
    return pd.Series(nasc_dict)


acoustic_data_df = acoustic_data["prc_nasc_df"]



# SQL(database_file, "drop", table_name="nasc_df")
# SQL(database_file, "validate", **kwargs)
# SQL(database_file, "create", table_name="nasc_df", primary_keys=["latitude", "longitude", "ping_time"], dataframe=nasc_data_df)
# SQL(database_file, "validate", **kwargs)
# SQL(database_file, "select", table_name="nasc_df")
# SQL(database_file, "insert", table_name="nasc_df", id_columns=["latitude", "longitude", "ping_time"], dataframe=nasc_data_df)
# SQL(database_file, "select", table_name="nasc_df")
# SQL(database_file, "insert", table_name="nasc_df", id_columns=["latitude", "longitude", "ping_time"], dataframe=nasc_data_df)
# SQL(database_file, "select", table_name="nasc_df")
# SQL(database_file, "insert", table_name="nasc_df", dataframe=nasc_data_df)
# SQL(database_file, "drop", table_name="nasc_df")
# SQL_DTYPES[type(dataframe["ping_time"][0]).__name__]

def process_acoustic_data(acoustic_data_df: pd.DataFrame, file_configuration: dict, 
                          echometrics: bool = True):

    # Integrate NASC (and compute the echometrics, if necessary)
    nasc_data_df = (
        acoustic_data_df.groupby(["longitude", "latitude", "ping_time"])
        .apply(lambda group: integrate_nasc(group, echometrics), include_groups=False)
        .reset_index()
    )
    # ---- Amend the dtypes if echometrics were computed
    if echometrics:
        nasc_data_df = (
            nasc_data_df
            .astype({"n_layers": int, "mean_Sv": float, "max_Sv": float, "nasc_db": float,
                    "center_of_mass": float, "dispersion": float, "evenness": float,
                    "aggregation": float, "occupied_area": float})
        )

    # Get the acoustics database file
    acoustics_db = file_configuration["database"]["acoustics"]

    # Insert the new data into the database and pull in the combined previous and new data combined
    full_nasc_df = sql_data_exchange(acoustics_db, dataframe=nasc_data_df, 
                                     table_name="nasc_df", 
                                     id_columns=["longitude", "latitude", "ping_time"],
                                     primary_keys=["longitude", "latitude", "ping_time"],
                                     output_type=pd.DataFrame)

    # Return the output
    return full_nasc_df

####################################################################################################
def reset_db_files(file_configuration: dict, table_exception: Optional[Union[str, List[str]]] = None):

    # Get all database files
    database_files = file_configuration["database"]

    # Iterate through all keys
    for _, db_file in database_files.items():
        # ---- Map the table names
        table_names = SQL(db_file, "map")
        # ---- Drop any noted exceptions
        if not isinstance(table_exception, list):
            table_exception = [table_exception]
        # ---- Drop exception table name
        if None not in table_exception:
            table_names = list(set(table_names) - set(table_exception))
        _ = [SQL(db_file, "drop", table_name=table) for table in table_names]
        # ---- Validate that all tables were removed        
        if set(table_names).intersection(set(SQL(table_names, "map"))):
            raise ValueError(
                f"Attempted reset of [{str(db_file)}] failed."
            )

SPATIAL_CONFIG_MAP = {
    "closest_haul": {
        "proximity": {
            "choices": ["distance", "time"],
        },
    },
    "global" : {},
    "griddify": {
        "bounds": {
            "longitude": {
                "types": [float]
            },
            "latitude": {
                "types": [float]
            },
            "northings": {
                "types": [float]
            },
            "eastings": {
                "types": [float]
            },
            "pairs": [("longitude", "latitude"), ("northings", "eastings")],
        },
        "grid_resolution": {
            "x_distance": {
                "types": float,
            },
            "y_distance": {
                "types": float,
            },
            "d_longitude": {
                "types": float,
            },
            "d_latitude": {
                "types": float,
            },
            "grid_size_x": {
                "types": int,
            },
            "grid_size_y": {
                "types": int,
            },
            "pairs": [("x_distance", "y_distance"), ("d_longitude", "d_latitude"), 
                      ("grid_size_x", "grid_size_y")],       
        },
    },
    "inpfc": {
        "stratum_names": {
                "types": [int, str]
            },
        "latitude_max": {
            "types": [float],
        },
    },
    "weighted_haul": {
        "proximity": {
            "choices": ["distance", "time"]
        },
    },
}



reset_db_files(file_configuration, table_exception = "files_read")
reset_db_files(file_configuration)

stamp = 20240714194248
stamp.astype(int)
int(stamp)
import re
from datetime import datetime

def infer_datetime_format(timestamp_str: Union[int, str]):
    patterns = {
        r"^\d{14}$": "%Y%m%d%H%M%S",             # YYYYMMDDHHMMSS
        r"^\d{8}$": "%Y%m%d",                     # YYYYMMDD
        r"^\d{6}$": "%H%M%S",                     # HHMMSS
        r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$": "%Y-%m-%d %H:%M:%S",  # YYYY-MM-DD HH:MM:SS
        r"^\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}$": "%Y/%m/%d %H:%M:%S",  # YYYY/MM/DD HH:MM:SS
        r"^\d{4}-\d{2}-\d{2}$": "%Y-%m-%d",       # YYYY-MM-DD
        r"^\d{4}/\d{2}/\d{2}$": "%Y/%m/%d"        # YYYY/MM/DD
    }
    
    for pattern, date_format in patterns.items():
        if re.match(pattern, timestamp_str):
            return date_format
    
    raise ValueError("Unknown timestamp format")

filter_dict = dict(species_filer=species_filter, trawl_filter=trawl_filter)

def biology_data_filter(biology_data: pd.DataFrame, filter_dict: dict):

    # Create dataframe copy
    data_copy = biology_data.copy()

    # Iterate through dictionary to apply filters (if present)
    for column, value in filter_dict.items():
        if column in data_copy.columns:
            data_copy = data_copy[data_copy[column] == value]

    # Return output
    return data_copy



df[(df['species_id'] == species_filter if 'species_id' in df.columns else True)]
df[(df["species_id"] == 17 if "species_id" in df.columns)]

(df[df["haul_num"] == 17 if "haul_num" in df.columns] else True)


from datetime import datetime

df = biology_output["trawl_info_df"]
df.loc[(df['species_id'] == species_filter if 'species_id' in df.columns else True), :]
df.index

biology_output["trawl_info_df"].reset_index().index
df = biology_output["catch_df"]
df = df.loc[0, :].to_frame().T
df.index
df.loc[(df['species_id'] == species_filter if 'species_id' in df.columns else True)]

def convert_datetime(timestamp: Union[int, str, pd.Series]):

    if isinstance(timestamp, pd.Series):
        test_timestamp = str(timestamp[0])
    else:
        test_timestamp = str(timestamp)

    # Approximate the datetime format
    datetime_format = infer_datetime_format(str(test_timestamp))

    #
    if isinstance(timestamp, pd.Series):
        return timestamp.apply(lambda x: datetime.strptime(x, datetime_format))
    else:
        return datetime.strptime(timestamp, datetime_format)
    
infer_datetime_format(stamp)
convert_datetime(stamp)
infer_datetime_format(202407)

# {'global': False, 'INPFC': True, 'closest_haul': False, 'weighted_haul': False}
file_configuration["geospatial"]["link_biology_acoustics"] = "INPFC"
file_configuration["geospatial"]
spatial_config = file_configuration["geospatial"]
###############

acoustic_data = self.input["acoustics"]
biology_data = self.input["biology"]

def load_spatial_data(acoustic_data: dict,
                      biology_data: dict,                      
                      file_configuration: dict,):
    
    # Extract spatial strata *only* if spatial information from the configuration settings
    # ---- Get (geo)spatial config
    spatial_config = file_configuration["geospatial"]
    # ---- Remove case sensitivity
    spatial_config = {key.lower(): value for key, value in spatial_config.items()}
    # ---- Extract the projection
    projection = spatial_config["projection"]
    # ---- Extract the biology-acoustics linking method options
    acoustics_biology_link = spatial_config["link_biology_acoustics"]

    # Validate the configuration
    validate_spatial_config(spatial_config)

    # Assign the spatial link constraints to the acoustic and biological data
    if acoustics_biology_link == "INPFC":
        apply_inpfc_definitions(acoustic_data, biology_data, spatial_config)



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
    
####################################################################################################  
# TEST: BIOLOGY FILE INGESTION CONFIGURATION
# NOTE: 
# ---- Run function: `load_validated_acoustic_data` using previously defined `file_configuration`
biology_data, file_configuration = load_biology_data(file_configuration)
biology_data
####################################################################################################
prc_nasc_df = acoustic_data["prc_nasc_df"]

def process_acoustic_data(acoustic_data_df: pd.DataFrame, file_configuration: dict, 
                          echometrics: bool = True):

    # Integrate NASC (and compute the echometrics, if necessary)
    nasc_data_df = (
        acoustic_data_df.groupby(["longitude", "latitude", "ping_time"])
        .apply(lambda group: integrate_nasc(group, echometrics), include_groups=False)
        .reset_index()
    )
    # ---- Amend the dtypes if echometrics were computed
    if echometrics:
        nasc_data_df = (
            nasc_data_df
            .astype({"n_layers": int, "mean_Sv": float, "max_Sv": float, "nasc_db": float,
                             "center_of_mass": float, "dispersion": float, "evenness": float,
                             "aggregation": float, "occupied_area": float})
        )

    # Get the name of the associated db file
    acoustics_db = file_configuration["database"]["acoustics"]
    # ---- Get current tables
    tables = SQL(acoustics_db, "inspect")
    
    # 
    if "nasc_df" not in tables:
        _ = SQL(acoustics_db, "insert", table_name="nasc_df", dataframe=nasc_data_df)
    else:
        # ---- 
        nasc_sql = SQL(acoustics_db, "select", table_name="nasc_df")
        # ----
        index_equiv = nasc_data_df[["longitude", "latitude", "ping_time"]].isin(nasc_sql)
        # ----
        bool_idx = index_equiv.apply(lambda x: np.all(x), axis=1)
        # ---- 
        _ = SQL(acoustics_db, "insert", table_name="nasc_df", dataframe=nasc_data_df.loc[~bool_idx])
        # ----
        nasc_data_df = pd.concat([nasc_sql, nasc_data_df], ignore_index=True)

    # Return the output
    return nasc_data_df


SQL(acoustics_db, command="drop", table_name="nasc_df")
SQL(acoustics_db, "inspect")

nasc_analysis = process_acoustic_data(acoustic_data["prc_nasc_df"], file_configuration)

SQL(acoustics_db, command="select", table_name="nasc_df")

TS_SLOPE = 20.0
TS_INTERCEPT = -68.0

# CONVERT TO TS
comb_lengths["ts"] = TS_SLOPE * np.log10(comb_lengths["length"]) + TS_INTERCEPT
# TO SIGMA_BS
comb_lengths["sigma_bs"] = 10 ** (comb_lengths["ts"] / 10)
# WEIGHTED MEAN SIGMA_BS
sigma_mean = np.average(comb_lengths["sigma_bs"], weights=comb_lengths["length_count"])

from typing import Optional
from echopop.utils import operations
from echopop.acoustics import ts_length_regression, to_linear, to_dB

__all__ = ["operations"]

# Meld bio datasets
length_datasets = biology_data["specimen_df"].meld(biology_data["length_df"], 
                                                   contrasts=["haul_num", "sex", "species_id", "length"])

# Create distribution
distrib_params = file_configuration["biology"]["length_distribution"]["bins"]

length_bins = np.linspace(**{key: value for key, value in zip(["start", "stop", "num"], distrib_params)}, dtype=float)
binwidth = np.diff(length_bins / 2.0).mean()
intervals = np.concatenate([length_bins[:1] - binwidth, length_bins + binwidth])
length_bins_df = pd.DataFrame({"bin": length_bins, "interval": pd.cut(length_bins, intervals)})
# 
length_datasets["length_bin"] = pd.cut(length_datasets["length"], bins=intervals, labels=length_bins_df["bin"])

stratify_key = file_configuration["geospatial"]["link_biology_acoustics"]

if stratify_key == "global":
    length_distribution = (
        length_datasets.pivot_table(columns=["sex"], index=["length_bin"], 
                                    values="length_count", aggfunc="sum", observed=False)
    )
    #
    length_distribution["total"] = length_distribution.sum(axis=1)

length_distribution.transpose()
SQL(biology_db, "drop", table_name="length_distribution")
# Get the name of the associated db file
biology_db = file_configuration["database"]["biology"]
# ---- Get current tables
tables = SQL(biology_db, "inspect")


if "length_distribution" not in tables:
    _ = SQL(biology_db, "insert", table_name="length_distribution", 
            dataframe=length_distribution.transpose())
    

SQL(biology_db, "select", table_name="length_distribution")
SQL(biology_db, "drop", table_name="length_distribution")
SQL(biology_db, "replace", table_name="length_distribution", dataframe=length_distribution.unstack().reset_index(name="count"))
length_distribution.unstack().reset_index(name="count")
mixed = SQL(biology_db, "select", table_name="length_distribution")
length_bins[:1]
from typing import Optional
from echopop.utils import operations
from echopop.acoustics import ts_length_regression, to_linear, to_dB

__all__ = ["operations"]

# Meld bio datasets
length_datasets = biology_data["specimen_df"].meld(biology_data["length_df"], 
                                                   contrasts=["haul_num", "species_id", "length"])

ts_length_parameters_spp = [
    spp
    for spp in file_configuration["acoustics"]["TS_length_regression_parameters"].values()
    if spp["number_code"] in np.unique(length_datasets.species_id).astype(int)
]

# ---- get species info
target_species = pd.DataFrame.from_dict(ts_length_parameters_spp)

ts_lengths_df = length_datasets.merge(
    target_species.drop("length_units", axis=1),
    left_on=["species_id"],
    right_on=["number_code"],
)
# ---- filter out other spp
length_datasets[length_datasets["species_id"].isin(target_species["number_code"])]

#
file_configuration["acoustics"]["TS_length_regression_parameters"][target_species["text_code"]]

def average_sigma_bs(length: Union[pd.DataFrame, float, int], TS_L_slope: Optional[float] = None, TS_L_intercept: Optional[float] = None, weighted: Optional[Union[float, int, str]] = None):

    # 
    if isinstance(length, pd.DataFrame):
        if "length" not in length.columns: 
            raise ValueError(
                "Column [`length`] missing from dataframe input `length`."
            )
        if "TS_L_slope" not in length.columns and TS_L_slope is None:
            raise ValueError(
                "Value [`TS_L_slope`] missing from dataframe input `length` and optional "
                "separate argument `TS_L_slope`."
            )
        if "TS_L_intercept" not in length.columns and TS_L_intercept is None:
            raise ValueError(
                "Value [`TS_L_intercept`] missing from dataframe input `length` and optional "
                "separate argument `TS_L_intercept`."
        )
    elif isinstance(length, float) or isinstance(length, int):
        if TS_L_slope is None:
            raise ValueError(
                "Argument [`TS_L_slope`] missing."
            )
        elif TS_L_slope is not None and not isinstance(TS_L_slope, float):
            raise TypeError(
                "Argument `TS_L_slope` must be type `float`."
        )
        if "TS_L_intercept" not in length.columns and TS_L_intercept is None:
            raise ValueError(
                "Argument [`TS_L_intercept`] missing."
        )
        elif TS_L_intercept is not None and not isinstance(TS_L_intercept, float):
            raise TypeError(
                "Argument `TS_L_intercept` must be type `float`."
        )

    #
    if TS_L_slope is None:
        TS_L_slope = length["TS_L_slope"]

    #
    if TS_L_intercept is None:
        TS_L_intercept = length["TS_L_intercept"]

    #
    if isinstance(length, pd.DataFrame):
        length_val = length["length"]

    ts_value = ts_length_regression(length_val, TS_L_slope, TS_L_intercept)
    sigma_bs_value = to_linear(ts_value)



    if isinstance(weighted, str):
        if weighted not in length.columns:
            raise ValueError(
                f"Argument [`weighted` (str)], '{weighted}', is not a column in argument `length` "
                f"(DataFrame)."
            )
        else: 
            return (sigma_bs_value * length[weighted]).sum() / length[weighted].sum()
    elif weighted is not None: 
        if weighted.size != sigma_bs_value.size:
            raise ValueError(
                f"Argument [`weighted` (float|int)] of size {weighted.size} does not match size of "
                f"argument [`length` (float|int)`] of size {sigma_bs_value.size}."
            )
        else:
            return (sigma_bs_value * weighted).sum() / weighted.sum()
    else:
        return sigma_bs_value.mean()

average_sigma_bs

ts_lengths_df.groupby(["haul_num"]).apply(average_sigma_bs).apply(lambda x: to_dB(x))
def integrate_nasc(prc_nasc_df: pd.DataFrame):

# Compute the number of layers
echometrics.update({
    "n_layers": acoustic_df["depth"][acoustic_df["NASC"] > 0.0].size
})

# Compute the index of aggregation
echometrics.update({
    "aggregation": 1 / echometrics["evenness"]
})

# Get the occupied area
echometrics.update({
    "occupied_area": (
        acoustic_df["dz"][acoustic_df["ABC"] > 0.0].sum() / acoustic_df["depth"].max()
    )
})




pd.read_fr
pd.read_sql(text(SQL_COMMANDS["select"].format(**kwargs)), con=connection)
engine = create_engine(f"sqlite:///{db_file}")
connection = engine.connect()
kwargs["dataframe"].to_sql(name=kwargs["table_name"], 
                                                  con=connection, 
                                                  if_exists="append", index=False)
connection.close()
engine.dispose()
SQL(db_file, "insert", table_name=table_name, columns="*", 
        filter_columns=insertion_filter,
        dataframe=df)

SQL(db_file, "select", table_name="files_read")
SQL(db_file, "select", table_name="catch_df")
SQL(db_file, "select", table_name="specimen_df")
SQL(db_file, "select", table_name="length_df")

def check_table_schema(connection, **kwargs):
    query = text(("PRAGMA table_info({table_name});").format(**kwargs))
    schema = connection.execute(query).fetchall()
    print("Table Schema:", schema)

check_table_schema(connection, table_name=table_name)

def insert_test_data(connection, table_name):
    test_data = pd.DataFrame({
        'trawl_partition': ['test'],
        'species_id': ['test'],
        'haul_weight': [0.0],
        'catch_percentage': [0.0],
        'haul_num': [1]
    })
    
    test_data.to_sql(name=table_name, con=connection, if_exists='append', index=False)
    print("Test data inserted.")

insert_test_data(connection, table_name)

kwargs = {}
command = "insert"
kwargs["table_name"] = "catch_df"
kwargs["dataframe"] = df
kwargs["filter_columns"] = insertion_filter
columns = "*"


re.compile(file_name_format)
pattern = file_name_format
pattern = pattern.replace('{DATE:YYYYMM}', r'(?P<DATE>\d{6})')
pattern = pattern.replace('{HAUL}', r'(?P<HAUL>\d+)')
pattern = pattern.replace('{FILE_ID}', r'(?P<FILE_ID>.+)')
regex = re.compile(pattern)
haul_values = []

file_name_format.search(file.name)
sub_df_lst = []
for file in subcsv_files:
    match = regex.search(file.name)
    if match:
        haul_value = match.group('HAUL')
        df = pd.read_csv(file, usecols=list(sub_config_map.keys()))
        df['HAUL'] = haul_value  # Append HAUL value as a new column
        sub_df_lst.append(df)
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

lat = 55.5000
lon = -134.2500
utm_code = int(utm_string_generator(lon, lat))
utm_proj = pyproj.Proj(f"epsg:{utm_code}")
utm_proj(lon, lat)
gpd.GeoDataFrame(geometry=gpd.points_from_xy(np.array([lon]), np.array([lat])), crs=projection).to_crs(utm_code)


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