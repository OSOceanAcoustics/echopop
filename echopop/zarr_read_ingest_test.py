import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union, Tuple, Optional
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
from echopop.live.live_core import LIVE_DATA_STRUCTURE, LIVE_FILE_FORMAT_MAP, LIVE_INPUT_FILE_CONFIG_MAP, SPATIAL_CONFIG_MAP
from echopop.live.live_data_loading import validate_data_directory
from echopop.live.sql_methods import SQL, SQL_COMMANDS, query_processed_files, format_sql_columns
from echopop.live import live_data_processing as eldp
from echopop.live import live_data_loading as eldl
from echopop.live.live_survey import LiveSurvey
from echopop.live.live_acoustics import preprocess_acoustic_data
from echopop.live.live_biology import preprocess_biology_data
from echopop.survey import Survey

survey_2019 = Survey("C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml", "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")
survey_2019.transect_analysis()
survey_2019.analysis["transect"]["biology"]["weight"]["weight_stratum_df"]
analysis_dict = survey_2019.analysis["transect"]

proportions_dict=analysis_dict["biology"]["proportions"]["number"]
length_weight_dict = analysis_dict["biology"]["weight"]
stratum_proportions_sexed["proportion_aged"] + stratum_proportions_sexed["proportion_unaged"]
####################################################################################################  
# TEST: YAML FILE CONFIGURATION
# ---- Define filepaths
self = LiveSurvey
live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"
# ---- Run function: `live_configuration`
file_configuration = self.config
files = biology_files

biology_output = initial_biology_output
file_configuration = self.config
table_name = "length_df"
df = filtered_biology_output[table_name]
database_file = biology_db
kwargs = dict(dataframe=df, table_name=table_name, id_columns=["id"], primary_keys=["id"], output_type=pd.DataFrame)

def process_biology_data(self):

    # Compute `sigma_bs` by sending it to the appropriate database table
    compute_sigma_bs(biology_unprocessed["specimen_df"], biology_unprocessed["length_df"], 
                     self.config)
    
    # Bin the length measurements of the biological data
    bin_length_data(biology_unprocessed, self.config["length_distribution"])

    # Compute the length-weight regression and add it to the SQL table
    length_weight_df = length_weight_regression(biology_unprocessed["specimen_df"], 
                                                self.config["length_distribution"],
                                                self.config)
    
    # Compute length-binned counts for the aggregated and individual-based measurements
    specimen_binned, specimen_binned_filtered, length_binned = (
        length_bin_counts(biology_unprocessed["length_df"], biology_unprocessed["specimen_df"], 
                          self.config)
    )

    # Compute the number proportions
    specimen_number_proportion, length_number_proportion, sex_number_proportions = (
        number_proportions(specimen_binned, specimen_binned_filtered, length_binned,
                           self.config)
    )
    
    # Compute the length-binned weights for the aggregated and individual-based measurements
    length_weight_binned, specimen_weight_binned = (
        length_bin_weights(biology_unprocessed["length_df"],
                           biology_unprocessed["specimen_df"],
                           length_weight_df,self.config)
    )

    # Calculate the average weights among male, female, and all fish
    fitted_weight_df = compute_average_weights(specimen_number_proportion, 
                                               length_number_proportion, 
                                               sex_number_proportions,
                                               length_weight_df,
                                               self.config["length_distribution"],
                                               self.config)

catch_data = self.input["biology"]["catch_df"]

# Get the spatial column name, if there is one
spatial_column = file_configuration["spatial_column"]
# ---- Append additional columns that will be used
contrast_columns = spatial_column + ["sex", "species_id"]

# Calculate grouped totals
# ---- Sum the net haul weights from station 1/unaged fish
catch_weights = catch_data.count_variable(
    contrasts=["species_id"] + spatial_column, 
    variable="haul_weight", fun="sum"
)
# ---- Rename resulting columns for both
catch_weights.rename(columns={"count": "total_weight"}, inplace=True)

# ---- Specimen
specimen_weights = specimen_weight_binned.sum().reset_index(name="total_weight")

specimen_weight_binned
# Calculate the sexed and total stratum weights for each sex among unaged fish
# ---- Sum the net haul weights from station 1/unaged fish
catch_weights = catch_data.count_variable(
    contrasts=["species_id"] + file_configuration["spatial_column"], 
    variable="haul_weight", fun="sum"
)
# ---- Rename resulting columns for both
catch_weights.rename(columns={"count": "total_weight"}, inplace=True)

# For the specimen data 
# ---- Sum the net haul weights from station 1/unaged fish
# ---- Specimen
specimen_weights_sex = (
    specimen_weight_binned
    .groupby(contrast_columns)["weight"]
    .sum()
)
# ---- Total (per stratum, if it exists)
specimen_weight_total = specimen_weights_sex.transpose().unstack(1).sum(axis=1)

# For the length (unaged) dataset
length_weights_sex = (
    length_weight_binned
    .groupby(contrast_columns)["weight_interp"]
    .sum()
)
# ---- Further reduce to the grand total (per stratum, if it exists)
length_weight_total = length_weights_sex.transpose().unstack(1).sum(axis=1)

# ---- Standardize the unaged sexed weights
length_weight_standardized = (
    (length_weights_sex / length_weight_total).unstack(0) 
    * catch_weights["total_weight"].to_numpy()
)

# Calculate the specimen weight proportions
# ---- Pivot weight bins
specimen_weight_binned_pvt = (
    specimen_weight_binned.pivot_table(
        columns=spatial_column,
        index=["length_bin", "species_id", "sex"],
        values="weight",
        observed = False
    )
)
# ---- Divide by the aged stratum weights (relative to only aged fish)
specimen_weight_proportions_pvt = (
    specimen_weight_binned_pvt / specimen_weight_total.to_numpy()
)
# ---- Pivot back to the desired format
specimen_weight_proportion = (
    specimen_weight_proportions_pvt
    .stack().reset_index(name="weight_proportion")
    .pivot_table(columns=stratum_column + ["species_id", "sex"], 
                 index="length_bin", values="weight_proportion")
)    
# ---- Calculate the internal (i.e. only aged fish) for each sex
within_specimen_sex_proportions = (
    specimen_weight_proportion.sum()
)

# Calculate the total strata weights
# ---- Index `catch_weights`
catch_weights_idx = catch_weights.set_index(stratum_column + ["species_id"])
# ---- Compute the spatially-stratified/grouped weights
spatial_weights = (
    pd.concat([specimen_weight_total.to_frame("total_weight"), catch_weights_idx])
    .pivot_table(
        columns=stratum_column, 
        aggfunc="sum", 
        values="total_weight", 
        observed=False
    )
)

# Calculate the weight proportions relative to the overall stratum weights
# ---- Aged
# -------- Reformat into dataframe and merge with total stratum weights
specimen_weights_binned_df = (
    specimen_weight_binned_pvt.stack()
    .to_frame("specimen_weight")
    .reset_index()
    .merge(spatial_weights.T.reset_index(), on=stratum_column)
)
# -------- Calculate proportions
specimen_weights_binned_df["weight_proportion_overall"] = (
    specimen_weights_binned_df["specimen_weight"] / specimen_weights_binned_df["total_weight"]
)
# -------- Consolidate to calculate the sexed proportions per stratum
specimen_weight_sex_proportions = specimen_weights_binned_df.groupby(stratum_column + ["species_id", "sex"])[
    "weight_proportion_overall"
].sum()
# ---- Unaged
# -------- Reformat into dataframe and merge with total stratum weights
length_weights_sex_standardized_df = (
    length_weight_standardized.stack()
    .to_frame("catch_weight")
    .reset_index()
    .merge(spatial_weights.T.reset_index(), on=stratum_column)
)
# -------- Calculate proportions
length_weights_sex_standardized_df["weight_proportion_overall"] = (
    length_weights_sex_standardized_df["catch_weight"]
    / length_weights_sex_standardized_df["total_weight"]
)
# -------- Back-calculate the sexed weight proportions relative to just unaged fish
# ------------ Aggregate proportions
length_total_sex_proportions = length_weights_sex_standardized_df.pivot_table(
    columns=["species_id", "sex"], index=stratum_column, values="weight_proportion_overall"
).transpose().unstack(["species_id"]).sum(axis=0)
# ------------ Re-compute the proportions
length_weight_sex_proportions = (
    length_weights_sex_standardized_df.pivot_table(
        index=["species_id", "sex"], columns=stratum_column, 
        values="weight_proportion_overall"
    )
    / length_total_sex_proportions.to_numpy()
)

# Compute the overall length-binned weight distributions among unaged fish
# ---- Extract the number proportions computed for unaged fish
length_number_proportions = length_number_proportion.copy()
# ---- Filter out values besides those computed for 'all' fish
length_number_proportions = length_number_proportions[length_number_proportions["sex"] == "all"]
# ---- Convert to a table
length_number_proportions_tbl = length_number_proportions.pivot_table(
    columns=stratum_column + ["species_id"],
    index=["length_bin"],
    values="proportion_number_length",
    aggfunc="sum",
    observed=False,
)
# ---- Extract the fitted weight values calculated for all fish
length_weight_all = length_weight_df[length_weight_df["sex"] == "all"]
# ---- Generate the fitted weight array
fitted_weights = length_weight_all.copy()
# ---- Get actual length bins in dataset
fitted_weights = fitted_weights[fitted_weights["length_bin"].isin(length_number_proportions["length_bin"])]
# ---- Apportion the averaged weights
length_apportioned_weights = length_number_proportions_tbl.T * fitted_weights["weight_fitted"].to_numpy()
# ---- Compute the average weight proportions per length bin per stratum
average_length_bin_weights = length_apportioned_weights.T / length_apportioned_weights.sum(axis=1)
# ---- Convert back to a DataFrame
average_length_bin_weights_df = average_length_bin_weights.unstack().reset_index(
    name="weight_proportion"
)

# Calculate the aged and unaged weight proportions
# ---- Aged
aged_proportions = specimen_weight_sex_proportions.unstack("sex").sum(axis=1)
# ---- Unaged
unaged_proportions = 1 - aged_proportions
# -------- Re-weight the unaged sexed proportions
unaged_weight_sex_proportions_overall = (
    (length_weight_sex_proportions * unaged_proportions.unstack().transpose()).astype(float).fillna(0.0)
)

unaged_proportions.unstack().transpose()
# Format the outputs
# ---- Aged: stratum-sex-age-length relative to aged and total weights
aged_overall_df = (
    specimen_weight_proportion.unstack()
    .reset_index(name="weight_proportions")
    .merge(
        specimen_weights_binned_df[
            stratum_column + ["length_bin", "sex", "species_id", "weight_proportion_overall"]
        ]
    )
)
# ---- Aged: stratum-sex relative to total weights
aged_sex_df =within_specimen_sex_proportions.reset_index(name="weight_proportion_aged").set_index(
        stratum_column + ["species_id", "sex"]
    )
# ---- Add the aged sex proportiosn relative to the overall survey
aged_sex_df["weight_proportion_overall_aged"] = specimen_weight_sex_proportions
# ---- Consolidate the aged and unaged sexed dataframes
# -------- Initialize the dataframe
aged_unaged_sex_proportions = aged_sex_df.reset_index().set_index(["species_id", "sex"] + stratum_column)
# --------- Add the within-unaged weight proportions
aged_unaged_sex_proportions["weight_proportion_unaged"] = (
    length_weight_sex_proportions.stack()
)
# --------- Add the overall-unaged weight proportions
aged_unaged_sex_proportions["weight_proportion_overall_unaged"] = (
    unaged_weight_sex_proportions_overall.stack()
)
# ---- Overall aged and unaged proportions
aged_unaged_proportions = aged_proportions.reset_index(name="aged_proportions")
# ---- Set index
aged_unaged_proportions.set_index(stratum_column + ["species_id"], inplace=True)
# -------- Add unaged proportions
aged_unaged_proportions["unaged_proportions"] = unaged_proportions#.reset_index()
# ---- Reset the index
aged_unaged_proportions = aged_unaged_proportions.reset_index()
####################################################################################################
# * Functionality for reading in processed acoustic data
# TODO: Expand data validator and limit cases to '*.zarr' (for now)
# TODO: Refactor "extra" components such as the validation steps, xarray-to-dataframe piping, etc.
# TODO: Documentation
file_settings = file_configuration["input_directories"]["acoustics"]
root_directory = file_configuration["data_root_dir"]


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



from echopop.live.live_core import SPATIAL_CONFIG_MAP

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

    # Create spatial dictionary that will be added as an `input`
    spatial_dict = {"link_method": acoustics_biology_link}

    # Assign the spatial link constraints to the acoustic and biological data
    if acoustics_biology_link == "INPFC":
        spatial_dict.update({"strata": create_inpfc_strata(spatial_config)})

    # Return the dictionary as an output
    return spatial_dict



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
        .apply(lambda group: integrate_nasc(group, echometrics))
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

biology_data = self.input["biology"]

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

def average_sigma_bs(length: Union[pd.DataFrame, float, int], 
                     TS_L_slope: Optional[float] = None, 
                     TS_L_intercept: Optional[float] = None, 
                     weighted: Optional[Union[float, int, str]] = None):

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

def parse_condition(condition):
    # Handle nested conditions and logical operators
    condition = condition.replace('&', ' AND ').replace('|', ' OR ')

    # Handle "IN" lists and replace square brackets with parentheses
    condition = re.sub(r'(\w+)\s*IN\s*\[(.*?)\]', lambda m: f"{m.group(1)} IN ({m.group(2)})", condition, flags=re.IGNORECASE)
    
    # Handle range conditions for BETWEEN, including floats
    condition = re.sub(r'(\d*\.\d+|\d+)\s*<=\s*(\w+)\s*<=\s*(\d*\.\d+|\d+)', 
                       lambda m: f"{m.group(2)} BETWEEN {m.group(1)} AND {m.group(3)}", condition)
    
    # Handle individual comparisons
    condition = re.sub(r'(\w+)\s*([<>!=]+)\s*(\d*\.\d+|\d+)', lambda m: f"{m.group(1)} {m.group(2)} {m.group(3)}", condition)
    condition = re.sub(r'(\w+)\s*([<>!=]+)\s*(\'[^\']*\')', lambda m: f"{m.group(1)} {m.group(2)} {m.group(3)}", condition)

    # Handle single equal sign
    condition = re.sub(r'(\w+)\s*=\s*(\d*\.\d+|\d+)', lambda m: f"{m.group(1)} = {m.group(2)}", condition)

    # Remove redundant spaces
    condition = re.sub(r'\s+', ' ', condition).strip()

    return condition

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