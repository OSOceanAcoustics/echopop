from echopop.survey import Survey
import copy
from pathlib import Path
from typing import List, Literal, Optional, Union
import numpy as np
import pandas as pd
import os
import re

from echopop.core import ECHOVIEW_EXPORT_MAP, REGION_EXPORT_MAP, BIODATA_HAUL_MAP, DATA_STRUCTURE, LAYER_NAME_MAP, NAME_CONFIG
from echopop.spatial.transect import export_transect_layers, export_transect_spacing
from echopop.utils.data_structure_utils import map_imported_datasets
from echopop.utils.validate_df import KSStrata, DATASET_DF_MODEL
from echopop.utils.operations import compile_patterns, extract_parts_and_labels, group_merge
from echopop.utils.load_nasc import validate_echoview_exports, validate_export_directories, get_transect_numbers, filter_export_regions, get_haul_transect_key, compile_patterns, consolidate_exports, construct_transect_region_key


survey = Survey( init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                 survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml" )

self = survey
index_variable: Union[str, List[str]] = ["transect_num", "interval"]
ingest_exports: Optional[Literal["echoview", "echopype"]] = "echoview"
region_class_column: str = "region_class"
transect_pattern: str = r"T(\d+)"
unique_region_id: str = "region_id"
verbose: bool = True
configuration_dict = self.config
input_dict = self.input
dataset_type = ["biological", "kriging", "stratification"]

# Re-initialize the input keys, if needed
# ---- Get name of the proposed input dictionary keys
input_keys = [
    (
        LAYER_NAME_MAP[key]["superlayer"][0]
        if LAYER_NAME_MAP[key]["superlayer"]
        else LAYER_NAME_MAP[key]["name"]
    )
    # for key in CONFIG_MAP.keys()
    for key in DATASET_DF_MODEL.keys()
    if key in list(dataset_type)
]

# ---- Map the complete datasets
imported_data = map_imported_datasets(input_dict)
# ---- Re-initialize if data already loaded to avoid duplication issues
if set(input_keys).issubset(imported_data):
    # ---- Reset the relevant keys
    input_dict.update({key: copy.deepcopy(DATA_STRUCTURE["input"][key]) for key in input_keys})

# Check whether data files defined from the configuration file exists
# ---- Generate flat JSON table comprising all configuration parameter names
flat_configuration_table = pd.json_normalize(configuration_dict).filter(regex="filename")

# Get the subset table if specific `dataset_type` is defined
if dataset_type:
    # ---- Get the outermost dictionary keys
    outer_keys = flat_configuration_table.columns.str.split(".").str[0]
    # ---- Get the associated column names
    matching_columns = flat_configuration_table.columns[outer_keys.isin(dataset_type)]
    # ---- Filter the columns
    flat_configuration_table = flat_configuration_table.filter(matching_columns)
# ---- Default to `CONFIG_MAP` keys otherwise
else:
    # dataset_type = list(CONFIG_MAP.keys())
    dataset_type = list(DATASET_DF_MODEL.keys())
# ---- Parse the flattened configuration table to identify data file names and paths
parsed_filenames = flat_configuration_table.values.flatten()
# ---- Evaluate whether either file is missing
data_existence = [
    (Path(configuration_dict["data_root_dir"]) / file).exists() for file in parsed_filenames
]

# Assign the existence status to each configuration file for error evaluation
# ---- Error evaluation and print message (if applicable)
if not all(data_existence):
    missing_data = parsed_filenames[~np.array(data_existence)]
    raise FileNotFoundError(f"The following data files do not exist: {missing_data}")

# Get the subset table if specific `dataset_type` is defined
if dataset_type:
    # ---- Get the outermost dictionary keys
    outer_keys = flat_configuration_table.columns.str.split(".").str[0]
    # ---- Get the associated column names
    matching_columns = flat_configuration_table.columns[outer_keys.isin(dataset_type)]
    # ---- Filter the columns
    flat_configuration_table = flat_configuration_table.filter(matching_columns)
# ---- Default to `CONFIG_MAP` keys otherwise
else:
    # dataset_type = list(CONFIG_MAP.keys())
    dataset_type = list(DATASET_DF_MODEL.keys())
# ---- Parse the flattened configuration table to identify data file names and paths
parsed_filenames = flat_configuration_table.values.flatten()
# ---- Evaluate whether either file is missing
data_existence = [
    (Path(configuration_dict["data_root_dir"]) / file).exists() for file in parsed_filenames
]

# Assign the existence status to each configuration file for error evaluation
# ---- Error evaluation and print message (if applicable)
if not all(data_existence):
    missing_data = parsed_filenames[~np.array(data_existence)]
    raise FileNotFoundError(f"The following data files do not exist: {missing_data}")


# Get the applicable `CONFIG_MAP` keys for the defined datasets
expected_datasets = set(DATASET_DF_MODEL.keys()).intersection(dataset_type)

dataset = "biological"

config_df = pd.DataFrame(configuration_dict[dataset])

datalayer = config_df.columns[0]

def extract_filename_sheetname(root_dir, df):
    if isinstance(df.iloc[0, 0], dict):  # Check if first element is a dict (nested case)
        # Flatten nested dictionaries and extract tuples
        return [
            (root_dir / d['filename'], d['sheetname'])
            for row in df.values
            for d in row
        ]
    else:
        # Handle case where `sheetname` might be a list
        return [
            (root_dir / df.at['filename', col], sheet) if isinstance(df.at['sheetname', col], list)
            else (root_dir / df.at['filename', col], df.at['sheetname', col])
            for col in df.columns
            for sheet in (df.at['sheetname', col] if isinstance(df.at['sheetname', col], list) 
                          else [df.at['sheetname', col]])
        ]

data_root_directory = Path(configuration_dict["data_root_dir"])
dataset = "stratification"
datalayer = "strata"

# Define validation settings from CONFIG_MAP
validation_settings = DATASET_DF_MODEL[dataset][datalayer]

# Define configuration settings w/ file + sheet names
config_settings = configuration_dict[dataset][datalayer]

# Create reference index of the dictionary path
config_map = [dataset, datalayer]

sheetname = config_settings["sh"]
sheet = [sheetname] if isinstance(sheetname, str) else sheetname
dataset_dict = {}
dataset_dict["file_name"] = np.repeat(config_settings["filename"], len(config_settings["sheetname"]))

data_root_directory / config_settings["filename"]
Path(str(dataset_dict["file_name"]) * len(config_settings["sheetname"]))


if isinstance(sheetname, list):
    dataset_dict["file_name"].append(dataset_dict["file_name"][0])
    


if dataset == "biological" 

configuration_dict[dataset][datalayer]
configuration_dict["biological"]["specimen"]

filename, sheetname = extract_filename_sheetname(data_root_directory, )


config_df = pd.DataFrame(configuration_dict).filter(expected_datasets).dropna(how="all")

dataset = "biological"
root_dir = Path(configuration_dict["data_root_dir"])

for dataset in config_df.columns.to_list():

    dataset_df = pd.DataFrame(config_df[dataset].dropna(how="all").to_dict())
    extract_filename_sheetname(root_dir, dataset_df)
    

[pd.DataFrame(config_df[dataset]) for dataset in ]

{key: extract_filename_sheetname(pd.DataFrame(configuration_dict[dataset])) for key in config_df.columns}

file_sheet_pairs = {key: extract_filename_sheetname(pd.DataFrame(configuration_dict[dataset])) for key in config_df.columns.to_list()}
{key: extract_filename_sheetname(pd.DataFrame(configuration_dict[dataset])) for key in config_df.columns.to_list()}
# Define validation settings from CONFIG_MAP
validation_settings = DATASET_DF_MODEL[dataset][datalayer]

# Create paired filename-sheetname tuples

# Define configuration settings w/ file + sheet names
config_settings = configuration_dict[dataset][datalayer]