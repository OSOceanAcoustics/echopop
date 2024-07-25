import yaml
import re

from pathlib import Path
from typing import Union, Tuple, Optional, List

import pandas as pd
import xarray as xr
import numpy as np

from .live_core import(
    LIVE_DATA_STRUCTURE,
    LIVE_FILE_FORMAT_MAP,
    LIVE_INPUT_FILE_CONFIG_MAP
)

from .sql_methods import SQL

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
        raise FileNotFoundError(
            f"The following configuration files do not exist: {missing_config}."
            )

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

def validate_data_directory(root_directory: str, file_settings: dict) -> List[Path]:

    # Get acoustic directory and initialization settings
    # ---- Create the full filepath 
    directory_path = Path(root_directory) / file_settings["directory"]
    # ---- Get the defined file extension
    file_extension = file_settings["extension"]

    # Validate filepath, columns, datatypes
    # ---- Error evaluation (if applicable)
    if not directory_path.exists():
        raise FileNotFoundError(
            f"The acoustic data directory [{directory_path}] does not exist."
        )
    
    # Validate that files even exist
    # ---- List available *.zarr files
    data_files = list(directory_path.glob(f"*{'.'+file_extension}"))
    # ---- Error evaluation (if applicable)
    if not data_files:
        raise FileNotFoundError(
            f"No `*.{file_extension}` files found in [{directory_path}]!"
        )
    
    # Return the output
    return data_files

def query_processed_files(root_directory: str, file_settings: dict, files: List[Path]) -> dict:

    # Get the database name 
    db_name = file_settings["database_name"]

    # Create filepath to the SQL database
    # ---- Create Path to SQL database file
    db_directory = Path(root_directory) / "database"
    # ---- Create the directory if it does not already exist
    db_directory.mkdir(parents=True, exist_ok=True)
    # ---- Complete path to the database file
    db_file = db_directory / db_name

    # Create a list of string-formatted Path names
    files_str = [str(file) for file in files]
    # ---- Create DataFrame
    current_files = pd.DataFrame(files_str, columns=["filepath"])

    # Check for the table `files_read`
    files_read_tbl = SQL(db_file, "validate", table_name="files_read")

    # Validate whether the table exists; if not, create the table and then insert
    if not files_read_tbl:
        # ---- Create table
        SQL(db_file, "create", table_name="files_read", dataframe=current_files, 
            primary_keys = ["filepath"])
        # ---- Populate table
        SQL(db_file, "insert", table_name="files_read", dataframe=current_files)
        # ---- Break early
        return files_str, db_file
    
    # Query already existing files
    previous_files = SQL(db_file, "select", table_name="files_read", output_type=str)
    # ---- Insert file list
    SQL(db_file, "insert", table_name="files_read", dataframe=current_files, id_columns="filepath")

    # Filter out previously processed files
    # ---- Apply filter by comparing sets and return the output
    return list(set(files_str) - set(previous_files)), db_file

def sql_data_exchange(database_file: Path, **kwargs):

    # Check whether the `table_name` table exists
    table_exists = SQL(database_file, "validate", **kwargs)

    # If empty and table does not exist
    if kwargs["dataframe"].empty and table_exists:
        return SQL(database_file, "select", **kwargs)

    # Create table if it does not exist and run the initial insertion
    if not table_exists:
        # ---- Create table
        SQL(database_file, "create", **kwargs)
        # ---- Ignore the `id_columns` argument, if present
        try:
            del kwargs["id_columns"]
        except KeyError:
            pass
        # ---- Insert into table        
        SQL(database_file, "insert", **kwargs)
        # ---- Return the initial dataframe
        return kwargs.get("dataframe")
    
    # Insert into the table
    SQL(database_file, "insert", **kwargs)
    
    # Select existing data frame the database and return the output
    return SQL(database_file, "select", **kwargs)

def read_acoustic_zarr(acoustic_files: Path) -> tuple:

    # Get the file-specific settings, datatypes, columns, etc.
    # ---- Get defined columns and datatypes from `LIVE_INPUT_FILE_CONFIG_MAP`
    acoustics_config_map = LIVE_INPUT_FILE_CONFIG_MAP["acoustics"]
    # ---- Create list of coordinate data variables
    specified_vars = list(acoustics_config_map["xarray_variables"].keys())
    # ---- Create set of coordinate variables
    specified_coords = list(acoustics_config_map["xarray_coordinates"].keys())      
    # ---- Concatenate into a full configuration map
    full_config_map = {**acoustics_config_map["xarray_coordinates"],
                        **acoustics_config_map["xarray_variables"]} 
    
    # Determine the file loading method for the `acoustic_files`
    if len(acoustic_files) > 1:
        zarr_data_ds = xr.open_mfdataset(acoustic_files, engine="zarr", chunks="auto", 
                                         data_vars=specified_vars, coords=specified_coords)
    else:
        zarr_data_ds = xr.open_dataset(acoustic_files[0], engine="zarr", chunks="auto")

    # Pre-process the Dataset, convert it to a DataFrame, and validate the structure
    # ---- Convert to a DataFrame
    zarr_data_df = zarr_data_ds.to_dataframe().reset_index()
    # ---- Check for any missing columns
    missing_columns = (
        [key for key in full_config_map.keys() if key not in zarr_data_df.columns]
    )
    # ---- Raise Error, if needed
    if missing_columns: 
        raise ValueError(
            f"The following columns are missing from at least one file: in "
            f"{', '.join(missing_columns)}!"
        )
    # ---- Select defined columns
    zarr_data_df_filtered = zarr_data_df[full_config_map.keys()].astype(full_config_map)

    # Gather some of the units
    data_units = {
        "longitude": zarr_data_ds.longitude.units,
        "latitude": zarr_data_ds.latitude.units,
        "frequency": zarr_data_ds.frequency_nominal.units,
    }

    # Return a Tuple
    return zarr_data_df_filtered, data_units

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

def preprocess_acoustic_data(prc_nasc_df: pd.DataFrame,
                             file_configuration: dict) -> pd.DataFrame:

    # Get acoustic processing settings
    acoustic_analysis_settings = file_configuration["acoustics"]
    # ---- Extract the fined acoustic frequency
    transmit_settings = acoustic_analysis_settings["transmit"]

    # Filter the dataset
    # ---- Configure `frequency_nominal`, if necessary
    prc_nasc_df["frequency_nominal"] = (
        configure_transmit_frequency(prc_nasc_df["frequency_nominal"],
                                     transmit_settings,
                                     acoustic_analysis_settings["dataset_units"]["frequency"])
    )
    # ---- Filter out any unused frequency coordinates
    prc_nasc_df_filtered = (
        prc_nasc_df[prc_nasc_df["frequency_nominal"] == transmit_settings["frequency"]]
    )

    # Remaining adjustments to the acoustic data prior to being passed to the `LiveSurvey` object
    # ---- Replace NASC `NaN` values with `0.0`
    prc_nasc_df_filtered.loc[:, "NASC"] = prc_nasc_df_filtered.loc[:, "NASC"].fillna(0.0)
    # ---- Drop the `frequency_nominal` column and return the output 
    return prc_nasc_df_filtered.drop(columns = ["frequency_nominal"])

def load_acoustic_data(file_configuration: dict) -> Tuple[pd.DataFrame]:

    # Get the acoustic file settings and root directory
    # ---- File settings
    file_settings = file_configuration["input_directories"]["acoustics"]
    # ---- Root directory
    root_directory = file_configuration["data_root_dir"]
    
    # Get and validate the acoustic data directory and files
    acoustic_files = validate_data_directory(root_directory, file_settings)

    # Query `acoustics.db` to process only new files (or create the db file in the first place)
    new_acoustic_files, file_configuration["database"]["acoustics"] = (
        query_processed_files(root_directory, file_settings, acoustic_files)  
    )  

    # Read in the acoustic data files
    if new_acoustic_files:
        # ! [REQUIRES DASK] ---- Read in the listed file
        prc_nasc_df, acoustic_data_units = read_acoustic_zarr(new_acoustic_files)
        # ---- Add the `acoustic_data_units` to the dictionary
        file_configuration["acoustics"]["dataset_units"] = acoustic_data_units
        # ---- Preprocess the acoustic dataset
        prc_nasc_df_processed = preprocess_acoustic_data(prc_nasc_df, file_configuration)
        # ---- Return output
        return prc_nasc_df_processed
    else:
        return None
    
def filter_filenames(directory_path: Path, filename_id: str, 
                     files: List[Path],
                     file_extension: str):

    # Drop the `{FIELD_ID}` tag identifier
    file_id_format = re.sub(r'\{FILE_ID:([^}]+)\}', r'\1', filename_id)
    # ---- Replace all other tags with `*` placeholders
    file_id_format = re.sub(r"\{[^{}]+\}", "*", file_id_format)
    # ---- Create Path object with the generalized format
    subfile_path_obj = directory_path.glob(f"{file_id_format}.{file_extension}")
    # ---- List all files that match this pattern
    subfile_str = [str(file) for file in list(subfile_path_obj)]

    # Convert list of proposed files from Path to String
    file_str = [str(file) for file in list(files)]
    
    # Find intersection with the proposed filenames and return the output
    return list(set(subfile_str).intersection(set(file_str)))

def compile_filename_format(file_name_format: str):

    # Create a copy of `file_name_format`
    regex_pattern = file_name_format
    
    # Iterate through the keys from `LIVE_FILE_FORMAT_MAP` to format a regex pattern
    for key, value in LIVE_FILE_FORMAT_MAP.items():
        regex_pattern = regex_pattern.replace(f"{{{key}}}", value["expression"])
    # ---- Replace the `FILE_ID` tag
    regex_pattern = re.sub(r'\{FILE_ID:(.+?)\}', r'(?P<FILE_ID>\1)', regex_pattern)

    # Compile the regex pattern and return the output
    return re.compile(regex_pattern)

def read_biology_csv(file: Path, pattern: re.Pattern, config_map: dict):

    # Read in the `*.csv` file
    df = pd.read_csv(file, usecols=list(config_map["dtypes"].keys()))

    # Validate the dataframe
    # ---- Check for any missing columns
    missing_columns = (
        [key for key in config_map["dtypes"].keys() if key not in df.columns]
    )
    # ---- Raise Error, if needed
    if missing_columns: 
        raise ValueError(
            f"The following columns are missing from [{file}]: {', '.join(missing_columns)}!"
        )
    # ---- Ensure the correct datatypes
    df_validated = df.astype(config_map["dtypes"])
    # ---- Replace column names and drop 
    df_validated = df_validated.rename(columns=config_map["names"])

    # Get the substring components that can be added to the DataFrame
    filename_substrings = re.findall(r'\{([^:}]+)(?::[^}]+)?}', pattern)
    # ---- Create sub-list of columns that can be added to the DataFrame
    valid_tags = list(set(["HAUL", "SPECIES_CODE"]).intersection(set(filename_substrings)))

    # Compile the filename regular expression
    compiled_regex = compile_filename_format(pattern)
    # ---- Create the `Match` object that will be used to parse the string
    match_obj = compiled_regex.search(file.name)

    # Iterate through the filename-derived tags and add them to the DataFrame
    for i in valid_tags: 
        matched_key = LIVE_FILE_FORMAT_MAP[i]
        df_validated[matched_key["name"]] = matched_key["dtype"](match_obj.group(i))

    # Return the resulting DataFrame
    return df_validated

def preprocess_biology_data(biology_output: dict, file_configuration: dict):
  
    # Get SQL database file
    biology_db = file_configuration["database"]["biology"]
    
    # Get contrasts used for filtering the dataset
    # ---- Species
    species_filter = file_configuration["species"]["number_code"]
    # ---- Trawl partition information
    trawl_filter = file_configuration["biology"]["catch"]["partition"]
    # ---- Create filter dictionary
    filter_dict = dict(species_id=species_filter, trawl_partition=trawl_filter)

    # Apply the filter
    filtered_biology_output = {
        key: biology_data_filter(df, filter_dict) 
        for key, df in biology_output.items() if isinstance(df, pd.DataFrame) and not df.empty
    }
    # ---- Swap this out if no new files are present
    if not filtered_biology_output:
        # ---- Get available tables
        table_list = list(set(SQL(biology_db, "map")) - set(["files_read"]))
        # ---- Plug into the dictionary
        filtered_biology_output.update({key: pd.DataFrame() for key in table_list})
    # ---- Initialize the results dictionary   
    results_dict = {key: pd.DataFrame() for key in filtered_biology_output.keys()}

    # Update the SQL database
    for table_name, df in filtered_biology_output.items():
        # ---- Get identifier columns
        key_columns = get_table_key_names(biology_db, filtered_biology_output, table_name)
        # ---- Create copy
        df = df.copy()
        # ---- Add an autoincrementing tag that will serve as a primary key and unique constraint
        df.loc[:, "id"] = "row" + df.index.astype(str) + "-" + "-".join(key_columns)        
        # ---- Insert the new data into the database & pull in the combined dataset
        table_df = sql_data_exchange(biology_db, 
                                     dataframe=df, 
                                     table_name=table_name, 
                                     id_columns=["id"],
                                     primary_keys=["id"],
                                     output_type=pd.DataFrame)
        # ---- Add to the outgoing dictionary (and drop SQL db identifier)
        results_dict.update({table_name: table_df.drop(columns="id")})
    
    # Return the output
    return results_dict

def get_table_key_names(db_file: Path, data_dict: dict, table_name: str) -> List[str]:

    # Get the data input column names
    if data_dict[table_name].empty:
        # ---- Inspect the table
        inspected_table = SQL(db_file, "inspect", table_name=table_name)
        # ---- Create a list of the data columns
        table_columns = list(inspected_table.keys())
    else:
        # ---- Get the DataFrame column names
        table_columns = data_dict[table_name].columns

    # Create a list of the primary keys
    key_columns = (
           set(table_columns)
           .intersection(["trawl_partition", "sex", "haul_num", "species_id", "longitude", 
                          "latitude"]) 
        )

    # Return a list of the output
    return list(key_columns)

def load_biology_data(file_configuration: dict):

    # Get the acoustic file settings and root directory
    # ---- File settings
    file_settings = file_configuration["input_directories"]["biology"]
    # ---- Root directory
    root_directory = file_configuration["data_root_dir"]

    # Get and validate the acoustic data directory and files
    biology_files = validate_data_directory(root_directory, file_settings)

    # Query `biology.db` to process only new files (or create the db file in the first place)
    # SQL(biology_db, "drop", table_name="files_read")
    new_biology_files, file_configuration["database"]["biology"] = (
        query_processed_files(root_directory, file_settings, biology_files)
    )

    # Get the file-specific settings, datatypes, columns, etc.
    # ---- Get defined columns and datatypes from `LIVE_INPUT_FILE_CONFIG_MAP`
    biology_config_map = LIVE_INPUT_FILE_CONFIG_MAP["biology"]
    # ---- Extract the expected file name ID's
    biology_file_ids = file_settings["file_name_formats"]
    # ---- Extract all of the file ids
    biology_config_ids = list(biology_file_ids.keys())
    # ---- Initialize the dictionary that will define this key in the `input` attribute
    biology_output = {f"{key}_df": pd.DataFrame() for key in biology_config_ids}
    # ---- Create filepath object
    directory_path = Path(file_configuration["data_root_dir"]) / file_settings["directory"]
    
    # Add SQL file to dict
    file_configuration["database"]["biology"] = (
        Path(file_configuration["data_root_dir"]) / "database" / file_settings["database_name"]         
    )

    # Iterate through the different biology datasets and read them in
    for dataset in list(biology_file_ids.keys()):
        # ---- Get dataset-specific file lists
        dataset_files = filter_filenames(directory_path, 
                                         file_settings["file_name_formats"][dataset], 
                                         new_biology_files, 
                                         file_settings["extension"])
        # ---- If there are dataset files available
        if dataset_files:
            # ---- Read in validated biology data
            dataframe_list = [read_biology_csv(Path(file), 
                                               file_settings["file_name_formats"][dataset], 
                                               biology_config_map[dataset]) 
                              for file in dataset_files]
            # ---- Concatenate the dataset
            dataframe_combined = pd.concat(dataframe_list, ignore_index=True)
            # ---- Lower-case sex
            if "sex" in dataframe_combined.columns: 
                dataframe_combined["sex"] = dataframe_combined["sex"].str.lower()
            # ---- Lower-case trawl partition type
            if "trawl_partition" in dataframe_combined.columns: 
                dataframe_combined["trawl_partition"] = dataframe_combined["trawl_partition"].str.lower()
            # ---- Reformat datetime column
            if "datetime" in dataframe_combined.columns:
                dataframe_combined["datetime"] = convert_datetime(dataframe_combined["datetime"])
            # ---- Add to the data dictionary
            biology_output[f"{dataset}_df"] = dataframe_combined

    # Pre-process and return the results
    return preprocess_biology_data(biology_output, file_configuration)

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

def validate_spatial_config(spatial_config: dict):

    # Check the link method
    # ---- Extract string-formatted method name
    link_method = spatial_config["link_biology_acoustics"].lower()
    # ---- Validate
    if link_method not in SPATIAL_CONFIG_MAP.keys():
        raise ValueError(
            f"Unexpected biology-acoustic linking parameter ([{link_method}]). Valid options "
            f"include: 'global', 'closest_haul', 'weighted_haul', 'griddify', and 'INPFC'."
        )
    
    # Verify that associated parameters are present in the configuration settings
    # ---- Get keys as a list
    config_keys = list(spatial_config.keys())
    # ---- Check for specific methods
    if link_method not in config_keys and link_method != "global":
        raise ValueError(
            f"No parameters provided for the biology-acoustic linking ([{link_method}])."
        )
    
    # Check key settings
    if link_method == "griddify": 
        validate_griddify_config(spatial_config, link_method)
    elif link_method == "inpfc": 
        validate_inpfc_config(spatial_config, link_method)
    elif link_method != "global": 
        validate_hauls_config(spatial_config, link_method)

def validate_hauls_config(spatial_config: dict, link_method: str):

    # Get the link method configuration map
    link_method_settings = SPATIAL_CONFIG_MAP[link_method]

    # Extract the defined settings
    input_method_settings = spatial_config[link_method]

    # Check for `proximity` 
    if "proximity" not in input_method_settings.keys():
        raise KeyError(
            "The following parameters are missing from the biology-acoustic linking method: "
            "'proximity'!"
        )
    
    # Evaluate valid options for `proximity`
    if input_method_settings["proximity"] not in link_method_settings["proximity"]["choices"]:
        raise KeyError(
            f"Value biology-acoustic linking method parameter `proximity` must be one of the : "
            f"following: {link_method_settings["proximity"]["choices"]}."
        )       
    
def validate_griddify_config(spatial_config: dict, link_method: str):

    # Get the link method configuration map
    link_method_settings = SPATIAL_CONFIG_MAP[link_method]

    # Extract the defined settings
    input_method_settings = spatial_config[link_method]

    # Check for the required keys
    key_diff = set(input_method_settings.keys()).difference(set(link_method_settings.keys()))
    # ---- Raise Error
    if key_diff:
        raise KeyError(
            f"The following parameters are missing from the biology-acoustic linking method: "
            f"{list(key_diff)}!"
        )    
    
    # Iterate through the keys to evaluate inputs
    for key in list(input_method_settings.keys()):
        # ---- Subset the input method config
        input = input_method_settings[key]
        # ---- Get the original config of the dtypes
        model = link_method_settings[key]
        # ---- Compare entries
        parameter_diff = set(input.keys()).difference(set(model.keys()))
        # ---- Raise Error
        if parameter_diff:
            raise KeyError(
                f"Unexpected parameter(s) ('{parameter_diff}') detected in '{link_method}' "
                f"configuration."
            )    
        # ---- Check if the appropriate coordinate pairs are present
        coordinate_pairs = [set(param).intersection(set(input.keys())) for param in model["pairs"]]
        # ---- Count the number of paired coordinates
        pair_counts = [len(param) for param in coordinate_pairs]
        # ---- If there are multiple pairs
        if (np.array(pair_counts) == 2).sum() != 1:
            raise ValueError(
                f"A single coordinate-pair is allowed (and required) within the '{key}' parameter "
                f"for the link method '{link_method}' defined via the following options: "
                f"{model["pairs"]}."
            )
        # ---- Check the datatypes
        for parameter in input.keys():
            # ---- Get the datatypes
            config_dtypes = model[parameter]["types"]
            # ---- Get input parameter
            input_parameter = input[parameter]
            # ---- If List
            if isinstance(config_dtypes, list):
                if not isinstance(input_parameter, list):
                    raise TypeError(
                        f"Biology-acoustic linking method argument '{parameter}' within '{key}' "
                        f"for method '{link_method}' must be contained within a list."
                    )
            else:
                input_parameter = [input_parameter]
                config_dtypes = [config_dtypes]
            # ---- Check correct datatypes
            if not np.all([type(value) in config_dtypes for value in input_parameter]):
                raise TypeError(
                    f"Biology-acoustic linking method argument '{parameter}' within '{key}' "
                    f"for method '{link_method}' must be one of the following types within a list: "
                    f"{config_dtypes}."
                )    

def validate_inpfc_config(spatial_config: dict, link_method: str):

    # Get the link method configuration map
    link_method_settings = SPATIAL_CONFIG_MAP[link_method]

    # Extract the defined settings
    input_method_settings = spatial_config[link_method]

    # Check for the required keys
    key_diff = set(input_method_settings.keys()).difference(set(link_method_settings.keys()))
    # ---- Raise Error
    if key_diff:
        raise KeyError(
            f"The following parameters are missing from the biology-acoustic linking method: "
            f"{list(key_diff)}!"
        )
    
    # Iterate through the keys to evaluate inputs
    for key in list(input_method_settings.keys()):
        # ---- Subset the input method config
        input = input_method_settings[key]
        # ---- Get the original config of the dtypes
        model = link_method_settings[key]["types"]
        # ---- Evaluate if a list 
        if not isinstance(input, list):
            raise TypeError(
                f"Biology-acoustic linking method argument '{key}' for method '{link_method}' must "
                f"be contained within a list."
            )
        # ---- Evaluate if it is a type within the list
        if not type(input[0]) in model:
            raise TypeError(
                f"Biology-acoustic linking method argument '{key}' for method '{link_method}' must "
                f"be one of the following types within a list: {model}."
            )    
        
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

