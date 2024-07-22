import yaml
import re

from pathlib import Path
from typing import Union, Tuple

import pandas as pd
import xarray as xr
import numpy as np

from .livecore import(
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

# TODO: Documentation
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
    
# TODO: Documentation
def read_biology_csv(file: Path, pattern: re.Pattern, config_settings: dict):

    # Read in the `*.csv` file
    df = pd.read_csv(file, usecols=list(config_settings["dtypes"].keys()))

    # Validate the dataframe
    # ---- Check for any missing columns
    missing_columns = (
        [key for key in config_settings["dtypes"].keys() if key not in df.columns]
    )
    # ---- Raise Error, if needed
    if missing_columns: 
        raise ValueError(
            f"The following columns are missing from [{file}]: {', '.join(missing_columns)}!"
        )
    # ---- Ensure the correct datatypes
    df_validated = df.astype(config_settings["dtypes"])
    # ---- Replace column names and drop 
    df_validated = df_validated.rename(columns=config_settings["names"])

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

# TODO: Documentation
# TODO: Refactor, break up cyclomatic complexity
def load_biology_data(file_configuration: dict, update_config: bool = True):

    # Get acoustic directory and initialization settings
    # ---- Files
    biology_file_settings = file_configuration["input_directories"]["biological"]
    # ---- General settings
    biology_analysis_settings = file_configuration["biology"]

    # Get the file-specific settings, datatypes, columns, etc.
    # ---- Get defined columns and datatypes from `LIVE_INPUT_FILE_CONFIG_MAP`
    biology_config_map = LIVE_INPUT_FILE_CONFIG_MAP["biology"]
        # ---- Extract the expected file name ID's
    biology_file_ids = biology_file_settings["file_name_formats"]
    # ---- Extract all of the file ids
    biology_config_ids = list(biology_file_ids.keys())
    # ---- Initialize the dictionary that will define this key in the `input` attribute
    biology_output = {f"{key}_df": pd.DataFrame() for key in biology_config_ids}
    # ---- Initialize the SQL dictionary
    sql_biology_output = {f"{key}_df": pd.DataFrame() for key in biology_config_ids}
    
    # Create full filepath
    biology_directory_path = (
        Path(file_configuration["data_root_dir"]) / biology_file_settings["directory"]
    )
    # ---- Directory check
    directory_existence = biology_directory_path.exists()
    # ---- Error evaluation (if applicable)
    if not directory_existence:
        raise FileNotFoundError(
            f"The acoustic data directory [{biology_directory_path}] does not exist."
        )
    # ---- Get the defined file extension
    file_extension = biology_file_settings["extension"]
    # ---- Create Path.glob generator object
    file_path_obj = biology_directory_path.glob(f"*{'.'+file_extension}")
    #---- Create list of `*.csv`` files
    csv_files = list(file_path_obj)
    # ---- Ensure files exist or raise error otherwise
    if len(csv_files) < 1:
        raise FileNotFoundError(
            f"No `*.csv` files found in [{biology_directory_path}]!"
        )
    else: 
        # ---- Create Path to SQL database file
        db_directory = Path(file_configuration["data_root_dir"]) / "database"
        # ---- Create the directory if it does not already exist
        db_directory.mkdir(parents=True, exist_ok=True)
        # ---- Complete path to `biology.db`
        db_file = db_directory / "biology.db"
        # ---- Query the external SQL database to see if the file tracking table exists
        tables = SQL(db_file, "inspect")
        # ---- Create a list of string-formatted Path names
        csv_files_str = [str(file) for file in csv_files]
        # ---- Create DataFrame
        current_files = pd.DataFrame(csv_files_str, columns=["filepath"])
        # ---- Create if it is missing and then advance `csv_files`
        if "files_read" not in tables:
            # ---- Insert into the SQL database file
            _ = SQL(db_file, "insert", table_name="files_read", columns="filepath",
                     dataframe=current_files)        
            # ---- Create empty list for later comparison
            new_files = []
        else:
            # ---- Pull already processed filenames
            previous_files = SQL(db_file, "select", table_name="files_read")
            # ---- Compare against the current filelist 
            new_files = (
                [file for file in csv_files_str if file not in set(previous_files["filepath"])]
            )  
            # ---- Create a DataFrame for the new files
            new_files_df = pd.DataFrame(new_files, columns=["filepath"])
            # ---- Insert into the SQL database file
            _ = SQL(db_file, "insert", table_name="files_read", dataframe=new_files_df) 

    # Iterate through each of the file ids and read in the data 
    for id in list(biology_file_ids.keys()): 
        # ---- Extract the specific config mapping for this tag/id
        sub_config_map = biology_config_map[id]
        # ---- Drop the `{FIELD_ID}` tag identifier
        file_id_format = re.sub(r'\{FILE_ID:([^}]+)\}', r'\1', biology_file_ids[id])
        # ---- Replace all other tags with `*` placeholders
        file_id_format = re.sub(r"\{[^{}]+\}", "*", file_id_format)
        # ---- Create Path object with the generalized format
        subfile_path_obj = biology_directory_path.glob(f"{file_id_format}.{file_extension}")
        # ---- List all files that match this pattern
        subcsv_files_str = [str(file) for file in list(subfile_path_obj)]
        # ---- Filter for only new files
        subset_files = set(subcsv_files_str).intersection(set(new_files))
        # ---- Pull from SQL database, if applicable
        if f"{id}_df" in tables:
            # ---- SELECT
            sql_df = SQL(db_file, "select", table_name=f"{id}_df", columns="*")
            # ---- Concatenate to the dictionary
            sql_biology_output[f"{id}_df"] = pd.concat([biology_output[f"{id}_df"], sql_df])
        # ---- Add data files not stored in SQL database
        if len(subset_files) > 0 or len(subset_files)== 0 and f"{id}_df" not in tables:
            if len(subset_files) > 0:
                file_list = subset_files
            else:
                file_list = subcsv_files_str
            # ---- Create a list of relevant dataframes
            sub_df_lst = [read_biology_csv(Path(file), biology_file_ids[id], sub_config_map) 
                            for file in file_list]
            # ---- Concatenate into a single DataFrame
            sub_df = pd.concat(sub_df_lst, ignore_index=True)
            # ---- Lower-case sex
            if "sex" in sub_df.columns:
                sub_df["sex"] = sub_df["sex"].str.lower()
            # ---- Concatenate to the dictionary DataFrame
            biology_output[f"{id}_df"] = pd.concat([biology_output[f"{id}_df"], sub_df])

    # Get contrasts used for filtering the dataset
    # ---- Species
    species_filter = file_configuration["species"]["number_code"]
    # ---- Trawl partition information
    trawl_filter = biology_analysis_settings["catch"]["partition"]
    # ---- Apply the filter
    filtered_biology_output = {
        key: df[
            (df['species_id'] == species_filter if 'species_id' in df.columns else True) &
            (df['trawl_partition'].str.lower() == trawl_filter if 'trawl_partition' in df.columns else True)
        ]
        for key, df in biology_output.items() if isinstance(df, pd.DataFrame) and not df.empty
    }

    # Update the SQL database
    for table_name, df in filtered_biology_output.items():
        # ---- Update        
        _ = SQL(db_file, "insert", table_name=table_name, columns="*", 
                dataframe=df)
        
    # Combine the two datasets 
    merged_output = {
        key: pd.concat([
            sql_biology_output.get(key, pd.DataFrame()), 
            filtered_biology_output.get(key, pd.DataFrame())
        ]).drop_duplicates().reset_index(drop=True)
        for key in set(sql_biology_output) | set(filtered_biology_output)
    }
    # ---- Return output
    if update_config:
        if file_configuration["database"]["biology"] is None: 
            file_configuration["database"]["biology"] = db_file
        return merged_output, file_configuration
    else:
        return merged_output

# TODO: Expand data validator and limit cases to '*.zarr' (for now)
# TODO: Refactor "extra" components such as the validation steps, xarray-to-dataframe piping, etc.
# TODO: Documentation
def load_acoustic_data(file_configuration: dict, update_config: bool = True) -> Tuple[pd.DataFrame, xr.Dataset]:
    # Get acoustic directory and initialization settings
    # ---- Files
    acoustic_file_settings = file_configuration["input_directories"]["acoustic"]
    # ---- General settings
    acoustic_analysis_settings = file_configuration["acoustics"]
    
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
    # ---- Initialize the dictionary that will define this key in the `input` attribute
    acoustics_output = {"prc_nasc_df": pd.DataFrame(), 
                        "nasc_df": pd.DataFrame()}
    # ---- Initialize the SQL dictionary
    # sql_acoustics_output = {"sv_df": pd.DataFrame()}

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
    # ---- Create Path.glob generator object (the case of a *.zarr file)
    file_path_obj = acoustic_directory_path.glob(f"*{'.'+file_extension}")
    # ---- Find all zarr files
    zarr_files = list(file_path_obj)
    # ---- Ensure files exist or raise error otherwise
    if len(zarr_files) < 1:
        raise FileNotFoundError(
            f"No `*.zarr` files found in [{acoustic_directory_path}]!"
        )
    else:
        # ---- Create Path to SQL database file
        db_directory = Path(file_configuration["data_root_dir"]) / "database"
        # ---- Create the directory if it does not already exist
        db_directory.mkdir(parents=True, exist_ok=True)
        # ---- Complete path to `biology.db`
        db_file = db_directory / "acoustics.db"
        # ---- Query the external SQL database to see if the file tracking table exists
        tables = SQL(db_file, "inspect")
        # ---- Create a list of string-formatted Path names
        zarr_files_str = [str(file) for file in zarr_files]
        # ---- Create DataFrame
        current_files = pd.DataFrame(zarr_files_str, columns=["filepath"])
        # ---- Create if it is missing and then advance `zarr_files`
        if "files_read" not in tables:
            # ---- Insert into the SQL database file
            _ = SQL(db_file, "insert", table_name="files_read", columns="filepath",
                    dataframe=current_files)        
            # ---- Create empty list for later comparison
            new_files = []
        else:
            # ---- Pull already processed filenames
            previous_files = SQL(db_file, "select", table_name="files_read")
            # ---- Compare against the current filelist 
            new_files = (
                [file for file in zarr_files_str if file not in set(previous_files["filepath"])]
            )  
            # ---- Create a DataFrame for the new files
            new_files_df = pd.DataFrame(new_files, columns=["filepath"])
            # ---- Insert into the SQL database file
            _ = SQL(db_file, "insert", table_name="files_read", dataframe=new_files_df) 

    # Find new files that have not yet been processed
    if not new_files: 
        subset_files = zarr_files
    else:
        subset_files = set(zarr_files).intersection(set(new_files))

    # Read in the `*.zarr` file(s)
    # ! [REQUIRES DASK] ---- Read in the listed file
    if len(subset_files) > 1:
        zarr_data_ds = xr.open_mfdataset(subset_files, engine="zarr", chunks="auto", 
                                         data_vars=specified_vars, coords=specified_coords)
    elif len(subset_files) == 1:
        zarr_data_ds = xr.open_dataset(subset_files[0], engine="zarr", chunks="auto")

    # Pre-process the Dataset, convert it to a DataFrame, and validate the structure
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
    zarr_data_df_filtered = zarr_data_df[full_config_map.keys()].astype(full_config_map)

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
    acoustics_output["prc_nasc_df"] = zarr_data_df_output.drop(columns = ["frequency_nominal"])
    # ---- Return output
    if update_config:
        if file_configuration["database"]["acoustics"] is None: 
            file_configuration["database"]["acoustics"] = db_file
        return acoustics_output, file_configuration
    else:
        return acoustics_output