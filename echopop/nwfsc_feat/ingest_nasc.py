from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Union

import pandas as pd
import numpy as np
import re
import functools
from ..ingest import read_csv_file
from ..core.echoview import ECHOVIEW_TO_ECHOPOP, ECHOVIEW_DATABASE_EXPORT_FILESET, ECHOVIEW_EXPORT_ROW_SORT

def map_transect_num(
    ev_export_paths: Dict[str, Generator], transect_pattern: str = r"T(\d+)"
) -> pd.DataFrame:
    """
    Extract transect numbers from Echoview export filenames without using loops.

    Parameters
    ----------
    ev_export_paths : Dict[str, Generator]
        Dictionary with keys "analysis", "cells", "intervals", "layers", with values being 
        `pathlib.Path.glob` generators pointing toward Echoview export *.csv files.
    transect_pattern : str
        Regex pattern to extract transect numbers from file names.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns for file_type, file_path, and transect_num.
    """
    # Compile the transect pattern regex
    compiled_pattern = re.compile(transect_pattern)
    
    # Extract all file info using a list comprehension and flatten the result
    records = [
        {"file_type": file_type, "file_path": path, "transect_num": float(compiled_pattern.search(str(path)).group(1))}
        for file_type, paths in ev_export_paths.items()
        for path in list(paths)
        if compiled_pattern.search(str(path))
    ]
    
    # Convert to DataFrame
    return pd.DataFrame(records)

def validate_transect_exports(transect_files_df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate that each transect number has all required export file types without using loops.
    
    For each transect number, there should be exactly four files:
    one each for "analysis", "cells", "intervals", and "layers".
    
    Parameters
    ----------
    transect_files_df : pd.DataFrame
        DataFrame with columns "file_type", "file_path", and "transect_num".
        
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only valid transect numbers that have all required file types.
    """
    # Create a pivot table to count file types by transect number
    type_count = pd.crosstab(
        transect_files_df["transect_num"], 
        transect_files_df["file_type"]
    )
    
    # Check that each transect has all required file types
    valid_mask = (type_count[list(ECHOVIEW_DATABASE_EXPORT_FILESET)] == 1).all(axis=1)
    valid_transects = type_count.index[valid_mask].tolist()
    
    # Filter the original DataFrame to include only valid transect numbers
    filtered_df = transect_files_df[transect_files_df["transect_num"].isin(valid_transects)].copy()
    
    # Log the removed transect numbers
    all_transects = transect_files_df["transect_num"].unique()
    removed_transects = np.setdiff1d(all_transects, valid_transects)
    if len(removed_transects) > 0:
        print(f"Removed transect numbers with incomplete file sets: {sorted(removed_transects)}")
        
    return filtered_df

def clean_echoview_cells_df(
    cells_df: pd.DataFrame,
    inplace: bool = False
) -> Optional[pd.DataFrame]:
    """
    Clean the various region class and name entries

    Parameters
    ----------
    cells_df : pd.DataFrame
        Cells export DataFrame
    inplace : bool, default False
        If True, perform operation in-place and return None

    Returns
    -------
    pd.DataFrame or None
        DataFrame with cleaned strings for the region class and name strings if inplace=False
        None if inplace=True
    """
    
    # Work on the original or a copy based on inplace parameter
    df = cells_df if inplace else cells_df.copy()    

    # Check column overlap
    region_cols = list(set({"region_class", "region_name"}).intersection(df.columns))

    if not region_cols:
        # No relevant columns found, return unchanged
        return None if inplace else df
    
    # Adjust "region_class" and/or "region_name" if needed
    df[region_cols] = (
        df[region_cols]
        .replace('"', "", regex=True)
        .replace(r"^\s*(.*?)\s*$", r"\1", regex=True)
    )

    # Set to lowercase
    if "region_class" in df.columns:
        df["region_class"] = df["region_class"].str.lower()

    # Return None for inplace=True, df otherwise
    return None if inplace else df


def impute_bad_coordinates(data: pd.DataFrame, column: str) -> None:
    """
    Impute bad or missing coordinates in a DataFrame.
    
    Parameters
    ----------
    data: pandas.DataFrame
        DataFrame containing georeferenced data with latitude and longitude coordinates
    column: str
        Name of the column containing coordinates to impute
        
    Returns
    -------
    None
        The function modifies the pandas.DataFrame in-place
    """
    # Find indices where values are "999.0" or `NaN`
    invalid_coords = data.index[(data[column] == 999.0) | (data[column].isna())].to_numpy()
    
    # Break from helper if there are no invalid coordinates
    if len(invalid_coords) == 0:
        return
        
    # Evaluate cases where there are multiple steps between indices (i.e. multiple groups)
    delta_invalid_coords = np.where(np.diff(invalid_coords) > 1)[0]
    
    # Digitize into discrete groups of "bad" data
    idx_groups = np.digitize(invalid_coords, invalid_coords[delta_invalid_coords], right=True) \
                 if len(delta_invalid_coords) > 0 else np.zeros(len(invalid_coords), dtype=int)
    
    # Iterate across the bad data groups
    for group in np.unique(idx_groups):
        # Get the matching group index
        in_group_idx = np.where(idx_groups == group)[0]
        # Parse the index
        in_group = invalid_coords[in_group_idx]
        
        # Case 1: Single or consecutive bad coordinates at head of transect
        if any(in_group == 0):
            # Rolling imputation
            for i in range(in_group.size - 1, -1, -1):
                data.loc[i, column] = 2 * data.loc[i + 1, column] - data.loc[i + 2, column]
                
        # Case 2: Single or consecutive bad coordinates at tail of transect
        elif any(in_group + 2 > len(data)):
            # Get starting index
            start_index = in_group.min()
            # Get ending index
            end_index = in_group.max()
            # Rolling imputation
            for i in range(start_index, end_index + 1, 1):
                data.loc[i, column] = 2 * data.loc[i - 1, column] - data.loc[i - 2, column]
                
        # Case 3: Single or consecutive bad coordinates in the middle of transect
        else:
            # Get starting index
            start_index = in_group.min() - 1
            # Get ending index
            end_index = in_group.max() + 1
            # Calculate the mean difference between "good/valid" coordinates
            step_mean = (data.loc[end_index, column] - data.loc[start_index, column]) / (
                in_group.size + 1
            )
            # Impute over "bad" points
            data.loc[in_group, column] = data.loc[start_index, column] + step_mean * np.arange(
                1, in_group.size + 1, 1
            )

def read_echoview_export(filename: Path, 
                         validator: Optional[Any] = None) -> pd.DataFrame:
    """
    Generic reader for Echoview export CSVs. Used for files like analysis, cells, layers, 
    intervals.

    Parameters
    ----------
    filename : pathlib.Path
        Full path to the NASC CSV file.
    validator : Any
        File-specific validator, if defined.

    Returns
    -------
    pd.DataFrame
        Cleaned and formatted data.
    """

    # Read the CSV file
    df = read_csv_file(filename)

    # Rename columns used by Echopop
    df.rename(columns=ECHOVIEW_TO_ECHOPOP, inplace=True)

    # TODO: Validation step would be here

    return df

def sort_echoview_export_df(
    export_df: pd.DataFrame,
    inplace: bool = False
) -> pd.DataFrame:
    """
    Sort the file rows and reset indices 
    
    Parameters
    ----------
    export_df : pd.DataFrame
        Echoview export DataFrame that can be either "analysis", "cells", "intervals", or "layers". 
    inplace : bool, default False
        If True, perform operation in-place and return None
        
    Returns
    -------
    pd.DataFrame or None
        Sorted and reindexed DataFrame if inplace=False, None otherwise 
    """

    # Work on the original or a copy based on inplace parameter
    df = export_df if inplace else export_df.copy()

    # Check columns against the appropriate sorting columns
    sort_cols = list(set(ECHOVIEW_EXPORT_ROW_SORT).intersection(df.columns))

    # Sort the columns
    df.sort_values(sort_cols, inplace=True)

    # Reindex
    df.reset_index(drop=True, inplace=True)
    
    # Return None for inplace=True, df otherwise
    return None if inplace else df

def update_transect_spacing(
    transect_data: pd.DataFrame, 
    default_transect_spacing: float,
    latitude_threshold: float = 60.0,
    inplace: bool = False
) -> Optional[pd.DataFrame]:
    """
    Calculate and update the maximum spacing between transects.
    
    Parameters
    ----------
    transect_data : pd.DataFrame
        DataFrame containing transect data with 'transect_num', 'latitude' columns
    default_transect_spacing : float
        Default spacing to use (in nautical miles)
    latitude_threshold : float, default 60.0
        Maximum latitude to consider for spacing calculations
    inplace : bool, default False
        If True, modify the input DataFrame in-place and return None
    
    Returns
    -------
    pd.DataFrame or None
        Updated DataFrame if inplace=False, None otherwise
    """
    # Work on a copy or the original based on inplace parameter
    df = transect_data if inplace else transect_data.copy()
    
    # Get unique transect numbers
    transect_number = np.unique(df["transect_num"])
    
    # Initialize max transect spacing column
    df["transect_spacing"] = default_transect_spacing
    
    # Iterate through the transects to determine the maximum spacing
    for i in range(len(transect_number)):
        if i >= 2:
            # ---- For 2 transects prior to the current transect
            lag_2_index = df.index[
                (df["transect_num"] == transect_number[i - 2])
                & (df["latitude"] < latitude_threshold)
            ]
            # ---- For 1 transect prior to the current transect
            lag_1_index = df.index[
                (df["transect_num"] == transect_number[i - 1])
            ]
            # ---- Current transect
            current_index = df.index[
                (df["transect_num"] == transect_number[i])
                & (df["latitude"] < latitude_threshold)
            ]
            
            # Check if we have data for all three transects
            if len(lag_2_index) > 0 and len(lag_1_index) > 0 and len(current_index) > 0:
                # ---- Calculate the mean transect latitude (lag-2)
                lag_2_latitude = df.loc[lag_2_index, "latitude"].mean()
                # ---- Calculate the mean transect latitude (current)
                current_latitude = df.loc[current_index, "latitude"].mean()
                # ---- Compute the difference in the latitudes of adjacent transects
                delta_latitude = np.abs(current_latitude - lag_2_latitude)
                # ---- Get latitude range for the lag-2 transect
                latitude_2_range = (
                    df.loc[lag_2_index, "latitude"].max()
                    - df.loc[lag_2_index, "latitude"].min()
                )
                # ---- Get latitude range for current transect
                latitude_range = (
                    df.loc[current_index, "latitude"].max()
                    - df.loc[current_index, "latitude"].min()
                )
                # ---- Assign maximum spacing
                if (
                    (delta_latitude <= 2.0 * default_transect_spacing * 1.1 / 30.0)
                    & (latitude_2_range < 1 / 6)
                    & (latitude_range < 1 / 6)
                ):
                    df.loc[lag_1_index, "transect_spacing"] = delta_latitude * 30.0
    
    # Return the updated dataframe or None if inplace
    return None if inplace else df    

def read_echoview_nasc(filename: Path,
                       transect_num: float,
                       validator: Optional[Any] = None) -> pd.DataFrame:
    """
    Generic reader for Echoview export CSVs. Used for files like analysis, cells, layers, 
    intervals.

    Parameters
    ----------
    filename : pathlib.Path
        Full path to the NASC CSV file.
    transect_num : float
        Transect number to use for filtering or labeling.

    Returns
    -------
    pd.DataFrame
        Cleaned and formatted DataFrame
    """

    # Read in the defined CSV file
    nasc_df = read_echoview_export(filename, validator)
    
    # Add transect number
    nasc_df["transect_num"] = transect_num
    
    # Fix latitude and longitude
    # ---- Latitude
    if "latitude" in nasc_df.columns:
        impute_bad_coordinates(nasc_df, "latitude")
    # ---- Longitude
    if "longitude" in nasc_df.columns:
        impute_bad_coordinates(nasc_df, "longitude")

    # Return the cleaned DataFrame
    return nasc_df


# def read_echoview_export(filename: str, transect_number: int, validator: MetaModel, 
#                          column_mapping: Optional[dict] = None) -> pd.DataFrame:
#     """
#     Read Echoview export file, validate it, and prepare it for further processing.
    
#     Parameters
#     ----------
#     filename: str
#         Export file path
#     transect_number: int
#         Transect number associated with the export file
#     validator: MetaModel
#         A ``pandera`` ``DataFrameModel``-class validator
#     column_mapping: Optional[dict]
#         Dictionary mapping Echoview column names to EchoPop column names
        
#     Returns
#     -------
#     pd.DataFrame
#         Validated and processed DataFrame
        
#     See Also
#     --------
#     echopop.load.ingest_echoview_exports
#         A more detailed description
#     """
#     # Read and prepare the CSV file
#     export_file = read_csv_file(filename)
    
#     # Validate
#     export_valid = validator.validate_df(export_file, filename)
    
#     # Assign transect number column
#     export_valid["transect_num"] = transect_number
    
#     # Assign the names if column_mapping is provided
#     if column_mapping:
#         export_valid.rename(columns=column_mapping, inplace=True)
    
#     # Get coordinate column names, if present
#     lon_col = next((col for col in export_valid.columns if "lon" in col.lower()), None)
#     lat_col = next((col for col in export_valid.columns if "lat" in col.lower()), None)
    
#     # Impute latitude
#     if lat_col is not None:
#         impute_bad_coordinates(export_valid, lat_col)
        
#     # Impute longitude
#     if lon_col is not None:
#         impute_bad_coordinates(export_valid, lon_col)
    
#     return export_valid


def echoview_nasc_to_df(
    filtered_df: pd.DataFrame 
) -> list[pd.DataFrame]:
    """
    Reads and returns Echoview NASC dataframes for each file in the input DataFrame.
    
    Parameters
    ----------
    filtered_df : pd.DataFrame
        DataFrame with columns "file_path" and "transect_num".
        
    Returns
    -------
    list[pd.DataFrame]
        List of parsed and validated DataFrames for each file.
    """
    return [read_echoview_nasc(row["file_path"], row["transect_num"]) 
            for _, row in filtered_df.iterrows()]

# # based on the current generate_dataframes
# def echoview_nasc_to_df(
#     files: list[Path], transect_num: list[int]
# ) -> Generator[pd.DataFrame, None, None]:
#     """
#     Generator that reads and returns Echoview NASC dataframes for given files and transects.

#     Parameters
#     ----------
#     files : list of Path
#         List of NASC CSV file paths.
#     transect_num : list of int
#         Transect numbers to include.

#     Yields
#     ------
#     pd.DataFrame
#         Parsed and validated DataFrame for each file/transect.
#     """
#     # don't need the yield structure
#     return [read_echoview_nasc(f, tnum) for f, tnum in zip(files, transect_num)]


# # based on the current export_transect_spacing
# def update_transect_spacing(
#     intervals_df: pd.DataFrame, default_transect_spacing: float
# ) -> pd.DataFrame:
#     """
#     Fill missing or invalid transect spacing values in interval data.

#     Parameters
#     ----------
#     intervals_df : pd.DataFrame
#         Interval data exported from Echoview.
#     default_transect_spacing : float
#         Fallback value for missing spacing.

#     Returns
#     -------
#     pd.DataFrame
#         Updated DataFrame with consistent spacing values.
#     """
#     pass


# # ONLY do data ingestion and organization, not writing out anything
# def merge_echoview_nasc(
#     nasc_path: Path, nasc_filename_pattern: str, default_transect_spacing: float
# ) -> pd.DataFrame:
#     """
#     Ingest and merge all Echoview NASC files (intervals, cells, layers).

#     Parameters
#     ----------
#     nasc_path : Path
#         Directory containing Echoview export files (*.csv).
#     default_transect_spacing : float
#         Default spacing to impute where missing.

#     Returns
#     -------
#     pd.DataFrame
#         Merged dataframe from intervals, cells, and layers.
#     """
    
#     # Get all echoview NASC files: analysis, cells, intervals, layers
#     # -- do not need to validate, since non-existing folders/files will error out automatically
#     # -- use .glob to get all NASC-related files
#     ev_nasc_files: dict = {
#         "analysis": nasc_path.glob("*\.csv"),  # PATTERN built from nasc_filename_pattern
#         "cells": nasc_path.glob("PATTERN"),
#         "intervals": nasc_path.glob("PATTERN"),
#         "layers": nasc_path.glob("PATTERN"),
#     }

#     # Get all transect numbers
#     # -- for each transect number, there should be 4 files
#     # -- store only transect numbers with a complete set (4) csv files
#     # -- raise warning for transect numbers with an incomplete set (<4) csv files
#     transect_num: list = get_transect_num(ev_nasc_files)

#     # Read and concat intervals, cells, and layers dataframes
#     # -- use current code in consolidate_exports
#     # -- but do not worry about validator at this time
    # df_intervals: pd.DataFrame = pd.concat(
    #     echoview_nasc_to_df(ev_nasc_files["intervals"], transect_num), ...
    # )
    # df_cells: pd.DataFrame = pd.concat(
    #     echoview_nasc_to_df(ev_nasc_files["cells"], transect_num), ...
    # )
    # df_layers: pd.DataFrame = pd.concat(
    #     echoview_nasc_to_df(ev_nasc_files["layers"], transect_num), ...
    # )

#     # Wrangle cells_df columns

#     # Update transect spacing
#     # TODO: what does update_transect_spacing do?
#     df_intervals = update_transect_spacing(df_intervals, default_transect_spacing)

#     # Explicitly merge the 3 dataframes
#     # -- do not need group_merge as a separate method
#     df_merged: pd.DataFrame
#     return df_merged


# def read_transect_region_file() -> pd.DataFrame:
#     """
#     Read external transect region definition file (if provided).

#     Parameters
#     ----------
#     filepath : str or Path
#         Path to the file defining transect-region mapping.

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame with region-classified transects.
#     """
#     pass


# # NOTE: can keep it but likely won't use in workflow,
# #       since years with this info the files are already made
# # same as the current
# def construct_transect_region_key(
#     df_merged: pd.DataFrame, region_class_mapping: dict
# ) -> pd.DataFrame:
#     """
#     Build key assigning transects to region classes using mappings.

#     Parameters
#     ----------
#     df_merged : pd.DataFrame
#         Output of merge_echoview_nasc.
#     region_class_mapping : dict
#         Dict from config defining region name -> region class mapping.

#     Returns
#     -------
#     pd.DataFrame
#         Mapping of transects to region classes.
#     """
#     pass


# # TODO: consider renaming this based on what the operation
# def load_transect_region_key(file_path: Union[str, Path]) -> pd.DataFrame:
#     """
#     Load mapping between transects and region classes.

#     Parameters
#     ----------
#     file_path : str or Path
#         path to the transect_region_key

#     Returns
#     -------
#     pd.DataFrame
#         Mapping of transects to region classes.
#     """
#     pass


# def compute_depth_layer_height(
#     df_cells: pd.DataFrame, region_names: Optional[List[str]] = None
# ) -> pd.DataFrame:
#     """
#     Compute mean depth and layer height from NASC cell data with selected regions.

#     Parameters
#     ----------
#     df_cells : pd.DataFrame
#         Cell-level NASC data.
#     region_names : list of str, optional
#         Regions to filter before computing mean depth and layer height.

#     Returns
#     -------
#     pd.DataFrame
#         Dataframe with mean depth, height per region or transect.
#     """
#     pass


# # Last section of the current ingest_echoview_exports()
# # TODO: in the current code you created interval_copy from the updated df_interval, can you not use df_merged to get the same info?
# def consolidate_echoview_nasc(
#     df_merged: pd.DataFrame,
#     region_names: List[str],
#     survey_identifier: str,
# ) -> pd.DataFrame:
#     """
#     Consolidate NASC data for selected regions into final output.

#     Parameters
#     ----------
#     df_merged : pd.DataFrame
#         Output from merge_echoview_nasc.
#     region_names : list of str
#         List of region names to include.

#     Returns
#     -------
#     pd.DataFrame
#         Final NASC dataframe (filtered and cleaned).
#     """
#     # survey_identifier controls the `filter_transect_intervals` component

#     pass
