from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Union

import pandas as pd
import numpy as np
import re
import functools
from ..ingest import read_csv_file
from ..core.echopop_columns import ECHOVIEW_TO_ECHOPOP

def map_transect_num(
    ev_export_paths: Dict[str, Generator], transect_pattern: str = r"T(\d+)" 
) -> Dict[float, Dict[str, Path]]:
    """
    Extract and map transect numbers to complete sets of Echoview export files.

    Parameters
    ----------
    ev_export_paths : Dict[str, Generator]
        Dictionary with keys "analysis", "cells", "intervals", "layers", with values being 
        `pathlib.Path.glob` generators pointing torward Echoview export *.csv files.
    transect_pattern : str
        Regex pattern to extract transect numbers from fil
        e names.

    Returns
    -------
    Dict[str, List[float, pathlib.Path]]:
        Dictionary with transect numbers as keys and dictionaries of file types ('analysis', 
        'cells', 'intervals', and 'layers') and filepaths as values.
    """

    # Compite the transect pattern regex
    compiled_pattern = re.compile(transect_pattern) 

    # Reduce dictitonary into a dictionary with list of tuples with transect number and filepath
    transect_file_mapping = functools.reduce(
        lambda acc, item: {
            **acc,
            item[0]: acc.get(item[0], []) + [
                (float(compiled_pattern.search(str(path)).group(1)), path)
                for path in list(item[1])
                if compiled_pattern.search(str(path))
            ]
        },
        ev_export_paths.items(),
        {}
    )
    
    # Return the transect_file_mapping
    return transect_file_mapping

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
        Cleaned and formatted data.
    """

    # Read in the defined CSV file
    nasc_df = read_echoview_export(filename, validator)
    
    # Add transect number
    nasc_df["transect_num"] = transect_num
    
    # Fix latitude and longitude
    # ---- Latitude
    impute_bad_coordinates(nasc_df, "latitude")
    # ---- Longitude
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
    file_transect_pairs: list[tuple[float, Path]]  # List of (transect_num, filepath) tuples
) -> list[pd.DataFrame]:
    """
    Reads and returns Echoview NASC dataframes for given file/transect pairs.

    Parameters
    ----------
    file_transect_pairs : list of tuple
        List of (transect_num: float, filepath: pathlib.Path) tuples.

    Returns
    -------
    list[pd.DataFrame]
        List of parsed and validated DataFrames for each file/transect.
    """
    # Note: We swap the order since the tuples are (transect_num, filepath)
    return [read_echoview_nasc(path, t_num) for t_num, path in file_transect_pairs]

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
#     df_intervals: pd.DataFrame = pd.concat(
#         echoview_nasc_to_df(ev_nasc_files["intervals"], transect_num), ...
#     )
#     df_cells: pd.DataFrame = pd.concat(
#         echoview_nasc_to_df(ev_nasc_files["cells"], transect_num), ...
#     )
#     df_layers: pd.DataFrame = pd.concat(
#         echoview_nasc_to_df(ev_nasc_files["layers"], transect_num), ...
#     )

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
