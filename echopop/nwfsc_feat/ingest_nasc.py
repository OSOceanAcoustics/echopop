from typing import List, Union, Generator, Optional
from pathlib import Path
import pandas as pd



def read_echoview_export(filename: Union[str, Path]) -> pd.DataFrame:
    """
    Generic reader for Echoview export CSVs. Used for files like analysis, cells, layers, intervals.
    
    Parameters
    ----------
    filename : str or Path
        Path to the Echoview export CSV file.

    Returns
    -------
    pd.DataFrame
        Parsed CSV as a DataFrame.
    """
    # Select only 1 type is sufficient since the list passed are from transects with a complete set  of 4 csv files
    pass


def get_transect_num(ev_nasc_files: dict[str, list[Path]]) -> list[int]:
    """
    Identify transect numbers with a complete set of required NASC export files.

    Parameters
    ----------
    ev_nasc_files : dict
        Dictionary with keys "analysis", "cells", "intervals", "layers" and values being lists of Path objects.

    Returns
    -------
    list of int
        List of valid transect numbers (integers).
    """
    pass


# based on the current read_echoview_export
def read_echoview_nasc(filename: Union[str, Path], transect_num: int) -> pd.DataFrame:
    """
    Read and process a single Echoview NASC CSV file for a given transect.

    Parameters
    ----------
    filename : str or Path
        Full path to the NASC CSV file.
    transect_num : int
        Transect number to use for filtering or labeling.

    Returns
    -------
    pd.DataFrame
        Cleaned and formatted data.
    """
    pass


# based on the current generate_dataframes
def echoview_nasc_to_df(files: list[Path], transect_num: list[int]) -> Generator[pd.DataFrame, None, None]:
    """
    Generator that reads and returns Echoview NASC dataframes for given files and transects.

    Parameters
    ----------
    files : list of Path
        List of NASC CSV file paths.
    transect_num : list of int
        Transect numbers to include.

    Yields
    ------
    pd.DataFrame
        Parsed and validated DataFrame for each file/transect.
    """
    # don't need the yield structure
    return [
        read_echoview_nasc(f, tnum)
        for f, tnum in zip(files, transect_num)
    ]
    

# based on the current export_transect_spacing
def update_transect_spacing(intervals_df: pd.DataFrame, default_transect_spacing: float) -> pd.DataFrame:
    """
    Fill missing or invalid transect spacing values in interval data.

    Parameters
    ----------
    intervals_df : pd.DataFrame
        Interval data exported from Echoview.
    default_transect_spacing : float
        Fallback value for missing spacing.

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with consistent spacing values.
    """
    pass


# ONLY do data ingestion and organization, not writing out anything
def merge_echoview_nasc(
    nasc_path: Union[str, Path],
    nasc_filename_pattern: str,
    default_transect_spacing: float
) -> pd.DataFrame:
    """
    Ingest and merge all Echoview NASC files (intervals, cells, layers).

    Parameters
    ----------
    nasc_path : str or Path
        Directory containing the NASC export files.
    nasc_filename_pattern : str
        Glob pattern to locate the NASC files.
    default_transect_spacing : float
        Default spacing to impute where missing.

    Returns
    -------
    pd.DataFrame
        Merged dataframe from intervals, cells, and layers.
    """
    # Get all echoview NASC files: analysis, cells, intervals, layers
    # -- do not need to validate, since non-existing folders/files will error out automatically
    # -- use .glob to get all NASC-related files
    ev_nasc_files: dict = {
        "analysis": nasc_path.glob("PATTERN"),  # PATTERN built from nasc_filename_pattern
        "cells": nasc_path.glob("PATTERN"),
        "intervals": nasc_path.glob("PATTERN"),
        "layers": nasc_path.glob("PATTERN"),
    }

    # Get all transect numbers
    # -- for each transect number, there should be 4 files
    # -- store only transect numbers with a complete set (4) csv files
    # -- raise warning for transect numbers with an incomplete set (<4) csv files
    transect_num: list = get_transect_num(ev_nasc_files)

    # Read and concat intervals, cells, and layers dataframes
    # -- use current code in consolidate_exports
    # -- but do not worry about validator at this time
    df_intervals: pd.DataFrame = pd.concat(echoview_nasc_to_df(ev_nasc_files["intervals"], transect_num), ...)
    df_cells: pd.DataFrame = pd.concat(echoview_nasc_to_df(ev_nasc_files["cells"], transect_num), ...)
    df_layers: pd.DataFrame = pd.concat(echoview_nasc_to_df(ev_nasc_files["layers"], transect_num), ...)

    # Wrangle cells_df columns

    # Update transect spacing
    # TODO: what does update_transect_spacing do?
    df_intervals = update_transect_spacing(df_intervals, default_transect_spacing)

    # Explicitly merge the 3 dataframes
    # -- do not need group_merge as a separate method
    df_merged: pd.DataFrame
    return df_merged


def read_transect_region_file() -> pd.DataFrame:
    """
    Read external transect region definition file (if provided).

    Parameters
    ----------
    filepath : str or Path
        Path to the file defining transect-region mapping.

    Returns
    -------
    pd.DataFrame
        DataFrame with region-classified transects.
    """
    pass


# NOTE: can keep it but likely won't use in workflow, 
#       since years with this info the files are already made
# same as the current
def construct_transect_region_key(
    df_merged: pd.DataFrame,
    region_class_mapping: dict
) -> pd.DataFrame:
    """
    Build key assigning transects to region classes using mappings.

    Parameters
    ----------
    df_merged : pd.DataFrame
        Output of merge_echoview_nasc.
    region_class_mapping : dict
        Dict from config defining region name -> region class mapping.

    Returns
    -------
    pd.DataFrame
        Mapping of transects to region classes.
    """
    pass


# TODO: consider renaming this based on what the operation
def load_transect_region_key(file_path: Union[str, Path]) -> pd.DataFrame:
    """
    Load mapping between transects and region classes.

    Parameters
    ----------
    file_path : str or Path
        path to the transect_region_key

    Returns
    -------
    pd.DataFrame
        Mapping of transects to region classes.
    """
    pass


def compute_depth_layer_height(
    df_cells: pd.DataFrame,
    region_names: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Compute mean depth and layer height from NASC cell data with selected regions.

    Parameters
    ----------
    df_cells : pd.DataFrame
        Cell-level NASC data.
    region_names : list of str, optional
        Regions to filter before computing mean depth and layer height.

    Returns
    -------
    pd.DataFrame
        Dataframe with mean depth, height per region or transect.
    """
    pass


# Last section of the current ingest_echoview_exports()
# TODO: in the current code you created interval_copy from the updated df_interval, can you not use df_merged to get the same info?
def consolidate_echoview_nasc(
    df_merged: pd.DataFrame,
    region_names: List[str],
    survey_identifier: str,
) -> pd.DataFrame:
    """
    Consolidate NASC data for selected regions into final output.

    Parameters
    ----------
    df_merged : pd.DataFrame
        Output from merge_echoview_nasc.
    region_names : list of str
        List of region names to include.

    Returns
    -------
    pd.DataFrame
        Final NASC dataframe (filtered and cleaned).
    """
    # survey_identifier constrols the `filter_transect_intervals` component

    pass
