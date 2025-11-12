from pathlib import Path
from typing import Dict, Union

import numpy as np
import pandas as pd


def load_single_stratum_sheet(
    strata_filepath: Path,
    sheet_name: str,
    column_name_map: Dict[str, str] = {},
) -> pd.DataFrame:
    """
    Load a single stratification sheet from an Excel file.

    Parameters
    ----------
    strata_filepath : Path
        Path to the Excel file containing stratification data
    sheet_name : str
        Name of the sheet to load
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names

    Returns
    -------
    pandas.DataFrame
        Processed stratification DataFrame
    """
    # Read Excel file into memory
    df = pd.read_excel(strata_filepath, sheet_name=sheet_name, index_col=None, header=0)

    # Force the column names to be lower case
    df.columns = df.columns.str.lower()

    # Rename columns if mapping is provided
    if column_name_map:
        df.rename(columns=column_name_map, inplace=True)

    return df


def load_strata(
    strata_filepath: Path,
    strata_sheet_map: Dict[str, str],
    column_name_map: Dict[str, str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Load stratification data from an Excel file with multiple sheets.

    Parameters
    ----------
    strata_filepath : Path
        Path to the Excel file containing stratification data
    strata_sheet_map : dict
        Dictionary mapping stratification types to sheet names
        (e.g., {"inpfc": "INPFC", "ks": "stratification1"})
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names
        (e.g., {"fraction_hake": "nasc_proportion", "haul": "haul_num"})

    Returns
    -------
    dict
        Dictionary containing stratification DataFrames keyed by stratification type

    Examples
    --------
    >>> sheet_map = {"inpfc": "INPFC", "ks": "stratification1"}
    >>> col_map = {"fraction_hake": "nasc_proportion", "haul": "haul_num"}
    >>> strata_data = load_stratification("strata_file.xlsx", sheet_map, col_map)
    """

    if not strata_filepath.exists():
        raise FileNotFoundError(f"Stratification file not found: {strata_filepath}")

    # Load each stratification sheet - dictionary comprehension instead of for loop
    strata_dict = {
        strata_type: load_single_stratum_sheet(strata_filepath, sheet_name, column_name_map)
        for strata_type, sheet_name in strata_sheet_map.items()
    }

    return strata_dict


def load_geostrata(
    geostrata_filepath: Union[str, Path],
    geostrata_sheet_map: Dict[str, str],
    column_name_map: Dict[str, str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Load geographic stratification data from an Excel file with multiple sheets.

    Parameters
    ----------
    geostrata_filepath : str or Path
        Path to the Excel file containing geographic stratification data
    geostrata_sheet_map : dict
        Dictionary mapping stratification types to sheet names
        (e.g., {"inpfc": "INPFC", "ks": "stratification1"})
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names
        (e.g., {"Latitude (upper limit)": "northlimit_latitude", "stratum": "stratum_num"})

    Returns
    -------
    dict
        Dictionary containing geographic stratification DataFrames keyed by stratification type,
        each with consolidated latitude intervals from INPFC and KS strata assignments
    """

    if not geostrata_filepath.exists():
        raise FileNotFoundError(f"Geographic stratification file not found: {geostrata_filepath}")

    # First load each sheet
    raw_geostrata_dict = {
        strata_type: load_single_stratum_sheet(geostrata_filepath, sheet_name, column_name_map)
        for strata_type, sheet_name in geostrata_sheet_map.items()
    }

    # Then process each dataframe with the extracted function
    processed_geostrata_dict = {
        strata_type: geostrata_bins(df) for strata_type, df in raw_geostrata_dict.items()
    }

    return processed_geostrata_dict


def geostrata_bins(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process a geographic stratification DataFrame by adding latitude intervals
    and renaming columns as needed.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with geographic stratification data

    Returns
    -------
    pd.DataFrame
        Processed DataFrame with latitude intervals added
    """
    result_df = df.copy()

    if "northlimit_latitude" in result_df.columns:
        # Sort by latitude
        result_df.sort_values(["northlimit_latitude"], inplace=True)

        # Create latitude bins
        latitude_bins = np.concatenate([[-90], result_df["northlimit_latitude"].unique(), [90]])

        # Add categorical intervals
        result_df["latitude_interval"] = pd.cut(result_df["northlimit_latitude"], latitude_bins)

    return result_df


def join_strata_by_haul(
    data: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    strata_df: Dict[str, pd.DataFrame],
    default_stratum: float = 0.0,
    stratum_name: str = "stratum_num",
) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Join stratification data by haul number or other matching column.

    Parameters
    ----------
    data : pd.DataFrame or Dict[str, pd.DataFrame]
        DataFrame or dictionary of DataFrames to join with strata
    strata_df : pd.DataFrame
        Specific stratification DataFrame with stratum-haul key information
    default_stratum : float
        Default stratum value when there are no matching/corresponding values
    stratum_name : str, default="stratum_num"
        Name of the column containing stratum information

    Returns
    -------
    Same type as input data with stratification added
    """

    # Get stratification columns (excluding join column)
    strata_cols = [col for col in strata_df.columns if col != "haul_num"]

    # Function to join a single DataFrame
    def join_single_df(df):
        if not isinstance(df, pd.DataFrame) or "haul_num" not in df.columns:
            return df
        if "haul_num" not in strata_df.columns:
            return df

        # Check if stratification columns already exist
        existing_cols = set(strata_cols).intersection(set(df.columns))
        if existing_cols:
            # Drop existing shared columns
            df = df.drop(columns=list(existing_cols))

        # Merge
        df_merged = df.merge(strata_df, on=["haul_num"], how="left")

        # Rename the stratum column name, if needed
        if stratum_name not in df_merged.columns:
            df_merged.rename(columns={"stratum_num": stratum_name}, inplace=True)

        # Replace missing strata with `default_stratum`
        df_merged[stratum_name] = df_merged[stratum_name].fillna(default_stratum)

        return df_merged

    # Apply based on input type
    if isinstance(data, pd.DataFrame):
        return join_single_df(data)
    elif isinstance(data, dict):
        return {key: join_single_df(df) for key, df in data.items()}
    else:
        raise TypeError("Input 'data' must be DataFrame or dictionary of DataFrames")


def join_geostrata_by_latitude(
    data: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    geostrata_df: pd.DataFrame,
    stratum_name: str = "stratum_num",
) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Join geographic stratification data by latitude intervals.

    Parameters
    ----------
    data : pd.DataFrame or Dict[str, pd.DataFrame]
        DataFrame or dictionary of DataFrames with latitude information
    geostrata_df : pd.DataFrame
        Geographic stratification DataFrame with latitude boundaries and stratum info
    stratum_name : str, default="stratum_num"
        Name of the column containing stratum information

    Returns
    -------
    Same type as input data with geostratification added
    """
    # Sort the geostrata DataFrame by latitude
    geostrata_df = geostrata_df.copy().sort_values("northlimit_latitude")

    # Create latitude bins
    latitude_bins = np.concatenate([[-90], geostrata_df["northlimit_latitude"].unique(), [90]])

    # Get geostrata columns (excluding northlimit_latitude)
    geostrata_cols = [col for col in geostrata_df.columns if col != "northlimit_latitude"]

    # Function to join a single DataFrame
    def join_single_df(df):
        if not isinstance(df, pd.DataFrame) or "latitude" not in df.columns:
            return df

        # Check if geostratification columns already exist
        existing_cols = set(geostrata_cols).intersection(set(df.columns))
        if existing_cols:
            # Drop existing geostratification columns if overwrite is True
            df = df.drop(columns=list(existing_cols))

        result = df.copy()

        # Create latitude intervals
        result[stratum_name] = pd.cut(
            result["latitude"],
            latitude_bins,
            right=False,
            labels=list(geostrata_df["stratum_num"]) + [1],
            ordered=False,
        )

        return result

    # Apply based on input type
    if isinstance(data, pd.DataFrame):
        return join_single_df(data)
    elif isinstance(data, dict):
        return {key: join_single_df(df) for key, df in data.items()}
    else:
        raise TypeError("Input 'data' must be DataFrame or dictionary of DataFrames")
