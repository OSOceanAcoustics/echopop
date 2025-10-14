from pathlib import Path
from typing import Dict, Union

import pandas as pd


def load_mesh_data(
    mesh_filepath: Union[str, Path], sheet_name: str, column_name_map: Dict[str, str] = {}
) -> pd.DataFrame:
    """
    Load mesh data from an Excel file.

    Parameters
    ----------
    mesh_filepath : str or Path
        Path to the Excel file containing mesh data
    sheet_name : str
        Name of the sheet containing the mesh data
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names

    Returns
    -------
    pd.DataFrame
        DataFrame containing the mesh data with centroid coordinates and fractions

    Examples
    --------
    >>> mesh_df = load_mesh_data("mesh_file.xlsx", "Grid")
    >>> mesh_df = load_mesh_data("mesh_file.xlsx", column_name_map={"centroid_latitude":
    "latitude"})
    """
    mesh_filepath = Path(mesh_filepath)

    if not mesh_filepath.exists():
        raise FileNotFoundError(f"Mesh data file not found: {mesh_filepath}")

    # Read Excel file into memory
    df = pd.read_excel(mesh_filepath, sheet_name=sheet_name, index_col=None, header=0)

    # Force the column names to be lower case
    df.columns = df.columns.str.lower()

    # Rename columns if mapping is provided
    if column_name_map:
        df.rename(columns=column_name_map, inplace=True)

    return df


def load_isobath_data(
    isobath_filepath: Union[str, Path], sheet_name: str, column_name_map: Dict[str, str] = {}
) -> pd.DataFrame:
    """
    Load isobath data from an Excel file.

    Parameters
    ----------
    isobath_filepath : str or Path
        Path to the Excel file containing mesh data
    sheet_name : str
        Name of the sheet containing the mesh data
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names

    Returns
    -------
    pd.DataFrame
        DataFrame containing the isobath data with longitude and latitude

    Examples
    --------
    >>> isobath_df = load_isobath_data("isobath_file.xlsx", "Sheet1")
    >>> isobath_df = load_isobath_data("isobath_file.xlsx", column_name_map={"latitude_200m":
    "latitude"})
    """
    mesh_filepath = Path(isobath_filepath)

    if not mesh_filepath.exists():
        raise FileNotFoundError(f"Mesh data file not found: {isobath_filepath}")

    # Read Excel file into memory
    df = pd.read_excel(isobath_filepath, sheet_name=sheet_name, index_col=None, header=0)

    # Force the column names to be lower case
    df.columns = df.columns.str.lower()

    # Rename columns if mapping is provided
    if column_name_map:
        df.rename(columns=column_name_map, inplace=True)

    return df
