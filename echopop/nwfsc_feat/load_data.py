import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd


def load_biological_data(
    biodata_filepath: Path,
    biodata_sheet_map: Dict[str, str],
    column_name_map: Dict[str, str] = None,
    subset_dict: Optional[Dict] = None,
    biodata_label_map: Optional[Dict[str, Dict]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Load biological data from a single Excel file with multiple sheets.
    Load biological data from a single Excel file with multiple sheets.

    Parameters
    ----------
    biodata_filepath : Path
        Path to the Excel file containing biological data
    biodata_sheet_map : dict
        Dictionary mapping dataset names to sheet names
        (e.g., {"specimen": "biodata_specimen", "length": "biodata_length", "catch":
        "biodata_catch"})
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names
        (e.g., {"frequency": "length_count", "haul": "haul_num"})
    subset_dict : dict, optional
        Subset dictionary containing ships and species_code for filtering
        Format: {"ships": {ship_id: {"survey": survey_id, "haul_offset": offset}}, "species_code":
        [codes]}
    biodata_label_map : dict, optional
        Dictionary mapping column names to value replacement dictionaries
        (e.g., {"sex": {1: "male", 2: "female", 3: "unsexed"}})

    Returns
    -------
    dict
        Dictionary containing processed biological DataFrames keyed by dataset name

    Examples
    --------
    >>> sheet_map = {"catch": "biodata_catch", "length": "biodata_length"}
    >>> subset = {"ships": {160: {"survey": 201906}}, "species_code": [22500]}
    >>> col_map = {"frequency": "length_count", "haul": "haul_num"}
    >>> label_map = {"sex": {1: "male", 2: "female", 3: "unsexed"}}
    >>> bio_data = load_biological_data("biodata.xlsx", sheet_map, col_map, subset, label_map)
    """

    if not biodata_filepath.exists():
        raise FileNotFoundError(f"Biological data file not found: {biodata_filepath}")

    # Load each biological dataset
    biodata_dict = {
        dataset_name: load_single_biological_sheet(
            biodata_filepath, sheet_name, column_name_map, subset_dict
        )
        for dataset_name, sheet_name in biodata_sheet_map.items()
    }

    # Apply label mappings if provided
    if biodata_label_map:
        # ---- For each column mapping in the label map
        for col, mapping in biodata_label_map.items():
            # ---- Apply to each dataframe that has that column
            for name, df in biodata_dict.items():
                if isinstance(df, pd.DataFrame) and col in df.columns:
                    df[col] = df[col].map(mapping).fillna(df[col])

    return biodata_dict


def apply_ship_survey_filters(
    df: pd.DataFrame,
    subset_dict: Optional[Dict] = None,
) -> pd.DataFrame:
    """
    Apply ship ID, survey ID, and haul offset filters to biological data.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to filter
    subset_dict : dict, optional
        Subset dictionary containing ships and species_code settings
        Format: {"ships": {ship_id: {"survey": survey_id, "haul_offset": offset}}, "species_code":
        [codes]}

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame
    """
    # Create copy
    df_filtered = df.copy()  # Fixed: df_initial -> df

    # Check if subset_dict exists
    if subset_dict is None:
        return df_filtered

    # Extract ship configuration from nested structure
    if "ships" in subset_dict:
        # ---- Get general ship ID's
        ship_ids = [*subset_dict["ships"].keys()]
        # ---- Get ship-specific configurations
        ship_config = subset_dict["ships"]
        # ---- Apply ship-based filter
        if "ship_id" in df_filtered.columns:
            df_filtered = df_filtered.loc[df_filtered["ship_id"].isin(ship_ids)]
        # ---- Collect survey IDs, if present
        survey_ids = [v["survey"] for v in ship_config.values() if "survey" in v]
        # ---- Apply survey-based filter
        if survey_ids and "survey" in df_filtered.columns:
            df_filtered = df_filtered.loc[df_filtered["survey"].isin(survey_ids)]
        # ---- Collect haul offsets, if any are present
        haul_offsets = {k: v["haul_offset"] for k, v in ship_config.items() if "haul_offset" in v}
        # ---- Apply haul number offsets, if defined
        if haul_offsets and "ship" in df_filtered.columns:
            df_filtered["haul_num"] = df_filtered["haul_num"] + df_filtered["ship"].map(
                haul_offsets
            ).fillna(0)

    # Filter species
    if "species_code" in subset_dict:
        df_filtered = df_filtered.loc[df_filtered["species_code"].isin(subset_dict["species_code"])]

    return df_filtered


def load_single_biological_sheet(
    biodata_filepath: Path,
    sheet_name: str,
    column_name_map: Dict = {},
    subset_dict: Optional[Dict] = None,
) -> pd.DataFrame:
    """
    Load and process a single biological data sheet.

    Parameters
    ----------
    biodata_filepath : Path
        Path to the Excel file
    sheet_name : str
        Name of the sheet to load
    column_name_map : dict, default={}
        Dictionary mapping original column names to new column names
    subset_dict : dict, optional
        Subset dictionary containing ships and species_code settings
        Format: {"ships": {ship_id: {"survey": survey_id, "haul_offset": offset}}, "species_code":
        [codes]}

    Returns
    -------
    pd.DataFrame
        Processed and validated biological data
    """
    # Read Excel file into memory
    df_initial = pd.read_excel(biodata_filepath, sheet_name=sheet_name, index_col=None, header=0)

    # Force the column names to be lower case
    df_initial.columns = df_initial.columns.str.lower()

    # Rename the columns
    df_initial.rename(columns=column_name_map, inplace=True)

    # Apply ship and survey filtering
    df_filtered = apply_ship_survey_filters(df_initial, subset_dict)

    return df_filtered


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
            # Drop existing stratification columns
            warnings.warn(
                f"Dropping existing stratification columns {existing_cols} from the dataframe."
            )
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


def load_kriging_variogram_params(
    geostatistic_params_filepath: Union[str, Path],
    sheet_name: str,
    column_name_map: Dict[str, str] = {},
) -> Tuple[Dict, Dict]:
    """
    Load kriging and variogram parameters from Excel file.

    Parameters
    ----------
    geostatistic_params_filepath : str or Path
        Path to the Excel file containing kriging parameters
    sheet_name : str
        Name of the sheet to load
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names

    Returns
    -------
    tuple
        Tuple containing (kriging_params_dict, variogram_params_dict)
    """

    if not geostatistic_params_filepath.exists():
        raise FileNotFoundError(
            f"Kriging parameter file not found: {geostatistic_params_filepath.as_posix()}"
        )

    # Read Excel file into memory and then transpose
    df_initial = pd.read_excel(
        geostatistic_params_filepath, sheet_name=sheet_name, index_col=None, header=None
    ).T

    # Take the values from the first row and redefine them as the column headers
    df_initial.columns = df_initial.iloc[0]
    df_initial = df_initial.drop(0)

    # Extract kriging parameters
    kriging_params = (
        df_initial.filter(regex="krig[.]")
        .rename(columns=lambda x: x.replace("krig.", ""))
        .rename(columns=column_name_map)
        .to_dict(orient="records")[0]
    )

    # Extract variogram parameters
    variogram_params = (
        df_initial.filter(regex="vario[.]")
        .rename(columns=lambda x: x.replace("vario.", ""))
        .rename(columns=column_name_map)
        .to_dict(orient="records")[0]
    )

    return kriging_params, variogram_params


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
