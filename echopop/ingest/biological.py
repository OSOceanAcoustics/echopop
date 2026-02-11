import copy
import itertools
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def load_single_biological_sheet(
    biodata_filepath: Path,
    sheet_name: str,
    column_name_map: Dict = {},
    subset_dict: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    Load and process a single biological data sheet.

    Parameters
    ----------
    biodata_filepath : pathlib.Path
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
    |pd.DataFrame|
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


def load_biological_data(
    biodata_filepath: Path,
    biodata_sheet_map: Dict[str, str],
    column_name_map: Dict[str, str] = None,
    subset_dict: Optional[Dict] = None,
    biodata_label_map: Optional[Dict[str, Dict]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Load biological data from a single Excel file with multiple sheets.

    Parameters
    ----------
    biodata_filepath : pathlib.Path
        Path to the Excel file containing biological data
    biodata_sheet_map : dict
        Dictionary mapping dataset names to sheet names
        (e.g., ``{"specimen": "biodata_specimen", "length": "biodata_length", "catch":
        "biodata_catch"}``)
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names
        (e.g., ``{"frequency": "length_count", "haul": "haul_num"}``)
    subset_dict : dict, optional
        Subset dictionary containing ships and species_code for filtering
        Format: ``{"ships": {ship_id: {"survey": survey_id, "haul_offset": offset}}, "species_code":
        [codes]}``
    biodata_label_map : dict, optional
        Dictionary mapping column names to value replacement dictionaries (e.g.,
        ``{"sex": {1: "male", 2: "female", 3: "unsexed"}}``)

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
    df : |pd.DataFrame|
        Input DataFrame to filter
    subset_dict : dict, optional
        Subset dictionary containing ships and species_code settings
        Format: ``{"ships": {ship_id: {"survey": survey_id, "haul_offset": offset}}, "species_code":
        [codes]}``
    Returns
    -------
    |pd.DataFrame|
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
        if "ship" in df_filtered.columns:
            df_filtered = df_filtered.loc[df_filtered["ship"].isin(ship_ids)]
        # ---- Collect survey IDs, if present
        survey_ids = list(
            itertools.chain.from_iterable(
                [
                    v["survey"] if isinstance(v["survey"], list) else [v["survey"]]
                    for v in ship_config.values()
                    if "survey" in v
                ]
            )
        )
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


def generate_composite_key(
    bio_data: Dict[str, pd.DataFrame],
    index_columns: List[str],
    adjust: Dict[str, np.number] = {},
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame, Dict[str, pd.DataFrame]]:

    # Validate that all columns exist in the underlying biodata
    valid_id_columns = {
        k: set(index_columns) <= set(v.columns) or "uid" in v.columns for k, v in bio_data.items()
    }
    # ---- Collect invalid entries
    invalid_entries = [k for k, v in valid_id_columns.items() if not v]
    # ---- Raise KeyError
    if invalid_entries:
        id_col_str = "', '".join(index_columns)
        entr_str = "', '".join(invalid_entries)
        raise KeyError(
            f"Both defined ID columns ('{id_col_str}') or preset 'uid' column missing from "
            f"the following biodata DataFrames: '{entr_str}'."
        )

    # Create copy
    bio_data_copy = copy.deepcopy(bio_data)

    # Apply adjustments
    for column, adjustment in adjust.items():
        # ---- Format new name
        new_name = column + "_adj"
        # ---- Apply adjustment
        bio_data_copy = {
            k: v.assign(
                **{
                    new_name: np.where(
                        v[column] + adjustment > 0, v[column] + adjustment, v[column]
                    )
                }
            )
            for k, v in bio_data_copy.items()
        }

    # Identify the constructor columns
    constructor_columns = [c if c not in adjust else c + "_adj" for c in index_columns]

    # Join the columns to create unique id column
    bio_data_copy = {
        k: v.assign(uid=v[constructor_columns].astype(str).agg("-".join, axis=1)).filter(
            list(bio_data[k].columns) + ["uid"]
        )
        for k, v in bio_data_copy.items()
    }

    # Get unique columns + uid for merging into downstream DataFrames
    unique_uid = pd.concat(
        [d.filter(index_columns + ["uid"]) for d in bio_data_copy.values()]
    ).drop_duplicates()

    # Return the unique IDs
    return unique_uid


def apply_composite_key(
    data: pd.DataFrame,
    composite_key: pd.DataFrame,
    haul_offset: Dict[np.number, np.number] = (),
) -> pd.DataFrame:

    # Apply adjustments to the original composite key columns
    composite_id = composite_key.copy()
    if len(haul_offset) == 2:
        composite_id["haul_num"] = np.where(
            composite_id["ship"] == haul_offset[0],
            composite_id["haul_num"] + haul_offset[1],
            composite_id["haul_num"],
        )

    # Drop "uid" from data
    if "uid" in data.columns:
        data.drop(columns=["uid"], inplace=True)

    # Left-merge
    shared_columns = [c for c in composite_id.columns if c in data.columns and c != "uid"]
    data_uid = data.merge(
        composite_id.filter(shared_columns + ["uid"]), on=shared_columns, how="left"
    )

    # Return
    return data_uid
