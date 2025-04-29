from typing import Union, Dict, Tuple
from pathlib import Path
import pandas as pd


# Break up load_data() and read_validated_data() into 3 separate functions to load biological, stratification, and kriging data
# Don't worry about validating input files for now


# TODO: read from the master spreadsheet
# TODO: combine in content of preprocess_biodata()
# TODO: filter to only output data from 1 species
def load_biological_data(
    root_path: Union[str, Path],
    file_path_dict: Dict,
    species_code: str
) -> Dict[pd.DataFrame]:
    """
    Load biological data from master biological data spreadsheet.

    Parameters
    ----------
    root_path : str or Path
        Path to the biological data file
    file_path_dict : dict
        Dictionary of paths to individual biological data files

    Returns
    -------
    df_length : pd.DataFrame
    """
    df_bio_dict: Dict = {
        "length": pd.DataFrame,
        "specimen": pd.DataFrame,
        "catch": pd.DataFrame,
    }
    return df_bio_dict


# TODO: combine in the current preprocess_spatial()
def load_stratification(
    root_path: Union[str, Path],
    file_path_dict: Dict,
) -> Dict[pd.DataFrame]:
    """
    Load stratification schemes from CSV

    Parameters
    ----------
    root_path : str or Path
        Path to stratification CSV
    file_path_dict : dict
        Dictionary of paths to individual stratification files

    Returns
    -------
    A dictionary of dataframes containing bio_strate and geo_strata info
    """

    # `bio_strata` is the current `strata`

    # In addition to the original operations in preprocess_spatial()
    # also do stratum renaming for inpfc originally in transect.py::save_transect_coordinates()
    # the two output dataframes should have the same column names
    #
    # BELOW CODE SNIPPET FROM save_transect_coordinates()
    # if stratum_def == "inpfc":
    #     stratum_rename = "stratum_num"
    # else:
    #     stratum_rename = stratum_col

    df_strata_dict: Dict = {
        "ks": pd.DataFrame,
        "inpfc": pd.DataFrame,
    }
    return df_strata_dict


# same as the current preprocess_biology_spatial()
def join_biological_stratification(
    df_bio_dict: Dict[pd.DataFrame], df_strata_dict: Dict[pd.DataFrame]
) -> Dict[pd.DataFrame]:
    return df_bio_dict


# same as the current preprocess_acoustic_spatial()
def join_acoustic_stratification(
    df_nasc: pd.DataFrame, df_strata_dict: Dict[pd.DataFrame]
) -> pd.DataFrame:
    return df_nasc


# same as the current preprocess_acoustic_biology_spatial()
def join_acoustic_all(
    df_nasc: pd.DataFrame,
    df_bio_dict: Dict[pd.DataFrame],
    df_strata_dict: Dict[pd.DataFrame],
) -> pd.DataFrame:
    return df_nasc


def consolidate_all_data(
    df_nasc: pd.DataFrame,
    df_bio_dict: Dict[pd.DataFrame],
    df_strata_dict: Dict[pd.DataFrame],
) -> pd.DataFrame:
    """
    Consolidate all input data.
    """
    df_bio_dict = join_biological_stratification(df_bio_dict, df_strata_dict)
    df_nasc = join_acoustic_stratification(df_nasc, df_strata_dict)
    df_nasc = join_acoustic_all(df_nasc, df_bio_dict, df_strata_dict)
    return df_nasc


# separate out components that are parameters
def load_kriging_variogram_params(
    root_path: Union[str, Path],
    file_path_dict: Dict,
    kriging_const: dict,
) -> Tuple[Dict, Dict]:
    # put what's in 
    kriging_params_dict: dict
    variogram_params_dict: dict
    return kriging_params_dict, variogram_params_dict
