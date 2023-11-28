from pathlib import Path
from typing import Set

import numpy as np
import pandas as pd

from ..utils.input_checks import check_column_names, check_existence_of_file


def _process_nasc_data(survey, nasc_var_types: dict) -> pd.DataFrame:
    """
    Loads in NASC data from the appropriate Excel file using the
    specified columns. Additionally, sets that data type of the
    columns in the DataFrame.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object
    nasc_var_types: dict
        A dictionary with string keys that are the NASC column names to
        grab and values are the types of those columns

    Returns
    -------
    df: pd.DataFrame
        A DataFrame filled with the requested NASC data
    """

    # select and check the appropriate nasc data file
    if survey.params["exclude_age1"]:
        # check existence of the file
        file_path = (
            survey.params["data_root_dir"] / survey.params["nasc_no_age1_filename"]
        )
        check_existence_of_file(file_path)

        df = pd.read_excel(
            file_path, sheet_name=survey.params["nasc_no_age1_sheetname"]
        )
    else:
        # check existence of the file
        file_path = (
            survey.params["data_root_dir"] / survey.params["nasc_all_ages_filename"]
        )
        check_existence_of_file(file_path)

        df = pd.read_excel(
            file_path, sheet_name=survey.params["nasc_all_ages_sheetname"]
        )

    check_column_names(df=df, expected_names=set(nasc_var_types.keys()), path_for_df=file_path)

    # obtaining those columns that are required
    df = df[nasc_var_types.keys()]

    # set data types of dataframe
    df = df.astype(nasc_var_types)

    if survey.params["survey_year"] < 2003:
        # TODO: it may be the case that we need to include lines 35-61 of
        #  EchoPro/general/load_files_parameters/get_NASC_data.m
        raise NotImplementedError(
            "Loading the NASC table for survey years less than 2003 has not been implemented!"
        )

    return df


def load_nasc_df(survey) -> pd.DataFrame:
    """
    Load VL interval-based NASC table.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object

    Returns
    -------
    Pandas Dataframe of NASC table.
    """

    # specify column names to grab and their corresponding type
    nasc_var_types = {
        "transect_num": int,
        "vessel_log_start": np.float64,
        "vessel_log_end": np.float64,
        "latitude": np.float64,
        "longitude": np.float64,
        "stratum_num": int,
        "transect_spacing": np.float64,
        "NASC": np.float64,
        "haul_num": int,
    }

    df = _process_nasc_data(survey, nasc_var_types)

    # set dataframe index
    df.set_index("transect_num", inplace=True)

    return df
