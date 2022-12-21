from pathlib import Path

import numpy as np
import pandas as pd

from ..utils.input_checks import check_column_names, check_existence_of_file

nasc_cols = {
    "transect_num",
    "vessel_log_start",
    "vessel_log_end",
    "latitude",
    "longitude",
    "stratum_num",
    "transect_spacing",
    "NASC",
    "haul_num",
}


def _check_nasc_df(nasc_df: pd.DataFrame, df_path: Path) -> None:
    """
    Ensures that the appropriate columns are
    contained in the NASC Dataframe.

    Parameters
    ----------
    nasc_df: pd.DataFrame
        The constructed NASC DataFrame
    df_path: Path
        The path to the Excel file used to construct the DataFrame
    """

    # TODO: should we add more in-depth checks here?

    check_column_names(df=nasc_df, expected_names=nasc_cols, path_for_df=df_path)


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

    _check_nasc_df(df, file_path)

    # obtaining those columns that are required
    df = df[
        [
            "transect_num",
            "vessel_log_start",
            "vessel_log_end",
            "latitude",
            "longitude",
            "stratum_num",
            "transect_spacing",
            "NASC",
            "haul_num",
        ]
    ].copy()

    # set data types of dataframe
    df = df.astype(
        {
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
    )

    if survey.params["survey_year"] < 2003:
        # TODO: it may be the case that we need to include lines 35-61 of
        #  EchoPro/general/load_files_parameters/get_NASC_data.m
        raise NotImplementedError(
            "Loading the NASC table for survey years less than 2003 has not been implemented!"
        )

    # set dataframe index
    df.set_index("transect_num", inplace=True)

    return df
