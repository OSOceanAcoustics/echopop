import numpy as np
import pandas as pd

from ..utils.input_checks_read import check_and_read


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
        df = check_and_read("NASC/no_age1", nasc_var_types, survey.params)
    else:
        df = check_and_read("NASC/all_ages", nasc_var_types, survey.params)

    if survey.params["survey_year"] < 2003:
        # TODO: it may be the case that we need to include lines 35-61 of
        #  EchoPro_matlab/general/load_files_parameters/get_NASC_data.m
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
