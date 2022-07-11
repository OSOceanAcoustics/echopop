import pandas as pd
import numpy as np


nasc_cols = {'Transect', 'VL start', 'VL end', 'Latitude', 'Longitude',
             'Stratum', 'Spacing', 'NASC', 'Assigned haul'}


def _check_nasc_df(nasc_df: pd.DataFrame):
    """
    Ensures that the appropriate columns are
    contained in the NASC Dataframe.

    TODO: should we add more in-depth checks here?
    """

    if not set(nasc_df.columns).intersection(nasc_cols):
        raise NameError("NASC dataframe does not contain all expected columns!")


def load_nasc_df(epro):
    """
    Load VL interval-based NASC table.

    Parameters
    ----------
    epro : EchoPro object
        An initialized EchoPro object

    Returns
    -------
    Pandas Dataframe of NASC table.
    """

    # select and check the appropriate nasc data file
    if epro.params['exclude_age1']:
        df = pd.read_excel(epro.params['data_root_dir'] + epro.params['nasc_no_age1_filename'],
                           sheet_name=epro.params['nasc_no_age1_sheetname'])
    else:
        df = pd.read_excel(epro.params['data_root_dir'] + epro.params['nasc_all_ages_filename'],
                           sheet_name=epro.params['nasc_all_ages_sheetname'])
    _check_nasc_df(df)

    # obtaining those columns that are required
    df = df[['Transect', 'VL start', 'VL end', 'Latitude', 'Longitude', 'Stratum', 'Spacing',
            'NASC', 'Assigned haul']].copy()

    # set data types of dataframe
    df = df.astype({'Transect': int, 'VL start': np.float64, 'VL end': np.float64,
                    'Latitude': np.float64, 'Longitude': np.float64, 'Stratum': int,
                    'Spacing': np.float64, 'NASC': np.float64, 'Assigned haul': int})

    # rename column TODO: in the future require Haul as the column name
    df.rename(columns={'Assigned haul': 'Haul'}, inplace=True)

    if epro.params['survey_year'] < 2003:
        # TODO: it may be the case that we need to include lines 35-61 of
        #  EchoPro/general/load_files_parameters/get_NASC_data.m
        raise NotImplementedError("Loading the NASC table for survey years less than 2003 has not been implemented!")

    # set dataframe index
    df.set_index('Transect', inplace=True)

    return df
