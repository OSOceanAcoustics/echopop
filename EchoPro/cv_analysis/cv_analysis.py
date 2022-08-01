import numpy as np
import pandas as pd
from geopy import distance
from ..numba_modules import compute_jolly_hampton
from typing import Tuple

"""
Contains functions for obtaining all data necessary 
for CV analysis and performs CV analysis for both 
data that is Kriged and data that is not Kriged.
"""


def get_transect_strata_info(lat_inpfc: Tuple[float],
                             biomass_table: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Computes transect and stratification information necessary
    for running the Jolly-Hampton algorithm.

    Parameters
    ----------
    lat_inpfc : Tuple[float]
        Bin values which represent the latitude bounds for
        each region within a survey (established by INPFC)
    biomass_table : pd.Dataframe
        DataFrame containing Longitude, Latitude, Spacing, and
        normalized_biomass_density columns

    Returns
    -------
    transect_info : pd.DataFrame
        Transect information needed by the JH algorithm
    strata_info : pd.DataFrame
        Stratum information needed by the JH algorithm
    """

    # compute transect values needed for distance calculation
    transect_info = pd.DataFrame(index=biomass_table.index.unique())
    transect_info["max_longitude"] = biomass_table['Longitude'].groupby(level=0).max()
    transect_info["min_longitude"] = biomass_table['Longitude'].groupby(level=0).min()
    transect_info["mean_latitude"] = biomass_table['Latitude'].groupby(level=0).mean()
    transect_info["mean_spacing"] = biomass_table['Spacing'].groupby(level=0).mean()

    # store the sum of the biomass for each transect
    transect_info["biomass"] = biomass_table['normalized_biomass_density'].groupby(level=0).sum()

    # compute the length of each transect
    transect_info["distance"] = transect_info.apply(
        lambda x: distance.distance(
            (x['mean_latitude'], x['min_longitude']),
            (x['mean_latitude'], x['max_longitude'])).nm,
        axis=1
    )

    # compute the area covered by each transect
    transect_info["area"] = transect_info.apply(lambda x: x["distance"] * x["mean_spacing"], axis=1)

    # bin the mean latitude using lat_inpfc, each bin represents a stratum
    strata_lat_stratum = pd.cut(transect_info['mean_latitude'],
                                lat_inpfc,
                                labels=range(len(lat_inpfc) - 1),
                                right=False).rename('lat_INPFC_stratum')

    strata_info = transect_info["area"].groupby(strata_lat_stratum).agg(['count', 'sum'])
    strata_info = strata_info.rename(columns={'count': 'num_transects', 'sum': 'total_transect_area'})

    # add this binning as a column
    transect_info['lat_INPFC_stratum'] = strata_lat_stratum

    # change index to lat_INPFC_bin
    transect_info = transect_info.reset_index().set_index(['lat_INPFC_stratum', 'Transect'])

    # create empty columns that will hold values in Jolly-Hampton algorithm
    strata_info["rhom"] = np.nan
    strata_info["biomass"] = np.nan
    strata_info["var_rhom"] = np.nan

    return transect_info, strata_info


def run_jolly_hampton(nr: int, lat_inpfc: Tuple[float],
                      biomass_table: pd.DataFrame, jh_fac: float,
                      seed: int = None) -> float:
    """
    Runs the Jolly-Hampton algorithm and computes
    the mean CV value over the given realizations.

    Parameters
    ----------
    nr : int
        The number of realizations to perform
    lat_inpfc : Tuple[float]
        Bin values which represent the latitude bounds for
        each region within a survey (established by INPFC)
    biomass_table : pd.Dataframe
        DataFrame containing Longitude, Latitude, Spacing, and
        normalized_biomass_density columns
    jh_fac : float
        Percent of transects to select in each stratum
        for the Jolly-Hampton algorithm
    seed : int
        Seed value for the random number generator

    Returns
    -------
    The mean Jolly-Hampton CV value.

    Notes
    -----
    The format of ``lat_inpfc`` should be such that it can be
    used by Pandas.cut.
    """

    transect_info, strata_info = get_transect_strata_info(lat_inpfc, biomass_table)

    # get numpy form of dataframe values, so we can use Numba
    transect_distances = transect_info['distance'].values.flatten()
    field = transect_info['biomass'].values.flatten()
    num_transects = strata_info['num_transects'].values.flatten()
    strata_nums = strata_info.index.values.to_numpy()
    total_transect_area = strata_info["total_transect_area"].values.flatten()

    # construct start and end indices for the strata
    start_end = [[np.sum(num_transects[:i]), np.sum(num_transects[:i]) + num_transects[i]] for i in
                 range(len(strata_nums))]
    s_e_ind = np.array(start_end)

    cv_jh_mean = compute_jolly_hampton(nr, jh_fac, num_transects, s_e_ind,
                                       transect_distances, field, total_transect_area,
                                       seed_val=seed)

    return cv_jh_mean



