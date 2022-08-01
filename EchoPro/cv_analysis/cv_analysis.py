import numpy as np
import pandas as pd
from geopy import distance
from ..numba_modules import compute_jolly_hampton


"""
Contains functions for obtaining all data necessary 
for CV analysis and performs CV analysis for both 
data that is Kriged and data that is not Kriged.
"""


def get_transect_strata_info(lat_inpfc, biomass_table):
    """
    Computes

    Parameters
    ----------
    lat_inpfc
    biomass_table

    Returns
    -------

    """

    transect_info = pd.DataFrame(index=biomass_table.index.unique())
    transect_info["max_longitude"] = biomass_table['Longitude'].groupby(level=0).max()
    transect_info["min_longitude"] = biomass_table['Longitude'].groupby(level=0).min()
    transect_info["mean_latitude"] = biomass_table['Latitude'].groupby(level=0).mean()
    transect_info["mean_spacing"] = biomass_table['Spacing'].groupby(level=0).mean()
    transect_info["biomass"] = biomass_table['normalized_biomass_density'].groupby(level=0).sum()

    transect_info["distance"] = transect_info.apply(
        lambda x: distance.distance(
            (x['mean_latitude'], x['min_longitude']),
            (x['mean_latitude'], x['max_longitude'])).nm,
        axis=1
    )

    transect_info["area"] = transect_info.apply(lambda x: x["distance"] * x["mean_spacing"], axis=1)

    # bin the mean latitude using lat_INPFC, each bin represents a stratum
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
    # transect_info
    strata_info["rhom"] = np.nan
    strata_info["biomass"] = np.nan
    strata_info["var_rhom"] = np.nan

    return transect_info, strata_info


def run_jolly_hampton(nr, lat_INPFC, biomass_table, jh_fac, seed=None):
    """
    Runs the Jolly Hampton algorithm to compute
    the CV value.

    Parameters
    ----------
    nr : int
        The number of realizations to perform
    lat_INPFC : list
        List of bins to be used by Pandas.cut
    biomass_table : Dataframe
        Table values used to compute the data needed
        for CV analysis
    jh_fac : float
        Percent of transects to select in each stratum
        for the Jolly Hampton Algorithm
    seed : int
        Seed value for the random number generator
    """

    transect_info, strata_info = get_transect_strata_info(lat_INPFC, biomass_table)

    # get numpy form of dataframe values so we can use Numba
    distance = transect_info['distance'].values.flatten()
    field = transect_info['biomass'].values.flatten()
    num_transects = strata_info['num_transects'].values.flatten()
    strata_nums = strata_info.index.values.to_numpy()
    total_transect_area = strata_info["total_transect_area"].values.flatten()

    # construct start and end indices for the strata
    start_end = [[np.sum(num_transects[:i]), np.sum(num_transects[:i]) + num_transects[i]] for i in
                 range(len(strata_nums))]
    s_e_ind = np.array(start_end)

    cv_jh_mean = compute_jolly_hampton(nr, jh_fac, num_transects, s_e_ind,
                                       distance, field, total_transect_area,
                                       seed_val=seed)

    return cv_jh_mean



