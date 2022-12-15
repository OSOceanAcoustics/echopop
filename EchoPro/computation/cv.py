import numpy as np
import pandas as pd
from geopy import distance
from .numba_functions import compute_jolly_hampton
from typing import Tuple

"""
Contains functions for obtaining all data necessary 
for CV analysis and performs CV analysis for both 
data that is Kriged and data that is not Kriged.
"""


def get_transect_strata_info_no_kriging(lat_inpfc: Tuple[float],
                                        biomass_table: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Computes transect and stratification information necessary
    for running the Jolly-Hampton algorithm for data that
    has not been Kriged.

    Parameters
    ----------
    lat_inpfc : Tuple[float]
        Bin values which represent the latitude bounds for
        each region within a survey (established by INPFC)
    biomass_table : pd.DataFrame
        DataFrame containing Longitude, Latitude, Spacing, and
        biomass_density_adult columns

    Returns
    -------
    transect_info : pd.DataFrame
        Transect information needed by the JH algorithm
    strata_info : pd.DataFrame
        Stratum information needed by the JH algorithm
    """

    # compute transect values needed for distance calculation
    transect_info = pd.DataFrame(index=biomass_table.index.unique())
    transect_info["max_longitude"] = biomass_table['longitude'].groupby(level=0).max()
    transect_info["min_longitude"] = biomass_table['longitude'].groupby(level=0).min()
    transect_info["mean_latitude"] = biomass_table['latitude'].groupby(level=0).mean()
    transect_info["mean_spacing"] = biomass_table['transect_spacing'].groupby(level=0).mean()

    # store the sum of the biomass for each transect
    transect_info["biomass"] = biomass_table['biomass_density_adult'].groupby(level=0).sum()

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

    return transect_info, strata_info


def get_transect_strata_info_kriged(lat_inpfc: Tuple[float],
                                    biomass_table: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Computes transect and stratification information necessary
    for running the Jolly-Hampton algorithm for data that
    has been Kriged.

    Parameters
    ----------
    lat_inpfc : Tuple[float]
        Bin values which represent the latitude bounds for
        each region within a survey (established by INPFC)
    biomass_table : pd.DataFrame
        DataFrame containing Latitude of centroid,
        Longitude of centroid, and biomass columns

    Returns
    -------
    transect_info : pd.DataFrame
        Transect information needed by the JH algorithm
    strata_info : pd.DataFrame
        Stratum information needed by the JH algorithm
    """

    # reduce biomass table to only essential columns
    reduced_table = biomass_table[["centroid_latitude",
                                   "centroid_longitude",
                                   "biomass"]].copy()

    # number of "virtual transects" within a latitude degree
    n_transect_per_lat = 5  # TODO: make this an input

    # latitude array with equal increment
    reduced_table["lat_eq_inc"] = np.round(
        reduced_table["centroid_latitude"] * n_transect_per_lat + 0.5) / n_transect_per_lat

    reduced_table.set_index("lat_eq_inc", inplace=True)

    # unique equal-spacing transects
    uniq_lat_eq_inc = np.unique(reduced_table.index.values)

    # compute transect values needed for distance calculation
    transect_info = pd.DataFrame(index=uniq_lat_eq_inc, dtype=np.float64)

    # store the sum of the biomass for each transect
    transect_info['biomass'] = reduced_table['biomass'].groupby(level='lat_eq_inc').sum()

    # store max and min of the longitude
    transect_info["max_longitude"] = reduced_table['centroid_longitude'].groupby(level=0).max()
    transect_info["min_longitude"] = reduced_table['centroid_longitude'].groupby(level=0).min()

    # compute and store the length (in nmi) of each transect
    transect_info["distance"] = transect_info.apply(
        lambda x: distance.distance(
            (x.name, x['min_longitude']),
            (x.name, x['max_longitude'])).nm,
        axis=1
    )

    # compute differences in unique latitudes
    uniq_lat_eq_inc_diff = np.diff(transect_info.index)
    mean_diff = np.mean(uniq_lat_eq_inc_diff)

    # compute the spacing between unique latitudes as nmi
    transect_info["spacing"] = np.nan
    transect_info["spacing"].iloc[1:-1] = uniq_lat_eq_inc_diff[:-1] * (60.0 / 2.0)
    transect_info["spacing"].iloc[0] = mean_diff * 60.0
    transect_info["spacing"].iloc[-1] = mean_diff * 60.0

    # compute the area (with units nmi^2) covered by each transect
    transect_info["area"] = transect_info.apply(lambda x: x["distance"] * x["spacing"], axis=1)

    # bin the mean latitude using lat_inpfc, each bin represents a stratum
    temp_df = pd.DataFrame(index=uniq_lat_eq_inc, dtype=np.float64)
    temp_df["unique_lat"] = transect_info.index.values
    strata_lat_stratum = pd.cut(temp_df["unique_lat"],
                                lat_inpfc, right=False,
                                labels=range(len(lat_inpfc) - 1)).rename('lat_INPFC_stratum')

    # include the last lat_inpfc right interval
    strata_lat_stratum.loc[lat_inpfc[-1]] = strata_lat_stratum.cat.categories.max()

    # create important strata information
    strata_info = transect_info["area"].groupby(strata_lat_stratum).agg(['count', 'sum'])
    strata_info = strata_info.rename(columns={'count': 'num_transects', 'sum': 'total_transect_area'})

    return transect_info, strata_info


def run_jolly_hampton(survey, nr: int, lat_inpfc: Tuple[float],
                      seed: int = None, kriged_data: bool = False) -> float:
    """
    Runs the Jolly-Hampton algorithm and computes
    the mean CV value over the given realizations.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object.
    nr : int
        The number of realizations to perform
    lat_inpfc : Tuple[float]
        Bin values which represent the latitude bounds for
        each region within a survey (established by INPFC)
    seed : int
        Seed value for the random number generator
    kriged_data : bool
        If True, perform CV analysis on Kriged data, otherwise
        perform CV analysis on data that has not been Kriged

    Returns
    -------
    The mean Jolly-Hampton CV value.

    Notes
    -----
    The format of ``lat_inpfc`` should be such that it can be
    used by Pandas.cut.
    """

    if kriged_data:
        transect_info, strata_info = get_transect_strata_info_kriged(lat_inpfc,
                                                                     survey.bio_calc.kriging_results_gdf)
    else:
        transect_info, strata_info = get_transect_strata_info_no_kriging(lat_inpfc,
                                                                         survey.bio_calc.transect_results_gdf)

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

    cv_jh_mean = compute_jolly_hampton(nr, survey.params["JH_fac"], num_transects, s_e_ind,
                                       transect_distances, field, total_transect_area,
                                       seed_val=seed)

    return cv_jh_mean



