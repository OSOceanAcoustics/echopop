import numpy as np
import numba as nb
import pandas as pd
import warnings
import sys
from geopy import distance
from numpy.random import default_rng


class CVAnalysis:
    """
    This class constructs all data necessary
    for CV analysis and performs CV analysis
    for both data that is kriged and data that
    is not kriged.

    Parameters
    ----------
    EPro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.EPro will also change this object.
    """

    def __init__(self, EPro = None):

        self.EPro = EPro

    def get_transect_strata_info(self, lat_INPFC, biomass_table):

        transect_info = pd.DataFrame(index=biomass_table.index.unique())
        transect_info["max_longitude"] = biomass_table['Longitude'].groupby(level=0).max()
        transect_info["min_longitude"] = biomass_table['Longitude'].groupby(level=0).min()
        transect_info["mean_latitude"] = biomass_table['Latitude'].groupby(level=0).mean()
        transect_info["mean_spacing"] = biomass_table['Spacing'].groupby(level=0).mean()
        transect_info["biomass"] = biomass_table['nwgt_total'].groupby(level=0).sum()

        transect_info["distance"] = transect_info.apply(
            lambda x: distance.distance(
                (x['mean_latitude'], x['min_longitude']),
                (x['mean_latitude'], x['max_longitude'])).nm,
            axis=1
        )

        transect_info["area"] = transect_info.apply(lambda x: x["distance"] * x["mean_spacing"], axis=1)

        # bin the mean latitude using lat_INPFC, each bin represents a stratum
        strata_lat_stratum = pd.cut(transect_info['mean_latitude'],
                                    lat_INPFC,
                                    labels=range(len(lat_INPFC) - 1),
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

    def compute_random_strata_info(self, transect_info, strata_info, rng):
        # TODO: we need to improve the speed of this computation!
        #  Maybe replace pandas with numpy here

        for stratum in strata_info.index:

            # randomly select transects in the stratum
            num_ind = round(self.EPro.params["JH_fac"] * strata_info.loc[stratum]['num_transects'])
            sel_ind = rng.choice(transect_info.loc[stratum].index, num_ind, replace=False)

            transect_stratum_info = transect_info.loc[(stratum, sel_ind), ['distance', 'biomass']]
            # transect_stratum_info = transect_info.loc[stratum].iloc[:8][['distance', 'biomass']]

            # transect-length weighting factor of the transects in the stratum
            wgt = transect_stratum_info['distance'] / transect_stratum_info['distance'].mean()

            # normalized biomass of the transects in the stratum
            rhom_trans_stratum = transect_stratum_info['biomass'] / transect_stratum_info['distance']

            # transect-length-normalized mean density of the stratum
            strata_info.loc[stratum, "rhom"] = (transect_stratum_info['biomass'] * transect_stratum_info[
                'distance']).sum() / transect_stratum_info[
                                                   'distance'].sum()

            # variance of the transect-length weighted biomass within the stratum
            if num_ind != 1:
                strata_info.loc[stratum, "var_rhom"] = (wgt ** 2 * (
                            rhom_trans_stratum - strata_info["rhom"].loc[stratum]) ** 2).sum() / (
                                                                   num_ind * (num_ind - 1))
            else:
                strata_info.loc[stratum, "var_rhom"] = (wgt ** 2 * (
                            rhom_trans_stratum - strata_info["rhom"].loc[stratum]) ** 2).sum() / (num_ind * num_ind)

            strata_info.loc[stratum, "biomass"] = transect_stratum_info['biomass'].sum()

        # area weighted variance of the "transect-length weighted biomass"
        CV = np.sqrt((strata_info["var_rhom"] * strata_info["total_transect_area"] ** 2).sum()) / (
                    strata_info["total_transect_area"] * strata_info["rhom"]).sum()

        return CV

    def run_jolly_hampton(self, nr, lat_INPFC, biomass_table):
        """
        Runs the Jolly Hampton algorithm to compute
        the CV

        Parameters
        ----------
        nr : int
            The number of realizations to perform
        lat_INPFC : list
            List of bins to be used by Pandas.cut
        biomass_table : Dataframe
            Table values used to compute the data needed
            for CV analysis.
        """

        rng = default_rng()  # TODO: make this an input of the function

        transect_info, strata_info = self.get_transect_strata_info(lat_INPFC, biomass_table)

        CV_JH_vals = np.full((nr, ), np.nan)

        # TODO: this can be parallelized ... if we want to/ need to
        # for i in range(nr):   # loop through realizations
        #
        #     CV_JH_vals[i] = self.compute_random_strata_info(transect_info, strata_info, rng)
        #
        #     print(CV_JH_vals[i])
        #     import sys
        #     sys.exit()
        #
        # CV_JH_mean = np.nanmean(CV_JH_vals)
        #
        # print(f"CV_JH_mean = {CV_JH_mean}")
        #
        # return CV_JH_mean

        return transect_info, strata_info, rng


