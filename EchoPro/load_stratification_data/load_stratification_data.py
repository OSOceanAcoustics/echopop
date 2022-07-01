import numpy as np
import pandas as pd
import xarray as xr


class LoadStrataData:
    """
    A Class that loads and processes all
    stratification data

    Parameters
    ----------
    stratification_index : int
        Index for the chosen stratification:
        0 = INPFC strata
        1 = KS (trawl)-based
        2-6 = geographically based but close to trawl-based stratification
        7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
        10 = one stratum for the whole survey
    root_dir : str
        The root directory that contains all stratification data
    filename_strata : str
        The name of the file used to construct strata_df/ds
    stratification_filename : str
        The name of the file used to construct geo_strata_df


    """

    def __init__(self, stratification_index, root_dir,
                 filename_strata, stratification_filename):

        self.stratification_index = stratification_index
        self.root_dir = root_dir
        self.filename_strata = filename_strata
        self.stratification_filename = stratification_filename

    def _load_stratafication_file(self):
        """
        Loads stratification file. Produces the variable self.strata_df,
        which is a Pandas dataframe representing stratification file.
        """
        # TODO: we should make stratification_index the sheet name

        # load stratification file
        if self.stratification_index != 1 and self.stratification_index != 0:
            raise NotImplementedError(f"stratification_index of {self.stratification_index} has not been implemented!")
        else:

            if self.stratification_index == 1:
                strata_df = pd.read_excel(self.root_dir + self.filename_strata, sheet_name='Base KS')
                strata_df = strata_df[['Year', 'Cluster name', 'Haul', 'wt']].copy()

                # set data types of dataframe
                strata_df = strata_df.astype({'Year': int, 'Cluster name': int,
                                              'Haul': int, 'wt': np.float64})

                strata_df.rename(columns={'Cluster name': 'strata'}, inplace=True)
                strata_df.set_index(['Haul', 'strata'], inplace=True)
                strata_df.sort_index(inplace=True)

            else:
                strata_df = pd.read_excel(self.root_dir + self.filename_strata, sheet_name='INPFC')
                strata_df = strata_df[['Year', 'INPFC', 'Haul', 'wt']].copy()

                # set data types of dataframe
                strata_df = strata_df.astype({'Year': int, 'INPFC': int,
                                              'Haul': int, 'wt': np.float64})

                strata_df.rename(columns={'INPFC': 'strata'}, inplace=True)
                strata_df.set_index(['Haul', 'strata'], inplace=True)
                strata_df.sort_index(inplace=True)

            return strata_df

    def _load_geographic_stratification(self):

        # TODO: we should make stratification_index the sheet name

        # load geographic stratification file
        if self.stratification_index != 1 and self.stratification_index != 0:
            raise NotImplementedError(f"stratification_index of {self.stratification_index} has not been implemented!")
        else:
            if self.stratification_index == 1:
                geo_strata_df = pd.read_excel(self.root_dir + self.stratification_filename, sheet_name='stratification1')
                geo_strata_df = geo_strata_df[['Strata index', 'Latitude (upper limit)']].copy()

                # set data types of dataframe
                geo_strata_df = geo_strata_df.astype({'Strata index': int,
                                                      'Latitude (upper limit)': np.float64})
            else:
                geo_strata_df = pd.read_excel(self.root_dir + self.stratification_filename, sheet_name='INPFC')
                geo_strata_df = geo_strata_df[['Strata index', 'Latitude (upper limit)']].copy()

                # set data types of dataframe
                geo_strata_df = geo_strata_df.astype({'Strata index': int,
                                                      'Latitude (upper limit)': np.float64})

        return geo_strata_df

    @staticmethod
    def _get_strata_ds(strata_df, specimen_df, length_df):
        """
        #TODO: fill in this documentation!

        Parameters
        ----------
        strata_df
        specimen_df
        length_df

        Returns
        -------

        """

        # TODO: do we need a DataSet or can all of this be put into a Dataframe?
        strata_df["length_average_haul"] = np.nan
        strata_df["TS_lin_haul"] = np.nan
        strata_df["sig_bs_haul"] = np.nan

        # select the indices that do not have nan in either Length or Weight
        spec_df = specimen_df[['Length', 'Weight']].copy()  # TODO: does it make sense to exclude based on wgt?
        spec_df = spec_df.dropna(how='any')

        for haul_num in spec_df.index.unique():

            all_len = spec_df.loc[haul_num]['Length']

            # TODO: replace with length_ds?
            if haul_num in length_df.index:
                all_len = np.concatenate([all_len,
                                          length_df.loc[haul_num]['Length']])

            # TODO: these two functions are specific to Hake, replace with input in the future
            TS0j = 20.0 * np.log10(all_len) - 68.0
            TSj = 10.0 * np.log10(np.nanmean(10.0 ** (TS0j / 10.0)))
            strata_df.loc[haul_num, "TS_lin_haul"] = TSj
            strata_df.loc[haul_num, "sig_bs_haul"] = 10.0 ** (TSj / 10.0)
            strata_df.loc[haul_num, "length_average_haul"] = np.nanmean(all_len)

        strata_ds = strata_df.to_xarray()

        strata_ds["sig_bs"] = strata_ds.sig_bs_haul.mean(dim="Haul", skipna=True)
        strata_ds["sig_b"] = 4.0 * np.pi * strata_ds["sig_bs"]

        return strata_df, strata_ds

    def get_strata_data(self, specimen_df, length_df):
        """
        Obtains the stratification data needed
        to conduct core routines of EchoPro.

        Returns
        -------
        # TODO: fill in this documentation!!
        """

        strata_df = self._load_stratafication_file()

        geo_strata_df = self._load_geographic_stratification()

        strata_df, strata_ds = self._get_strata_ds(strata_df, specimen_df, length_df)

        return strata_df, strata_ds, geo_strata_df
