import numpy as np
import pandas as pd
import xarray as xr
import warnings
import sys
from ..load_biological_data import LoadBioData


class LoadStrataData:
    """
    A Class that loads and processes all
    stratification data

    Parameters
    ----------
    EPro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.EPro will also change this object.
    """

    def __init__(self, EPro = None):

        self.EPro = EPro

    def load_stratafication_file(self, stratification_index):
        """
        Loads stratification file. Produces the variable self.strata_df,
        which is a Pandas dataframe representing stratification file.

        Parameters
        ----------
        stratification_index : int  # TODO: this looks to mirror KS_stratification, it seems to be unnecessary to use this
                                    # TODO: Ask Chu if it is ok to remove this.
            Index for the chosen stratification
            0 = INPFC strata
            1 = KS (trawl)-based
            2-6 = geographically based but close to trawl-based stratification
            7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
            10 = one stratum for the whole survey
        """
        # TODO: we should make stratification_index the sheet name

        # load stratification file
        if stratification_index != 1 and stratification_index != 0:
            raise NotImplementedError(f"stratification_index of {stratification_index} has not been implemented!")
        else:

            if stratification_index == 1:
                strata_df = pd.read_excel(self.EPro.params['data_root_dir'] +
                                                    self.EPro.params['filename_strata'],
                                                    sheet_name='Base KS')
                strata_df = strata_df[['Year', 'Cluster name', 'Haul', 'wt']].copy()

                # set data types of dataframe
                strata_df = strata_df.astype({'Year': int,
                                                                  'Cluster name': int,
                                                                  'Haul': int,
                                                                  'wt': np.float64})

                strata_df.rename(columns={'Cluster name': 'strata'}, inplace=True)
                strata_df.set_index(['Haul', 'strata'], inplace=True)
                strata_df.sort_index(inplace=True)

            else:
                strata_df = pd.read_excel(self.EPro.params['data_root_dir'] +
                                                    self.EPro.params['filename_strata'],
                                                    sheet_name='INPFC')
                strata_df = strata_df[['Year', 'INPFC', 'Haul', 'wt']].copy()

                # set data types of dataframe
                strata_df = strata_df.astype({'Year': int,
                                                                  'INPFC': int,
                                                                  'Haul': int,
                                                                  'wt': np.float64})

                strata_df.rename(columns={'INPFC': 'strata'}, inplace=True)
                strata_df.set_index(['Haul', 'strata'], inplace=True)
                strata_df.sort_index(inplace=True)

            self.EPro.strata_df = strata_df

    def load_geographic_stratification(self, stratification_index):

        # TODO: we should make stratification_index the sheet name

        # load geographic stratification file
        if stratification_index != 1 and stratification_index != 0:
            raise NotImplementedError(f"stratification_index of {stratification_index} has not been implemented!")
        else:
            if stratification_index == 1:
                geo_strata_df = pd.read_excel(self.EPro.params['data_root_dir'] +
                                          self.EPro.params['stratification_filename'],
                                          sheet_name='stratification1')
                geo_strata_df = geo_strata_df[['Strata index', 'Latitude (upper limit)']].copy()

                # set data types of dataframe
                geo_strata_df = geo_strata_df.astype({'Strata index': int,
                                                      'Latitude (upper limit)': np.float64})
            else:
                geo_strata_df = pd.read_excel(self.EPro.params['data_root_dir'] +
                                          self.EPro.params['stratification_filename'],
                                          sheet_name='INPFC')
                geo_strata_df = geo_strata_df[['Strata index', 'Latitude (upper limit)']].copy()

                # set data types of dataframe
                geo_strata_df = geo_strata_df.astype({'Strata index': int,
                                                      'Latitude (upper limit)': np.float64})

        self.geo_strata_df = geo_strata_df

    def __check_for_empty_strata(self):
        """
        A function that fills in empty strata when necessary.
        """

        if self.EPro.bio_strata.empty:

            raise NotImplementedError("Filling in of empty strata has not been implemented!")

        # TODO: Do we have to fill-in trawl information for stratum number < first non-empty stratum?
        #  Matlab get_historical_strata_data lines 144-169
        #  These lines are used when no data is present. Originally created to account for data from
        #  Fisheries data (outside of FEAT group). It is known that this creates some error/uncertainty
        #  of the biomass estimate.

        # TODO: Do we have to interpolate trawl information to stranum (no tralws) by using the parameters
        #  of the adjacient strata (stratum)
        #  Matlab get_historical_strata_data lines 172-194
        #  Say we have strata 1 and 3, we create strata 2 using the average of 1 and 3. This is for
        #  other fishery data not the FEAT data

        # TODO: Do we have to extrpolate trawl information to stranum (no tralws) > last non-empty stratum
        #  (situation in using the Observer data or A-SHOP data)
        #  Matlab get_historical_strata_data lines 195-216
        #  Say we want 6 strata (INPFC), but we have data for strata 2, 4, 5. Then this code copies
        #  stratum 2 to stratum 1, and average strata 2 and 4 to get stratum 3, then copy stratum 5 to stratum 6.

    # def __get_bio_strata(self, stratification_index, KS_stratification, transect_reduction_fraction):
    #     """
    #     A function that obtains a subset of the strata data corresponding
    #     to the biological data. This information is stored in the
    #     Pandas Dataframe ``self.bio_strata``.
    #     """
    #
    #     if stratification_index == 1:
    #
    #         if self.EPro.params['CAN_strata_num0']:
    #             raise NotImplementedError("Filled CAN_strata_num0 value has not been implemented!")
    #
    #         if self.EPro.params['exclude_age1']:
    #
    #             # TODO: Do we need to have a condition to check for para.proc.age1_haul?
    #             # TODO: According to Chu, this seems to be unused.
    #             # raise NotImplementedError(f"age1_haul")
    #
    #             if transect_reduction_fraction != 0.0:
    #                 raise NotImplementedError("transect reduction has not been implemented!")
    #             else:
    #
    #                 # select only stratum not equal to zero
    #                 strata_mask = self.EPro.strata_df['strata'] != 0
    #                 self.EPro.bio_strata = self.EPro.strata_df[strata_mask][
    #                     ['strata', 'wt']].reset_index().set_index('Cluster name')
    #     else:
    #
    #         if self.EPro.params['CAN_strata_num0']:
    #             raise NotImplementedError("Filled CAN_strata_num0 value has not been implemented!")
    #
    #         if self.EPro.params['exclude_age1']:
    #
    #             # TODO: Do we need to have a condition to check for para.proc.age1_haul?
    #             # TODO: According to Chu, this seems to be unused.
    #             # raise NotImplementedError(f"age1_haul")
    #
    #             if transect_reduction_fraction != 0.0:
    #                 raise NotImplementedError("transect reduction has not been implemented!")
    #             else:
    #
    #                 # select only stratum not equal to zero
    #                 strata_mask = self.EPro.strata_df['strata'] != 0
    #                 self.EPro.bio_strata = self.EPro.strata_df[strata_mask][['strata', 'wt']].reset_index().set_index('strata')
    #
    #     # TODO: investigate this further! Currently no techniques for filling in the strata have been implemented.
    #     self.__check_for_empty_strata()

    def __get_expanded_data(self, ind):
        """
        Collect the expanded length data from length_ds for
        male, female, and unsexed genders.

        Parameters
        ----------
        ind : int
            index of haul to select data from
        """
        # Collect expanded length data from length_ds
        if ind in self.EPro.length_ds.Haul:
            male_val_len = np.repeat(self.EPro.length_ds.Length.values,
                                     np.nan_to_num(
                                         self.EPro.length_ds.Frequency.sel(Haul=ind, Sex=1).values).astype(dtype=int),
                                     axis=0)

            female_val_len = np.repeat(self.EPro.length_ds.Length.values,
                                       np.nan_to_num(
                                           self.EPro.length_ds.Frequency.sel(Haul=ind, Sex=2).values).astype(dtype=int),
                                       axis=0)

            unsexed_ind = np.logical_and(self.EPro.length_ds.Sex.values != 1, self.EPro.length_ds.Sex.values != 2)
            if np.sum(unsexed_ind) > 1:
                raise NotImplementedError("Unsexed values greater than 1 have not been implemented!")
            else:
                unsexed_val_len = np.repeat(self.EPro.length_ds.Length.values,
                                            np.nan_to_num(
                                                self.EPro.length_ds.Frequency.sel(Haul=ind,
                                                                                  Sex=3).values).astype(dtype=int),
                                            axis=0)
        else:
            male_val_len = []
            female_val_len = []
            unsexed_val_len = []

        return male_val_len, female_val_len, unsexed_val_len

    def __get_len_data_specimen(self, gen, ind):
        """
        Collect length data from specimen_df for a specific gender.

        Parameters
        ----------
        gen : float
            Gender
        ind : int
            index of haul to select data from
        """
        selected_ind = (self.EPro.specimen_df['Sex'] == gen) & pd.notnull(self.EPro.specimen_df['Age'])
        if ind in selected_ind[selected_ind].index:
            val_specimen = self.EPro.specimen_df[selected_ind].loc[ind]['Length']

            if isinstance(val_specimen, np.float64):
                val_specimen = [val_specimen]
            else:
                val_specimen = val_specimen.values
        else:
            val_specimen = []

        return val_specimen

    def __get_bio_len_age_haul(self):

        # construct hake-haul number array
        self.EPro.params['bio_hake_trawl_num'] = np.union1d(self.EPro.length_ds.Haul.values,
                                                       self.EPro.specimen_df.index.unique().values)

        # initialize parameters with a numpy array
        len_bio_hake_len_bin = len(self.EPro.params['bio_hake_len_bin'])
        len_bio_hake_trawl_num = len(self.EPro.params['bio_hake_trawl_num'])

        # TODO: uncomment below if we want unsexed to be included
        self.EPro.params['bio_len_haul_U'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.EPro.params['bio_aged_len_haul_U'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))

        self.EPro.params['bio_len_haul_M'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.EPro.params['bio_len_haul_F'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.EPro.params['bio_aged_len_haul_M'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.EPro.params['bio_aged_len_haul_F'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))

        j = 0
        for i in self.EPro.params['bio_hake_trawl_num']:

            male_val_len, female_val_len, unsexed_val_len = self.__get_expanded_data(i)

            # Collect length data from specimen_df for males
            male_val_specimen = self.__get_len_data_specimen(1.0, i)

            # Collect length data from specimen_df for females
            female_val_specimen = self.__get_len_data_specimen(2.0, i)

            unsexed_ind = np.logical_and(self.EPro.length_ds.Sex.values != 1, self.EPro.length_ds.Sex.values != 2)
            if np.sum(unsexed_ind) > 1:
                raise NotImplementedError("Unsexed values greater than 1 have not been implemented!")
            else:
                # Collect length data from specimen_df for unsexed
                unsexed_val_specimen = self.__get_len_data_specimen(3.0, i)

            # store length data from length_ds and specimen_df
            self.EPro.params['bio_len_haul_M'][:, j] = LoadBioData.get_bin_counts(
                np.concatenate([male_val_len, male_val_specimen]), self.EPro.params['bio_hake_len_bin'])

            self.EPro.params['bio_len_haul_F'][:, j] = LoadBioData.get_bin_counts(
                np.concatenate([female_val_len, female_val_specimen]), self.EPro.params['bio_hake_len_bin'])

            self.EPro.params['bio_len_haul_U'][:, j] = LoadBioData.get_bin_counts(
                np.concatenate([unsexed_val_len, unsexed_val_specimen]), self.EPro.params['bio_hake_len_bin'])

            self.EPro.params['bio_aged_len_haul_M'][:, j] = LoadBioData.get_bin_counts(male_val_specimen,
                                                                                    self.EPro.params['bio_hake_len_bin'])
            self.EPro.params['bio_aged_len_haul_F'][:, j] = LoadBioData.get_bin_counts(female_val_specimen,
                                                                                    self.EPro.params['bio_hake_len_bin'])

            self.EPro.params['bio_aged_len_haul_U'][:, j] = LoadBioData.get_bin_counts(unsexed_val_specimen,
                                                                                    self.EPro.params['bio_hake_len_bin'])
            j += 1

        self.EPro.params['bio_len_haul_ALL'] = self.EPro.params['bio_len_haul_M'] + \
                                               self.EPro.params['bio_len_haul_F'] + \
                                               self.EPro.params['bio_len_haul_U']

        self.EPro.params['bio_aged_len_haul_ALL'] = self.EPro.params['bio_aged_len_haul_M'] + self.EPro.params[
            'bio_aged_len_haul_F'] + self.EPro.params['bio_aged_len_haul_U']

    def __remove_empty_bio_len_age_haul(self):
        """
        Deletes columns in bio_len_age_haul arrays that are empty.
        """

        aged_del_ind = self.EPro.params['bio_aged_len_haul_ALL'].sum(axis=0)

    @staticmethod
    def get_bin_ind(input_data: np.array, centered_bins: np.array):
        """
        This function manually computes bin counts given ``input_data``. This
        function is computing the histogram of ``input_data`` using
        bins that are centered, rather than bins that are on the edge.
        The first value is between negative infinity and the first bin
        center plus the bin width divided by two. The last value is
        between the second to last bin center plus the bin width
        divided by two to infinity.


        Parameters
        ----------
        input_data: numpy array
            The data to create a histogram of.
        centered_bins: numpy array
            An array that specifies the bin centers.

        Returns
        -------
        hist_ind: numpy array
            The index values of input_data corresponding to the histogram

        """

        bin_diff = np.diff(centered_bins) / 2.0

        hist_ind = [np.argwhere(input_data <= centered_bins[0] + bin_diff[0]).flatten()]

        for i in range(len(centered_bins) - 2):
            # get values greater than lower bound
            g_lb = centered_bins[i] + bin_diff[i] < input_data

            # get values less than or equal to the upper bound
            le_ub = input_data <= centered_bins[i + 1] + bin_diff[i + 1]

            hist_ind.append(np.argwhere(g_lb & le_ub).flatten())

        hist_ind.append(np.argwhere(input_data > centered_bins[-2] + bin_diff[-1]).flatten())

        return hist_ind

    def get_age_key_das(self, spec_w_strata, bins_len, bins_age):
        """
        Constructs the DataArrays corresponding to the keys
        associated with age.

        Parameters
        ----------
        spec_w_strata : Pandas Dataframe
            Species Dataframe with the strata as an index and
            columns "Length", "Age", "Weight".
        bins_len : Numpy array
            1D numpy array specifies the bin centers for lengths.
        bins_age : Numpy array
            1D numpy array specifies the bin centers for ages.

        Returns
        -------
        The Xarray DataArrays: age_len_key_da, age_len_key_wgt_da,
        and age_len_key_norm_da with coordinates
        (strata, bins_len, bins_age).

        """

        unique_strata = np.sort(spec_w_strata.index.unique())

        age_len_key = np.empty((unique_strata.shape[0],
                                bins_len.shape[0], bins_age.shape[0]), dtype=int)

        age_len_key_wgt = np.empty((unique_strata.shape[0],
                                    bins_len.shape[0], bins_age.shape[0]), dtype=np.float64)

        age_len_key_norm = np.empty((unique_strata.shape[0],
                                     bins_len.shape[0], bins_age.shape[0]), dtype=np.float64)

        stratum_i = 0
        for stratum in unique_strata:

            input_arr_len = spec_w_strata.loc[stratum].Length.values
            input_arr_age = spec_w_strata.loc[stratum].Age.values
            input_arr_wgt = spec_w_strata.loc[stratum].Weight.values

            age_bins_ind = self.get_bin_ind(input_arr_age, bins_age)

            i = 0
            for arr_ind in age_bins_ind:
                temp = self.get_bin_ind(input_arr_len[arr_ind], bins_len)
                age_len_key[stratum_i, :, i] = np.array([i.shape[0] for i in temp])

                age_len_key_wgt[stratum_i, :, i] = np.array([np.sum(input_arr_wgt[arr_ind[i]]) for i in temp])
                i += 1

            # normalized age_len_ky
            age_len_key_norm[stratum_i, :, :] = age_len_key[stratum_i, :, :] / np.sum(age_len_key[stratum_i, :, :])

            stratum_i += 1

        age_len_key_da = xr.DataArray(data=age_len_key,
                                      coords={'strata': unique_strata,
                                              'len_bins': bins_len, 'age_bins': bins_age})

        age_len_key_wgt_da = xr.DataArray(data=age_len_key_wgt,
                                          coords={'strata': unique_strata,
                                                  'len_bins': bins_len, 'age_bins': bins_age})

        age_len_key_norm_da = xr.DataArray(data=age_len_key_norm,
                                           coords={'strata': unique_strata,
                                                   'len_bins': bins_len, 'age_bins': bins_age})

        return age_len_key_da, age_len_key_wgt_da, age_len_key_norm_da

    @staticmethod
    def get_length_val_reg_vals(len_name, val_name, df: pd.DataFrame = None):
        """
        Obtains the regression values for a provided length and value.

        Parameters
        ----------
        len_name : str
            Name of column in df corresponding to the length
        val_name : str
            Name of column in df corresponding to the value
        df : Pandas Dataframe
            Dataframe containing the length and value columns


        Returns
        -------
        reg_w0 and reg_p used in the calculation
        reg_w0 * bio_hake_len_bin ** reg_p
        """

        # select the indices that do not have nan in either Length or Weight
        len_wgt_nonull = np.logical_and(df[len_name].notnull(), df[val_name].notnull())

        df_no_null = df.loc[len_wgt_nonull]

        L = df_no_null[len_name].values
        V = df_no_null[val_name].values

        # length-value regression for all trawls (male & female)
        x = np.log10(L)
        y = np.log10(V)

        p = np.polyfit(x, y, 1)  # linear regression

        reg_w0 = 10.0 ** p[1]
        reg_p = p[0]

        return reg_w0, reg_p

    def generate_length_val_key(self, bio_hake_len_bin, reg_w0, reg_p,
                                len_name, val_name, df: pd.DataFrame = None):
        """
        process length weight data (all hake trawls) to obtain
        (1) length-weight regression or length-weight-key
        (2) length-age keys

        Parameters
        ----------
        df: Pandas Dataframe
            Pandas Dataframe describing the length weight data

        """

        # select the indices that do not have nan in either Length or Weight
        len_wgt_nonull = np.logical_and(df[len_name].notnull(), df[val_name].notnull())

        df_no_null = df.loc[len_wgt_nonull]

        L = df_no_null[len_name].values
        V = df_no_null[val_name].values

        # total number of fish individuals at length specified by bio_hake_len_bin
        len_nALL = self.get_bin_ind(L, bio_hake_len_bin)
        len_nALL = np.array([i.shape[0] for i in len_nALL])

        # normalized length-key
        norm_len_key_ALL = len_nALL / sum(len_nALL)

        if (reg_p is None) or (reg_w0 is None):

            # length-value regression for all trawls (male & female)
            x = np.log10(L)
            y = np.log10(V)

            p = np.polyfit(x, y, 1)  # linear regression

            reg_w0 = 10.0 ** p[1]
            reg_p = p[0]

        # value at length or length-value-key
        # length-weight-key per fish over entire survey region (an array)
        len_val_ALL = reg_w0 * bio_hake_len_bin ** reg_p

        # create length-value structured relations
        for i in range(len(len_nALL)):

            # bins with less than 5 samples will be replaced by the regression curve
            if len_nALL[i] >= 5:
                ind = (bio_hake_len_bin[i] - 1 < L) & (
                            L <= bio_hake_len_bin[i] + 1)
                len_val_ALL[i] = np.mean(V[ind])

        return len_val_ALL, len_nALL, norm_len_key_ALL

    def get_weight_key_das(self, spec_w_strata, bins_len, reg_w0, reg_p, len_name='Length',
                           val_name='Weight'):

        unique_strata = np.sort(spec_w_strata.index.unique())

        len_wgt_key = np.empty((unique_strata.shape[0], bins_len.shape[0]))
        len_wgt_key[:, :] = np.nan

        len_key = np.empty((unique_strata.shape[0], bins_len.shape[0]))
        len_key[:, :] = np.nan

        len_key_norm = np.empty((unique_strata.shape[0], bins_len.shape[0]))
        len_key_norm[:, :] = np.nan

        stratum_i = 0
        for stratum in unique_strata:
            l_w_k, l_k, l_k_n = self.generate_length_val_key(bins_len, reg_w0, reg_p,
                                                             len_name=len_name,
                                                             val_name=val_name,
                                                             df=spec_w_strata.loc[stratum])

            len_wgt_key[stratum_i, :] = l_w_k
            len_key[stratum_i, :] = l_k
            len_key_norm[stratum_i, :] = l_k_n

            stratum_i += 1

        len_wgt_key_da = xr.DataArray(data=len_wgt_key,
                                      coords={'strata': unique_strata, 'len_bins': bins_len})

        len_key_da = xr.DataArray(data=len_key,
                                  coords={'strata': unique_strata, 'len_bins': bins_len})

        len_key_norm_da = xr.DataArray(data=len_key_norm,
                                       coords={'strata': unique_strata, 'len_bins': bins_len})

        return len_wgt_key_da, len_key_da, len_key_norm_da

    def get_len_keys_strata(self, len_strata, bins_len):

        input_data = len_strata['Length'].values
        len_ind = self.get_bin_ind(input_data, bins_len)
        len_key = np.array([i.shape[0] for i in len_ind])
        len_key_n = len_key / np.sum(len_key)

        input_data = len_strata[len_strata['Sex'] == 1]['Length'].values
        len_ind = self.get_bin_ind(input_data, bins_len)
        len_key = np.array([i.shape[0] for i in len_ind])
        len_key_M = len_key / np.sum(len_key)

        input_data = len_strata[len_strata['Sex'] == 2]['Length'].values
        len_ind = self.get_bin_ind(input_data, bins_len)
        len_key = np.array([i.shape[0] for i in len_ind])
        len_key_F = len_key / np.sum(len_key)

        return len_key_n, len_key_M, len_key_F

    def get_age_related_sums(self, spec_strata_M, spec_strata_F, bins_len, bins_age):

        input_df = pd.concat([spec_strata_M, spec_strata_F], axis=0)

        # age_len_key_da, _, age_len_key_norm_da = self.get_age_key_das(input_df, bins_len, bins_age)
        # age_len_key_M_da, _, age_len_key_norm_M_da = self.get_age_key_das(spec_strata_M, bins_len, bins_age)
        # age_len_key_F_da, _, age_len_key_norm_F_da = self.get_age_key_das(spec_strata_F, bins_len, bins_age)

        _, _, age_len_key_norm_da = self.get_age_key_das(input_df, bins_len, bins_age)
        _, _, age_len_key_norm_M_da = self.get_age_key_das(spec_strata_M, bins_len, bins_age)
        _, _, age_len_key_norm_F_da = self.get_age_key_das(spec_strata_F, bins_len, bins_age)

        len_age_key_sum = age_len_key_norm_da.isel(strata=0).sum(dim='age_bins').values

        len_age_key_M_sum = age_len_key_norm_M_da.isel(strata=0).sum(dim='age_bins').values

        len_age_key_F_sum = age_len_key_norm_F_da.isel(strata=0).sum(dim='age_bins').values

        return len_age_key_sum, len_age_key_M_sum, len_age_key_F_sum

    def get_biomass_constants(self, spec_w_strata, length_explode_df, bins_len, bins_age):
        """
        Obtains the constants associated with each stratum,
        which are used in the biomass density calculation
        """

        spec_strata_ind = spec_w_strata.index.unique()
        len_strata_ind = length_explode_df.index.unique()
        strata_ind = spec_strata_ind.intersection(len_strata_ind).values

        total_N = np.empty(strata_ind.shape[0], dtype=np.float64)

        spec_M_prop = np.empty(strata_ind.shape[0], dtype=np.float64)
        spec_F_prop = np.empty(strata_ind.shape[0], dtype=np.float64)

        len_M_prop = np.empty(strata_ind.shape[0], dtype=np.float64)
        len_F_prop = np.empty(strata_ind.shape[0], dtype=np.float64)

        fac2_ALL = np.empty(strata_ind.shape[0], dtype=np.float64)
        fac1_ALL = np.empty(strata_ind.shape[0], dtype=np.float64)

        fac1_M = np.empty(strata_ind.shape[0], dtype=np.float64)
        fac1_F = np.empty(strata_ind.shape[0], dtype=np.float64)

        fac2_M = np.empty(strata_ind.shape[0], dtype=np.float64)
        fac2_F = np.empty(strata_ind.shape[0], dtype=np.float64)

        len_wgt_prod = np.empty(strata_ind.shape[0], dtype=np.float64)
        len_wgt_M_prod = np.empty(strata_ind.shape[0], dtype=np.float64)
        len_wgt_F_prod= np.empty(strata_ind.shape[0], dtype=np.float64)

        len_weight_ALL, _, _ = self.generate_length_val_key(bins_len, reg_w0=None, reg_p=None,
                                                            len_name='Length',
                                                            val_name='Weight', df=spec_w_strata)

        ind = 0
        for stratum in strata_ind:

            spec_strata = spec_w_strata.loc[stratum]
            spec_strata_M = spec_strata[spec_strata['Sex'] == 1]
            spec_strata_F = spec_strata[spec_strata['Sex'] == 2]
            len_strata = length_explode_df.loc[stratum]
            len_strata_M = len_strata[len_strata['Sex'] == 1]
            len_strata_F = len_strata[len_strata['Sex'] == 2]

            len_age_key_sum, len_age_key_M_sum, len_age_key_F_sum \
                = self.get_age_related_sums(spec_strata_M, spec_strata_F, bins_len, bins_age)

            # total number of sexed fish at stations 1 and 2
            # total_N[ind] = spec_strata_M.shape[0] + spec_strata_F.shape[0] + len_strata_M.shape[0] + len_strata_F.shape[0]

            total_N[ind] = spec_strata_M.shape[0] + spec_strata_F.shape[0] + len_strata.shape[0]

            # total_N[ind] = spec_strata.shape[0] + len_strata.shape[0] - 1  # TODO: Correct this!

            # proportion of males/females in station 2
            spec_M_prop[ind] = spec_strata_M.shape[0] / total_N[ind]
            spec_F_prop[ind] = spec_strata_F.shape[0] / total_N[ind]

            # proportion of males/females in station 1
            len_M_prop[ind] = len_strata_M.shape[0] / total_N[ind]
            len_F_prop[ind] = len_strata_F.shape[0] / total_N[ind]

            # total proportion of sexed fish in station 2
            fac2_ALL[ind] = spec_M_prop[ind] + spec_F_prop[ind]

            # total proportion of sexed fish in station 1
            fac1_ALL[ind] = 1.0 - fac2_ALL[ind]

            fac1_M[ind] = fac1_ALL[ind] * 1.0
            fac1_F[ind] = fac1_ALL[ind] * 1.0  # TODO: something looks wrong here, ask Chu

            fac1_M[ind] = fac1_M[ind] / (fac1_M[ind] + spec_M_prop[ind])
            fac1_F[ind] = fac1_F[ind] / (fac1_F[ind] + spec_F_prop[ind])

            fac2_M[ind] = spec_M_prop[ind] / (fac1_M[ind] + spec_M_prop[ind])
            fac2_F[ind] = spec_F_prop[ind] / (fac1_F[ind] + spec_F_prop[ind])

            # average of total proportion of sexed fish pertaining to station 1
            # fac1_ALL = fac1_ALL / (fac1_ALL + fac2_ALL) # TODO: do we need to calculate this?

            # average of total proportion of sexed fish pertaining to station 2
            # fac2_ALL = fac2_ALL / (fac1_ALL + fac2_ALL) # TODO: do we need to calculate this?

            len_key_n, len_key_M, len_key_F = self.get_len_keys_strata(len_strata, bins_len)

            len_wgt_prod[ind] = np.dot(fac1_ALL[ind] * len_key_n + fac2_ALL[ind] * len_age_key_sum, len_weight_ALL)
            len_wgt_M_prod[ind] = np.dot(fac1_M[ind] * len_key_M + fac2_M[ind] * len_age_key_M_sum, len_weight_ALL)
            len_wgt_F_prod[ind] = np.dot(fac1_F[ind] * len_key_F + fac2_F[ind] * len_age_key_F_sum, len_weight_ALL)

            ind += 1

        total_N_da = xr.DataArray(data=total_N, coords={'strata': strata_ind})

        spec_M_prop_da = xr.DataArray(data=spec_M_prop, coords={'strata': strata_ind})
        spec_F_prop_da = xr.DataArray(data=spec_F_prop, coords={'strata': strata_ind})

        len_M_prop_da = xr.DataArray(data=len_M_prop, coords={'strata': strata_ind})
        len_F_prop_da = xr.DataArray(data=len_F_prop, coords={'strata': strata_ind})

        len_wgt_prod_da = xr.DataArray(data=len_wgt_prod, coords={'strata': strata_ind})
        len_wgt_M_prod_da = xr.DataArray(data=len_wgt_M_prod, coords={'strata': strata_ind})
        len_wgt_F_prod_da = xr.DataArray(data=len_wgt_F_prod, coords={'strata': strata_ind})


        # bio_calc_const = xr.Dataset({'total_N': total_N_da, 'spec_M_prop': spec_M_prop_da,
        #                              'spec_F_prop': spec_F_prop_da, 'len_M_prop': len_M_prop_da,
        #                              'len_F_prop_da': len_F_prop_da, 'fac1_ALL': fac1_ALL,
        #                              'fac2_ALL': fac2_ALL, 'fac1_M': fac1_M, 'fac1_F': fac1_F,
        #                              'fac2_M': fac2_M, 'fac2_F': fac2_F})

        bio_calc_const = xr.Dataset({'spec_M_prop': spec_M_prop_da, 'spec_F_prop': spec_F_prop_da,
                                     'len_M_prop': len_M_prop_da, 'len_F_prop': len_F_prop_da,
                                     'len_wgt_prod': len_wgt_prod_da, 'len_wgt_M_prod': len_wgt_M_prod_da,
                                     'len_wgt_F_prod': len_wgt_F_prod_da, 'total_N': total_N_da})

        return bio_calc_const

    # from EchoPro.load_stratification_data import LoadStrataData
    # strata_class = LoadStrataData(epro_2019)
    # # get the bins for the lengths
    # bins_len = epro_2019.params['bio_hake_len_bin']
    # # get the bins for the ages
    # bins_age = epro_2019.params['bio_hake_age_bin']
    #
    # bc = strata_class.get_biomass_constants(spec_w_strata, length_explode_df, bins_len, bins_age)
    # bc
    #
    # nntk_male = bio_dense_df.apply(lambda x: np.round(
    #     x.n_A * (bc.len_M_prop.sel(strata=x.Stratum).values + bc.spec_M_prop.sel(strata=x.Stratum).values)), axis=1)
    # nntk_female = bio_dense_df.apply(lambda x: np.round(
    #     x.n_A * (bc.len_F_prop.sel(strata=x.Stratum).values + bc.spec_F_prop.sel(strata=x.Stratum).values)), axis=1)
    #
    # bio_dense_df['nntk_male'] = nntk_male
    # bio_dense_df['nntk_female'] = nntk_female
    #
    # nWgt_male_int = bio_dense_df.apply(lambda x: x.nntk_male * bc.len_wgt_M_prod.sel(strata=x.Stratum).values, axis=1)
    # nWgt_female_int = bio_dense_df.apply(lambda x: x.nntk_female * bc.len_wgt_F_prod.sel(strata=x.Stratum).values,
    #                                      axis=1)
    #
    # bio_dense_df['nWgt_male'] = nWgt_male_int
    # bio_dense_df['nWgt_female'] = nWgt_female_int
    #
    # nWgt_unsexed_int = bio_dense_df.apply(
    #     lambda x: (x.n_A - x.nntk_male - x.nntk_female) * bc.len_wgt_prod.sel(strata=x.Stratum).values, axis=1)
    # bio_dense_df['nWgt_unsexed'] = nWgt_unsexed_int
    #
    # bio_dense_df['nWgt_total'] = bio_dense_df['nWgt_male'] + bio_dense_df['nWgt_female'] + bio_dense_df['nWgt_unsexed']
    #
    # age_len_key_da, age_len_key_wgt_da, age_len_key_norm_da = strata_class.get_age_key_das(spec_w_strata,
    #                                                                                        bins_len, bins_age)
    #
    # # TODO: it would probably be better to do an average of station 1 and 2 here... (Chu doesn't do this)
    # age_len_key_wgt_norm_da = age_len_key_wgt_da / age_len_key_wgt_da.sum(dim=['len_bins', 'age_bins'])
    #
    # # each stratum's multiplier once normalized weight has been calculated
    # age2_wgt_proportion_da = 1.0 - age_len_key_wgt_norm_da.isel(age_bins=0).sum(
    #     dim='len_bins') / age_len_key_wgt_norm_da.sum(dim=['len_bins', 'age_bins'])
    #
    # nWgt_total_2_prop = bio_dense_df.apply(lambda x: x.nWgt_total * age2_wgt_proportion_da.sel(strata=x.Stratum).values,
    #                                        axis=1)
    #
    # bio_dense_df['nWgt_total_2_prop'] = nWgt_total_2_prop

    def get_strata_data(self, stratification_index, KS_stratification,
                        transect_reduction_fraction: float = 0.0):
        """
        Get quantities and keys associated with strata for:
        1. trawl information
        2. Length - key
        3. Length-weight key
        4. Length-age - key for abundance & biomass

        Parameters
        ----------
        stratification_index : int  # TODO: this looks to mirror KS_stratification, it seems to be unnecessary to use this
                                    # TODO: Ask Chu if it is ok to remove this.
            Index for the chosen stratification
            0 = INPFC strata
            1 = KS (trawl)-based
            2-6 = geographically based but close to trawl-based stratification
            7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
            10 = one stratum for the whole survey
        KS_stratification : int
            Specifies the type of stratification to be used.
            0 = Pre-Stratification or customized stratification (geographically defined)
            1 = Post-Stratification (KS-based or trawl-based)
        transect_reduction_fraction : float
            Reduction fraction for transect TODO: should this be 5 or 0.05?
        """

        # TODO: is it necessary to have this? It looks like it is not necessary.
        # if para.proc.stratification_index == 0 & para.proc.KS_stratification == 0 # INPFC with
        # if dat0(1, 1) == 2015
        #     # remove age1 hauls but not those used in transect_region_hual files, to be consistent
        #     with K - S stratification
        #     para.proc.age1_haul = [4 9 10 13 15 18 23 30 35 39 41 48 60 69 227];

        if self.EPro.params['bio_data_type'] != 3:

            print("Do we need to set stratum_id or just use strata_df? Look into this!")
            # stratum_id = self.strata_df['Cluster name'] or self.strata_df['INPFC']

        else:
            # TODO: If this statement is triggered, then this is probably not the right place to
            #  put this code because is could potential be ran in bootstrapping.
            raise NotImplementedError("Loading the stratification file has not been "
                                      + f"implemented for bio_data_type = {self.EPro.params['bio_data_type']}")

        # TODO: these lines look unnecessary for the Python version
        # data.bio_acoust.haul_wgt_tbl = dat(:, 2: 3);
        # data.bio.haul_strata = unique(dat(:, 1));

        # TODO: maybe an option like this would be useful, currently all age 0 are removed
        #  we would need to modify self.__get_bio_strata() too
        # rm_strata_ind_0: True  # If true, removes strata index 0 from strata file e.g. Cluster name=0 or
        #                        # INPFC=0 in 2019/US&CAN strata 2019_final.xlsx
        # removes strata index 0 -- corresponds to age 0
        # if self.EPro.params['rm_strata_ind_0']:
        #     if stratification_index == 1:
        #         self.EPro.strata_df = self.EPro.strata_df[self.EPro.strata_df['Cluster name'] != 0]
        #     else:
        #         self.EPro.strata_df = self.EPro.strata_df[self.EPro.strata_df['INPFC'] != 0]

        # self.__get_bio_strata(stratification_index, KS_stratification, transect_reduction_fraction)
        #
        # self.__get_bio_len_age_haul()

        # self.__remove_empty_bio_len_age_haul()

        return