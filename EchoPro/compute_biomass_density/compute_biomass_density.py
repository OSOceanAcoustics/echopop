import numpy as np
import pandas as pd
import xarray as xr


class ComputeBiomassDensity:
    """
    A class that computes the biomass density
    of an animal population based off of NASC
    values and other associated data.

    Parameters
    ----------
    epro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.epro will also change this object.
    """

    def __init__(self, epro = None):

        self.epro = epro

        self.bio_hake_len_bin = epro.params['bio_hake_len_bin']
        self.bio_hake_age_bin = epro.params['bio_hake_age_bin']

    @staticmethod
    def _get_bin_ind(input_data: np.array, centered_bins: np.array):
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

    def _generate_length_val_key(self, len_name: str, val_name: str, df: pd.DataFrame = None):
        """
        Generates a length-value key by
        1. Binning the lengths
        2. Fitting a linear regression model to the length and values points
        3. Set bins with greater than 5 samples equal to the mean of
        the values corresponding to that bin
        4. Set bins with less than 5 samples equal to regression curve
        calculated in step 2

        Parameters
        ----------
        len_name : str 
            Name of length column in df 
        val_name : str 
            Name of value column in df
        df: pd.DataFrame
            Pandas Dataframe containing the length and value data

        Returns
        -------
        len_val_key np.ndarray 
            1D array representing the length-value key 
        """

        # select the indices that do not have nan in either Length or Weight
        len_wgt_nonull = np.logical_and(df[len_name].notnull(), df[val_name].notnull())
        df_no_null = df.loc[len_wgt_nonull]

        # get numpy version of those entries that are not NaN
        L = df_no_null[len_name].values
        V = df_no_null[val_name].values

        # binned length indices
        len_bin_ind = self._get_bin_ind(L, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([i.shape[0] for i in len_bin_ind])

        # length-value regression
        p = np.polyfit(np.log10(L), np.log10(V), 1)

        # obtain linear regression parameters
        reg_w0 = 10.0 ** p[1]
        reg_p = p[0]

        # value at length or length-value-key.
        # length-value-key per fish over entire survey region (an array)
        len_val_reg = reg_w0 * self.bio_hake_len_bin ** reg_p

        # set length-value key to the mean of the values in the bin
        len_val_key = np.array([np.mean(V[ind]) if ind.size > 0 else 0.0 for ind in len_bin_ind])

        # replace those bins with less than 5 samples with the regression value
        less_five_ind = np.argwhere(len_bin_cnt < 5).flatten()
        len_val_key[less_five_ind] = len_val_reg[less_five_ind]

        return len_val_key

    def get_norm_len_key_station_1(self, df: pd.DataFrame):
        """
        Computes the normalized length key for
        data obtained from station 1 i.e. data
        that tells you how many fish are of a
        particular length.

        Parameters
        ----------
        df : pd.DataFrame
            Data from station 1 that has columns Length and
            Frequency.

        Returns
        -------
        A numpy array of the normalized length key i.e. the
        count of each bin divided by the total of all bin counts
        """

        # get numpy arrays of the Length and Frequency columns
        length_arr = df.Length.values
        freq_arr = df.Frequency.values

        # binned length indices
        len_ind = self._get_bin_ind(length_arr, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([np.sum(freq_arr[i]) for i in len_ind])

        return len_bin_cnt / np.sum(len_bin_cnt)

    def get_norm_len_key_station_2(self, df: pd.DataFrame):
        """
        Computes the normalized length key for
        data obtained from station 2 i.e. data
        that does not have a frequency associated
        with it.

        Parameters
        ----------
        df : pd.DataFrame
            Data from station 2 that has a Length column

        Returns
        -------
        A numpy array of the normalized length key i.e. the
        count of each bin divided by the total of all bin counts
        """

        # numpy array of Length column
        length_arr = df.Length.values

        # binned length indices
        len_bin_ind = self._get_bin_ind(length_arr, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([i.shape[0] for i in len_bin_ind])

        return len_bin_cnt / np.sum(len_bin_cnt)

    def _get_age_weight_keys(self, df, age_name, len_name):

        # TODO: document/comment this!

        input_arr_len = df[len_name].values
        input_arr_age = df[age_name].values
        input_arr_wgt = df.Weight.values

        age_bins_ind = self._get_bin_ind(input_arr_age, self.bio_hake_age_bin)
        len_bin_ind = self._get_bin_ind(input_arr_len[age_bins_ind[0]], self.bio_hake_len_bin)

        return np.array([np.sum(input_arr_wgt[age_bins_ind[0][i]]) for i in len_bin_ind]).sum() / input_arr_wgt.sum()

    def _compute_proportions(self, spec_strata_M, spec_strata_F,
                             len_strata, len_strata_M, len_strata_F):

        # TODO: comment/document this!

        # total number of sexed fish at stations 1 and 2
        total_N = spec_strata_M.shape[0] + spec_strata_F.shape[0] + len_strata.Frequency.sum()
        # total_N = spec_strata.shape[0] + (len_strata.Frequency.sum())  # TODO: This is what it should be

        # proportion of males/females in station 2
        spec_M_prop = spec_strata_M.shape[0] / total_N
        spec_F_prop = spec_strata_F.shape[0] / total_N

        # proportion of males/females in station 1
        len_M_prop = len_strata_M.Frequency.sum() / total_N
        len_F_prop = len_strata_F.Frequency.sum() / total_N

        # total proportion of sexed fish in station 2
        fac2_ALL = spec_M_prop + spec_F_prop

        # total proportion of sexed fish in station 1
        fac1_ALL = 1.0 - fac2_ALL

        fac1_M = fac1_ALL * 1.0
        fac1_F = fac1_ALL * 1.0  # TODO: something looks wrong here, ask Chu

        fac1_M = fac1_M / (fac1_M + spec_M_prop)
        fac1_F = fac1_F / (fac1_F + spec_F_prop)

        fac2_M = spec_M_prop / (fac1_M + spec_M_prop)
        fac2_F = spec_F_prop / (fac1_F + spec_F_prop)

        # average of total proportion of sexed fish pertaining to station 1
        fac1_ALL = fac1_ALL / (fac1_ALL + fac2_ALL)  # TODO: do we need to calculate this?

        # average of total proportion of sexed fish pertaining to station 2
        fac2_ALL = fac2_ALL / (fac1_ALL + fac2_ALL)  # TODO: do we need to calculate this?

        return spec_M_prop + len_M_prop, spec_F_prop + len_F_prop, fac1_M, fac1_F, fac2_M, fac2_F, fac1_ALL, fac2_ALL

    def _fill_len_wgt_prod(self, bio_const_df: pd.DataFrame,
                           stratum: int, spec_drop_df: pd.DataFrame,
                           length_drop_df: pd.DataFrame, len_weight_key: np.array):
        """


        Parameters
        ----------
        bio_const_df : pd.DataFrame
        stratum : int
        spec_drop_df : pd.DataFrame
        length_drop_df : pd.DataFrame
        len_weight_key : np.array

        Returns
        -------

        """

        # TODO: document and comment!

        # get specimen in the stratum and split sort into males and females
        spec_stratum = spec_drop_df.loc[stratum]
        spec_strata_m = spec_stratum[spec_stratum['Sex'] == 1]
        spec_strata_f = spec_stratum[spec_stratum['Sex'] == 2]

        # get lengths in the stratum and split sort into males and females
        len_strata = length_drop_df.loc[stratum]
        len_strata_m = len_strata[len_strata['Sex'] == 1]
        len_strata_f = len_strata[len_strata['Sex'] == 2]

        # get the normalized length keys for station 1
        len_key_all_s1 = self.get_norm_len_key_station_1(len_strata)
        len_key_m_s1 = self.get_norm_len_key_station_1(len_strata_m)
        len_key_f_s1 = self.get_norm_len_key_station_1(len_strata_f)

        # get the normalized length keys for station 2
        # len_key_all_s2 = self.get_norm_len_key_station_2(spec_stratum)  # TODO: this is what should be done!
        spec_stratum_mf = pd.concat([spec_strata_m, spec_strata_f], axis=0)
        len_key_all_s2 = self.get_norm_len_key_station_2(spec_stratum_mf)
        len_key_m_s2 = self.get_norm_len_key_station_2(spec_strata_m)
        len_key_f_s2 = self.get_norm_len_key_station_2(spec_strata_f)

        # TODO: stopped here!

        m_prop, f_prop, fac1_m, fac1_f, fac2_m, \
        fac2_f, fac1_all, fac2_all = self._compute_proportions(spec_strata_m, spec_strata_f,
                                                               len_strata, len_strata_m, len_strata_f)

        bio_const_df.M_prop.loc[stratum] = m_prop
        bio_const_df.F_prop.loc[stratum] = f_prop

        bio_const_df.len_wgt_prod.loc[stratum] = np.dot(fac1_all * len_key_all_s1 + fac2_all * len_key_all_s2,
                                                        len_weight_key)
        bio_const_df.len_wgt_M_prod.loc[stratum] = np.dot(fac1_m * len_key_m_s1 + fac2_m * len_key_m_s2,
                                                          len_weight_key)
        bio_const_df.len_wgt_F_prod.loc[stratum] = np.dot(fac1_f * len_key_f_s1 + fac2_f * len_key_f_s2,
                                                          len_weight_key)

        return bio_const_df

    def _get_biomass_constants(self):
        """
        Obtains the constants associated with each stratum,
        which are used in the biomass density calculation.
        Specifically, we obtain the following constants for
        each stratum:
        * M_prop -- proportion of males
        * F_prop -- proportion of females
        * len_wgt_prod -- product of the length-key and
        the weight-key for the whole population
        * len_wgt_M_prod -- product of the length-key and
        the weight-key for the male population
        * len_wgt_F_prod -- product of the length-key and
        the weight-key for the female population

        Returns
        -------
        bio_const_df : pd.Dataframe
            Dataframe with index of stratum and columns
            corresponding to the constants specified
            above.
        """

        # determine the strata that are in both specimen_df and length_df
        spec_strata_ind = self.epro.specimen_df.index.unique()
        len_strata_ind = self.epro.length_df.index.unique()
        strata_ind = spec_strata_ind.intersection(len_strata_ind).values

        # obtain the length-weight key for all specimen data
        len_weight_ALL = self._generate_length_val_key(len_name='Length',
                                                       val_name='Weight',
                                                       df=self.epro.specimen_df)

        # select the indices that do not have nan in either Length or Weight
        spec_drop = self.epro.specimen_df.dropna(how='any')

        # select the indices that do not have nan in either Length or Weight
        length_drop_df = self.epro.length_df.dropna(how='any')

        # initialize dataframe that will hold all important calculated constants
        bio_const_df = pd.DataFrame(columns=['M_prop', 'F_prop', 'len_wgt_prod',
                                             'len_wgt_M_prod', 'len_wgt_F_prod'],
                                    index=strata_ind, dtype=np.float64)

        # for each stratum compute the necessary constants
        for stratum in strata_ind:
            bio_const_df = self._fill_len_wgt_prod(bio_const_df, stratum, spec_drop,
                                                        length_drop_df, len_weight_ALL)

        return bio_const_df

    def _add_stratum_column(self):
        """
        Adds the "stratum" column to self.epro.strata_df
        and self.epro.length_df. Additionally, this
        function will set the index to "stratum".
        """

        # get df relating the haul to the stratum
        strata_haul_df = self.epro.strata_df.reset_index()[['Haul', 'stratum']].set_index('Haul')

        # add stratum column to strata_df and set it as the index
        self.epro.specimen_df['stratum'] = strata_haul_df.loc[self.epro.specimen_df.index]
        self.epro.specimen_df.set_index('stratum', inplace=True)

        # add stratum column to length_df and set it as the index
        self.epro.length_df['stratum'] = strata_haul_df.loc[self.epro.length_df.index]
        self.epro.length_df.set_index('stratum', inplace=True)

    def get_final_biomass_table(self):

        # TODO: clean up this code!! Is it necessary to create bio_dense_df?

        # add stratum column to length and specimen df and set it as the index
        self._add_stratum_column()

        bc_df = self._get_biomass_constants()

        # TODO: should we include the below code that is commented out?
        # # calculates the interval for the area calculation
        # interval = (nasc_df['VL start'].iloc[1:].values - nasc_df['VL start'].iloc[:-1].values)
        # last_interval = nasc_df['VL end'].iloc[-1] - nasc_df['VL start'].iloc[-1]
        #
        # interval = np.concatenate([interval, np.array([last_interval])])
        #
        # median_interval = np.median(interval)
        #
        # # remove outliers at the end of the transect
        # ind_outliers = np.argwhere(np.abs(interval - median_interval) > 0.05).flatten()
        # interval[ind_outliers] = nasc_df['VL end'].values[ind_outliers] - nasc_df['VL start'].values[ind_outliers]

        bio_dense_df = self.epro.nasc_df[['Stratum', 'NASC', 'Haul']].copy()
        # bio_dense_df['interval'] = interval

        wgt_vals = self.epro.strata_df.reset_index().set_index('Haul')['wt']
        wgt_vals_ind = wgt_vals.index

        mix_sa_ratio = self.epro.nasc_df.apply(lambda x: wgt_vals[x.Haul] if x.Haul in wgt_vals_ind else 0.0, axis=1)
        self.epro.nasc_df['mix_sa_ratio'] = mix_sa_ratio

        bio_dense_df['n_A'] = self.epro.nasc_df.apply(
            lambda x: np.round((x.mix_sa_ratio * x.NASC) / float(self.epro.strata_sig_b.loc[x.Stratum])),
            axis=1)

        # bio_dense_df['A'] = bio_dense_df['interval'] * nasc_df['Spacing']
        # bio_dense_df['N_A'] = bio_dense_df['n_A'] * bio_dense_df['A']

        # expand the bio constants dataframe so that it corresponds to bio_dense_df
        bc_expanded_df = bc_df.loc[bio_dense_df.Stratum.values]

        bio_dense_df['nntk_male'] = np.round(bio_dense_df.n_A.values * bc_expanded_df.M_prop.values)
        bio_dense_df['nntk_female'] = np.round(bio_dense_df.n_A.values * bc_expanded_df.F_prop.values)

        bio_dense_df['nWgt_male'] = bio_dense_df.nntk_male.values * bc_expanded_df.len_wgt_M_prod.values
        bio_dense_df['nWgt_female'] = bio_dense_df.nntk_female.values * bc_expanded_df.len_wgt_F_prod.values

        bio_dense_df['nWgt_unsexed'] = (bio_dense_df.n_A.values - bio_dense_df.nntk_male.values -
                                        bio_dense_df.nntk_female.values) * bc_expanded_df.len_wgt_prod.values

        bio_dense_df['nWgt_total'] = bio_dense_df['nWgt_male'] + bio_dense_df['nWgt_female'] + bio_dense_df[
            'nWgt_unsexed']

        # TODO: This is necessary to match the Matlab output
        #  in theory this should be done when we load the df,
        #  however, this changes the results slightly.
        spec_drop = self.epro.specimen_df.dropna(how='any')

        # age_len_key_da, age_len_key_wgt_da, age_len_key_norm_da = self.get_age_key_das(spec_drop)

        # TODO: it would probably be better to do an average of station 1 and 2 here... (Chu doesn't do this)
        # age_len_key_wgt_norm_da = age_len_key_wgt_da / age_len_key_wgt_da.sum(dim=['len_bins', 'age_bins'])

        # bran_cool = age_len_key_wgt_norm_da.isel(age_bins=0).sum(dim='len_bins') / age_len_key_wgt_norm_da.sum(dim=['len_bins', 'age_bins'])
        # print(f"bran_cool = {bran_cool}")
        # print(age_len_key_wgt_da.isel(age_bins=0).sum(dim='len_bins'))
        # print(" ")

        # each stratum's multiplier once normalized weight has been calculated
        stratum_ind = spec_drop.index.unique()
        age2_wgt_proportion_df = pd.DataFrame(columns=['val'], index=stratum_ind)
        for i in stratum_ind:
            # print(spec_drop.loc[i].Weight.sum())
            out = self._get_age_weight_keys(spec_drop.loc[i], age_name="Age", len_name="Length")
            age2_wgt_proportion_df.loc[i].val = 1.0 - out

        # each stratum's multiplier once normalized weight has been calculated
        # age2_wgt_proportion_da = 1.0 - age_len_key_wgt_norm_da.isel(age_bins=0).sum(
        #     dim='len_bins') / age_len_key_wgt_norm_da.sum(dim=['len_bins', 'age_bins'])

        # nWgt_total_2_prop = bio_dense_df.apply(
        #     lambda x: x.nWgt_total * float(age2_wgt_proportion_da.sel(strata=x.Stratum)),
        #     axis=1)

        nWgt_total_2_prop = bio_dense_df.apply(
            lambda x: x.nWgt_total * float(age2_wgt_proportion_df.loc[x.Stratum]), axis=1)

        # minimal columns to do Jolly Hampton CV on data that has not been kriged
        self.epro.final_biomass_table = self.epro.nasc_df[['Latitude', 'Longitude', 'Stratum', 'Spacing']].copy()
        self.epro.final_biomass_table["normalized_biomass_density"] = nWgt_total_2_prop

        print("We are using our own biomass density calculation!")
