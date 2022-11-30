import numpy as np
import pandas as pd
import geopandas as gpd
from typing import Tuple, List, Optional, Union


class ComputeBiomassDensity:
    """
    A class that computes the biomass density
    of an animal population based off of NASC
    values and other associated data.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        self.survey will also change this object.
    """

    def __init__(self, survey=None):

        self.survey = survey

        self.bio_hake_len_bin = survey.params['bio_hake_len_bin']
        self.bio_hake_age_bin = survey.params['bio_hake_age_bin']

        # initialize class variables that will be set downstream
        self.length_df = None
        self.strata_df = None
        self.specimen_df = None
        self.nasc_df = None
        self.final_biomass_table = None
        self.krig_results_gdf = None
        self.bio_const_df = None   # biomass constants for each stratum
        self.age2_wgt_prop_df = None
        self.strata_sig_b = None

    def _get_strata_sig_b(self) -> None:
        """
        Computes the backscattering cross-section (sigma_b),
        using the strata, specimen, and length dataframes.
        These values are then stored in self.strata_sig_b
        as a Pandas series with index "stratum_num".
        """

        # TODO: the target strength functions are specific to Hake, replace with input in the future

        # initialize sig_bs_haul column in strata_df
        self.strata_df["sig_bs_haul"] = np.nan

        # select the indices that do not have nan in either length or weight
        spec_df = self.specimen_df[['length', 'weight']].copy()
        spec_df = spec_df.dropna(how='any')

        for haul_num in spec_df.index.unique():

            # lengths from specimen file associated with index haul_num
            spec_len = spec_df.loc[haul_num]['length']

            if haul_num in self.length_df.index:

                # add lengths from length file associated with index haul_num
                length_len = self.length_df.loc[haul_num]['length'].values
                length_count = self.length_df.loc[haul_num]['length_count'].values

                # empirical relation for target strength
                TS0j_length = 20.0 * np.log10(length_len) - 68.0

                # sum of target strengths
                sum_TS0j_length = np.nansum((10.0 ** (TS0j_length / 10.0)) * length_count)

                # total number of values used to calculate sum_TS0j_length
                num_length = np.nansum(length_count)

            else:

                # sum of target strengths
                sum_TS0j_length = 0.0

                # total number of values used to calculate sum_TS0j_length
                num_length = 0.0

            # empirical relation for target strength
            TS0j_spec = 20.0 * np.log10(spec_len) - 68.0

            # sum of target strengths
            sum_TS0j_spec = np.nansum(10.0 ** (TS0j_spec / 10.0))

            # mean differential backscattering cross-section for each haul
            self.strata_df.loc[haul_num,
                                      "sig_bs_haul"] = (sum_TS0j_spec + sum_TS0j_length)/(num_length + TS0j_spec.size)

        # mean backscattering cross-section for each stratum
        self.strata_sig_b = 4.0 * np.pi * self.strata_df['sig_bs_haul'].groupby('stratum_num').mean()

    def _fill_missing_strata_indices(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        When selecting a subset of the transects, it is possible that some
        strata do not have important constants defined. This function fills in these
        missing values with artificial data. This is done as follows for
        each missing stratum:
        - If the value is only known for 1 stratum, then all missing stratum will be
        filled with the value of the known stratum, otherwise the below items
        will determine how the value is filled:
            - If the missing stratum index is less than the minimum known stratum index,
            then the missing stratum value will be set to the value of the minimum known
            stratum index.
            - If the missing stratum index is greater than the maximum known stratum index,
            then the missing stratum value will be set to the value of the maximum known
            stratum index.
            - If the missing stratum index is between the minimum and maximum known stratum
            indices, then the missing value will be set to the average of the values provided
            by the two closest known stratum indices.

        Parameters
        ----------
        df: pd.DataFrame
            A Dataframe with stratum as its index

        Returns
        -------
        df: pd.DataFrame
            ``df`` with all missing strata data filled in
        """

        # all strata indices that are missing
        missing_strata = set(self.nasc_df["stratum_num"].unique()) - set(df.index.values)

        # if there are no missing strata then do nothing
        if len(missing_strata) == 0:
            return df

        max_known_strata = df.index.max()
        min_known_strata = df.index.min()

        # if there is only 1 stratum with a value, fill all missing
        # strata with the only known value
        fill_w_one_val = False
        if len(df.index) == 1:
            fill_w_one_val = True

        for strata in missing_strata:

            if fill_w_one_val:

                # fill all missing strata with the only known value
                df.loc[strata] = df.loc[max_known_strata]

            else:

                if strata > max_known_strata:

                    # fill value with the value at the largest known stratum
                    df.loc[strata] = df.loc[max_known_strata]

                elif strata < min_known_strata:

                    # fill value with the value at the smallest known stratum
                    df.loc[strata] = df.loc[min_known_strata]

                else:

                    # get the two indices of strata_sig_b.index that are closest to strata
                    strata_ind_diff = np.abs(strata - df.index)
                    idx = np.argpartition(strata_ind_diff, 2)

                    # get values at two smallest indices
                    val_1 = df.iloc[idx[0]]
                    val_2 = df.iloc[idx[1]]

                    # average two closest values
                    average_vals = (val_1 + val_2) / 2.0

                    # fill in value with the average value of the two closest known strata
                    df.loc[strata] = average_vals

        return df

    @staticmethod
    def _get_bin_ind(input_data: np.ndarray,
                     centered_bins: np.ndarray) -> List[np.ndarray]:
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
        input_data: np.ndarray
            The data to create a histogram of.
        centered_bins: np.ndarray
            An array that specifies the bin centers.

        Returns
        -------
        hist_ind: list
            The index values of input_data corresponding to the histogram

        """

        # get the distance between bin centers
        bin_diff = np.diff(centered_bins) / 2.0

        # fill the first bin
        hist_ind = [np.argwhere(input_data <= centered_bins[0] + bin_diff[0]).flatten()]

        for i in range(len(centered_bins) - 2):
            # get values greater than lower bound
            g_lb = centered_bins[i] + bin_diff[i] < input_data

            # get values less than or equal to the upper bound
            le_ub = input_data <= centered_bins[i + 1] + bin_diff[i + 1]

            # fill bin
            hist_ind.append(np.argwhere(g_lb & le_ub).flatten())

        # fill in the last bin
        hist_ind.append(np.argwhere(input_data > centered_bins[-2] + bin_diff[-1]).flatten())

        return hist_ind

    def _add_stratum_column(self) -> None:
        """
        Adds the ``stratum_num`` column to self.strata_df
        and self.length_df. Additionally, this
        function will set the index to ``stratum_num``.
        """

        # get df relating the haul to the stratum
        strata_haul_df = self.strata_df.reset_index()[['haul_num', 'stratum_num']].set_index('haul_num')

        # add stratum_num column to specimen_df and set it as the index
        self.specimen_df['stratum_num'] = strata_haul_df.loc[self.specimen_df.index]
        self.specimen_df.set_index('stratum_num', inplace=True)

        # add stratum_num column to length_df and set it as the index
        self.length_df['stratum_num'] = strata_haul_df.loc[self.length_df.index]
        self.length_df.set_index('stratum_num', inplace=True)

    def _generate_length_val_key(self, len_name: str, val_name: str,
                                 df: pd.DataFrame = None) -> np.ndarray:
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
        len_val_key : np.ndarray
            1D array representing the length-value key 
        """

        # select the indices that do not have nan in either length or value
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

    def _get_norm_len_key_station_1(self, df: pd.DataFrame) -> np.ndarray:
        """
        Computes the normalized length key for
        data obtained from station 1 i.e. data
        that tells you how many fish are of a
        particular length.

        Parameters
        ----------
        df : pd.DataFrame
            Data from station 1 that has columns ``length`` and ``length_count``.

        Returns
        -------
        A numpy array of the normalized length key i.e. the
        count of each bin divided by the total of all bin counts
        """

        # get numpy arrays of the length and length_count columns
        length_arr = df.length.values
        length_count_arr = df.length_count.values

        # binned length indices
        len_ind = self._get_bin_ind(length_arr, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([np.sum(length_count_arr[i]) for i in len_ind])

        return len_bin_cnt / np.sum(len_bin_cnt)

    def _get_norm_len_key_station_2(self, df: pd.DataFrame) -> np.ndarray:
        """
        Computes the normalized length key for
        data obtained from station 2 i.e. data
        that does not have a frequency associated
        with it.

        Parameters
        ----------
        df : pd.DataFrame
            Data from fish measurement station 2 that has a ``length`` column

        Returns
        -------
        A numpy array of the normalized length key i.e. the
        count of each bin divided by the total of all bin counts
        """

        # numpy array of length column
        length_arr = df.length.values

        # binned length indices
        len_bin_ind = self._get_bin_ind(length_arr, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([i.shape[0] for i in len_bin_ind])

        return len_bin_cnt / np.sum(len_bin_cnt)

    @staticmethod
    def _compute_proportions(spec_strata_m: pd.DataFrame, spec_strata_f: pd.DataFrame,
                             len_strata: pd.DataFrame, len_strata_m: pd.DataFrame,
                             len_strata_f: pd.DataFrame) -> Tuple[list, list, list, list]:
        """
        Computes proportions needed for biomass density
        calculation.

        Parameters
        ----------
        spec_strata_m : pd.Dataframe
            Subset of specimen_df corresponding to the stratum and males
        spec_strata_f : pd.Dataframe
            Subset of specimen_df corresponding to the stratum and females
        len_strata : pd.Dataframe
            Subset of length_df corresponding to the stratum
        len_strata_m : pd.Dataframe
            Subset of length_df corresponding to the stratum and males
        len_strata_f : pd.Dataframe
            Subset of length_df corresponding to the stratum and females

        Returns
        -------
        gender_prop : list
            List of male and female proportions for both stations, respectively
        fac1 : list
            List of average fraction of sexed fish from station 1
        fac2 : list
            List of average fraction of sexed fish from station 2
        tot_prop : list
            List of average of total proportion of sexed fish
        """

        # total number of sexed fish at stations 1 and 2
        total_n = spec_strata_m.shape[0] + spec_strata_f.shape[0] + len_strata.length_count.sum()
        # total_n = spec_strata.shape[0] + (len_strata.length_count.sum())  # TODO: This is what it should be

        # proportion of males/females in station 2
        spec_m_prop = spec_strata_m.shape[0] / total_n
        spec_f_prop = spec_strata_f.shape[0] / total_n

        # proportion of males/females in station 1
        len_m_prop = len_strata_m.length_count.sum() / total_n
        len_f_prop = len_strata_f.length_count.sum() / total_n

        # total proportion of sexed fish in station 2
        tot_prop2 = spec_m_prop + spec_f_prop

        # total proportion of sexed fish in station 1
        tot_prop1 = 1.0 - tot_prop2

        # fraction of males and females in station 1
        fac1_m = tot_prop1 * 1.0  # TODO: something looks wrong here
        fac1_f = tot_prop1 * 1.0  # TODO: something looks wrong here

        # average fraction of males and females from station 1?
        fac1_m = fac1_m / (fac1_m + spec_m_prop)
        fac1_f = fac1_f / (fac1_f + spec_f_prop)

        # average fraction of males and females from station 2?
        fac2_m = spec_m_prop / (fac1_m + spec_m_prop)
        fac2_f = spec_f_prop / (fac1_f + spec_f_prop)

        # average of total proportion of sexed fish pertaining to station 1
        tot_prop1 = tot_prop1 / (tot_prop1 + tot_prop2)  # TODO: do we need to calculate this? Denominator will be 1

        # average of total proportion of sexed fish pertaining to station 2
        tot_prop2 = tot_prop2 / (tot_prop1 + tot_prop2)  # TODO: do we need to calculate this? Denominator will be 1

        # list of male and female proportions for both stations
        gender_prop = [spec_m_prop + len_m_prop, spec_f_prop + len_f_prop]

        # list of average fraction of sexed fish from station 1
        fac1 = [fac1_m, fac1_f]

        # list of average fraction of sexed fish from station 2
        fac2 = [fac2_m, fac2_f]

        # list of average of total proportion of sexed fish
        tot_prop = [tot_prop1, tot_prop2]

        return gender_prop, fac1, fac2, tot_prop

    def _fill_len_wgt_prod(self, bio_const_df: pd.DataFrame,
                           stratum: int, spec_drop_df: pd.DataFrame,
                           length_drop_df: pd.DataFrame,
                           len_weight_key: np.array) -> pd.DataFrame:
        """
        Fills in the biomass constant dataframe for a
        particular stratum.

        Parameters
        ----------
        bio_const_df : pd.DataFrame
            Biomass constant dataframe to fill
        stratum : int
            Stratum to fill
        spec_drop_df : pd.DataFrame
            specimen_df with NaN values dropped
        length_drop_df : pd.DataFrame
            length_df with NaN values dropped
        len_weight_key : np.array
            length-weight key for all specimen data

        Returns
        -------
        bio_const_df : pd.DataFrame
            Biomass constant dataframe with stratum filled in
        """

        # get specimen in the stratum and split into males and females
        spec_stratum = spec_drop_df.loc[stratum]
        spec_strata_m = spec_stratum[spec_stratum['sex'] == 1]
        spec_strata_f = spec_stratum[spec_stratum['sex'] == 2]

        # get lengths in the stratum and split into males and females
        len_strata = length_drop_df.loc[stratum]
        len_strata_m = len_strata[len_strata['sex'] == 1]
        len_strata_f = len_strata[len_strata['sex'] == 2]

        # get the normalized length keys for station 1
        len_key_all_s1 = self._get_norm_len_key_station_1(len_strata)
        len_key_m_s1 = self._get_norm_len_key_station_1(len_strata_m)
        len_key_f_s1 = self._get_norm_len_key_station_1(len_strata_f)

        # get the normalized length keys for station 2
        # len_key_all_s2 = self._get_norm_len_key_station_2(spec_stratum)  # TODO: this is what should be done!
        spec_stratum_mf = pd.concat([spec_strata_m, spec_strata_f], axis=0)
        len_key_all_s2 = self._get_norm_len_key_station_2(spec_stratum_mf)
        len_key_m_s2 = self._get_norm_len_key_station_2(spec_strata_m)
        len_key_f_s2 = self._get_norm_len_key_station_2(spec_strata_f)

        gender_prop, fac1, fac2, tot_prop = self._compute_proportions(spec_strata_m, spec_strata_f,
                                                                      len_strata, len_strata_m,
                                                                      len_strata_f)

        # fill df with bio constants needed for biomass density calc
        bio_const_df.M_prop.loc[stratum] = gender_prop[0]
        bio_const_df.F_prop.loc[stratum] = gender_prop[1]
        bio_const_df.len_wgt_prod.loc[stratum] = np.dot(tot_prop[0] * len_key_all_s1 + tot_prop[1] * len_key_all_s2,
                                                        len_weight_key)
        bio_const_df.len_wgt_M_prod.loc[stratum] = np.dot(fac1[0] * len_key_m_s1 + fac2[0] * len_key_m_s2,
                                                          len_weight_key)
        bio_const_df.len_wgt_F_prod.loc[stratum] = np.dot(fac1[1] * len_key_f_s1 + fac2[1] * len_key_f_s2,
                                                          len_weight_key)

        return bio_const_df

    def _get_biomass_constants(self) -> None:
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

        Notes
        -----
        The following class variable is created in this function:
        bio_const_df : pd.Dataframe
            Dataframe with index of stratum and columns
            corresponding to the constants specified
            above.
        """

        # determine the strata that are in both specimen_df and length_df
        spec_strata_ind = self.specimen_df.index.unique()
        len_strata_ind = self.length_df.index.unique()
        strata_ind = spec_strata_ind.intersection(len_strata_ind).values

        # obtain the length-weight key for all specimen data
        len_weight_spec = self._generate_length_val_key(len_name='length',
                                                        val_name='weight',
                                                        df=self.specimen_df)

        # select the indices that do not have nan in either Length or Weight
        spec_drop = self.specimen_df.dropna(how='any')

        # select the indices that do not have nan in either Length or Weight
        length_drop_df = self.length_df.dropna(how='any')

        # initialize dataframe that will hold all important calculated constants
        bio_const_df = pd.DataFrame(columns=['M_prop', 'F_prop', 'len_wgt_prod',
                                             'len_wgt_M_prod', 'len_wgt_F_prod'],
                                    index=strata_ind, dtype=np.float64)

        # for each stratum compute the necessary constants
        for stratum in strata_ind:
            bio_const_df = self._fill_len_wgt_prod(bio_const_df, stratum, spec_drop,
                                                   length_drop_df, len_weight_spec)

        self.bio_const_df = bio_const_df

    @staticmethod
    def _get_interval(nasc_df: pd.DataFrame) -> np.ndarray:
        """
        Calculates the interval needed for the area calculation.

        Parameters
        ----------
        nasc_df : pd.DataFrame
            NASC df with ``vessel_log_start`` and ``vessel_log_end`` columns

        Returns
        -------
        interval : np.ndarray
            Intervals needed for area calculation
        """

        # calculate interval for all values except for the last interval
        interval = (nasc_df['vessel_log_start'].iloc[1:].values - nasc_df['vessel_log_start'].iloc[:-1].values)

        # calculate last interval
        last_interval = nasc_df['vessel_log_end'].iloc[-1] - nasc_df['vessel_log_start'].iloc[-1]

        # combines all intervals
        interval = np.concatenate([interval, np.array([last_interval])])

        # get median interval, so we can use it to remove outliers
        median_interval = np.median(interval)

        # remove outliers at the end of the transect
        ind_outliers = np.argwhere(np.abs(interval - median_interval) > 0.05).flatten()
        interval[ind_outliers] = nasc_df['vessel_log_end'].values[ind_outliers] - nasc_df['vessel_log_start'].values[ind_outliers]

        return interval

    def _get_tot_norm_wgt(self, n_A: pd.Series) -> np.ndarray:
        """
        Calculates the total normalized weight
        for each NASC value.

        Parameters
        ----------
        n_A : pd.Series
            Series representing the nautical areal density

        Returns
        -------
        The total normalized weight
        """

        # expand the bio constants dataframe so that it corresponds to nasc_df
        bc_expanded_df = self.bio_const_df.loc[self.nasc_df.stratum_num.values]

        # compute the normalized biomass density for males and females
        nntk_male = np.round(n_A.values * bc_expanded_df.M_prop.values)
        nntk_female = np.round(n_A.values * bc_expanded_df.F_prop.values)

        # compute the normalized weight for males, females, and unsexed
        nwgt_male = nntk_male * bc_expanded_df.len_wgt_M_prod.values
        nwgt_female = nntk_female * bc_expanded_df.len_wgt_F_prod.values
        nwgt_unsexed = (n_A.values - nntk_male - nntk_female) * bc_expanded_df.len_wgt_prod.values

        # TODO: the returned value is called normalized, but it isn't between [0, 1]! Different name or bug?

        # compute the total normalized weight
        return nwgt_male + nwgt_female + nwgt_unsexed

    def _get_age_weight_key(self, df: Union[pd.DataFrame, pd.Series]) -> float:
        """
        Computes the normalized weight of animals
        in the first age bin.

        Parameters
        ----------
        df : pd.DataFrame or pd.Series
            species_df with NaNs dropped for a particular stratum

        Returns
        -------
        A float value corresponding to the normalized
        weight for the first age bin.

        Notes
        -----
        The input ``df`` is often a DataFrame, however, when a subset of
        data is selected it can become a Series.
        """

        # account for the case when df is a Series
        if isinstance(df, pd.Series):
            # get numpy arrays of length, age, and weight
            input_arr_len = np.array([df.length])
            input_arr_age = np.array([df.age])
            input_arr_wgt = np.array([df.weight])

        else:
            # get numpy arrays of length, age, and weight
            input_arr_len = df.length.values
            input_arr_age = df.age.values
            input_arr_wgt = df.weight.values

        # bin the ages
        age_bins_ind = self._get_bin_ind(input_arr_age, self.bio_hake_age_bin)

        # bin those lengths that correspond to the lengths in the first age bin
        len_bin_ind = self._get_bin_ind(input_arr_len[age_bins_ind[0]], self.bio_hake_len_bin)

        # normalized weight for the first age bin
        return np.array([np.sum(input_arr_wgt[age_bins_ind[0][i]]) for i in len_bin_ind]).sum() / input_arr_wgt.sum()

    def _get_age2_wgt_prop(self) -> None:
        """
        Obtains the multiplier for each stratum to be applied to the total
        normalized weight. The weight corresponds to the age 2 weight proportion.
        """

        # TODO: This is necessary to match the Matlab output
        #  in theory this should be done when we load the df,
        #  however, this changes the results slightly.
        spec_drop = self.specimen_df.dropna(how='any')

        # each stratum's multiplier once normalized weight has been calculated
        stratum_ind = spec_drop.index.unique()
        self.age2_wgt_prop_df = pd.DataFrame(columns=['val'], index=stratum_ind, dtype=np.float64)
        for i in stratum_ind:
            self.age2_wgt_prop_df.loc[i].val = 1.0 - self._get_age_weight_key(spec_drop.loc[i])

    def _construct_biomass_table(self, norm_bio_dense: np.array) -> None:
        """
        Constructs self.final_biomass_table, which
        contains the normalized biomass density.

        Parameters
        ----------
        norm_bio_dense : np.array
            Numpy array of normalized biomass density
        """

        # minimal columns to do Jolly Hampton CV on data that has not been kriged
        final_df = self.nasc_df[['latitude', 'longitude', 'stratum_num', 'transect_spacing']].copy()
        final_df["normalized_biomass_density"] = norm_bio_dense

        # TODO: should we include the below values in the final biomass table?
        # calculates the interval for the area calculation
        # final_df["interval"] = self._get_interval(self.nasc_df)

        # calculate the area corresponding to the NASC value
        # final_df["Area"] = interval * self.nasc_df['transect_spacing']

        # calculate the total number of fish in a given area
        # final_df["N_A"] = n_A * A

        # construct GeoPandas DataFrame to simplify downstream processes
        self.final_biomass_table = gpd.GeoDataFrame(final_df,
                                                    geometry=gpd.points_from_xy(final_df.longitude,
                                                                                final_df.latitude))

    def set_class_variables(self, selected_transects: Optional[List] = None) -> None:
        """
        Set class variables corresponding to the Dataframes from ``survey``,
        which hold all necessary data for the biomass calculation.

        Parameters
        ----------
        selected_transects : list or None
            The subset of transects used in the biomass calculation

        Notes
        -----
        If ``selected_transects`` is not ``None``, then all ``survey`` Dataframes will
        be copied, else all ``survey`` Dataframes except ``nasc_df`` will be copied.
        """

        if selected_transects is not None:

            transect_vs_haul = self.survey.haul_to_transect_mapping_df["transect_num"].dropna().astype(int).reset_index().set_index("transect_num")

            # TODO: do a check that all hauls are mapped to a transect

            sel_transects = transect_vs_haul.index.unique().intersection(selected_transects).values
            sel_hauls = transect_vs_haul.loc[sel_transects]["haul_num"].unique()

            sel_hauls_length = self.survey.length_df.index.intersection(sel_hauls).unique()
            sel_haul_strata = self.survey.strata_df.index.get_level_values("haul_num").intersection(sel_hauls).unique()
            sel_haul_specimen = self.survey.specimen_df.index.intersection(sel_hauls).unique()

            self.length_df = self.survey.length_df.loc[sel_hauls_length].copy()
            self.strata_df = self.survey.strata_df.loc[sel_haul_strata].copy()
            self.specimen_df = self.survey.specimen_df.loc[sel_haul_specimen].copy()

            # select nasc data based on haul_num, so we do not select a stratum that is not in length/specimen data
            self.nasc_df = self.survey.nasc_df.loc[sel_transects]

        else:
            self.length_df = self.survey.length_df.copy()
            self.strata_df = self.survey.strata_df.copy()
            self.specimen_df = self.survey.specimen_df.copy()
            self.nasc_df = self.survey.nasc_df

    def get_final_biomass_table(self, selected_transects: Optional[List] = None) -> None:
        """
        Orchestrates the calculation of the normalized
        biomass density and creation of
        self.final_biomass_table, which contains
        the normalized biomass density and associated
        useful variables.

        Parameters
        ----------
        selected_transects : list or None
            The subset of transects used in the biomass calculation
        """

        self.set_class_variables(selected_transects)

        # get the backscattering cross-section for each stratum
        self._get_strata_sig_b()

        # add stratum_num column to length and specimen df and set it as the index
        self._add_stratum_column()

        self._get_biomass_constants()

        self._get_age2_wgt_prop()

        # fill in missing strata constants
        self.strata_sig_b = self._fill_missing_strata_indices(df=self.strata_sig_b.copy())
        self.bio_const_df = self._fill_missing_strata_indices(df=self.bio_const_df.copy())
        self.age2_wgt_prop_df = self._fill_missing_strata_indices(df=self.age2_wgt_prop_df.copy())

        # calculate proportion coefficient for mixed species
        wgt_vals = self.strata_df.reset_index().set_index('haul_num')['fraction_hake']
        wgt_vals_ind = wgt_vals.index
        mix_sa_ratio = self.nasc_df.apply(lambda x: wgt_vals[x.haul_num] if x.haul_num in wgt_vals_ind else 0.0, axis=1)

        # calculate the nautical areal density
        n_A = np.round((mix_sa_ratio*self.nasc_df.NASC) /
                       self.strata_sig_b.loc[self.nasc_df.stratum_num].values)

        # total normalized weight for each n_A value
        nwgt_total = self._get_tot_norm_wgt(n_A)

        # obtain normalized biomass density
        # TODO: the computed value is called normalized, but it isn't between [0, 1]! Different name or bug?
        norm_bio_dense = nwgt_total * self.age2_wgt_prop_df.loc[self.nasc_df.stratum_num.values].values.flatten()

        self._construct_biomass_table(norm_bio_dense)
