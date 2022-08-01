import numpy as np
import pandas as pd
from typing import Tuple, List
import EchoPro.survey as Survey


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

    def __init__(self, survey: Survey = None):

        self.survey = survey

        self.bio_hake_len_bin = survey.params['bio_hake_len_bin']
        self.bio_hake_age_bin = survey.params['bio_hake_age_bin']

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
        Adds the "stratum" column to self.survey.strata_df
        and self.survey.length_df. Additionally, this
        function will set the index to "stratum".
        """

        # get df relating the haul to the stratum
        strata_haul_df = self.survey.strata_df.reset_index()[['Haul', 'stratum']].set_index('Haul')

        # add stratum column to strata_df and set it as the index
        self.survey.specimen_df['stratum'] = strata_haul_df.loc[self.survey.specimen_df.index]
        self.survey.specimen_df.set_index('stratum', inplace=True)

        # add stratum column to length_df and set it as the index
        self.survey.length_df['stratum'] = strata_haul_df.loc[self.survey.length_df.index]
        self.survey.length_df.set_index('stratum', inplace=True)

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

    def _get_norm_len_key_station_1(self, df: pd.DataFrame) -> np.ndarray:
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

    def _get_norm_len_key_station_2(self, df: pd.DataFrame) -> np.ndarray:
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
        total_n = spec_strata_m.shape[0] + spec_strata_f.shape[0] + len_strata.Frequency.sum()
        # total_n = spec_strata.shape[0] + (len_strata.Frequency.sum())  # TODO: This is what it should be

        # proportion of males/females in station 2
        spec_m_prop = spec_strata_m.shape[0] / total_n
        spec_f_prop = spec_strata_f.shape[0] / total_n

        # proportion of males/females in station 1
        len_m_prop = len_strata_m.Frequency.sum() / total_n
        len_f_prop = len_strata_f.Frequency.sum() / total_n

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
        spec_strata_m = spec_stratum[spec_stratum['Sex'] == 1]
        spec_strata_f = spec_stratum[spec_stratum['Sex'] == 2]

        # get lengths in the stratum and split into males and females
        len_strata = length_drop_df.loc[stratum]
        len_strata_m = len_strata[len_strata['Sex'] == 1]
        len_strata_f = len_strata[len_strata['Sex'] == 2]

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

    def _get_biomass_constants(self) -> pd.DataFrame:
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
        spec_strata_ind = self.survey.specimen_df.index.unique()
        len_strata_ind = self.survey.length_df.index.unique()
        strata_ind = spec_strata_ind.intersection(len_strata_ind).values

        # obtain the length-weight key for all specimen data
        len_weight_spec = self._generate_length_val_key(len_name='Length',
                                                        val_name='Weight',
                                                        df=self.survey.specimen_df)

        # select the indices that do not have nan in either Length or Weight
        spec_drop = self.survey.specimen_df.dropna(how='any')

        # select the indices that do not have nan in either Length or Weight
        length_drop_df = self.survey.length_df.dropna(how='any')

        # initialize dataframe that will hold all important calculated constants
        bio_const_df = pd.DataFrame(columns=['M_prop', 'F_prop', 'len_wgt_prod',
                                             'len_wgt_M_prod', 'len_wgt_F_prod'],
                                    index=strata_ind, dtype=np.float64)

        # for each stratum compute the necessary constants
        for stratum in strata_ind:
            bio_const_df = self._fill_len_wgt_prod(bio_const_df, stratum, spec_drop,
                                                   length_drop_df, len_weight_spec)

        return bio_const_df

    @staticmethod
    def _get_interval(nasc_df: pd.DataFrame) -> np.ndarray:
        """
        Calculates the interval needed for the area calculation.

        Parameters
        ----------
        nasc_df : pd.DataFrame
            NASC df with 'VL start' and 'VL end' columns

        Returns
        -------
        interval : np.ndarray
            Intervals needed for area calculation
        """

        # calculate interval for all values except for the last interval
        interval = (nasc_df['VL start'].iloc[1:].values - nasc_df['VL start'].iloc[:-1].values)

        # calculate last interval
        last_interval = nasc_df['VL end'].iloc[-1] - nasc_df['VL start'].iloc[-1]

        # combines all intervals
        interval = np.concatenate([interval, np.array([last_interval])])

        # get median interval, so we can use it to remove outliers
        median_interval = np.median(interval)

        # remove outliers at the end of the transect
        ind_outliers = np.argwhere(np.abs(interval - median_interval) > 0.05).flatten()
        interval[ind_outliers] = nasc_df['VL end'].values[ind_outliers] - nasc_df['VL start'].values[ind_outliers]

        return interval

    def _get_tot_norm_wgt(self, bc_df: pd.DataFrame, n_A: pd.Series) -> np.ndarray:
        """
        Calculates the total normalized weight
        for each NASC value.

        Parameters
        ----------
        bc_df : pd.Dataframe
            Dataframe of biomass constants for each stratum
        n_A : pd.Series
            Series representing the nautical areal density

        Returns
        -------
        The total normalized weight
        """

        # expand the bio constants dataframe so that it corresponds to nasc_df
        bc_expanded_df = bc_df.loc[self.survey.nasc_df.Stratum.values]

        # compute the normalized biomass density for males and females
        nntk_male = np.round(n_A.values * bc_expanded_df.M_prop.values)
        nntk_female = np.round(n_A.values * bc_expanded_df.F_prop.values)

        # compute the normalized weight for males, females, and unsexed
        nwgt_male = nntk_male * bc_expanded_df.len_wgt_M_prod.values
        nwgt_female = nntk_female * bc_expanded_df.len_wgt_F_prod.values
        nwgt_unsexed = (n_A.values - nntk_male - nntk_female) * bc_expanded_df.len_wgt_prod.values

        # compute the total normalized weight
        return nwgt_male + nwgt_female + nwgt_unsexed

    def _get_age_weight_key(self, df: pd.DataFrame) -> float:
        """
        Computes the normalized weight of animals
        in the first age bin.

        Parameters
        ----------
        df : pd.DataFrame
            species_df with NaNs dropped for a particular stratum

        Returns
        -------
        A float value corresponding to the normalized
        weight for the first age bin.
        """

        # get numpy arrays of length, age, and weight
        input_arr_len = df.Length.values
        input_arr_age = df.Age.values
        input_arr_wgt = df.Weight.values

        # bin the ages
        age_bins_ind = self._get_bin_ind(input_arr_age, self.bio_hake_age_bin)

        # bin those lengths that correspond to the lengths in the first age bin
        len_bin_ind = self._get_bin_ind(input_arr_len[age_bins_ind[0]], self.bio_hake_len_bin)

        # normalized weight for the first age bin
        return np.array([np.sum(input_arr_wgt[age_bins_ind[0][i]]) for i in len_bin_ind]).sum() / input_arr_wgt.sum()

    def _get_normalized_biomass_density(self, nwgt_total: np.ndarray) -> np.ndarray:
        """
        Computes and returns the normalized
        biomass density.

        Parameters
        ----------
        nwgt_total : np.array
            Total normalized weight for each n_A value

        Returns
        -------
        Numpy array of normalized biomass density
        """

        # TODO: This is necessary to match the Matlab output
        #  in theory this should be done when we load the df,
        #  however, this changes the results slightly.
        spec_drop = self.survey.specimen_df.dropna(how='any')

        # each stratum's multiplier once normalized weight has been calculated
        stratum_ind = spec_drop.index.unique()
        age2_wgt_prop_df = pd.DataFrame(columns=['val'], index=stratum_ind, dtype=np.float64)
        for i in stratum_ind:
            age2_wgt_prop_df.loc[i].val = 1.0 - self._get_age_weight_key(spec_drop.loc[i])

        # normalized biomass density
        return nwgt_total * age2_wgt_prop_df.loc[self.survey.nasc_df.Stratum.values].values.flatten()

    def _construct_biomass_table(self, norm_bio_dense: np.array) -> None:
        """
        Constructs self.survey.final_biomass_table, which
        contains the normalized biomass density.

        Parameters
        ----------
        norm_bio_dense : np.array
            Numpy array of normalized biomass density
        """

        # minimal columns to do Jolly Hampton CV on data that has not been kriged
        self.survey.final_biomass_table = self.survey.nasc_df[['Latitude', 'Longitude', 'Stratum', 'Spacing']].copy()
        self.survey.final_biomass_table["normalized_biomass_density"] = norm_bio_dense

        # TODO: should we include the below values in the final biomass table?
        # calculates the interval for the area calculation
        # self.survey.final_biomass_table["interval"] = self._get_interval(self.survey.nasc_df)

        # calculate the area corresponding to the NASC value
        # self.survey.final_biomass_table["Area"] = interval * self.survey.nasc_df['Spacing']

        # calculate the total number of fish in a given area
        # self.survey.final_biomass_table["N_A"] = n_A * A

    def get_final_biomass_table(self) -> None:
        """
        Orchestrates the calculation of the normalized
        biomass density and creation of
        self.survey.final_biomass_table, which contains
        the normalized biomass density and associated
        useful variables.
        """

        # add stratum column to length and specimen df and set it as the index
        self._add_stratum_column()

        bc_df = self._get_biomass_constants()

        # calculate proportion coefficient for mixed species
        wgt_vals = self.survey.strata_df.reset_index().set_index('Haul')['wt']
        wgt_vals_ind = wgt_vals.index
        mix_sa_ratio = self.survey.nasc_df.apply(lambda x: wgt_vals[x.Haul] if x.Haul in wgt_vals_ind else 0.0, axis=1)

        # calculate the nautical areal density
        n_A = np.round((mix_sa_ratio*self.survey.nasc_df.NASC) /
                       self.survey.strata_sig_b.loc[self.survey.nasc_df.Stratum].values)

        # total normalized weight for each n_A value
        nwgt_total = self._get_tot_norm_wgt(bc_df, n_A)

        norm_bio_dense = self._get_normalized_biomass_density(nwgt_total)

        self._construct_biomass_table(norm_bio_dense)
