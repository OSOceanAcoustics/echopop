from typing import List, Optional, Tuple, Union

import geopandas as gpd
import numpy as np
import pandas as pd

from ..utils.binning import get_bin_ind


class ComputeTransectVariables:
    """
    A class that computes several variables (e.g. abundance, biomass)
    describing the animal population at the provided transect points
    using NASC values and other associated data.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        self.survey will also change this object.
    """

    def __init__(self, survey=None):

        self.survey = survey

        self.bio_hake_len_bin = survey.params["bio_hake_len_bin"]
        self.bio_hake_age_bin = survey.params["bio_hake_age_bin"]

        # initialize class variables that will be set downstream
        self.length_df = None
        self.strata_df = None
        self.specimen_df = None
        self.nasc_df = None
        self.transect_results_gdf = None
        self.transect_results_male_gdf = None
        self.transect_results_female_gdf = None
        self.kriging_results_gdf = None
        self.kriging_results_male_gdf = None
        self.kriging_results_female_gdf = None
        self.bio_param_df = None  # biomass parameters for each stratum
        self.weight_fraction_adult_df = None
        self.weight_fraction_all_ages_df = None
        self.weight_fraction_all_ages_female_df = None
        self.weight_fraction_all_ages_male_df = None
        self.num_fraction_adult_df = None
        self.strata_sig_b = None
        self.transect_bin_abundance_male_df = None
        self.transect_bin_abundance_female_df = None
        self.transect_bin_abundance_df = None
        self.transect_bin_biomass_male_df = None
        self.transect_bin_biomass_female_df = None
        self.transect_bin_biomass_df = None
        self.kriging_bin_abundance_male_df = None
        self.kriging_bin_abundance_female_df = None
        self.kriging_bin_abundance_df = None
        self.kriging_bin_biomass_male_df = None
        self.kriging_bin_biomass_female_df = None
        self.kriging_bin_biomass_df = None
        self.all_strata = None
        self.missing_strata = None
        self.percentage_transects_selected = None
        self.sel_tran_strata_choice = dict()
        self.stratum_choices = dict()
        self.strata_sig_b_df = None
        self.specimen_all_df = None
        self.bin_ds = None
        self.jollyhampton_cv = None

    def _get_strata_sig_b(self) -> None:
        """
        Computes the backscattering cross-section (sigma_b),
        using the strata, specimen, and length dataframes.
        These values are then stored in self.strata_sig_b
        as a Pandas series with index "stratum_num".
        """

        # TODO: the target strength functions are specific to Hake, replace with input in the future

        # initialize sig_bs_haul column in strata_sig_b_df
        self.strata_sig_b_df["sig_bs_haul"] = np.nan

        # select the indices that do not have nan in either length or weight
        spec_df = self.specimen_df[["length", "weight"]].copy()
        spec_df = spec_df.dropna(how="any")

        for haul_num in spec_df.index.unique():

            # lengths from specimen file associated with index haul_num
            spec_len = spec_df.loc[haul_num]["length"]

            if haul_num in self.length_df.index:

                # add lengths from length file associated with index haul_num
                length_len = self.length_df.loc[haul_num]["length"].values
                length_count = self.length_df.loc[haul_num]["length_count"].values

                # empirical relation for target strength
                TS0j_length = 20.0 * np.log10(length_len) - 68.0

                # sum of target strengths
                sum_TS0j_length = np.nansum(
                    (10.0 ** (TS0j_length / 10.0)) * length_count
                )

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
            self.strata_sig_b_df.loc[haul_num, "sig_bs_haul"] = (
                sum_TS0j_spec + sum_TS0j_length
            ) / (num_length + TS0j_spec.size)

        # mean backscattering cross-section for each stratum
        self.strata_sig_b = (
            4.0
            * np.pi
            * self.strata_sig_b_df["sig_bs_haul"].groupby("stratum_num").mean()
        )

        # include placeholder for strata without values
        self.missing_strata = []
        for stratum in self.all_strata:
            if stratum not in self.strata_sig_b.index:
                self.missing_strata.append(stratum)
                self.strata_sig_b[stratum] = np.nan

    def set_strata_for_missing_strata(self) -> None:
        """
        Constructs and sets the dictionary ``self.sel_tran_strata_choice``
        with keys as values in ``self.missing_strata`` and values
        as a list with the first and second element corresponding
        to the stratum that are closest and less than or greater
        than the missing stratum, respectively.
        """

        # construct array of all known strata
        known_strata_arr = np.array(
            [i for i in self.all_strata if i not in self.missing_strata]
        )

        # determine the strata that should replace the missing strata
        for m_strat in self.missing_strata:

            # get bool array of values less than m_strat
            less_than_m_strat = known_strata_arr < m_strat
            greater_than_m_strat = known_strata_arr > m_strat

            # assign stratum values
            if not any(less_than_m_strat):
                new_stratum_g = min(known_strata_arr[greater_than_m_strat])
                new_stratum_l = None

            elif not any(greater_than_m_strat):
                new_stratum_l = max(known_strata_arr[less_than_m_strat])
                new_stratum_g = None

            else:
                new_stratum_g = min(known_strata_arr[greater_than_m_strat])
                new_stratum_l = max(known_strata_arr[less_than_m_strat])

            # store stratum values for m_strat
            self.sel_tran_strata_choice[m_strat] = [new_stratum_l, new_stratum_g]

    def set_stratum_choice(self) -> None:
        """
        Constructs and set the dictionary ``self.stratum_choices``,
        which is a dictionary the specifies what strata or stratum
        should be used to select data (e.g. specimen, length DataFrames).
        This routine is necessary to ensure that transect selection can
        be done easily.
        """

        # fill in strata or stratum that should be used for all strata
        for stratum in self.all_strata:

            if stratum not in self.missing_strata:
                stratum_choice = stratum
            else:
                stratum_choice = [
                    s for s in self.sel_tran_strata_choice[stratum] if s is not None
                ]

            # assign stratum choices for stratum value
            self.stratum_choices[stratum] = stratum_choice

    def _fill_missing_strata_sig_b(self) -> None:
        """
        Fills in missing strata data for the DataFrame
        ``self.strata_sig_b`` using criteria established
        in EchoPro Matlab.
        """

        for m_strat in self.missing_strata:

            # obtain less and greater than stratum with respect to missing stratum
            less_than_m_strat, greater_than_m_strat = self.sel_tran_strata_choice[
                m_strat
            ]

            if (less_than_m_strat is None) and (greater_than_m_strat is not None):
                # replace missing value with value at next filled stratum greater than m_strat
                self.strata_sig_b.loc[m_strat] = self.strata_sig_b.loc[
                    greater_than_m_strat
                ]

            elif (less_than_m_strat is not None) and (greater_than_m_strat is None):
                # replace missing value with value at next filled stratum less than m_strat
                self.strata_sig_b.loc[m_strat] = self.strata_sig_b.loc[
                    less_than_m_strat
                ]

            else:
                # replace missing value with average of two closest filled strata
                self.strata_sig_b.loc[m_strat] = (
                    self.strata_sig_b.loc[less_than_m_strat]
                    + self.strata_sig_b.loc[greater_than_m_strat]
                ) / 2.0

    def _add_stratum_column(self) -> None:
        """
        Adds the ``stratum_num`` column to self.strata_df
        and self.length_df. Additionally, this
        function will set the index to ``stratum_num``.
        """

        # get df relating the haul to the stratum
        strata_haul_df = self.strata_df.reset_index()[
            ["haul_num", "stratum_num"]
        ].set_index("haul_num")

        # add stratum_num column to specimen_df and set it as the index
        self.specimen_df["stratum_num"] = strata_haul_df.loc[self.specimen_df.index]
        self.specimen_df.set_index("stratum_num", inplace=True)

        if self.percentage_transects_selected is not None:
            # TODO: this is necessary to mimic the Matlab code (may be able to optimize this)
            self.specimen_all_df["stratum_num"] = strata_haul_df.loc[
                self.specimen_all_df.index
            ]
            self.specimen_all_df.set_index("stratum_num", inplace=True)

        # add stratum_num column to length_df and set it as the index
        self.length_df["stratum_num"] = strata_haul_df.loc[self.length_df.index]
        self.length_df.set_index("stratum_num", inplace=True)

    def _generate_length_val_conversion(
        self, len_name: str, val_name: str, df: pd.DataFrame = None
    ) -> np.ndarray:
        """
        Generates a length-to-value conversion by
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
            1D array representing the length-to-value conversion
        """

        # select the indices that do not have nan in either length or value
        len_wgt_nonull = np.logical_and(df[len_name].notnull(), df[val_name].notnull())
        df_no_null = df.loc[len_wgt_nonull]

        # get numpy version of those entries that are not NaN
        L = df_no_null[len_name].values
        V = df_no_null[val_name].values

        # binned length indices
        len_bin_ind = get_bin_ind(L, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([i.shape[0] for i in len_bin_ind])

        # length-value regression
        p = np.polyfit(np.log10(L), np.log10(V), 1)

        # obtain linear regression parameters
        reg_w0 = 10.0 ** p[1]
        reg_p = p[0]

        # value at length or length-value-key.
        # length-value-key per fish over entire survey region (an array)
        len_val_reg = reg_w0 * self.bio_hake_len_bin**reg_p

        # set length-value key to the mean of the values in the bin
        len_val_key = np.array(
            [np.mean(V[ind]) if ind.size > 0 else 0.0 for ind in len_bin_ind]
        )

        # replace those bins with less than 5 samples with the regression value
        less_five_ind = np.argwhere(len_bin_cnt < 5).flatten()
        len_val_key[less_five_ind] = len_val_reg[less_five_ind]

        return len_val_key

    def _get_distribution_lengths_station_1(self, df: pd.DataFrame) -> np.ndarray:
        """
        Computes the length distribution from
        data obtained from station 1 i.e. data
        that tells you how many fish are of a
        particular length.

        Parameters
        ----------
        df : pd.DataFrame
            Data from station 1 that has columns ``length`` and ``length_count``.

        Returns
        -------
        A numpy array of the length distribution i.e. the
        count of each bin divided by the total of all bin counts
        """

        # get numpy arrays of the length and length_count columns
        length_arr = df.length.values
        length_count_arr = df.length_count.values

        # binned length indices
        len_ind = get_bin_ind(length_arr, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([np.sum(length_count_arr[i]) for i in len_ind])

        return len_bin_cnt / np.sum(len_bin_cnt)

    def _get_distribution_lengths_station_2(self, df: pd.DataFrame) -> np.ndarray:
        """
        Computes the length distribution from
        data obtained from station 2 i.e. data
        that does not have a frequency associated
        with it.

        Parameters
        ----------
        df : pd.DataFrame
            Data from fish measurement station 2 that has a ``length`` column

        Returns
        -------
        A numpy array of the length distribution i.e. the
        count of each bin divided by the total of all bin counts
        """

        # numpy array of length column
        length_arr = df.length.values

        # binned length indices
        len_bin_ind = get_bin_ind(length_arr, self.bio_hake_len_bin)

        # total number of lengths in a bin
        len_bin_cnt = np.array([i.shape[0] for i in len_bin_ind])

        return len_bin_cnt / np.sum(len_bin_cnt)

    @staticmethod
    def _compute_proportions(
        spec_strata_m: pd.DataFrame,
        spec_strata_f: pd.DataFrame,
        len_strata: pd.DataFrame,
        len_strata_m: pd.DataFrame,
        len_strata_f: pd.DataFrame,
    ) -> Tuple[list, list, list, list]:
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
        total_n = (
            spec_strata_m.shape[0]
            + spec_strata_f.shape[0]
            + len_strata.length_count.sum()
        )
        # total_n = spec_strata.shape[0] + (len_strata.length_count.sum())  # TODO: This is what it should be  # noqa

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
        tot_prop1 = tot_prop1 / (
            tot_prop1 + tot_prop2
        )  # TODO: do we need to calculate this? Denominator will be 1

        # average of total proportion of sexed fish pertaining to station 2
        tot_prop2 = tot_prop2 / (
            tot_prop1 + tot_prop2
        )  # TODO: do we need to calculate this? Denominator will be 1

        # list of male and female proportions for both stations
        gender_prop = [spec_m_prop + len_m_prop, spec_f_prop + len_f_prop]

        # list of average fraction of sexed fish from station 1
        fac1 = [fac1_m, fac1_f]

        # list of average fraction of sexed fish from station 2
        fac2 = [fac2_m, fac2_f]

        # list of average of total proportion of sexed fish
        tot_prop = [tot_prop1, tot_prop2]

        return gender_prop, fac1, fac2, tot_prop

    def _fill_averaged_weight(
        self,
        bio_param_df: pd.DataFrame,
        stratum: int,
        spec_drop_df: pd.DataFrame,
        length_drop_df: pd.DataFrame,
        length_to_weight_conversion: np.array,
    ) -> pd.DataFrame:
        """
        Fills in the biomass parameter dataframe for a
        particular stratum.

        Parameters
        ----------
        bio_param_df : pd.DataFrame
            Biomass parameter dataframe to fill
        stratum : int
            Stratum to fill
        spec_drop_df : pd.DataFrame
            specimen_df with NaN values dropped, corresponding to ``stratum``
        length_drop_df : pd.DataFrame
            length_df with NaN values dropped, corresponding to ``stratum``
        length_to_weight_conversion : np.array
            length-to-weight conversion (i.e. an array that contains the corresponding
            weight of the length bins) for all specimen data

        Returns
        -------
        bio_param_df : pd.DataFrame
            Biomass parameter dataframe with stratum filled in
        """

        # get specimen in the stratum and split into males and females
        spec_stratum = spec_drop_df
        spec_strata_m = spec_stratum[spec_stratum["sex"] == 1]
        spec_strata_f = spec_stratum[spec_stratum["sex"] == 2]

        # get lengths in the stratum and split into males and females
        len_strata = length_drop_df
        len_strata_m = len_strata[len_strata["sex"] == 1]
        len_strata_f = len_strata[len_strata["sex"] == 2]

        # get the distribution lengths for station 1
        distribution_length_s1 = self._get_distribution_lengths_station_1(len_strata)
        distribution_length_m_s1 = self._get_distribution_lengths_station_1(
            len_strata_m
        )
        distribution_length_f_s1 = self._get_distribution_lengths_station_1(
            len_strata_f
        )

        # get the distribution lengths for station 2
        # TODO: this is what should be done!
        # distribution_length_s2 = self._get_distribution_lengths_station_2(spec_stratum)
        spec_stratum_mf = pd.concat([spec_strata_m, spec_strata_f], axis=0)
        distribution_length_s2 = self._get_distribution_lengths_station_2(
            spec_stratum_mf
        )
        distribution_length_m_s2 = self._get_distribution_lengths_station_2(
            spec_strata_m
        )
        distribution_length_f_s2 = self._get_distribution_lengths_station_2(
            spec_strata_f
        )

        gender_prop, fac1, fac2, tot_prop = self._compute_proportions(
            spec_strata_m, spec_strata_f, len_strata, len_strata_m, len_strata_f
        )

        # fill df with bio parameters needed for biomass density calc
        bio_param_df.M_prop.loc[stratum] = gender_prop[0]
        bio_param_df.F_prop.loc[stratum] = gender_prop[1]
        bio_param_df.averaged_weight.loc[stratum] = np.dot(
            tot_prop[0] * distribution_length_s1 + tot_prop[1] * distribution_length_s2,
            length_to_weight_conversion,
        )
        bio_param_df.averaged_weight_M.loc[stratum] = np.dot(
            fac1[0] * distribution_length_m_s1 + fac2[0] * distribution_length_m_s2,
            length_to_weight_conversion,
        )
        bio_param_df.averaged_weight_F.loc[stratum] = np.dot(
            fac1[1] * distribution_length_f_s1 + fac2[1] * distribution_length_f_s2,
            length_to_weight_conversion,
        )

        return bio_param_df

    def _get_biomass_parameters(self) -> None:
        """
        Obtains the parameters associated with each stratum,
        which are used in the biomass density calculation.
        Specifically, we obtain the following parameters for
        each stratum:
        * M_prop -- proportion of males
        * F_prop -- proportion of females
        * averaged_weight-- the averaged weight for the whole population
        * averaged_weight_M -- averaged_weight for the male population
        * averaged_weight_F -- averaged_weight for the female population

        Notes
        -----
        The following class variable is created in this function:
        bio_param_df : pd.Dataframe
            Dataframe with index of stratum and columns
            corresponding to the parameters specified
            above.
        """

        # obtain the length-to-weight conversion for all specimen data
        # TODO: to match Matlab version of selection of transects we use
        #  self.specimen_all_df instead of self.specimen_df
        length_to_weight_conversion_spec = self._generate_length_val_conversion(
            len_name="length", val_name="weight", df=self.specimen_all_df
        )

        # select the indices that do not have nan in either Length or Weight
        spec_drop = self.specimen_df.dropna(how="any")

        # select the indices that do not have nan in either Length or Weight
        length_drop_df = self.length_df.dropna(how="any")

        # initialize dataframe that will hold all important calculated parameters
        bio_param_df = pd.DataFrame(
            columns=[
                "M_prop",
                "F_prop",
                "averaged_weight",
                "averaged_weight_M",
                "averaged_weight_F",
            ],
            index=self.all_strata,
            dtype=np.float64,
        )

        # for each stratum compute the necessary parameters
        for stratum in self.all_strata:

            # appropriately select data for stratum
            spec_in = spec_drop.loc[self.stratum_choices[stratum]]
            length_in = length_drop_df.loc[self.stratum_choices[stratum]]

            bio_param_df = self._fill_averaged_weight(
                bio_param_df,
                stratum,
                spec_in,
                length_in,
                length_to_weight_conversion_spec,
            )

        self.bio_param_df = bio_param_df

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
        interval = (
            nasc_df["vessel_log_start"].iloc[1:].values
            - nasc_df["vessel_log_start"].iloc[:-1].values
        )

        # calculate last interval
        last_interval = (
            nasc_df["vessel_log_end"].iloc[-1] - nasc_df["vessel_log_start"].iloc[-1]
        )

        # combines all intervals
        interval = np.concatenate([interval, np.array([last_interval])])

        # get median interval, so we can use it to remove outliers
        median_interval = np.median(interval)

        # remove outliers at the end of the transect
        ind_outliers = np.argwhere(np.abs(interval - median_interval) > 0.05).flatten()
        interval[ind_outliers] = (
            nasc_df["vessel_log_end"].values[ind_outliers]
            - nasc_df["vessel_log_start"].values[ind_outliers]
        )

        return interval

    def _get_age_weight_num_proportions(
        self, df: Union[pd.DataFrame, pd.Series]
    ) -> Tuple[float, float]:
        """
        Computes the proportion of animals in a provided age bin for
        both the weight and number of animals.

        Parameters
        ----------
        df : pd.DataFrame or pd.Series
            species_df with NaNs dropped for a particular stratum

        Returns
        -------
        age_len_prop: float
            The length proportion for a particular age
        age_wgt_prop
            The weight proportion for a particular age

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
        age_bins_ind = get_bin_ind(input_arr_age, self.bio_hake_age_bin)

        # bin those lengths that correspond to the lengths in the first age bin
        len_bin_ind = get_bin_ind(input_arr_len[age_bins_ind[0]], self.bio_hake_len_bin)

        # the length proportion of animals for a given age
        age_len_prop = len(input_arr_len[age_bins_ind[0]]) / len(input_arr_len)

        # the weight proportion of animals for a given age
        age_wgt_prop = (
            np.array(
                [np.sum(input_arr_wgt[age_bins_ind[0][i]]) for i in len_bin_ind]
            ).sum()
            / input_arr_wgt.sum()
        )

        # length and weight for the first age bin
        return age_len_prop, age_wgt_prop

    def _get_all_age_weight_proportions(
        self, df: Union[pd.DataFrame, pd.Series], age_bin_ind: int
    ) -> float:
        """
        Computes the proportion of animals in a provided age bin for
        the weight of animals.

        Parameters
        ----------
        df : pd.DataFrame or pd.Series
            species_df with NaNs dropped for a particular stratum
        age_bin_ind: int
            The age bin index to calculate information for

        Returns
        -------
        age_wgt_prop
            The weight proportion for a particular age

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
        age_bins_ind = get_bin_ind(input_arr_age, self.bio_hake_age_bin)

        # bin those lengths that correspond to the lengths in the given age bin
        len_bin_ind = get_bin_ind(
            input_arr_len[age_bins_ind[age_bin_ind]], self.bio_hake_len_bin
        )

        if self.survey.params["exclude_age1"] is True:

            # return 0.0, since no data should be included here
            if age_bin_ind == 0:
                return 0.0

            # bin those lengths that correspond to the lengths in the first age bin
            len_bin_ind_0 = get_bin_ind(
                input_arr_len[age_bins_ind[0]], self.bio_hake_len_bin
            )

            # get the weight of the first age bin
            wgt_age_0 = np.array(
                [np.sum(input_arr_wgt[age_bins_ind[0][i]]) for i in len_bin_ind_0]
            ).sum()

            # get total weight minus the first age bin weight
            denominator_wgt = input_arr_wgt.sum() - wgt_age_0
        else:

            # get total weight
            denominator_wgt = input_arr_wgt.sum()

        # the weight of animals in a given age bin
        numerator_wgt = np.array(
            [np.sum(input_arr_wgt[age_bins_ind[age_bin_ind][i]]) for i in len_bin_ind]
        ).sum()

        # the weight proportion of animals for a given age
        if denominator_wgt != 0.0:
            age_wgt_prop = numerator_wgt / denominator_wgt
        else:
            age_wgt_prop = numerator_wgt

        # weight for the given age bin
        return age_wgt_prop

    def _get_weight_num_fraction_adult(self) -> None:
        """
        Obtains the multipliers for each stratum to be applied to the total
        areal biomass density and abundance. The values correspond to the
        age 2 fractions.
        """

        # TODO: This is necessary to match the Matlab output
        #  in theory this should be done when we load the df,
        #  however, this changes the results slightly.
        spec_drop = self.specimen_df.dropna(how="any")

        # each stratum's multiplier once areal biomass density has been calculated
        self.weight_fraction_adult_df = pd.DataFrame(
            columns=["val"], index=self.all_strata, dtype=np.float64
        )
        self.num_fraction_adult_df = pd.DataFrame(
            columns=["val"], index=self.all_strata, dtype=np.float64
        )

        for stratum in self.all_strata:

            # select specimen data
            spec_in = spec_drop.loc[self.stratum_choices[stratum]]

            age_len_prop, age_wgt_prop = self._get_age_weight_num_proportions(spec_in)

            self.weight_fraction_adult_df.loc[stratum].val = abs(1.0 - age_wgt_prop)
            self.num_fraction_adult_df.loc[stratum].val = abs(1.0 - age_len_prop)

    def _get_weight_fraction_all_ages(self) -> None:
        """
        Obtains the multipliers for each stratum to be applied to the total
        areal biomass density and abundance. The values correspond to all age bins.
        """

        # TODO: This is necessary to match the Matlab output
        #  in theory this should be done when we load the df,
        #  however, this changes the results slightly.
        spec_drop = self.specimen_df.dropna(how="any")

        # obtain the male and female entries of spec_drop
        spec_drop_M = spec_drop[spec_drop["sex"] == 1]
        spec_drop_F = spec_drop[spec_drop["sex"] == 2]

        # the number of age bins
        bin_length = len(self.bio_hake_age_bin)

        # each stratum's multiplier once areal biomass density has been calculated
        self.weight_fraction_all_ages_df = pd.DataFrame(
            columns=["age_bin_" + str(i + 1) for i in range(bin_length)],
            index=self.all_strata,
            dtype=np.float64,
        )
        self.weight_fraction_all_ages_male_df = pd.DataFrame(
            columns=["age_bin_" + str(i + 1) for i in range(bin_length)],
            index=self.all_strata,
            dtype=np.float64,
        )
        self.weight_fraction_all_ages_female_df = pd.DataFrame(
            columns=["age_bin_" + str(i + 1) for i in range(bin_length)],
            index=self.all_strata,
            dtype=np.float64,
        )

        for stratum in self.all_strata:

            # select specimen data
            spec_in = spec_drop.loc[self.stratum_choices[stratum]]
            spec_in_M = spec_drop_M.loc[self.stratum_choices[stratum]]
            spec_in_F = spec_drop_F.loc[self.stratum_choices[stratum]]

            # obtain the weight fraction for all age bins and a given stratum
            for j in range(bin_length):
                age_wgt_prop = self._get_all_age_weight_proportions(spec_in, j)
                self.weight_fraction_all_ages_df.loc[stratum][
                    "age_bin_" + str(j + 1)
                ] = age_wgt_prop

                age_wgt_prop_M = self._get_all_age_weight_proportions(spec_in_M, j)
                self.weight_fraction_all_ages_male_df.loc[stratum][
                    "age_bin_" + str(j + 1)
                ] = age_wgt_prop_M

                age_wgt_prop_F = self._get_all_age_weight_proportions(spec_in_F, j)
                self.weight_fraction_all_ages_female_df.loc[stratum][
                    "age_bin_" + str(j + 1)
                ] = age_wgt_prop_F

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

        # set flag to determine if all transects have been selected
        all_transects_selected = True

        # list containing all unique transects present in NASC data
        all_transects = list(self.survey.nasc_df.index.unique())

        if selected_transects is not None:

            # sort list of transects so they can be compared
            selected_transects.sort()
            all_transects.sort()

            # compare transects
            all_transects_selected = selected_transects == all_transects

        if (selected_transects is not None) and (not all_transects_selected):

            # calculate the percentage of transects selected
            self.percentage_transects_selected = len(selected_transects) / len(
                all_transects
            )

            # get mapping between transects and hauls
            # note: the index will have duplicates
            transect_vs_haul = (
                self.survey.haul_to_transect_mapping_df
                .reset_index()
                .set_index("transect_num")
            )

            # get transects and hauls based off of mapping and selected transects
            sel_transects = (
                transect_vs_haul.index.unique().intersection(selected_transects).values
            )
            sel_hauls = transect_vs_haul.loc[sel_transects]["haul_num"].unique()

            # obtain the hauls to use for the length, strata, and specimen DataFrames
            sel_hauls_length = self.survey.length_df.index.intersection(
                sel_hauls
            ).unique()
            sel_haul_strata = (
                self.survey.strata_df.index.get_level_values("haul_num")
                .intersection(sel_hauls)
                .unique()
            )
            sel_haul_specimen = self.survey.specimen_df.index.intersection(
                sel_hauls
            ).unique()

            # select a subset of length, strata, and specimen data
            self.length_df = self.survey.length_df.loc[sel_hauls_length].copy()
            self.strata_sig_b_df = self.survey.strata_df.loc[sel_haul_strata].copy()
            self.specimen_df = self.survey.specimen_df.loc[sel_haul_specimen].copy()

            # select nasc data using the user provided selected transects
            self.nasc_df = self.survey.nasc_df.loc[selected_transects].copy()

            # set strata and specimen DataFrames that contain the full set of Data
            # TODO: set variables containing all data to match Matlab output
            self.strata_df = self.survey.strata_df.copy()
            self.specimen_all_df = self.survey.specimen_df.copy()

        else:
            self.percentage_transects_selected = None
            self.length_df = self.survey.length_df.copy()
            self.strata_df = self.survey.strata_df.copy()
            self.specimen_df = self.survey.specimen_df.copy()
            self.nasc_df = self.survey.nasc_df
            self.strata_sig_b_df = self.strata_df
            self.specimen_all_df = self.specimen_df

    def _set_numerical_density(self, bc_expanded_df: pd.DataFrame) -> None:
        """
        Calculates numerical density variables (such as numerical density
        for males and females) and then assigns them to the appropriate
        DataFrames.

        Parameters
        ----------
        bc_expanded_df: pd.DataFrame
            An expanded bio parameters dataframe

        Notes
        -----
        This function does not return anything, instead, the created
        variables are directly added to class variable DataFrames.
        """

        # calculate the areal numerical density
        self.transect_results_gdf["numerical_density"] = np.round(
            (self.mix_sa_ratio * self.nasc_df.NASC)
            / self.strata_sig_b.loc[self.nasc_df.stratum_num].values
        )

        # compute the areal numerical density for males and females
        self.transect_results_male_gdf["numerical_density"] = np.round(
            self.transect_results_gdf["numerical_density"].values
            * bc_expanded_df.M_prop.values
        )
        self.transect_results_female_gdf["numerical_density"] = np.round(
            self.transect_results_gdf["numerical_density"].values
            * bc_expanded_df.F_prop.values
        )

        # compute areal numerical density for adults
        self.transect_results_gdf["numerical_density_adult"] = (
            self.transect_results_gdf["numerical_density"]
            * self.num_fraction_adult_df.loc[self.nasc_df.stratum_num].values.flatten()
        )

    def _set_biomass_density(self, bc_expanded_df: pd.DataFrame) -> None:
        """
        Calculates total areal biomass density variables for each NASC value
        and then assigns them to the appropriate DataFrames.

        Parameters
        ----------
        bc_expanded_df: pd.DataFrame
            An expanded bio parameters dataframe

        Notes
        -----
        This function does not return anything, instead, the created
        variables are directly added to the class variable DataFrames.
        """

        # compute the areal biomass density for males, females, and unsexed
        self.transect_results_male_gdf["biomass_density"] = (
            self.transect_results_male_gdf["numerical_density"]
            * bc_expanded_df.averaged_weight_M.values
        )
        self.transect_results_female_gdf["biomass_density"] = (
            self.transect_results_female_gdf["numerical_density"]
            * bc_expanded_df.averaged_weight_F.values
        )
        biomass_density_unsexed = (
            self.transect_results_gdf["numerical_density"]
            - self.transect_results_male_gdf["numerical_density"]
            - self.transect_results_female_gdf["numerical_density"]
        ) * bc_expanded_df.averaged_weight.values

        # compute the total areal biomass density
        self.transect_results_gdf["biomass_density"] = (
            self.transect_results_male_gdf["biomass_density"]
            + self.transect_results_female_gdf["biomass_density"]
            + biomass_density_unsexed
        )

        # compute the total biomass density for adults
        self.transect_results_gdf["biomass_density_adult"] = (
            self.transect_results_gdf["biomass_density"]
            * self.weight_fraction_adult_df.loc[
                self.nasc_df.stratum_num
            ].values.flatten()
        )

    def _set_abundance(self, bc_expanded_df: pd.DataFrame) -> None:
        """
        Calculates abundance variables for each NASC value and then assigns
        them to the appropriate DataFrames.

        Parameters
        ----------
        bc_expanded_df: pd.DataFrame
            An expanded bio parameters dataframe

        Notes
        -----
        This function does not return anything, instead, the created
        variables are directly added to the class variable DataFrames.
        """

        # calculate the abundance in a given area
        self.transect_results_gdf["abundance"] = (
            self.mix_sa_ratio
            * self.nasc_df.NASC
            * self.transect_results_gdf["interval_area_nmi2"]
        ) / self.strata_sig_b.loc[self.nasc_df.stratum_num].values

        # Account for removed transects
        # TODO: this is done in the Matlab code (might be worth investigating)
        if self.percentage_transects_selected is not None:
            self.transect_results_gdf["abundance"] = (
                self.transect_results_gdf["abundance"]
                / self.percentage_transects_selected
            )

        # compute the abundance of males and females
        self.transect_results_male_gdf["abundance"] = (
            self.transect_results_gdf["abundance"] * bc_expanded_df.M_prop.values
        )
        self.transect_results_female_gdf["abundance"] = (
            self.transect_results_gdf["abundance"] * bc_expanded_df.F_prop.values
        )

        # create variable to improve readability
        fraction_adult_stratum_df = self.num_fraction_adult_df.loc[
            self.nasc_df.stratum_num
        ].values.flatten()

        # obtain the abundance for adults
        self.transect_results_male_gdf["abundance_adult"] = (
            self.transect_results_male_gdf["abundance"] * fraction_adult_stratum_df
        )
        self.transect_results_female_gdf["abundance_adult"] = (
            self.transect_results_female_gdf["abundance"] * fraction_adult_stratum_df
        )
        self.transect_results_gdf["abundance_adult"] = (
            self.transect_results_gdf["abundance"] * fraction_adult_stratum_df
        )

    def _set_biomass(self, bc_expanded_df: pd.DataFrame) -> None:
        """
        Calculates biomass variables for each NASC value and then assigns
        them to the appropriate DataFrames.

        Parameters
        ----------
        bc_expanded_df: pd.DataFrame
            An expanded bio parameters dataframe

        Notes
        -----
        This function does not return anything, instead, the created
        variables are directly added to the class variable DataFrames.

        All biomass values are calculated using abundance, instead
        of using the biomass density.
        """

        # calculate the biomass for males, females, and unsexed
        self.transect_results_female_gdf["biomass"] = (
            self.transect_results_female_gdf["abundance"]
            * bc_expanded_df.averaged_weight_F.values
        )
        self.transect_results_male_gdf["biomass"] = (
            self.transect_results_male_gdf["abundance"]
            * bc_expanded_df.averaged_weight_M.values
        )
        biomass_unsexed = (
            self.transect_results_gdf["abundance"]
            - self.transect_results_female_gdf["abundance"]
            - self.transect_results_male_gdf["abundance"]
        ) * bc_expanded_df.averaged_weight.values

        # compute the total biomass for each NASC value
        self.transect_results_gdf["biomass"] = (
            biomass_unsexed
            + self.transect_results_male_gdf["biomass"]
            + self.transect_results_female_gdf["biomass"]
        )

        # create variable to improve readability
        fraction_adult_stratum_df = self.weight_fraction_adult_df.loc[
            self.nasc_df.stratum_num
        ].values.flatten()

        # obtain the biomass for adults
        self.transect_results_female_gdf["biomass_adult"] = (
            self.transect_results_female_gdf["biomass"] * fraction_adult_stratum_df
        )
        self.transect_results_male_gdf["biomass_adult"] = (
            self.transect_results_male_gdf["biomass"] * fraction_adult_stratum_df
        )
        self.transect_results_gdf["biomass_adult"] = (
            self.transect_results_gdf["biomass"] * fraction_adult_stratum_df
        )

    @staticmethod
    def _compute_biomass_all_ages(
        weight_fraction_all_ages_df: pd.DataFrame,
        results_gdf: gpd.GeoDataFrame,
    ) -> None:
        """
        Compute the biomass for each age bin based off the provided input.
        Additionally, add computed variables to the input ``results_gdf``.

        Parameters
        ----------
        weight_fraction_all_ages_df: pd.DataFrame
            A DataFrame containing the weight fraction at each age
            bin for all strata (corresponds to ``biomass_column``)
        results_gdf: gpd.GeoDataFrame
            A GeoDataFrame containing the columns ``biomass_adult, stratum_num``
            and where the biomass at each age bin should be stored
        """

        # obtain stratum column from input gdf
        stratum_vals = results_gdf["stratum_num"]

        for bin_str in weight_fraction_all_ages_df:
            # expand the weight fraction for all ages
            expanded_age_bin = (
                weight_fraction_all_ages_df[bin_str].loc[stratum_vals].values
            )

            results_gdf["biomass_" + bin_str] = (
                expanded_age_bin * results_gdf["biomass_adult"]
            )

    def set_adult_NASC(self) -> None:
        """
        Computes NASC values corresponding to the adult animal population
        and assigns them to `self.transect_results_gdf``
        """

        # get the normalized length-age distribution
        len_age_dist_all_norm = (
            self.bin_ds.len_age_dist_all
            / self.bin_ds.len_age_dist_all.sum(dim=["len_bin", "age_bin"])
        )

        # create adult NASC proportion coefficient
        nasc_fraction_adult_df = pd.DataFrame(
            columns=["val"], index=len_age_dist_all_norm.stratum_num, dtype=np.float64
        )

        # for stratum in self.all_strata:
        for stratum in len_age_dist_all_norm.stratum_num.values:
            sig_bs_aged_ave = np.sum(
                self.survey.params["sig_b_coef"]
                * np.matmul(
                    (self.survey.params["bio_hake_len_bin"] ** 2),
                    len_age_dist_all_norm.sel(stratum_num=stratum).values,
                )
            )

            temp = self.survey.params["sig_b_coef"] * np.matmul(
                (self.survey.params["bio_hake_len_bin"] ** 2),
                len_age_dist_all_norm.sel(stratum_num=stratum).isel(age_bin=0).values,
            )

            age1_nasc_proportion = temp / sig_bs_aged_ave

            nasc_fraction_adult_df.loc[stratum] = abs(1.0 - age1_nasc_proportion)

        # Identify and populate strata missing in nasc_fraction_adult_df
        # (and therefore in self.bin_ds.len_age_dist_all)
        # based on surrounding strata as pre-allocated in stratum_choices
        # TODO: This scheme actually reflects the simpler implementation from v0.1.0-alpha.
        #   Revisit the Matlab code and related missing-stratum handling in this module
        #   https://github.com/uw-echospace/EchoPro_matlab/blob/26b939c9c7c0ccbf8828d402e249e1fdf6ed2b5a/general/load_files_parameters/get_historical_strata_data.m#L731-L795 # noqa
        missing_stratum_num = (
            set(self.nasc_df.stratum_num.values) - set(len_age_dist_all_norm.stratum_num.values)
        )
        for stratum in missing_stratum_num:
            nasc_fraction_adult_df.loc[stratum] = (
                nasc_fraction_adult_df.loc[self.stratum_choices[stratum]]['val'].mean()
            )
        # Storing nasc_fraction_adult_df is not required but will be very useful
        # in the interim for ongoing development and testing
        self.nasc_fraction_adult_df = nasc_fraction_adult_df

        fraction_adult_stratum_df = nasc_fraction_adult_df.loc[
            self.nasc_df.stratum_num
        ].values.flatten()

        # Calculate NASC for adults and assign values to results gdf
        self.transect_results_gdf["NASC_adult"] = self.nasc_df["NASC"] * fraction_adult_stratum_df

    def _construct_results_gdf(self) -> None:
        """
        Constructs self.transect_results_gdf, which contains the
        variables such as biomass density, biomass, abundance, and
        abundance density (numerical density).
        """

        # initialize GeoDataFrames that will hold final results, using nasc_df variables
        temp_df = self.nasc_df[
            ["latitude", "longitude", "stratum_num", "transect_spacing"]
        ].copy(deep=True)
        self.transect_results_gdf = gpd.GeoDataFrame(
            temp_df, geometry=gpd.points_from_xy(temp_df.longitude, temp_df.latitude)
        )
        self.transect_results_male_gdf = self.transect_results_gdf.copy(deep=True)
        self.transect_results_female_gdf = self.transect_results_gdf.copy(deep=True)

        # calculate proportion coefficient for mixed species
        # TODO: note we use all strata_df data every time to match Matlab output
        wgt_vals = self.strata_df.reset_index().set_index("haul_num")["fraction_hake"]
        wgt_vals_ind = wgt_vals.index
        self.mix_sa_ratio = self.nasc_df.apply(
            lambda x: wgt_vals[x.haul_num] if x.haul_num in wgt_vals_ind else 0.0,
            axis=1,
        )
        self.mix_sa_ratio.name = "hake_mix_coefficient"

        # expand the bio parameters dataframe so that it corresponds to nasc_df
        bc_expanded_df = self.bio_param_df.loc[self.nasc_df.stratum_num]

        # calculate and assign numerical density values
        self._set_numerical_density(bc_expanded_df)

        # calculate and assign biomass density values
        self._set_biomass_density(bc_expanded_df)

        # calculate the area corresponding to the NASC value
        self.transect_results_gdf["interval"] = self._get_interval(self.nasc_df)
        self.transect_results_gdf["interval_area_nmi2"] = (
            self.transect_results_gdf["interval"] * self.nasc_df["transect_spacing"]
        )

        # calculate and assign abundance values
        self._set_abundance(bc_expanded_df)

        # calculate and assign biomass values
        self._set_biomass(bc_expanded_df)

        # calculate and add male biomass for all ages to Transect results
        self._compute_biomass_all_ages(
            self.weight_fraction_all_ages_male_df,
            self.transect_results_male_gdf,
        )

        # calculate and add female biomass for all ages to Transect results
        self._compute_biomass_all_ages(
            self.weight_fraction_all_ages_female_df,
            self.transect_results_female_gdf,
        )

        # calculate and add biomass for all ages to Transect results
        self._compute_biomass_all_ages(
            self.weight_fraction_all_ages_df,
            self.transect_results_gdf,
        )

    def get_transect_results_gdf(
        self, selected_transects: Optional[List] = None
    ) -> None:
        """
        Orchestrates the construction of ``self.transect_results_gdf``,
        ``self.transect_results_male_gdf``, and ``self.transect_results_female_gdf``,
        which are GeoDataFrames that contain variables over the transect
        points (e.g. abundance, biomass).

        Parameters
        ----------
        selected_transects : list or None
            The subset of transects used in the calculations
        """

        # store the unique strata values, so they can be used later
        self.all_strata = (
            self.survey.strata_df.index.get_level_values(1).unique().values
        )

        # remove strata index 0 (always done in Matlab version)
        self.all_strata = np.delete(
            self.all_strata, np.argwhere(self.all_strata == 0)[0, 0]
        )

        self.set_class_variables(selected_transects)

        # get the backscattering cross-section for each stratum
        self._get_strata_sig_b()

        # add stratum_num column to length and specimen df and set it as the index
        self._add_stratum_column()

        # identify missing strata, assign strata to missing stratum, fill missing data
        self.set_strata_for_missing_strata()
        self.set_stratum_choice()
        self._fill_missing_strata_sig_b()

        self._get_biomass_parameters()

        self._get_weight_num_fraction_adult()

        self._get_weight_fraction_all_ages()

        self._construct_results_gdf()
