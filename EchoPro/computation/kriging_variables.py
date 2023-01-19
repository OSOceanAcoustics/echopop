from typing import List, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr


class ComputeKrigingVariables:
    # TODO: see if we want to change the name of this class
    # TODO: this class may not correctly account for bootstrapping!
    """
    A class that computes key variables corresponding to each
    Kriging mesh point, using output produced by Kriging.

    Parameters
    ----------
    krig : ComputeKrigingVariables
        An initialized ComputeKrigingVariables object. Note that any change to
        self.krig will also change this object.
    """

    def __init__(self, krig=None):

        self.krig = krig

    @staticmethod
    def _get_bin_ind(input_data: np.ndarray, bin_edges: np.ndarray) -> List[np.ndarray]:
        """
        This function manually finds those indices in ``input_data``
        that are in each bin.

        Parameters
        ----------
        input_data: np.ndarray
            The data to bin
        bin_edges: np.ndarray
            An array that specifies the bin edges.

        Returns
        -------
        hist_ind: list
            The index values of ``input_data`` corresponding to the histogram

        Notes
        -----
        The construction of the bin counts differs from the method
        produced in `biomass_desnity.py`. This is because the Matlab
        version of EchoPro is inconsistent in how it is binning.
        """

        # initialize list that will hold the indices
        hist_ind = []

        for i in range(len(bin_edges) - 1):

            # get values greater or equal than lower bound
            g_lb = bin_edges[i] <= input_data

            # get values less than the upper bound
            le_ub = input_data < bin_edges[i + 1]

            # fill bin
            hist_ind.append(np.argwhere(g_lb & le_ub).flatten())

        # fill in the last bin
        hist_ind.append(np.argwhere(input_data >= bin_edges[-1]).flatten())

        return hist_ind

    @staticmethod
    def _get_bin_ind_age(
        input_data: np.ndarray, age_bins: np.ndarray
    ) -> List[np.ndarray]:
        """
        This function manually finds the indices of ``input_data``
        for each age bin. An age bin here is described as all values
        equal to the provided age


        Parameters
        ----------
        input_data: np.ndarray
            The data to bin
        age_bins: np.ndarray
            An array that specifies the age for each bin

        Returns
        -------
        hist_ind: list
            The index values of ``input_data`` corresponding to the histogram

        Notes
        -----
        The construction of the bin counts differs from the method
        produced in `biomass_desnity.py`. This is because the Matlab
        version of EchoPro is inconsistent in how it is binning ages.
        """

        # initialize list that will hold indices
        hist_ind = []

        for age in age_bins:

            # determine indices with a value equal to age
            hist_ind.append(np.argwhere(input_data == age).flatten())

        return hist_ind

    @staticmethod
    def _initialize_ds(
        stratum_ind: np.ndarray, len_bin: np.ndarray, age_bin: np.ndarray
    ) -> xr.Dataset:
        """
        Initializes the parameter Dataset.

        Parameters
        ----------
        stratum_ind: np.ndarray
            A one dimensional array specifying all strata values
        len_bin: np.ndarray
            A one dimensional array specifying the length bin values
        age_bin: np.ndarray
            A one dimensional array specifying all age bin values

        Returns
        -------
        ds: xr.Dataset
            The parameter Dataset with all variables initialized
        """

        # initialize Dataset that will hold length age distributions
        ds = xr.Dataset(
            data_vars={
                "num_M": ("stratum", np.zeros(len(stratum_ind))),
                "num_F": ("stratum", np.zeros(len(stratum_ind))),
                "total_weight": ("stratum", np.zeros(len(stratum_ind))),
                "len_age_weight_prop_all": ("stratum", np.zeros(len(stratum_ind))),
                "len_age_weight_prop_M": ("stratum", np.zeros(len(stratum_ind))),
                "len_age_weight_prop_F": ("stratum", np.zeros(len(stratum_ind))),
                "aged_proportion": ("stratum", np.zeros(len(stratum_ind))),
                "unaged_proportion": ("stratum", np.zeros(len(stratum_ind))),
                "unaged_M_wgt_proportion": ("stratum", np.zeros(len(stratum_ind))),
                "unaged_F_wgt_proportion": ("stratum", np.zeros(len(stratum_ind))),
                "weight_len_all_normalized": (
                    ["stratum", "len_bin"],
                    np.zeros((len(stratum_ind), len(len_bin))),
                ),
                "len_age_dist_all": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_dist_M": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_dist_F": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_weight_dist_all": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_weight_dist_M": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_weight_dist_F": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_weight_dist_all_normalized": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_weight_dist_M_normalized": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
                "len_age_weight_dist_F_normalized": (
                    ["stratum", "len_bin", "age_bin"],
                    np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
                ),
            },
            coords={
                "stratum": ("stratum", stratum_ind),
                "len_bin": ("len_bin", len_bin),
                "age_bin": ("age_bin", age_bin),
            },
        )

        return ds

    def _get_len_wgt_distributions(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Computes the length-weight distributions for each gender
        and all animals using the specimen data

        Returns
        -------
        len_wgt_M: np.ndarray
            The length-weight distribution for males
        len_wgt_F: np.ndarray
            The length-weight distribution for females
        len_wgt_all: np.ndarray
            The length-weight distribution for all genders

        Notes
        -----
        These calculations use all specimen data, rather than the
        specimen data with NA values dropped. This is necessary to
        match the Matlab output.
        """

        len_wgt_M = self.krig.survey.bio_calc._generate_length_val_conversion(
            len_name="length",
            val_name="weight",
            df=self.krig.survey.bio_calc.specimen_df[
                self.krig.survey.bio_calc.specimen_df["sex"] == 1
            ],
        )

        len_wgt_F = self.krig.survey.bio_calc._generate_length_val_conversion(
            len_name="length",
            val_name="weight",
            df=self.krig.survey.bio_calc.specimen_df[
                self.krig.survey.bio_calc.specimen_df["sex"] == 2
            ],
        )

        len_wgt_all = self.krig.survey.bio_calc._generate_length_val_conversion(
            len_name="length",
            val_name="weight",
            df=self.krig.survey.bio_calc.specimen_df,
        )

        return len_wgt_M, len_wgt_F, len_wgt_all

    def _set_age_distribution_data(
        self, stratum: int, ds: xr.Dataset, df_M: pd.DataFrame, df_F: pd.DataFrame
    ) -> None:
        """
        Computes distributions for each age using the input DataFrames.
        Additionally, assigns these computed quantities to the Dataset ``ds``.

        Parameters
        ----------
        stratum: int
            The stratum corresponding to the input data (used when
            assigning values to ``ds``)
        ds: xr.Dataset
            The Dataset where computed quantities should be assigned to
        df_M: pd.DataFrame
            A DataFrame specifying length, age, and weight measurements
            of a male animal within a stratum
        df_F: pd.DataFrame
            A DataFrame specifying length, age, and weight measurements
            of a female animal within a stratum
        """

        # account for the case when df is a Series
        if isinstance(df_M, pd.Series):
            # get numpy arrays of length, age, and weight
            input_arr_len_M = np.array([df_M.length])
            input_arr_age_M = np.array([df_M.age])
            input_arr_wgt_M = np.array([df_M.weight])

            input_arr_len_F = np.array([df_F.length])
            input_arr_age_F = np.array([df_F.age])
            input_arr_wgt_F = np.array([df_F.weight])

        else:
            # get numpy arrays of length, age, and weight
            input_arr_len_M = df_M.length.values
            input_arr_age_M = df_M.age.values
            input_arr_wgt_M = df_M.weight.values

            input_arr_len_F = df_F.length.values
            input_arr_age_F = df_F.age.values
            input_arr_wgt_F = df_F.weight.values

        # bin the ages
        # TODO: binning is occurring differently than in biomass_denisty.py! It may be
        #  better and more consistent to use the function self.krig.survey.bio_calc._get_bin_ind
        age_bins_ind_M = self._get_bin_ind_age(
            input_arr_age_M, self.krig.survey.bio_calc.bio_hake_age_bin
        )
        age_bins_ind_F = self._get_bin_ind_age(
            input_arr_age_F, self.krig.survey.bio_calc.bio_hake_age_bin
        )

        # round input lengths
        # TODO: this is necessary to match the Matlab output and may not be necessary!
        input_arr_len_M = np.round(input_arr_len_M)
        input_arr_len_F = np.round(input_arr_len_F)

        # compute distributions for each age bin
        for age_bin in range(len(age_bins_ind_M)):

            # bin those lengths that correspond to the lengths in the given age bin
            # TODO: binning is occurring differently than in biomass_denisty.py!
            len_bin_ind_M = self._get_bin_ind(
                input_arr_len_M[age_bins_ind_M[age_bin]],
                self.krig.survey.bio_calc.bio_hake_len_bin,
            )
            len_bin_ind_F = self._get_bin_ind(
                input_arr_len_F[age_bins_ind_F[age_bin]],
                self.krig.survey.bio_calc.bio_hake_len_bin,
            )

            # get the distribution of weight for a particular age bin
            ds.sel(stratum=stratum).len_age_weight_dist_M[:, age_bin] = np.array(
                [
                    np.sum(input_arr_wgt_M[age_bins_ind_M[age_bin]][i])
                    for i in len_bin_ind_M
                ]
            )

            ds.sel(stratum=stratum).len_age_weight_dist_F[:, age_bin] = np.array(
                [
                    np.sum(input_arr_wgt_F[age_bins_ind_F[age_bin]][i])
                    for i in len_bin_ind_F
                ]
            )

            ds.sel(stratum=stratum).len_age_weight_dist_all[:, age_bin] = (
                ds.sel(stratum=stratum).len_age_weight_dist_M[:, age_bin]
                + ds.sel(stratum=stratum).len_age_weight_dist_F[:, age_bin]
            )

            # get the distribution of lengths for a particular age bin
            ds.sel(stratum=stratum).len_age_dist_M[:, age_bin] = np.array(
                [len(i) for i in len_bin_ind_M]
            )

            ds.sel(stratum=stratum).len_age_dist_F[:, age_bin] = np.array(
                [len(i) for i in len_bin_ind_F]
            )

            ds.sel(stratum=stratum).len_age_dist_all[:, age_bin] = (
                ds.sel(stratum=stratum).len_age_dist_M[:, age_bin]
                + ds.sel(stratum=stratum).len_age_dist_F[:, age_bin]
            )

        # obtain normalized distributions
        ds.sel(stratum=stratum).len_age_weight_dist_all_normalized[
            :, :
        ] = ds.len_age_weight_dist_all.sel(stratum=stratum) / np.nansum(
            ds.len_age_weight_dist_all.sel(stratum=stratum)
        )

        ds.sel(stratum=stratum).len_age_weight_dist_M_normalized[
            :, :
        ] = ds.len_age_weight_dist_M.sel(stratum=stratum) / np.nansum(
            ds.len_age_weight_dist_M.sel(stratum=stratum)
        )

        ds.sel(stratum=stratum).len_age_weight_dist_F_normalized[
            :, :
        ] = ds.len_age_weight_dist_F.sel(stratum=stratum) / np.nansum(
            ds.len_age_weight_dist_F.sel(stratum=stratum)
        )

    def _set_total_weight(
        self, haul_nums: np.ndarray, stratum: int, ds: xr.Dataset
    ) -> float:
        """
        Computes the total weight within the stratum using data from
        both stations. Additionally, assigns the variable ``total_weight``
        to ``ds``.

        Parameters
        ----------
        haul_nums: np.ndarray
            A one dimensional array of haul numbers within the stratum
            under consideration
        stratum: int
            The stratum corresponding to the input data (used when
            assigning values to ``ds``)
        ds: xr.Dataset
            The Dataset where computed quantities should be assigned to

        Returns
        -------
        wgt_station_1: float
            The total weight at station 1 for the stratum
        """

        # calculate the total weight within the stratum for station 1
        wgt_station_1 = np.array(
            [
                self.krig.survey.catch_df.loc[j][
                    "haul_weight"
                ].sum()  # TODO: make catch_df a bio_calc variable?
                for j in haul_nums
                if (j in self.krig.survey.catch_df.index)
                and (j in self.krig.survey.length_df.index)
            ]
        ).sum()

        # calculate the total weight within the stratum for station 2
        wgt_station_2 = self.krig.survey.bio_calc.specimen_df.loc[stratum][
            "weight"
        ].sum()

        # the total weight within the stratum
        ds.total_weight.loc[stratum] = wgt_station_1 + wgt_station_2

        return wgt_station_1

    @staticmethod
    def _get_length_based_wgt_interp(
        haul_nums: np.ndarray,
        df: pd.DataFrame,
        len_bin: np.ndarray,
        len_wgt: np.ndarray,
    ) -> float:
        """
        Obtains the weight of animals within a stratum using ``length_df``
        data and interpolation.

        Parameters
        ----------
        haul_nums: np.ndarray
            A one dimensional array of haul numbers within the stratum
            under consideration
        df: pd.DataFrame
            A DataFrame containing data from ``length_df``
        len_bin: np.ndarray
            Length bin values
        len_wgt: np.ndarray
            The length-weight distribution corresponding to ``df``

        Returns
        -------
        final_wgt: float
            The total weight of animals within the stratum
        """

        # initialize weight value
        final_wgt = 0.0

        # obtain the hauls that are within df
        len_haul = [round(df.loc[j]["length"]) for j in haul_nums if j in df.index]

        if len_haul:

            # get the number of lengths for each haul
            len_haul_M_counts = np.concatenate(
                [df.loc[j]["length_count"].values for j in haul_nums if j in df.index]
            )

            # calculate the weight of animals within the stratum
            len_haul_M = np.concatenate(len_haul)
            final_wgt = (
                np.interp(len_haul_M, len_bin, len_wgt) * len_haul_M_counts
            ).sum()

        return final_wgt

    def _get_length_df_based_wgt(
        self,
        stratum: int,
        ds: xr.Dataset,
        haul_nums: np.ndarray,
        length_df_M: pd.DataFrame,
        length_df_F: pd.DataFrame,
        len_bin: np.ndarray,
        len_wgt_M: np.ndarray,
        len_wgt_F: np.ndarray,
        len_wgt_all: np.ndarray,
    ) -> Tuple[float, float]:
        """
        Computes the weight of a particular stratum based off of the
        length data. Additionally, assigns values to the ``ds`` variable
        ``weight_len_all_normalized``.

        Parameters
        ----------
        stratum: int
            The stratum corresponding to the input data (used when
            assigning values to ``ds``)
        ds: xr.Dataset
            The Dataset where computed quantities should be assigned to
        haul_nums: np.ndarray
            A one dimensional array of haul numbers within the stratum
            under consideration
        length_df_M: pd.DataFrame
            Male data obtained from ``length_df``
        length_df_F: pd.DataFrame
            Female data obtained from ``length_df``
        len_bin: np.ndarray
            Length bin values
        len_wgt_M: np.ndarray
            The length-weight distribution for males
        len_wgt_F: np.ndarray
            The length-weight distribution for females
        len_wgt_all: np.ndarray
            The length-weight distribution for all genders

        Returns
        -------
        male_wgt: float
            The total weight of males within the stratum
        female_wgt: float
            The total weight of females within the stratum
        """

        # calculate the weight of male and female animals within the stratum
        male_wgt = self._get_length_based_wgt_interp(
            haul_nums, length_df_M, len_bin, len_wgt_M
        )
        female_wgt = self._get_length_based_wgt_interp(
            haul_nums, length_df_F, len_bin, len_wgt_F
        )

        # obtain all haul numbers within the stratum and in length_df
        hauls_in_all = [j for j in haul_nums if j in self.krig.survey.length_df.index]

        if hauls_in_all:

            # get normalized length distribution of length_df data
            len_dist_station1_normalized = (
                self.krig.survey.bio_calc._get_distribution_lengths_station_1(
                    self.krig.survey.length_df.loc[hauls_in_all]
                )
            )

            # obtain the weight per unit length distribution
            weight_len_all = len_wgt_all * len_dist_station1_normalized

            # normalized weight per unit length distribution
            ds.sel(stratum=stratum).weight_len_all_normalized[:] = (
                weight_len_all / weight_len_all.sum()
            )

        return female_wgt, male_wgt

    def _set_proportion_parameters(
        self,
        stratum: int,
        ds: xr.Dataset,
        male_wgt: float,
        female_wgt: float,
        wgt_station_1: float,
    ) -> None:
        """
        Calculates and assigns proportion parameters to ``ds`` for a stratum.

        Parameters
        ----------
        stratum: int
            The stratum corresponding to the input data (used when
            assigning values to ``ds``)
        ds: xr.Dataset
            The Dataset where computed quantities should be assigned to
        male_wgt: float
            The total weight of males within the stratum
        female_wgt: float
            The total weight of females within the stratum
        wgt_station_1: float
            The total weight at station 1 for the stratum
        """

        # calculate and assign the len_age_weight proportions
        ds.len_age_weight_prop_all.loc[stratum] = np.nansum(
            ds.len_age_weight_dist_all.sel(stratum=stratum)
        ) / ds.total_weight.sel(stratum=stratum)

        ds.len_age_weight_prop_M.loc[stratum] = np.nansum(
            ds.len_age_weight_dist_M.sel(stratum=stratum)
        ) / ds.total_weight.sel(stratum=stratum)

        ds.len_age_weight_prop_F.loc[stratum] = np.nansum(
            ds.len_age_weight_dist_F.sel(stratum=stratum)
        ) / ds.total_weight.sel(stratum=stratum)

        # calculate and assign aged proportions
        ds.aged_proportion.loc[stratum] = (
            ds.len_age_weight_prop_M.loc[stratum]
            + ds.len_age_weight_prop_F.loc[stratum]
        )

        # calculate and assign the unaged proportions
        ds.unaged_proportion.loc[stratum] = 1.0 - ds.aged_proportion.loc[stratum]

        # obtain the normalized weight of station 1 for males and females
        if (male_wgt != 0.0) and (female_wgt != 0.0):
            nM_wgt1 = wgt_station_1 * male_wgt / (male_wgt + female_wgt)
            nF_wgt1 = wgt_station_1 * female_wgt / (male_wgt + female_wgt)
        else:
            nM_wgt1 = 0.0
            nF_wgt1 = 0.0

        # obtain length and gender based weight proportion
        Len_M_wgt_proportion = nM_wgt1 / ds.total_weight.sel(stratum=stratum).values
        Len_F_wgt_proportion = nF_wgt1 / ds.total_weight.sel(stratum=stratum).values

        # obtain the proportion of males and females
        if (Len_M_wgt_proportion == 0.0) and (Len_F_wgt_proportion == 0.0):
            M_proportion = 0.5
            F_proportion = 0.5
        else:
            M_proportion = Len_M_wgt_proportion / (
                Len_M_wgt_proportion + Len_F_wgt_proportion
            )
            F_proportion = Len_F_wgt_proportion / (
                Len_M_wgt_proportion + Len_F_wgt_proportion
            )

        # calculate and assign the unaged weight proportion of males and females
        ds.unaged_M_wgt_proportion.loc[stratum] = (
            ds.unaged_proportion.loc[stratum].values * M_proportion
        )
        ds.unaged_F_wgt_proportion.loc[stratum] = (
            ds.unaged_proportion.loc[stratum].values * F_proportion
        )

    def generate_parameter_ds(self) -> xr.Dataset:
        """
        Creates a Dataset containing parameters that are
        necessary to compute Kriging result variables at
        each mesh point.

        Returns
        -------
        ds: xr.Dataset
            A Dataset containing useful parameters
        """

        # obtain specimen DataFrames without NA values
        # TODO: This is necessary to match the Matlab output
        #  in theory this should be done when we load the df,
        #  however, this changes the results slightly.
        spec_drop = self.krig.survey.bio_calc.specimen_df.dropna(how="any")
        spec_drop_M = spec_drop[spec_drop["sex"] == 1]
        spec_drop_F = spec_drop[spec_drop["sex"] == 2]

        # obtain gender based length DataFrames
        # TODO: does not use bio_calc.length_df value (it would give incorrect
        #  answers as that df drops values)
        length_df_M = self.krig.survey.length_df[self.krig.survey.length_df["sex"] == 1]
        length_df_F = self.krig.survey.length_df[self.krig.survey.length_df["sex"] == 2]

        # get all unique stratum values
        stratum_ind = spec_drop.index.unique()

        # get age and length bins (created to reduce clutter)
        len_bin = self.krig.survey.params["bio_hake_len_bin"]
        age_bin = self.krig.survey.params["bio_hake_age_bin"]

        # initialize Dataset that will hold all parameters
        ds = self._initialize_ds(stratum_ind, len_bin, age_bin)

        # obtain a mapping of hauls to strata
        haul_vs_stratum = self.krig.survey.bio_calc.strata_df.reset_index()[
            ["haul_num", "stratum_num"]
        ]

        # get length-weight distributions (includes all ages in quantity)
        len_wgt_M, len_wgt_F, len_wgt_all = self._get_len_wgt_distributions()

        for i in stratum_ind:

            # obtain haul numbers that are in the stratum i
            haul_nums = haul_vs_stratum[haul_vs_stratum["stratum_num"] == i][
                "haul_num"
            ].values

            # calculate and set the number of animals in the stratum based on gender
            ds.num_M.loc[i] = len(spec_drop_M.loc[i])
            ds.num_F.loc[i] = len(spec_drop_F.loc[i])

            # set age distribution related data
            self._set_age_distribution_data(
                stratum=i,
                ds=ds,
                df_M=spec_drop_M.loc[i],
                df_F=spec_drop_F.loc[i],
            )

            # assign the total weight in both station to ds and obtain the weight in station 1
            wgt_station_1 = self._set_total_weight(haul_nums, i, ds)

            # get length_df based weight for males and females
            female_wgt, male_wgt = self._get_length_df_based_wgt(
                i,
                ds,
                haul_nums,
                length_df_M,
                length_df_F,
                len_bin,
                len_wgt_M,
                len_wgt_F,
                len_wgt_all,
            )

            # calculate and assign proportion parameters to ds
            self._set_proportion_parameters(i, ds, male_wgt, female_wgt, wgt_station_1)

        return ds

    def _set_gender_biomass(self, ds: xr.Dataset) -> None:
        """
        Calculates the biomass for males and females at
        each mesh point. Additionally, adds the corresponding
        variables to ``kriging_results_gdf``.

        Parameters
        ----------
        ds: xr.Dataset
            A Dataset containing all parameters necessary to compute
            desired variables
        """

        # obtain stratum column from Kriging results
        stratum_vals = self.krig.survey.bio_calc.kriging_results_gdf[
            "stratum_num"
        ].values

        # create variables to improve readability
        biomass_density_adult_mean = self.krig.survey.bio_calc.kriging_results_gdf[
            "biomass_density_adult_mean"
        ]
        cell_area_nmi2 = self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]

        # calculate the aged biomass for males and females
        dist_weight_M_sum = (
            (ds.len_age_weight_dist_M_normalized * ds.len_age_weight_prop_M)
            .sum(dim=["len_bin", "age_bin"])
            .sel(stratum=stratum_vals)
            .values
        )
        dist_weight_F_sum = (
            (ds.len_age_weight_dist_F_normalized * ds.len_age_weight_prop_F)
            .sum(dim=["len_bin", "age_bin"])
            .sel(stratum=stratum_vals)
            .values
        )

        biomass_aged_M = biomass_density_adult_mean * dist_weight_M_sum * cell_area_nmi2
        biomass_aged_F = biomass_density_adult_mean * dist_weight_F_sum * cell_area_nmi2

        # calculate the unaged biomass for males and females
        unaged_M_wgt_prop_expan = ds.unaged_M_wgt_proportion.sel(
            stratum=stratum_vals
        ).values
        unaged_F_wgt_prop_expan = ds.unaged_F_wgt_proportion.sel(
            stratum=stratum_vals
        ).values

        biomass_unaged_M = (
            biomass_density_adult_mean * unaged_M_wgt_prop_expan * cell_area_nmi2
        )
        biomass_unaged_F = (
            biomass_density_adult_mean * unaged_F_wgt_prop_expan * cell_area_nmi2
        )

        # calculate and assign the total biomass of males and females
        self.krig.survey.bio_calc.kriging_results_male_gdf["biomass_adult"] = (
            biomass_aged_M + biomass_unaged_M
        )
        self.krig.survey.bio_calc.kriging_results_female_gdf["biomass_adult"] = (
            biomass_aged_F + biomass_unaged_F
        )

    def _set_abundance(self) -> None:
        """
        Calculates the abundance for males, females, and all sexes at
        each mesh point. Additionally, adds the corresponding
        variables to ``kriging_results_gdf``.
        """

        # expand the bio parameters dataframe so that it corresponds to mesh points
        averaged_weight_expanded = (
            self.krig.survey.bio_calc.bio_param_df.averaged_weight.loc[
                self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
            ]
        )

        # calculate and add abundance to Kriging results
        self.krig.survey.bio_calc.kriging_results_gdf["abundance_adult"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_adult"]
            / averaged_weight_expanded.values
        )

        self.krig.survey.bio_calc.kriging_results_male_gdf["abundance_adult"] = (
            self.krig.survey.bio_calc.kriging_results_male_gdf["biomass_adult"]
            / averaged_weight_expanded.values
        )

        self.krig.survey.bio_calc.kriging_results_female_gdf["abundance_adult"] = (
            self.krig.survey.bio_calc.kriging_results_female_gdf["biomass_adult"]
            / averaged_weight_expanded.values
        )

    def _set_biomass_cell_CV(self):
        """
        Compute the coefficient of Variation (CV) of biomass
        at each grid cell using Kriging output. Additionally,
        assigns the created variable to ``kriging_results_gdf``.
        """

        # create variables to improve readability
        biomass_density_adult = self.krig.survey.bio_calc.transect_results_gdf[
            "biomass_density_adult"
        ].values
        biomass_density_adult_mean = self.krig.survey.bio_calc.kriging_results_gdf[
            "biomass_density_adult_mean"
        ]
        cell_area_nmi2 = self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        biomass_density_adult_var = self.krig.survey.bio_calc.kriging_results_gdf[
            "biomass_density_adult_var"
        ]
        kriging_A0 = self.krig.survey.params["kriging_A0"]

        C0 = np.std(biomass_density_adult, ddof=1) ** 2
        Bn = np.nansum(biomass_density_adult_mean * cell_area_nmi2) * 1e-9
        self.krig.survey.bio_calc.kriging_results_gdf["biomass_adult_cell_CV"] = (
            kriging_A0
            * np.sqrt(biomass_density_adult_var * C0)
            * 1e-9
            / Bn
            * np.sqrt(len(biomass_density_adult_var))
        )

    def _compute_biomass_all_ages(
        self,
        weight_fraction_all_ages_df: pd.DataFrame,
        biomass_column: pd.DataFrame,
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
        biomass_column: pd.DataFrame
            A DataFrame corresponding to a biomass column from the Kriging results
        results_gdf: gpd.GeoDataFrame
            A GeoDataFrame where the biomass at each age bin should be stored
        """

        # obtain stratum column from Kriging results
        stratum_vals = results_gdf["stratum_num"].values

        for bin_str in weight_fraction_all_ages_df:

            # expand the weight fraction for all ages
            expanded_age_bin = (
                weight_fraction_all_ages_df[bin_str].loc[stratum_vals].values
            )

            results_gdf["biomass_" + bin_str] = expanded_age_bin * biomass_column

    def set_variables(self, ds: xr.Dataset) -> None:
        """
        Calculates variables over Kriging mesh points that are useful
        for analysis (e.g. abundance, NASC, CV). Additionally, assigns
        these variables to the appropriate Kriging results GeoDataFrame.

        Parameters
        ----------
        ds: xr.Dataset
            A Dataset containing all parameters necessary to compute
            desired variables
        """

        # calculate and add the male and female biomass to Kriging results
        self._set_gender_biomass(ds)

        # calculate and add abundance variables to Kriging results
        self._set_abundance()

        # add sig_b values to Kriging results
        self.krig.survey.bio_calc.kriging_results_gdf[
            "sig_b"
        ] = self.krig.survey.bio_calc.strata_sig_b.loc[
            self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
        ].values

        # calculate and add NASC to Kriging results
        self.krig.survey.bio_calc.kriging_results_gdf["NASC"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["abundance_adult"]
            * self.krig.survey.bio_calc.kriging_results_gdf["sig_b"]
        )

        # calculate and add biomass CV at each grid cell to Kriging results
        self._set_biomass_cell_CV()

        # calculate and add male biomass for all ages to Kriging results
        self._compute_biomass_all_ages(
            self.krig.survey.bio_calc.weight_fraction_all_ages_male_df,
            self.krig.survey.bio_calc.kriging_results_male_gdf["biomass_adult"],
            self.krig.survey.bio_calc.kriging_results_male_gdf,
        )

        # calculate and add female biomass for all ages to Kriging results
        self._compute_biomass_all_ages(
            self.krig.survey.bio_calc.weight_fraction_all_ages_female_df,
            self.krig.survey.bio_calc.kriging_results_female_gdf["biomass_adult"],
            self.krig.survey.bio_calc.kriging_results_female_gdf,
        )

        # calculate and add biomass for all ages to Kriging results
        self._compute_biomass_all_ages(
            self.krig.survey.bio_calc.weight_fraction_all_ages_df,
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_adult"],
            self.krig.survey.bio_calc.kriging_results_gdf,
        )
