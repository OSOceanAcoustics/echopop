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
    def _get_bin_ind(
        input_data: np.ndarray, centered_bins: np.ndarray
    ) -> list:  # List[np.ndarray]:
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

        # fill the first bin
        hist_ind = []  # np.argwhere(input_data < centered_bins[0]).flatten()]

        for i in range(len(centered_bins) - 1):
            # get values greater than lower bound
            g_lb = centered_bins[i] <= input_data

            # get values less than or equal to the upper bound
            le_ub = input_data < centered_bins[i + 1]

            # fill bin
            hist_ind.append(np.argwhere(g_lb & le_ub).flatten())

        # fill in the last bin
        hist_ind.append(np.argwhere(input_data >= centered_bins[-1]).flatten())

        return hist_ind

    @staticmethod
    def _get_bin_ind_age(
        input_data: np.ndarray, centered_bins: np.ndarray
    ) -> list:  # List[np.ndarray]:
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

        # fill the first bin
        hist_ind = []  # np.argwhere(input_data < centered_bins[0]).flatten()]

        for i in centered_bins:

            cond = input_data == i

            # fill bin
            hist_ind.append(np.argwhere(cond).flatten())

        return hist_ind

    def _set_len_age_distribution_data(self, stratum, ds, df_M, df_F):

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
        # age_bins_ind_M = self.krig.survey.bio_calc._get_bin_ind(
        #     input_arr_age_M, self.krig.survey.bio_calc.bio_hake_age_bin
        # )
        # age_bins_ind_F = self.krig.survey.bio_calc._get_bin_ind(
        #     input_arr_age_F, self.krig.survey.bio_calc.bio_hake_age_bin
        # )

        # TODO: binning is occurring differently than in biomass_denisty.py!
        age_bins_ind_M = self._get_bin_ind_age(
            input_arr_age_M, self.krig.survey.bio_calc.bio_hake_age_bin
        )
        age_bins_ind_F = self._get_bin_ind_age(
            input_arr_age_F, self.krig.survey.bio_calc.bio_hake_age_bin
        )

        input_arr_len_M = np.round(input_arr_len_M)
        input_arr_len_F = np.round(input_arr_len_F)

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

            # get the distribution of weight for males using the length, age, weight data
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

    def _generate_len_age_distributions(self):

        # TODO: This is necessary to match the Matlab output
        #  in theory this should be done when we load the df,
        #  however, this changes the results slightly.
        spec_drop = self.krig.survey.bio_calc.specimen_df.dropna(how="any")

        spec_drop_M = spec_drop[spec_drop["sex"] == 1]
        spec_drop_F = spec_drop[spec_drop["sex"] == 2]

        # TODO: does not use bio_calc.length_df value (it would give incorrect answers)
        length_df_M = self.krig.survey.length_df[self.krig.survey.length_df["sex"] == 1]
        length_df_F = self.krig.survey.length_df[self.krig.survey.length_df["sex"] == 2]

        stratum_ind = spec_drop.index.unique()

        len_bin = self.krig.survey.params["bio_hake_len_bin"]
        age_bin = self.krig.survey.params["bio_hake_age_bin"]

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

        # a mapping of hauls to strata
        haul_vs_stratum = self.krig.survey.bio_calc.strata_df.reset_index()[
            ["haul_num", "stratum_num"]
        ]

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

        for i in stratum_ind:

            haul_nums = haul_vs_stratum[haul_vs_stratum["stratum_num"] == i][
                "haul_num"
            ].values

            ds.num_M.loc[i] = len(spec_drop_M.loc[i])
            ds.num_F.loc[i] = len(spec_drop_F.loc[i])

            self._set_len_age_distribution_data(
                stratum=i,
                ds=ds,
                df_M=spec_drop_M.loc[i],
                df_F=spec_drop_F.loc[i],
            )

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

            wgt_station_2 = self.krig.survey.bio_calc.specimen_df.loc[i]["weight"].sum()

            ds.total_weight.loc[i] = (
                wgt_station_1 + wgt_station_2
            )  # TODO this is off for stratum 5

            ds.len_age_weight_prop_all.loc[i] = np.nansum(
                ds.len_age_weight_dist_all.sel(stratum=i)
            ) / ds.total_weight.sel(stratum=i)

            ds.len_age_weight_prop_M.loc[i] = np.nansum(
                ds.len_age_weight_dist_M.sel(stratum=i)
            ) / ds.total_weight.sel(stratum=i)

            ds.len_age_weight_prop_F.loc[i] = np.nansum(
                ds.len_age_weight_dist_F.sel(stratum=i)
            ) / ds.total_weight.sel(stratum=i)

            ds.sel(stratum=i).len_age_weight_dist_all_normalized[
                :, :
            ] = ds.len_age_weight_dist_all.sel(stratum=i) / np.nansum(
                ds.len_age_weight_dist_all.sel(stratum=i)
            )

            ds.sel(stratum=i).len_age_weight_dist_M_normalized[
                :, :
            ] = ds.len_age_weight_dist_M.sel(stratum=i) / np.nansum(
                ds.len_age_weight_dist_M.sel(stratum=i)
            )

            ds.sel(stratum=i).len_age_weight_dist_F_normalized[
                :, :
            ] = ds.len_age_weight_dist_F.sel(stratum=i) / np.nansum(
                ds.len_age_weight_dist_F.sel(stratum=i)
            )

            ds.aged_proportion.loc[i] = (
                ds.len_age_weight_prop_M.loc[i] + ds.len_age_weight_prop_F.loc[i]
            )

            ds.unaged_proportion.loc[i] = 1.0 - ds.aged_proportion.loc[i]

            len_haul_M = [
                round(length_df_M.loc[j]["length"])
                for j in haul_nums
                if j in length_df_M.index
            ]

            if len_haul_M:
                len_haul_M_counts = np.concatenate(
                    [
                        length_df_M.loc[j]["length_count"].values
                        for j in haul_nums
                        if j in length_df_M.index
                    ]
                )
                len_haul_M = np.concatenate(len_haul_M)

                male_wgt = (
                    np.interp(len_haul_M, len_bin, len_wgt_M) * len_haul_M_counts
                ).sum()
            else:
                male_wgt = 0.0

            len_haul_F = [
                round(length_df_F.loc[j]["length"])
                for j in haul_nums
                if j in length_df_M.index
            ]
            if len_haul_F:
                len_haul_F_counts = np.concatenate(
                    [
                        length_df_F.loc[j]["length_count"].values
                        for j in haul_nums
                        if j in length_df_F.index
                    ]
                )
                len_haul_F = np.concatenate(len_haul_F)

                female_wgt = (
                    np.interp(len_haul_F, len_bin, len_wgt_F) * len_haul_F_counts
                ).sum()
            else:
                female_wgt = 0.0

            hauls_in_all = [
                j for j in haul_nums if j in self.krig.survey.length_df.index
            ]

            if hauls_in_all:

                len_dist_station1_normalized = (
                    self.krig.survey.bio_calc._get_distribution_lengths_station_1(
                        self.krig.survey.length_df.loc[hauls_in_all]
                    )
                )

                weight_len_all = len_wgt_all * len_dist_station1_normalized

                # normalized weight per unit length distribution
                ds.sel(stratum=i).weight_len_all_normalized[:] = (
                    weight_len_all / weight_len_all.sum()
                )

            else:

                ds.sel(stratum=i).weight_len_all_normalized[:] = 0.0

            # print(" ")

            if (male_wgt != 0.0) and (female_wgt != 0.0):
                nM_wgt1 = wgt_station_1 * male_wgt / (male_wgt + female_wgt)
                nF_wgt1 = wgt_station_1 * female_wgt / (male_wgt + female_wgt)
            else:
                nM_wgt1 = 0.0
                nF_wgt1 = 0.0

            Len_M_wgt_proportion = nM_wgt1 / ds.total_weight.sel(stratum=i).values
            Len_F_wgt_proportion = nF_wgt1 / ds.total_weight.sel(stratum=i).values

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

            ds.unaged_M_wgt_proportion.loc[i] = (
                ds.unaged_proportion.loc[i].values * M_proportion
            )
            ds.unaged_F_wgt_proportion.loc[i] = (
                ds.unaged_proportion.loc[i].values * F_proportion
            )

        return ds

    def set_variables(self, ds):

        ds["dist_weight_M_sum"] = (
            ds.len_age_weight_dist_M_normalized * ds.len_age_weight_prop_M
        ).sum(dim=["len_bin", "age_bin"])
        ds["dist_weight_F_sum"] = (
            ds.len_age_weight_dist_F_normalized * ds.len_age_weight_prop_F
        ).sum(dim=["len_bin", "age_bin"])

        # TODO: create variable for self.krig.survey.bio_calc.
        #  kriging_results_gdf["stratum_num"].values

        # aged_prop_mesh = ds.aged_proportion.sel(
        #     stratum=self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
        # ).values
        # # len_age_wgt_norm_mesh = ds.len_age_weight_dist_all_normalized.sel(
        # # stratum=krig_results["stratum_num"].values).values
        # unaged_prop_mesh = ds.unaged_proportion.sel(
        #     stratum=self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
        # ).values

        # TODO: in Chu's current calculation there is a normalized term that is summed
        #  and included in the below expression as it is normalized this sum will always be one
        # Wgt_len_age_ALL_ii = (
        #     self.krig.survey.bio_calc.kriging_results_gdf["biomass_density_adult_mean"]
        #     * aged_prop_mesh
        #     * 1.0
        #     * self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        # )
        # Wgt_len_ALL_ii = (
        #     self.krig.survey.bio_calc.kriging_results_gdf["biomass_density_adult_mean"]
        #     * unaged_prop_mesh
        #     * 1.0
        #     * self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        # )
        #
        # # TODO: this may prove that we do not need to create the
        # #  vars Wgt_len_age_ALL_ii and Wgt_len_ALL_ii
        # self.krig.survey.bio_calc.kriging_results_gdf["test"] = (
        #     Wgt_len_age_ALL_ii + Wgt_len_ALL_ii
        # )

        # expand the bio parameters dataframe so that it corresponds to nasc_df
        averaged_weight_expanded = (
            self.krig.survey.bio_calc.bio_param_df.averaged_weight.loc[
                self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
            ]
        )

        # self.krig.survey.bio_calc.kriging_results_gdf["abundance"] = (
        #     Wgt_len_age_ALL_ii + Wgt_len_ALL_ii
        # ) / averaged_weight_expanded.values

        dist_weight_M_sum_expanded = (
            ds["dist_weight_M_sum"]
            .sel(
                stratum=self.krig.survey.bio_calc.kriging_results_gdf[
                    "stratum_num"
                ].values
            )
            .values
        )

        dist_weight_F_sum_expanded = (
            ds["dist_weight_F_sum"]
            .sel(
                stratum=self.krig.survey.bio_calc.kriging_results_gdf[
                    "stratum_num"
                ].values
            )
            .values
        )

        Wgt_len_age_M_ii = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_density_adult_mean"]
            * dist_weight_M_sum_expanded
            * self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        )

        Wgt_len_age_F_ii = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_density_adult_mean"]
            * dist_weight_F_sum_expanded
            * self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        )

        unaged_M_wgt_prop_expan = ds.unaged_M_wgt_proportion.sel(
            stratum=self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
        ).values

        unaged_F_wgt_prop_expan = ds.unaged_F_wgt_proportion.sel(
            stratum=self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
        ).values

        Wgt_len_M_ii = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_density_adult_mean"]
            * unaged_M_wgt_prop_expan
            * self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        )

        Wgt_len_F_ii = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_density_adult_mean"]
            * unaged_F_wgt_prop_expan
            * self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        )

        self.krig.survey.bio_calc.kriging_results_gdf["biomass_male"] = (
            Wgt_len_age_M_ii + Wgt_len_M_ii
        )

        self.krig.survey.bio_calc.kriging_results_gdf["biomass_female"] = (
            Wgt_len_age_F_ii + Wgt_len_F_ii
        )

        # calculate and add abundance to Kriging results
        self.krig.survey.bio_calc.kriging_results_gdf["abundance"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_adult"]
            / averaged_weight_expanded.values
        )

        self.krig.survey.bio_calc.kriging_results_gdf["abundance_male"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_male"]
            / averaged_weight_expanded.values
        )

        self.krig.survey.bio_calc.kriging_results_gdf["abundance_female"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_female"]
            / averaged_weight_expanded.values
        )

        self.krig.survey.bio_calc.kriging_results_gdf[
            "sig_b"
        ] = self.krig.survey.bio_calc.strata_sig_b.loc[
            self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
        ].values

        self.krig.survey.bio_calc.kriging_results_gdf["NASC"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["abundance"]
            * self.krig.survey.bio_calc.kriging_results_gdf["sig_b"]
        )
