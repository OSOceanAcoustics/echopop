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
    Krig : ComputeKrigingVariables
        An initialized ComputeKrigingVariables object. Note that any change to
        self.krig will also change this object.
    """

    def __init__(self, krig=None):

        self.krig = krig

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
        age_bins_ind_M = self.krig.survey.bio_calc._get_bin_ind(
            input_arr_age_M, self.krig.survey.bio_calc.bio_hake_age_bin
        )
        age_bins_ind_F = self.krig.survey.bio_calc._get_bin_ind(
            input_arr_age_F, self.krig.survey.bio_calc.bio_hake_age_bin
        )

        for age_bin in range(len(age_bins_ind_M)):
            # bin those lengths that correspond to the lengths in the first age bin
            len_bin_ind_M = self.krig.survey.bio_calc._get_bin_ind(
                input_arr_len_M[age_bins_ind_M[age_bin]],
                self.krig.survey.bio_calc.bio_hake_len_bin,
            )

            len_bin_ind_F = self.krig.survey.bio_calc._get_bin_ind(
                input_arr_len_F[age_bins_ind_F[age_bin]],
                self.krig.survey.bio_calc.bio_hake_len_bin,
            )

            ds.sel(stratum=stratum).len_age_weight_dist_M[:, age_bin] = np.array(
                [
                    np.sum(input_arr_wgt_M[age_bins_ind_M[age_bin][i]])
                    for i in len_bin_ind_M
                ]
            )

            ds.sel(stratum=stratum).len_age_weight_dist_F[:, age_bin] = np.array(
                [
                    np.sum(input_arr_wgt_F[age_bins_ind_F[age_bin][i]])
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
        haul_vs_stratum = self.krig.survey.strata_df.reset_index()[
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

        # len_wgt_all = self.krig.survey.bio_calc._generate_length_val_conversion(
        #     len_name="length",
        #     val_name="weight",
        #     df=self.krig.survey.bio_calc.specimen_df[
        #         self.krig.survey.bio_calc.specimen_df["sex"] == 1
        #         ] + self.krig.survey.bio_calc.specimen_df[
        #         self.krig.survey.bio_calc.specimen_df["sex"] == 2
        #         ],
        # )

        # print(len_wgt_M)
        # print(len_wgt_F)

        for i in stratum_ind:

            ds.num_M.loc[i] = len(spec_drop_M.loc[i])
            ds.num_F.loc[i] = len(spec_drop_F.loc[i])

            self._set_len_age_distribution_data(
                stratum=i,
                ds=ds,
                df_M=spec_drop_M.loc[i],
                df_F=spec_drop_F.loc[i],
            )

            haul_nums = haul_vs_stratum[haul_vs_stratum["stratum_num"] == i][
                "haul_num"
            ].values

            wgt_station_1 = np.array(
                [
                    self.krig.survey.catch_df.loc[j]["haul_weight"].sum()
                    for j in haul_nums
                ]
            ).sum()

            wgt_station_2 = self.krig.survey.bio_calc.specimen_df.loc[i]["weight"].sum()

            ds.total_weight.loc[i] = wgt_station_1 + wgt_station_2

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

            print(f"male_wgt = {male_wgt}")
            print(f"female_wgt = {female_wgt}")

        ds["dist_weight_M"] = (
            ds.len_age_weight_dist_M_normalized * ds.len_age_weight_prop_M
        )
        ds["dist_weight_F"] = (
            ds.len_age_weight_dist_F_normalized * ds.len_age_weight_prop_F
        )

        return ds
