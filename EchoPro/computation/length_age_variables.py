import pandas as pd
import xarray as xr


def _compute_len_age_abundance(
    abundance_df: pd.DataFrame, ds: xr.Dataset, sex: str, kriging_vals: bool
):
    """
    Computes the abundance at each length and age bin for a specified
    gender using the input abundance and parameter Dataset.

    Parameters
    ----------
    abundance_df: pd.DataFrame

    ds: xr.Dataset

    sex: str

    kriging_vals: bool

    Returns
    -------

    """

    # obtain only those strata that are defined in abundance_df
    defined_stratum = abundance_df.index.unique().values

    # compute the total number of animals in both stations
    total_N = ds.num_M + ds.num_F + ds.station_1_N

    # compute the proportion of the sex using station 2 data
    Len_Age_sex_proportion = ds[f"num_{sex}"] / total_N

    # compute the proportion of the sex using station 1 data
    # TODO: add this in to get the unaged bin
    Len_sex_proportion = ds[f"station_1_N_{sex}"] / total_N

    # get the normalized distribution of data in station 2
    Len_Age_key_sex_norm = ds[f"len_age_dist_{sex}"] / ds[f"len_age_dist_{sex}"].sum(
        ["len_bin", "age_bin"]
    )

    # sum together all abundance values in each stratum
    # TODO: we should probably rename ds coordinate to stratum_num
    if kriging_vals:
        abundance_sum_stratum = (
            abundance_df.groupby(level=0).sum()["abundance_adult"].to_xarray()
        )
    else:
        abundance_sum_stratum = (
            abundance_df.groupby(level=0).sum()["abundance"].to_xarray()
        )

    # get the abundance for the sex for each stratum (station 1)
    N_len = abundance_sum_stratum * Len_sex_proportion.sel(stratum_num=defined_stratum)

    # get the abundance for the sex for each stratum (station 2)
    N_len_age = abundance_sum_stratum * Len_Age_sex_proportion.sel(
        stratum_num=defined_stratum
    )

    Len_Age_Matrix_Acoust_unaged = (
        N_len
        * ds[f"len_dist_station1_normalized_{sex}"].sel(stratum_num=defined_stratum)
    ).sum("stratum_num")

    # get the abundance for the sex at each length and age bin
    Len_Age_Matrix_Acoust = (
        N_len_age * Len_Age_key_sex_norm.sel(stratum_num=defined_stratum)
    ).sum("stratum_num")

    # TODO: add this in when we implement the Kriging version
    # # redistribute the age 1 data if Kriging values are being used
    # if self.survey.params["exclude_age1"] and kriging_vals:
    #     self._redistribute_age1_data(Len_Age_Matrix_Acoust)

    return Len_Age_Matrix_Acoust, Len_Age_Matrix_Acoust_unaged


def _compute_len_age_biomass(biomass_df, ds, kriging_vals: bool):
    # TODO: document!

    # obtain only those strata that are defined in biomass_df
    defined_stratum = biomass_df.index.unique().values

    # obtain the total biomass for each stratum
    if kriging_vals:
        biomass_sum_stratum = (
            biomass_df.groupby(level=0).sum()["biomass_adult"].to_xarray()
        )
    else:
        biomass_sum_stratum = biomass_df.groupby(level=0).sum()["biomass"].to_xarray()

    # get the abundance for the sex at each length and age bin
    Len_Age_Matrix_biomass = (
        biomass_sum_stratum
        * ds.len_age_weight_dist_all_normalized.sel(stratum_num=defined_stratum)
    ).sum("stratum_num")

    # TODO: add this in when we implement the Kriging version
    # redistribute the age 1 data if Kriging values are being used
    # if self.survey.params["exclude_age1"] and kriging_vals:
    #     self._redistribute_age1_data(Len_Age_Matrix_biomass)
    # TODO: need to redistribute the data in a special way for biomass

    return Len_Age_Matrix_biomass


def get_len_age_abundance(gdf, ds, kriging_vals: bool):

    # TODO: document!

    if kriging_vals:
        abundance_df = gdf[["abundance_adult", "stratum_num"]]
    else:
        abundance_df = gdf[["abundance", "stratum_num"]]

    abundance_df = abundance_df.reset_index(drop=True).set_index("stratum_num")

    len_age_abundance_list = []
    for sex in ["M", "F"]:
        aged_da, unaged_da = _compute_len_age_abundance(
            abundance_df, ds, sex=sex, kriging_vals=False
        )

        aged_df = aged_da.to_pandas()

        final_df = pd.concat([aged_df, unaged_da.to_pandas()], axis=1)

        # create and assign column names
        final_df.columns = ["age_bin_" + str(column) for column in aged_df.columns] + [
            "Un-aged"
        ]

        # remove row header produced by xarray
        final_df.index.name = ""

        # create and assign index names
        final_df.index = ["len_bin_" + str(row_ind) for row_ind in final_df.index]

        len_age_abundance_list.append(final_df)

    # create and add the length age abundance df for both genders to list
    len_age_abundance_list += [len_age_abundance_list[0] + len_age_abundance_list[1]]

    return len_age_abundance_list


def get_len_age_biomass(gdf_all, gdf_male, gdf_female, ds, kriging_vals: bool):
    # TODO: document

    len_age_biomass_list = []
    for gdf in [gdf_male, gdf_female, gdf_all]:

        biomass_df = gdf[["biomass", "stratum_num"]]
        biomass_df = biomass_df.reset_index(drop=True).set_index("stratum_num")

        # obtain the biomass at the length and age bins
        final_df = _compute_len_age_biomass(
            biomass_df, ds, kriging_vals=kriging_vals
        ).to_pandas()

        # remove column and row header produced by xarray
        final_df.columns.name = ""
        final_df.index.name = ""

        # create and assign column names
        final_df.columns = ["age_bin_" + str(column) for column in final_df.columns]

        # create and assign index names
        final_df.index = ["len_bin_" + str(row_ind) for row_ind in final_df.index]

        len_age_biomass_list.append(final_df)

    return len_age_biomass_list
