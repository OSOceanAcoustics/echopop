import pandas as pd
import xarray as xr


def _get_len_age_abundance(
    abundance_df: pd.DataFrame, ds: xr.Dataset, sex: str, kriging_vals: bool
):
    # TODO: document!

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
            abundance_df.groupby(level=0)
            .sum()["abundance_adult"]
            .to_xarray()
            .rename({"stratum_num": "stratum"})
        )
    else:
        abundance_sum_stratum = (
            abundance_df.groupby(level=0)
            .sum()["abundance"]
            .to_xarray()
            .rename({"stratum_num": "stratum"})
        )

    # get the abundance for the sex for each stratum (station 1)
    N_len = abundance_sum_stratum * Len_sex_proportion.sel(stratum=defined_stratum)

    # get the abundance for the sex for each stratum (station 2)
    N_len_age = abundance_sum_stratum * Len_Age_sex_proportion.sel(
        stratum=defined_stratum
    )

    Len_Age_Matrix_Acoust_unaged = (
        N_len * ds[f"len_dist_station1_normalized_{sex}"].sel(stratum=defined_stratum)
    ).sum("stratum")

    # get the abundance for the sex at each length and age bin
    Len_Age_Matrix_Acoust = (
        N_len_age * Len_Age_key_sex_norm.sel(stratum=defined_stratum)
    ).sum("stratum")

    # TODO: add this in when we implement the Kriging version
    # # redistribute the age 1 data if Kriging values are being used
    # if self.survey.params["exclude_age1"] and kriging_vals:
    #     self._redistribute_age1_data(Len_Age_Matrix_Acoust)

    return Len_Age_Matrix_Acoust, Len_Age_Matrix_Acoust_unaged


def _get_len_age_biomass(biomass_df, ds, kriging_vals: bool):
    # TODO: document!

    # obtain only those strata that are defined in biomass_df
    defined_stratum = biomass_df.index.unique().values

    # obtain the total biomass for each stratum
    if kriging_vals:
        biomass_sum_stratum = (
            biomass_df.groupby(level=0)
            .sum()["biomass_adult"]
            .to_xarray()
            .rename({"stratum_num": "stratum"})
        )
    else:
        biomass_sum_stratum = (
            biomass_df.groupby(level=0)
            .sum()["biomass"]
            .to_xarray()
            .rename({"stratum_num": "stratum"})
        )

    # get the abundance for the sex at each length and age bin
    Len_Age_Matrix_biomass = (
        biomass_sum_stratum
        * ds.len_age_weight_dist_all_normalized.sel(stratum=defined_stratum)
    ).sum("stratum")

    # TODO: add this in when we implement the Kriging version
    # redistribute the age 1 data if Kriging values are being used
    # if self.survey.params["exclude_age1"] and kriging_vals:
    #     self._redistribute_age1_data(Len_Age_Matrix_biomass)
    # TODO: need to redistribute the data in a special way for biomass

    return Len_Age_Matrix_biomass


def compute_length_age_variables(
    transect_results_gdf, ds, compute_transect: bool, compute_kriging: bool
):

    abundance_df = transect_results_gdf[["abundance", "stratum_num"]]

    # TODO: use this for Kriging
    # abundance_df = kriging_results_gdf[["abundance_adult", "stratum_num"]]

    temp_M, temp_unaged_M = _get_len_age_abundance(
        abundance_df, ds, sex="M", kriging_vals=False
    )
    temp_F, temp_unaged_F = _get_len_age_abundance(
        abundance_df, ds, sex="F", kriging_vals=False
    )

    temp_F.to_pandas()

    temp_total = temp_M + temp_F

    print(temp_total)
