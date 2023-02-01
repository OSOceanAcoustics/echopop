from typing import List, Tuple

import geopandas as gpd
import pandas as pd
import xarray as xr


def _compute_len_age_abundance(
    abundance_df: pd.DataFrame, ds: xr.Dataset, sex: str, kriging_vals: bool
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Computes the abundance at each length and age bin for a specified
    gender using the input abundance DataFrame and parameter Dataset.

    Parameters
    ----------
    abundance_df: pd.DataFrame
        A DataFrame with index ``stratum_num`` and column with abundance
        values. If DataFrame corresponds to Kriging data then the column
        should be named ``abundance_adult``, else it should be named
        ``abundance``.
    ds: xr.Dataset
        A Dataset produced by the module ``parameters_dataset.py``, which
        contains all parameters necessary for computation
    sex: {'M', 'F'}
        A string specifying the gender for the abundance computation
    kriging_vals: bool
        If True, the abundance data was produced by Kriging, else
        it is Transect based

    Returns
    -------
    Len_Age_Matrix_Acoust: pd.DataFrame
        A DataFrame representing the abundance at each length and age
        bin for the specified gender
    Len_Age_Matrix_Acoust_unaged: pd.DataFrame
        A DataFrame representing the abundance for the unaged data
        at each length bin
    """

    # obtain only those strata that are defined in abundance_df
    defined_stratum = abundance_df.index.unique().values

    # compute the total number of animals in both stations
    total_N = ds.num_M + ds.num_F + ds.station_1_N

    # compute the proportion of the sex using station 2 data
    Len_Age_sex_proportion = ds[f"num_{sex}"] / total_N

    # compute the proportion of the sex using station 1 data
    Len_sex_proportion = ds[f"station_1_N_{sex}"] / total_N

    # get the normalized distribution of data in station 2
    Len_Age_key_sex_norm = ds[f"len_age_dist_{sex}"] / ds[f"len_age_dist_{sex}"].sum(
        ["len_bin", "age_bin"]
    )

    # sum together all abundance values in each stratum
    abundance_var_name = "abundance_adult" if kriging_vals else "abundance"
    abundance_sum_stratum = (
        abundance_df.groupby(level=0).sum()[abundance_var_name].to_xarray()
    )

    # get the abundance for the sex for each stratum (station 1)
    N_len = abundance_sum_stratum * Len_sex_proportion.sel(stratum_num=defined_stratum)

    # get the abundance for the sex for each stratum (station 2)
    N_len_age = abundance_sum_stratum * Len_Age_sex_proportion.sel(
        stratum_num=defined_stratum
    )

    # compute the abundance for the unaged data at each length bin
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
    # if exclude_age1 and kriging_vals:
    #     self._redistribute_age1_data(Len_Age_Matrix_Acoust)

    return Len_Age_Matrix_Acoust.to_pandas(), Len_Age_Matrix_Acoust_unaged.to_pandas()


def _compute_len_age_biomass(
    biomass_df: pd.DataFrame, ds: xr.Dataset, kriging_vals: bool
) -> pd.DataFrame:
    """
    Computes the biomass at each length and age bin using the
    input biomass DataFrame and parameter Dataset.

    Parameters
    ----------
    biomass_df: pd.DataFrame
        A DataFrame with index ``stratum_num`` and column with biomass
        values. If DataFrame corresponds to Kriging data then the column
        should be named ``biomass_adult``, else it should be named
        ``biomass``.
    ds: xr.Dataset
        A Dataset produced by the module ``parameters_dataset.py``, which
        contains all parameters necessary for computation
    kriging_vals: bool
        If True, the biomass data was produced by Kriging, else
        it is Transect based

    Returns
    -------
    Len_Age_Matrix_biomass: pd.DataFrame
        A DataFrame representing the biomass at each length and age
        bin for provided biomass data
    """

    # obtain only those strata that are defined in biomass_df
    defined_stratum = biomass_df.index.unique().values

    # sum together all biomass values in each stratum
    biomass_var_name = "biomass_adult" if kriging_vals else "biomass"
    biomass_sum_stratum = (
        biomass_df.groupby(level=0).sum()[biomass_var_name].to_xarray()
    )

    # get the biomass for the sex at each length and age bin
    Len_Age_Matrix_biomass = (
        biomass_sum_stratum
        * ds.len_age_weight_dist_all_normalized.sel(stratum_num=defined_stratum)
    ).sum("stratum_num")

    # TODO: add this in when we implement the Kriging version
    # TODO: need to redistribute the data in a special way for biomass
    # redistribute the age 1 data if Kriging values are being used
    # if exclude_age1 and kriging_vals:
    #     self._redistribute_age1_data(Len_Age_Matrix_biomass)

    return Len_Age_Matrix_biomass.to_pandas()


def get_len_age_abundance(
    gdf: gpd.GeoDataFrame, ds: xr.Dataset, kriging_vals: bool
) -> List[pd.DataFrame]:
    """
    Obtains and initiates the computation of the abundance at
    each length and age bin for male, female, and all
    gender data.

    Parameters
    ----------
    gdf: gpd.GeoDataFrame
        A GeoDataFrame with column ``stratum_num`` and a column with abundance
        values. If the GeoDataFrame corresponds to Kriging data then the column
        abundance column should be named ``abundance_adult``, else it should
        be named ``abundance``.
    ds: xr.Dataset
        A Dataset produced by the module ``parameters_dataset.py``, which
        contains all parameters necessary for computation
    kriging_vals: bool
        If True, the abundance data was produced by Kriging, else
        it is Transect based

    Returns
    -------
    len_age_abundance_list: list of pd.DataFrame
        A list where each element is a DataFrame containing the abundance at
        each length and age bin for male, female, and all genders (in
        that order)
    """

    # obtain abundance data
    if kriging_vals:
        abundance_df = gdf[["abundance_adult", "stratum_num"]]
    else:
        abundance_df = gdf[["abundance", "stratum_num"]]

    # make stratum_num the index of the DataFrame
    abundance_df = abundance_df.reset_index(drop=True).set_index("stratum_num")

    # compute, format, and store the abundance data for male and females
    len_age_abundance_list = []
    for sex in ["M", "F"]:

        # compute the abundance at each length and age bin for a gender
        aged_df, unaged_df = _compute_len_age_abundance(
            abundance_df, ds, sex=sex, kriging_vals=kriging_vals
        )

        # combine the aged and unaged abundance data
        final_df = pd.concat([aged_df, unaged_df], axis=1)

        # create and assign column names
        final_df.columns = ["age_bin_" + str(column) for column in aged_df.columns] + [
            "Un-aged"
        ]

        # remove row header produced by xarray
        final_df.index.name = ""

        # create and assign index names
        final_df.index = ["len_bin_" + str(row_ind) for row_ind in final_df.index]

        # store DataFrame
        len_age_abundance_list.append(final_df)

    # create and add the length-age abundance df for both genders to list
    len_age_abundance_list += [len_age_abundance_list[0] + len_age_abundance_list[1]]

    return len_age_abundance_list


def get_len_age_biomass(
    gdf_all: gpd.GeoDataFrame,
    gdf_male: gpd.GeoDataFrame,
    gdf_female: gpd.GeoDataFrame,
    ds: xr.Dataset,
    kriging_vals: bool,
) -> List[pd.DataFrame]:
    """
    Obtains and initiates the computation of the biomass at
    each length and age bin for male, female, and all
    gender data.

    Parameters
    ----------
    gdf_all: gpd.GeoDataFrame
        A GeoDataFrame with column ``stratum_num`` and column corresponding
        to the biomass produced by including all genders. If GeoDataFrame
        corresponds to Kriging data then the column should be named
        ``biomass_adult``, else it should be named ``biomass``.
    gdf_male: gpd.GeoDataFrame
        A GeoDataFrame with column ``stratum_num`` and column corresponding
        to the biomass produced by including only males. If GeoDataFrame
        corresponds to Kriging data then the column should be named
        ``biomass_adult``, else it should be named ``biomass``.
    gdf_female: gpd.GeoDataFrame
        A GeoDataFrame with column ``stratum_num`` and column corresponding
        to the biomass produced by including only females. If GeoDataFrame
        corresponds to Kriging data then the column should be named
        ``biomass_adult``, else it should be named ``biomass``.
    ds: xr.Dataset
        A Dataset produced by the module ``parameters_dataset.py``, which
        contains all parameters necessary for computation
    kriging_vals: bool
        If True, the biomass data was produced by Kriging, else
        it is Transect based

    Returns
    -------
    len_age_biomass_list: list of pd.DataFrame
        A list where each element is a DataFrame containing the biomass at
        each length and age bin for male, female, and all genders (in
        that order)
    """

    # compute, format, and store the biomass data for male, females, and all genders
    len_age_biomass_list = []
    for gdf in [gdf_male, gdf_female, gdf_all]:

        # obtain the biomass DataFrame
        # TODO: may need to account for Kriging here
        biomass_df = gdf[["biomass", "stratum_num"]]

        # make stratum_num the index of the DataFrame
        biomass_df = biomass_df.reset_index(drop=True).set_index("stratum_num")

        # obtain the biomass at the length and age bins
        final_df = _compute_len_age_biomass(biomass_df, ds, kriging_vals=kriging_vals)

        # remove column and row header produced by xarray
        final_df.columns.name = ""
        final_df.index.name = ""

        # create and assign column names
        final_df.columns = ["age_bin_" + str(column) for column in final_df.columns]

        # create and assign index names
        final_df.index = ["len_bin_" + str(row_ind) for row_ind in final_df.index]

        # store biomass data
        len_age_biomass_list.append(final_df)

    return len_age_biomass_list
