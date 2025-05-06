from typing import Dict

import pandas as pd
import xarray as xr


def assemble_proportions(
    dict_df_number_proportion: Dict[pd.DataFrame],
    dict_df_weight_proportion: Dict[pd.DataFrame],
) -> xr.Dataset:
    """
    Assemble xr.Dataset of number and weight proportions from dictionaries of dataframes.

    Parameters
    ----------
    dict_df_number_proportion : dict
        Dictionary containing multiple dataframes with aged and unaged number proportions
    dict_df_weight_proportion : dict
        Dictionary containing multiple dataframes with aged and unaged weight proportions

    Returns
    -------
    xr.Dataset
        Dataset containing proportions across stratum, sex, age_bin, and length bin.

        # TODO: sketch of ds_proportions
        - dimensions: stratum, sex, length_bin, age_bin
        - variables:
          - abundance_unaged: (stratum, sex, length_bin)
          - abundance_aaged: (stratum, sex, length_bin, age_bin)
          - biomass_aged: (stratum, sex, length_bin, age_bin)
          - biomass_unaged: (stratum, sex, length_bin)


    """
    ds_proportions: xr.Dataset
    return ds_proportions


def apportion_biomass(
    df_nasc: pd.DataFrame,
    ds_proportions: xr.Dataset,
) -> xr.Dataset:
    """
    Apportion biomass across stratum, sex, age_bin, and length bin.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the apportioned biomass across age and length bins.
        Each
    """

    ds_kriged_apportioned: xr.Dataset
    # dimensions: x, y, sex, length_bin, age_bin
    # coorindates: lon, lat, sex, length_bin, age_bin
    # variables:
    # -- stratum: (lat, lon)
    # -- biomass_aged: (lat, lon, stratum, sex, length_bin, age_bin)
    # -- biomass_unaged: (lat, lon, sex, length_bin)

    return ds_kriged_apportioned


# The current `impute_kriged_values` function
def fill_missing_aged_from_unaged(
    ds_kriged_apportioned: xr.Dataset,
    ds_proportions: xr.Dataset,
) -> xr.Dataset:
    """
    Fill missing length bins in the aged dataset using unaged data.
    """
    pass


# The current section in biology.py that starts with comment:
# "# Additional reapportionment if age-1 fish are excluded"
def reallocate_age1(
    ds_kriged_apportioned: xr.Dataset,
    ds_proportions: xr.Dataset,
) -> xr.Dataset:
    """
    Reallocate age-1 biomass to age-2+ fish.
    """
    pass
