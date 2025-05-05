from typing import Literal, Tuple, Dict

import numpy as np
import pandas as pd
import xarray as xr


# from the current fit_length_weight_relationship()
def length_weight_regression(
    df_specimen: pd.DataFrame,  # from df_bio_dict from load_data.py
    length_bins: np.array,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # NOTE: I think you can remove the fitted relationship from the output
    #       unless it is used elsewhere in the codebase
    # NOTE: does not need `species` as an input argument anymore
    #       since the loaded biological data will only contain 1 species
    #       from load_data.load_biological_data()
    # TODO: assemble df_length_bins from `length_bins` in the function if needed
    df_length_weight: pd.DataFrame
    df_regression: pd.DataFrame
    return df_length_weight, df_regression


# from the current quantize_number_counts()
def fish_count(
    df_specimen: pd.DataFrame,  # from df_bio_dict from load_data.py
    df_length: pd.DataFrame,  # from df_bio_dict from load_data.py
    aged: bool,
    sexed: bool,
) -> pd.DataFrame:
    # NOTE: make the function only do 1 thing for each call, not all of it, so output is simpler
    pass


# TODO: keeping inputs/outputs are dataframes for now,
#       think about changing to xarray dataarray later
#       When we are ready to use xarray, consider the following:
#       -- dimensions: stratum, sex, age_bin, length_bin
#       -- types of proportions
#          -- Age proportions across all strata
#          -- Sex proportions (male/female/unsexed) by stratum
#          -- Length proportions for aged fish by stratum and sex
#          -- Length proportions for unaged fish by stratum and sex
#       -- Some of the above can be derived from others
#       -- Some of the above can be used to form *_overall
# from the current number_proportions()
def number_proportions(
    df_aged_counts,
    df_unaged_counts,
) -> Dict[pd.DataFrame]:
    """
    Calculate number proportions for aged and unaged fish.

    Parameters
    ----------
    df_aged : pd.DataFrame
        Counts of aged fish by stratum, sex, length bin, and age bin
    df_unaged : pd.DataFrame
        Counts of unaged fish by stratum, sex, and length bin

    Returns
    -------
    xr.Dataset
        Dataset containing proportions with dimensions (stratum, sex, length, age)
    """
    # NOTE: Do only 1 species in one call
    pass


# TODO: keeping inputs/outputs are dataframes for now,
#       think about changing to xarray dataarray later
# The current `quantize_weights` function
def weight_distributions(
    df_specimen: pd.DataFrame,
    df_length: pd.DataFrame,
    df_length_weight: pd.DataFrame,  # from get_length_weight_regression() above
    aged: bool,
) -> pd.DataFrame:
    # NOTE: I cannot figure out what's going on in here
    # but I think the output is better represented as a N-D array
    # you can use xarray.Dataset to organize this type of data
    # it will also help with the slicing you have to do in the downstream weight_proportions
    if aged:
        df_sex_length_age: pd.DataFrame  # if xr.DataArray - dimensions: stratum, sex, length, age
        return df_sex_length_age  # the current specimen_table_sexed
    else:
        df_sex_length: pd.DataFrame  # if xr.DataArray - dimensions: stratum, sex, length, age
        return df_sex_length  # the current length_table_sexed


# TODO: keeping inputs/outputs are dataframes for now,
#       think about changing to xarray dataarray later
#       When we are ready to use xarray, output can be a dataarray:
#       -- dimensions: stratum, sex
#       -- values in the dataarray are the averaged weights for each sex and stratum combination
# The current `fit_length_weights` function
def stratum_averaged_weight() -> pd.DataFrame:
    """
    Calculate the weight proportion for each length-bin across all and sexed fish

    Parameters
    ----------
    proportions_dict: dict
        Dictionary containing multiple dataframes with aged and unaged quantized/binned proportions
    length_weight_dict: dict
        Dictionary containing length-weight regression terms and fitted values
    stratum_col: str
        Name of stratum column
    """
    df_fitted_weight: pd.DataFrame
    return df_fitted_weight


# NOTE: keeping inputs/outputs are dataframes for now,
#       think about changing to xarray dataarray later
#       When we are ready to use xarray, consider the following:
#       -- dimensions: stratum, sex, length, age
#       -- types of proportions:
#          -- Aged fish by stratum/sex/age/length
#          -- Unaged fish by stratum/sex/length
#          -- Combined aged+unaged proportions by sex
#          -- Overall aged vs unaged proportions
#       -- Some of the above can be derived from others
# The current weight_proportions()
def weight_proportions(
    df_catch: pd.DataFrame,
    df_number_proportion_dict: dict,  # number proportions: do we actually need this?
    df_weight_proportion_dict: dict,  # weight proportions
    df_length_weight: pd.DataFrame,  # length-weight regression
) -> Dict[pd.DataFrame]:
    pass
