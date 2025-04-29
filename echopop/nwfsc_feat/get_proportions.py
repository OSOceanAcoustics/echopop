from typing import Literal, Tuple
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
    sexed: bool
) -> pd.DataFrame:
    # NOTE: make the function only do 1 thing for each call, not all of it, so output is simpler
    pass


# from the current number_proportions()
# TODO: write out the math will make the function logic cleaner
def number_proportion(
    df_aged, df_unaged, df_weight,  # only include the necessary dfs -- I lost track
    proportion_type: Literal["age", "sex", "unaged_length", "unaged_length"]
) -> pd.DataFrame:
    # NOTE: Do only 1 species in one call
    pass


# from the current quantize_weights()
def weight_distributions(
    df_specimen: pd.DataFrame,
    df_length: pd.DataFrame,
    df_length_weight: pd.DataFrame,  # from get_length_weight_regression() above
    aged: bool
) -> xr.DataArray:
    # NOTE: I cannot figure out what's going on in here
    # but I think the output is better represented as a N-D array
    # you can use xarray.Dataset to organize this type of data
    # it will also help with the slicing you have to do in the downstream weight_proportions
    if aged:
        da_sex_length_age: xr.DataArray  # dimensions: stratum, sex, length, age
        return da_sex_length_age  # the current specimen_table_sexed
    else:
        da_sex_length: xr.DataArray  # dimensions: stratum, sex, length, age
        return da_sex_length  # the current length_table_sexed



# from the current fit_length_weights()
# TODO: write out the math will make the function logic cleaner
def stratum_averaged_weight() -> pd.DataFrame:
    pass 


# from the current weight_proportions()
def weight_proportion():
    pass