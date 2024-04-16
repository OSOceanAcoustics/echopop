from functools import reduce
from typing import Callable, List, Union

import numpy as np
import pandas as pd

from ..utils.monkey_patch_dataframe import patch_method_to_DataFrame


# !!! TODO: Assess whether `bin_variable` is necessary
# ==== Currently works with dataframes as a whole, but
# ==== it is probably more efficient/sensible to refit
# ==== this function to create an array of bin intervals that
# ==== can then be used to define a dataframe column. It may
# ==== also make more sense to simply remove this function entirely
# ==== since it's not particularly onerous to directly use
# ==== `pd.cut(...)`. If this function were to stay, it'd probably
# ==== be useful to enable multiple `bin_values` and `bin_variable`
# ==== inputs so that the function doesn't need to be repeated for
# ==== multiple bin variables like age and length
@patch_method_to_DataFrame(pd.DataFrame)
def bin_variable(dataframe: pd.DataFrame, bin_values: np.ndarray, bin_variable: str):
    """
    Discretizes the target variable into user-defined bins

    Parameters
    ----------
    dataframe: pd.DataFrame
        A DataFrame object containing length and age variables required
        for analyses handling binned variables
    bin_values: np.ndarray
        An array of discrete bin values used to construct the intervals for
        the variable of interest
    bin_variable: str
        The variable that will be binned into discrete intervals
    Notes
    -----
    This will add a column to defined dataframes that groups length and age values
    into each bin that is explicitly defined in the function arguments
    """

    return dataframe.assign(  # input dataframe
        **{f"{bin_variable}_bin": lambda x: pd.cut(x[bin_variable], bin_values)}
    )  # assign bin


@patch_method_to_DataFrame(pd.DataFrame)
def bin_stats(
    dataframe: pd.DataFrame,
    bin_variable: str,
    bin_values: np.ndarray,
    contrasts: Union[str, List[str]] = [],
    variables: Union[str, List[str]] = ["length", "weight"],
    functions: Union[str, List[str]] = ["mean", "size"],
):
    """
    Quantizes dataset given user-defined intervals/bins

    Parameters
    ----------
    dataframe: pd.DataFrame
        A DataFrame object containing length and age variables required
        for analyses handling binned variables
    bin_values: np.ndarray
        An array of discrete bin values used to construct the intervals for
        the variable of interest
    bin_variable: str
        The variable that will be binned into discrete intervals
    contrasts: str or List[str]
        Additional variables to group data by, such as sex or species
    variables: str or List[str]
        Data columns to quantize
    functions: str or List[str] or List[Callable]
        Summary statistics or other functions applied to quantized variables
    """
    # Ensure that the contrasts, if there are any, are a list in case input is just a str
    con_lst = [contrasts] if isinstance(contrasts, str) else contrasts

    # Ensure variables are a list in case input is just a str
    var_lst = [variables] if isinstance(variables, str) else variables

    # Ensure functions are contained within a list
    fun_lst = [functions] if isinstance(functions, (str, Callable)) else functions

    # Rename functions as needed and if they exist
    FUNCTION_ALIASES = {
        "mean": "mean",
        "size": "n",
    }

    # Check against an allowed list
    invalid_functions = set(functions) - FUNCTION_ALIASES.keys()

    if invalid_functions:
        raise ValueError(
            f"""Invalid aggregation functions provided: {invalid_functions}. Only
            {FUNCTION_ALIASES.keys()} are allowed."""
        )

    # Construction the aggregation dictionary to apply the summary statistics toward
    aggregation_dict = {
        variable: [(f"{FUNCTION_ALIASES.get(agg, agg)}_{variable}", agg) for agg in fun_lst]
        for variable in var_lst
    }

    return (
        dataframe.bin_variable(  # input dataframe
            bin_values, bin_variable
        )  # discretize variable into bins )
        .groupby(
            [f"{bin_variable}_bin"] + con_lst, observed=False
        )  # group by these variables/contrasts
        .agg(aggregation_dict)  # apply specified functions
        .replace(np.nan, 0)  # replace NaN w/ 0's
        .droplevel(level=0, axis=1)  # drop the column indices
        .reset_index()  # reset the row indices
    )


@patch_method_to_DataFrame(pd.DataFrame)
def count_variable(
    dataframe: pd.DataFrame, contrasts: Union[str, List[str]], variable: str, fun: str
):
    """
    Quantizes dataset given user-defined intervals/bins

    Parameters
    ----------
    dataframe: pd.DataFrame
        A DataFrame object containing length and age variables required
        for analyses handling binned variables
    contrasts: str or List[str]
        Additional variables to group data by, such as sex or species
    variable: str
        Data column to bin and count
    """
    return (
        dataframe.reset_index(drop=True)  # input dataframe
        .groupby(contrasts, observed=False)
        .agg({variable: [("count", fun)]})
        .replace(np.nan, 0)
        .droplevel(level=0, axis=1)
        .reset_index()
        .sort_values(contrasts)
    )


@patch_method_to_DataFrame(pd.DataFrame)
def meld(specimen_dataframe: pd.DataFrame, length_dataframe: pd.DataFrame):
    """
    Concatenates the specimen and length dataframes using a shared format

    Parameters
    ----------
    specimen_dataframe: pd.DataFrame
        A DataFrame object containing data from the specimen dataset
    length_dataframe: pd.DataFrame
        A DataFrame object containing data from the length dataset
    """
    # Reorganize the specimen dataframe so it matches the format
    # of the length dataframe w/ length counts
    specimen_stacked = (
        specimen_dataframe.copy()
        .groupby(
            ["stratum_num", "species_id", "sex", "group", "station", "length", "length_bin"],
            observed=False,
        )[["length"]]
        .apply(lambda x: len(x), include_groups=True)
        .reset_index(name="length_count")
    )

    # Concatenate the data frames and return
    return pd.concat([specimen_stacked, length_dataframe], join="inner").reset_index(drop=True)


@patch_method_to_DataFrame(pd.DataFrame)
def stretch(
    dataframe,
    variable,
    variable_contrast="sex",
    index_variables=["transect_num", "latitude", "longitude", "stratum_num"],
    sep="_",
    suffix="\\w+",
):
    """
    Melts dataframe into a parseable format

    Parameters
    ----------
    dataframe: pd.DataFrame
        A DataFrame object containing pertinent biological data
    variable: str
        Data variable name
    variable_contrast: str
        The name of the column that will be used to index the data variable
    index_variables: str or List
        A list or string of additional indexing/metadata variables that will be joined to the
        data-of-interest
    sep: str
        A character indicating the separation of the variable names in the wide format, to be
        stripped from the names in the long format
    suffix: str
        A regular expression capturing the wanted suffixes
    """
    # Ensure variables are a list in case input is just a str
    idx_lst = [index_variables] if isinstance(index_variables, str) else index_variables

    # Prepare the dataframe for pivoting from wide to long
    # ---- Filter out the dataframe columns with the target index variables
    dataframe_reduced = (
        # Select the index variables
        dataframe.filter(items=idx_lst)
        # Join with the target variable
        .join(dataframe.filter(regex=variable)).drop_duplicates()
    )

    # Pivot from wide to long
    return pd.wide_to_long(
        df=dataframe_reduced,
        stubnames=variable,
        i=idx_lst,
        j=variable_contrast,
        sep=sep,
        suffix=suffix,
    ).reset_index()


@patch_method_to_DataFrame(pd.DataFrame)
def group_merge(dataframe, dataframes_to_add, inner_on, outer_on, how="outer", drop_na=True):
    """
    Group merge operation

    Parameters
    ----------
    dataframe: pd.DataFrame
        Parent dataframe
    dataframes_to_add: list
        Dataframes that will be added to the parent dataframe
    inner_on: str or list
        Index/column name that is used for an inner merge of dataframes
        within `dataframes_to_add`
    outer_on: str or List
        Index/column name that is used for an outer merge of both the parent
        and to-be-added dataframes
    how: str
        Merge method
    drop_na: boolean
        Flag determining whether np.nan values are removed

    Notes
    ----------
    This is a lazy method for recursively merging multiple dataframes at once
    """
    # Ensure that both the 'dataframes' and 'on' arguments are lists
    # dataframes
    df_lst = [dataframes_to_add] if isinstance(dataframes_to_add, str) else dataframes_to_add

    # ---- on
    inner_on_lst = [inner_on] if isinstance(inner_on, str) else inner_on
    outer_on_lst = [outer_on] if isinstance(outer_on, str) else outer_on

    # Merge the dataframes that will be joined to the original dataframe
    filtered_dataframes = [df for df in df_lst if any(col in df.columns for col in inner_on_lst)]

    # Excluded dataframes
    excluded_dataframes = [
        df for df in df_lst if any(col not in df.columns for col in inner_on_lst)
    ]

    # Begin inner merge
    innermost_frames_to_add = reduce(
        lambda left, right: pd.merge(left, right, on=inner_on_lst + outer_on_lst, how=how),
        filtered_dataframes,
    )
    inner_frames_to_add = reduce(
        lambda left, right: pd.merge(left, right, on=outer_on_lst, how=how),
        excluded_dataframes + [innermost_frames_to_add],
    )

    # Find union of column names that will be used to join everything
    union_lst = dataframe.filter(items=inner_frames_to_add.columns).columns.tolist()

    # Merge together and remove NaN values depending on argument 'drop_na'
    if drop_na:
        merged_frame = dataframe.dropna().merge(
            inner_frames_to_add.dropna(), on=union_lst, how="outer"
        )
    else:
        merged_frame = dataframe.merge(inner_frames_to_add, on=union_lst, how="outer")

    # Return output
    return merged_frame
