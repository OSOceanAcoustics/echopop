from typing import Any, Dict, List, Optional, Union
import pandas as pd
import numpy as np
import functools
from itertools import product
from . import utils

def compute_binned_counts(
    data: pd.DataFrame,
    groupby_cols: List[str],
    count_col: str,
    agg_func: str = "size",
    exclude_filters: Dict[str, Any] = None
) -> pd.DataFrame:
    """
    Compute binned counts with grouping and optional filtering.
    
    Parameters
    ----------
    data : pd.DataFrame
        Input DataFrame
    groupby_cols : list
        Columns to group by
    count_col : str  
        Column to aggregate
    agg_func : str, default "size"
        Aggregation function to apply: "size", "sum", "count", etc.
    exclude_filters : dict, optional
        Column-value pairs to exclude. Format: {column: value_to_exclude}
        
    Returns
    -------
    pd.DataFrame
        Grouped counts
    """
    df = data.copy()
    
    # Apply exclusion filters
    if exclude_filters:
        for col, exclude_val in exclude_filters.items():
            if col in df.columns:
                df = df.loc[df[col] != exclude_val]
    
    # Apply aggregation
    if agg_func == "size":
        return (
            df.groupby(groupby_cols, observed=False)[count_col].size()
            .to_frame("count").reset_index()
        )
    else:
        return (
            df.groupby(groupby_cols, observed=False)[count_col]
            .agg(agg_func)
            .to_frame("count").reset_index()
        )

def number_proportions(
    *dataframes: pd.DataFrame, 
    group_columns: List[str] = ["stratum_num"],
    column_aliases: Optional[List[str]] = None,
    exclude_filters: Optional[Union[Dict[str, Any], List[Optional[Dict[str, Any]]]]] = None
) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Calculate number proportions from one or more count DataFrames.
    
    This function computes both within-group proportions and overall proportions
    across all provided dataframes. It can handle different DataFrame structures,
    as long as they all have a 'count' column and share the grouping columns.
    
    Parameters
    ----------
    *dataframes : pd.DataFrame
        One or more DataFrames containing count data. Each must have a 'count' column.
    group_columns : List[str], default ["stratum_num"]
        Columns to group by for calculating totals.
    column_aliases : List[str], optional
        Custom names for the dataframes, used in column naming and dictionary keys.
        If not provided, will use default names like "df_0", "df_1", etc.
    exclude_filters : Optional[Union[Dict[str, Any], List[Optional[Dict[str, Any]]]]], default None
        Filters to exclude rows from dataframes:
        - If Dict: Apply the same filter to all dataframes (current behavior)
        - If List: Apply each filter to its corresponding dataframe
          Use None for any dataframe that shouldn't be filtered
        
    Returns
    -------
    Union[pd.DataFrame, Dict[str, pd.DataFrame]]
        If only one DataFrame is provided, returns that DataFrame with added proportion columns.
        If multiple DataFrames are provided, returns a dictionary with keys based on column_aliases
        or default names ("df_0", "df_1", etc.)
        
    Notes
    -----
    The function calculates two types of proportions:
    - "proportion": Within-group proportion (count divided by total count within the same DataFrame/group)
    - "proportion_overall": Overall proportion (count divided by the sum of counts across all DataFrames)
    
    Examples
    --------
    >>> # Single DataFrame
    >>> result = number_proportions(aged_counts_df)  # Uses default "stratum_num" grouping
    >>> 
    >>> # Multiple DataFrames with aliases
    >>> result = number_proportions(
    ...     aged_counts_df, unaged_counts_df,
    ...     column_aliases=["aged", "unaged"]
    ... )
    """

    # Convert args to list for consistent processing
    df_list = list(dataframes)

    # Apply filters if provided
    if exclude_filters is not None:
        # Case 1: Dictionary - apply to all dataframes
        if isinstance(exclude_filters, dict):
            df_list = [utils.apply_filters(df, exclude_filter=exclude_filters) for df in df_list]        
        # Case 2: List - apply each filter to corresponding dataframe
        elif isinstance(exclude_filters, list):
            df_list = [
                utils.apply_filters(df, exclude_filter=filt) if filt is not None else df
                for df, filt in zip(df_list, exclude_filters)
            ]
    
    if not df_list:
        raise ValueError("At least one DataFrame must be provided")

    # Ensure all inputs have "count" column
    for i, df in enumerate(df_list):
        if "count" not in df.columns:
            raise ValueError(f"DataFrame {i} does not have a 'count' column")

    # Initialize the DataFrame with unique group combinations
    count_total = pd.DataFrame(
        product(*(
            pd.concat([df[col] for df in df_list if col in df.columns]).unique()
            for col in group_columns
        )),
        columns=group_columns
    )

    # Create dynamic column names based on number of DataFrames
    if column_aliases and len(column_aliases) >= len(df_list):
        total_cols = [f"total_{i}" for i in column_aliases[:len(df_list)]]
    else:
        total_cols = [f"total_df_{i}" for i in range(len(df_list))]

    # Map totals for each DataFrame and generate merge operations
    merge_operations = map(
        lambda i: count_total.merge(
            df_list[i].groupby(group_columns)["count"].sum().reset_index().rename(
                columns={"count": total_cols[i]}
            ),
            on=group_columns,
            how="left"
        ),
        range(len(df_list))
    )

    # Reduce merge operations
    count_total = functools.reduce(
        lambda left, right: pd.merge(
            left, right, on=group_columns, how="left"
        ),
        merge_operations,
        count_total
    ).fillna(0)

    # Calculate grand total across all DataFrames
    count_total["total_overall"] = count_total[total_cols].sum(axis=1)

    # Set index
    count_total_idx = count_total.set_index(group_columns)
    
    # Create an output dictionary to store processed DataFrames
    output_dict = {}

    # Calculate proportions for each DataFrame based on total_cols
    for i, total_col in enumerate(total_cols):
        # Get the current dataframe
        df = df_list[i].copy().set_index(group_columns)

        # Reindex, if needed, for correct broadcasting
        count_total_ridx = count_total_idx.reindex(df.index.values)
        
        # Get the alias for this dataframe (for column naming)
        alias = column_aliases[i] if column_aliases and i < len(column_aliases) else f"df_{i}"
        
        # Calculate within-group proportion
        df["proportion"] = (df["count"] / count_total_ridx[total_col]).astype(float).fillna(0.)
        
        # Calculate overall proportion
        df["proportion_overall"] = (
            df["count"] / count_total_ridx["total_overall"]
        ).astype(float).fillna(0.)
        
        # Add to output dictionary
        output_dict[f"{alias}"] = df.reset_index()

    # Return either the dictionary (if multiple DataFrames)
    # ---- or just the single DataFrame (if only one)
    result = list(output_dict.values())[0] if len(df_list) == 1 else output_dict
    return result
    

def apply_weight_interpolation(
    target_df: pd.DataFrame,
    interpolators: Dict,
    dependent_var: str,
    independent_var: str,    
    count_col: str,
    contrast_vars: Optional[Union[str, List[str]]] = None
) -> pd.DataFrame:
    """
    Apply weight interpolation to a DataFrame based on interpolators.
    
    Uses the provided interpolators to estimate weights for each row in the target
    DataFrame, optionally multiplying by a count column (e.g., for expanding 
    interpolated weights to a length frequency count).
    
    Parameters
    ----------
    target_df : pd.DataFrame
        DataFrame to apply interpolation to
    interpolators : Dict
        Dictionary of interpolator functions from group_interpolator_creator
    dependent_var : str
        Column name of the dependent variable to be created (e.g., 'weight')
    independent_var : str
        Column name of the independent variable (e.g., 'length')
    count_col : str
        Column containing count values to multiply interpolated values by
    contrast_vars : str, List[str], or None, default None
        Column(s) used as contrast variables in the interpolators.
        If None or empty list, uses the global interpolator.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added dependent_var column containing interpolated values
        (multiplied by count_col values if present)
    
    Notes
    -----
    When using contrast variables, rows with contrast values not found in the 
    interpolators will receive NaN values for the dependent variable.
    """
    # Create copy to avoid modifying original
    result_df = target_df.copy()
    
    # Check if we're using a global interpolator
    if '_global_' in interpolators:
        # Apply global interpolation
        global_interp = interpolators['_global_']
        if global_interp is not None:
            result_df[dependent_var] = result_df[independent_var].apply(
                lambda x: global_interp(x) if not pd.isna(x) else np.nan
            )
        else:
            result_df[dependent_var] = np.nan
    else:
        # Convert single string to list for consistent handling
        if isinstance(contrast_vars, str):
            contrast_vars = [contrast_vars]
        elif contrast_vars is None:
            contrast_vars = []
        
        # Define interpolation function using contrast variables
        def get_interpolated_value(row):
            # Create key tuple from contrast variables
            if len(contrast_vars) == 1:
                key = row[contrast_vars[0]]
            else:
                key = tuple(row[var] for var in contrast_vars)
            
            # Get independent variable value
            ind_val = row[independent_var]
            
            # Apply interpolator if available for this key
            if key in interpolators and interpolators[key] is not None and not pd.isna(ind_val):
                return interpolators[key](ind_val)
            
            return np.nan
        
        # Apply interpolation
        result_df[dependent_var] = result_df.apply(get_interpolated_value, axis=1)
    
    # Multiply by count if count_col exists and is in the dataframe
    if count_col in result_df.columns:
        result_df[dependent_var] *= result_df[count_col]
    
    return result_df


def binned_weights(
    length_dataset: pd.DataFrame,
    interpolate: bool = True,
    table_index: List[str] = ["length_bin"],
    table_cols: List[str] = [],
    length_weight_dataset: Optional[pd.DataFrame] = None,
    include_filter: Optional[Dict[str, Any]] = None,
    contrast_vars: Optional[Union[str, List[str]]] = None,  
) -> pd.DataFrame:
    """
    Process length-weight data with optional interpolation to create a weight table.
    
    Creates a pivot table of weights by length bins and optionally other stratification
    variables. Can use either direct weights (when interpolate=False) or interpolated 
    weights based on length-weight relationships.
    
    Parameters
    ----------
    length_dataset : pd.DataFrame
        Dataset with length measurements, must contain 'length_bin' and 'length_count'
        columns when using interpolation, or direct 'weight' values when not interpolating
    interpolate : bool, default True
        Whether to use interpolation for weights. If True, interpolates weights based on 
        length-weight relationships. If False, uses existing weight values in length_dataset.
    table_index : List[str], default ["length_bin"]
        Columns to use as index in the pivot table
    table_cols : List[str], default []
        Variable(s) to stratify the final pivot table by.
        These will be included as columns in the pivot table.
    length_weight_dataset : pd.DataFrame, optional
        Dataset with length-weight relationships. Required when interpolate=True.
        Must contain 'length_bin' and 'weight_fitted' columns for interpolation.
        Not used when interpolate=False.
    include_filter : Dict[str, Any], optional
        Filter to apply to both datasets (e.g., to include only certain sexes)
    contrast_vars : str, List[str], or None
        Variable(s) to use for contrast in interpolation.
        If None or empty list, a global interpolator is used.
    
    Returns
    -------
    pd.DataFrame
        Pivot table with weights stratified by specified variables
    
    Raises
    ------
    ValueError
        If interpolate=True but length_weight_dataset is None
    
    Notes
    -----
    The function expects 'length_bin' objects to have a 'mid' property that 
    returns the midpoint of each length interval.
    
    When interpolate=False, rows with missing weights will be dropped.
    
    Examples
    --------
    >>> # Using interpolation with sex-specific length-weight relationships
    >>> weights_by_sex = binned_weights(
    ...     length_dataset=length_freq_df,
    ...     length_weight_dataset=length_weight_model,
    ...     interpolate=True,
    ...     table_cols=["stratum_num", "sex"],
    ...     contrast_vars="sex"
    ... )
    
    >>> # Using direct weights (no interpolation) - no need for length_weight_dataset
    >>> specimen_weights = binned_weights(
    ...     length_dataset=specimen_df,
    ...     interpolate=False,
    ...     table_cols=["stratum_num", "sex", "age_bin"],
    ...     include_filter={"sex": ["female", "male"]}
    ... )
    """
    
    # Validation check
    if interpolate and length_weight_dataset is None:
        raise ValueError("length_weight_dataset must be provided when interpolate=True")
    
    # Apply filters if provided
    if include_filter:
        length_dataset = utils.apply_filters(length_dataset, include_filter)
    
    # Working copy of the dataset
    result_dataset = length_dataset.copy()
    
    if interpolate:
        # Apply filters if provided
        if include_filter:
            # ---- This is applied here since `length_weight_dataset` is optional
            length_weight_dataset = utils.apply_filters(length_weight_dataset, include_filter)
        
        # Extract length from the interval categories
        length_weight_dataset.loc[:, "length"] = length_weight_dataset.loc[:, "length_bin"].apply(
            lambda x: x.mid
        ).astype(float)
        
        # Create interpolators
        interpolators = utils.group_interpolator_creator(
            grouped_data=length_weight_dataset,
            independent_var="length",
            dependent_var="weight_fitted",
            contrast_vars=contrast_vars
        )
        
        # Apply interpolation
        result_dataset = apply_weight_interpolation(
            target_df=result_dataset,
            interpolators=interpolators,
            dependent_var="weight",
            independent_var="length",
            count_col="length_count",
            contrast_vars=contrast_vars
        )
    else:
        # For non-interpolation case, just use existing weight data
        # Filter out rows with missing weights
        result_dataset = result_dataset.dropna(subset=["weight"])

    # Prepare columns for pivot table
    contrast_list = []
    if contrast_vars:
        contrast_list = [contrast_vars] if isinstance(contrast_vars, str) else list(contrast_vars)

    cols_list = []
    if len(table_cols) > 0:
        cols_list = [table_cols] if isinstance(table_cols, str) else list(table_cols)
    
    # Combine and filter pivot columns
    pivot_columns = list(set(contrast_list + cols_list))
    pivot_columns = [col for col in pivot_columns if col in result_dataset.columns]
    
    # Create the pivot table
    if pivot_columns:
        # With stratification and/or contrast variables as columns
        weight_table = result_dataset.pivot_table(
            index=table_index,
            columns=pivot_columns,
            values="weight",
            aggfunc="sum",
            observed=False
        ).astype(float).fillna(0.0)
    else:
        # Without columns (single series)
        weight_table = result_dataset.pivot_table(
            index=table_index,
            values="weight",
            aggfunc="sum",
            observed=False
        ).astype(float).fillna(0.0)
    
    return weight_table

def calculate_adjusted_proportions(
    group_keys: List[str],
    aggregate_table: pd.DataFrame,
    sex_proportions_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate adjusted proportions across multiple groups.
    
    Parameters
    ----------
    group_keys : List[str]
        List of group keys/names (e.g., ["aged", "unaged"])
    aggregate_table : pd.DataFrame
        Table with aggregate proportions by group
    sex_proportions_table : pd.DataFrame
        Table with sex-specific proportions by group
        
    Returns
    -------
    pd.DataFrame
        DataFrame with adjusted proportions with multi-index (group, sex)
        
    Notes
    -----
    This function calculates adjusted proportions for multiple datasets.
    The first group is used as a reference, and other groups are adjusted
    based on their relative proportions.
    
    Examples
    --------
    >>> group_keys = ["aged", "unaged"]
    >>> adjusted = calculate_adjusted_proportions(
    ...     group_keys,
    ...     aggregate_table,
    ...     sex_proportions_table
    ... )
    """
    adjusted_props = {}
    first_group = group_keys[0]  # Use first group as reference
    
    # For all other groups, calculate relative to the first group
    for group in group_keys[1:]:
        adjusted_props[group] = aggregate_table.loc[group] / (
            aggregate_table.loc[group] + sex_proportions_table.loc[first_group]
        )
    
    # Calculate first group proportion using all other adjusted proportions
    if len(group_keys) > 1:
        adjusted_props[first_group] = sex_proportions_table.loc[first_group] / (
            sex_proportions_table.loc[first_group] + 
            pd.concat([adjusted_props[g] for g in group_keys[1:]], axis=0)
        )
    else:
        # If only one group, it gets 100% of the proportions
        adjusted_props[first_group] = sex_proportions_table.loc[first_group]
    
    # Combine into a single DataFrame with proper multi-index
    adjusted_proportions = pd.concat(
        [adjusted_props[group] for group in group_keys],
        keys=group_keys,
        names=["group", "sex"]
    )
    
    return adjusted_proportions


def calculate_grouped_weights(
    binned_weight_table_pvt: pd.DataFrame,
    length_proportions_pvt: pd.DataFrame,
    length_proportions_pvt_all: pd.DataFrame,
    aggregate_table: pd.DataFrame,
    adjusted_proportions: pd.DataFrame,
    group_keys: List[str],
) -> pd.DataFrame:
    """
    Calculate average weights by sex category across groups.
    
    Parameters
    ----------
    binned_weight_table_pvt : pd.DataFrame
        Pivot table with weight values by sex and length bin
    length_proportions_pvt : pd.DataFrame
        Proportions by group, sex, and length bin
    length_proportions_pvt_all : pd.DataFrame
        Aggregated proportions by group and length bin
    aggregate_table : pd.DataFrame
        Table with aggregate proportions by group
    adjusted_proportions : pd.DataFrame
        Adjusted proportions for each group and sex
    group_keys : List[str]
        List of group keys/names (e.g., ["aged", "unaged"])
        
    Returns
    -------
    pd.DataFrame
        DataFrame with average weights by sex category
        
    Notes
    -----
    The function calculates weighted average weights for different sex categories:
    - For "all", it uses aggregate proportions across all groups
    - For individual sexes, it uses sex-specific adjusted proportions
    
    Examples
    --------
    >>> weights = calculate_grouped_weights(
    ...     binned_weights,
    ...     length_proportions_group,
    ...     length_proportions_all,
    ...     aggregate_proportions,
    ...     adjusted_proportions,
    ...     ["aged", "unaged"]
    ... )
    """
    weights_dict = {}

    # Define the sex categories
    sex_categories = ["all", "female", "male"]
    
    # Weight for all sexes combined
    if "all" in sex_categories and "all" in binned_weight_table_pvt.index:
        weight_all = binned_weight_table_pvt.loc["all"]["weight_fitted"].dot(
            sum(
                length_proportions_pvt_all.loc[group] * aggregate_table.loc[group]
                for group in group_keys
            )
        )
        weights_dict["all"] = weight_all
    
    # Weight for each sex individually
    for sex in [s for s in sex_categories if s != "all"]:
        if sex in binned_weight_table_pvt.index.get_level_values(0):
            weight_base = binned_weight_table_pvt.loc["all" if "all" in binned_weight_table_pvt.index else sex]["weight_fitted"]
            weight_sex = weight_base.dot(
                sum(
                    length_proportions_pvt.loc[(group, sex)] * adjusted_proportions.loc[(group, sex)]
                    for group in group_keys if (group, sex) in length_proportions_pvt.index
                )
            )
            weights_dict[sex] = weight_sex
    
    # Combine into final DataFrame
    return pd.DataFrame(weights_dict)

def calculate_within_group_proportions(
    proportions_dict: Dict[str, pd.DataFrame],
    group_cols: List[str],
) -> pd.DataFrame:
    """
    Calculate within-group proportions for each DataFrame.
    
    Parameters
    ----------
    proportions_dict : Dict[str, pd.DataFrame]
        Dictionary of DataFrames with proportion data, keyed by group name
    group_cols : List[str]
        Columns to group by for calculating within-group proportions (e.g., ["stratum_num", "sex"])
    
    Returns
    -------
    pd.DataFrame
        Combined DataFrame with within-group proportions and original proportions
        
    Notes
    -----
    This function assumes each DataFrame has a 'proportion' column and a 'length_bin' column.
    The within-group proportion is calculated as the proportion divided by the sum of proportions
    within each grouping level.
    
    Examples
    --------
    >>> props_dict = {'aged': aged_df, 'unaged': unaged_df}
    >>> within_props = calculate_within_group_proportions(
    ...     props_dict,
    ...     ["stratum_num", "sex"]
    ... )
    """
    series = [
        (
            df.assign(
                within_group_proportion = df["proportion"] /
                df.groupby(group_cols)["proportion"].transform("sum")
            )
            .groupby(group_cols + ["length_bin"], observed=False)
            [["within_group_proportion", "proportion"]]
            .sum().reset_index().assign(group=key)
        )
        for key, df in proportions_dict.items()
    ]
    return pd.concat(series, axis=0)

def stratum_averaged_weight(
    proportions_dict: Dict[str, pd.DataFrame],
    binned_weight_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate stratum-specific average weights across multiple datasets with different proportions.
    
    This function combines length-at-age distributions from different groups (e.g., aged and unaged fish)
    to produce weighted average weights by sex and stratum. It accounts for different proportions
    of individuals in each group and correctly weights the contribution of each group to the final
    average.
    
    Parameters
    ----------
    proportions_dict : Dict[str, pd.DataFrame]
        Dictionary of DataFrames with proportion data, keyed by group name (e.g., 'aged', 'unaged').
        Each DataFrame must contain the following columns:
          - 'stratum_num': Stratum identifier
          - 'sex': Sex category (e.g., 'female', 'male')
          - 'length_bin': Length bin category
          - 'proportion': Proportion within each group/stratum
          - 'proportion_overall': Overall proportion across all strata
        
    binned_weight_table : pd.DataFrame
        DataFrame with fitted weight values by sex and length bin.
        Must contain the following columns:
          - 'sex': Sex category (e.g., 'female', 'male', 'all')
          - 'length_bin': Length bin category matching those in proportions_dict
          - 'weight_fitted': Weight values for each combination
        
    Returns
    -------
    pd.DataFrame
        DataFrame with average weights by stratum and sex.
        Index: stratum_num
        Columns: Sex categories ('all', 'female', 'male')
        
    Notes
    -----
    - This function assumes that the length bins in proportions_dict and binned_weight_table match
    - The function requires an 'all' sex category in binned_weight_table for calculating combined weights
    - Missing strata or sex categories will be excluded from the final results
    
    Examples
    --------
    >>> # Create proportion dictionaries
    >>> aged_data = pd.DataFrame({
    ...     'stratum_num': [1, 1, 2, 2],
    ...     'sex': ['female', 'male', 'female', 'male'],
    ...     'length_bin': ['(10, 20]', '(20, 30]', '(10, 20]', '(20, 30]'],
    ...     'proportion': [0.6, 0.4, 0.7, 0.3],
    ...     'proportion_overall': [0.3, 0.2, 0.35, 0.15]
    ... })
    >>> unaged_data = pd.DataFrame({
    ...     'stratum_num': [1, 1, 2, 2],
    ...     'sex': ['female', 'male', 'female', 'male'],
    ...     'length_bin': ['(10, 20]', '(20, 30]', '(10, 20]', '(20, 30]'],
    ...     'proportion': [0.5, 0.5, 0.6, 0.4],
    ...     'proportion_overall': [0.25, 0.25, 0.3, 0.2]
    ... })
    >>> props_dict = {'aged': aged_data, 'unaged': unaged_data}
    >>>
    >>> # Create weight table
    >>> weight_table = pd.DataFrame({
    ...     'sex': ['female', 'male', 'all', 'female', 'male', 'all'],
    ...     'length_bin': ['(10, 20]', '(10, 20]', '(10, 20]', '(20, 30]', '(20, 30]', '(20, 30]'],
    ...     'weight_fitted': [0.5, 0.4, 0.45, 1.2, 1.0, 1.1]
    ... })
    >>>
    >>> # Calculate average weights
    >>> weights = stratum_averaged_weight(props_dict, weight_table)
    >>> print(weights)
              all    female    male
    1      0.7750    0.8500   0.7000
    2      0.7850    0.8550   0.7100
    """
    # Get the grouping keys from the number proportions
    group_keys = list(proportions_dict.keys())

    # Create aggregate table with summed proportions per group (e.g. strata)
    aggregate_proportions = utils.create_grouped_table(
        proportions_dict, 
        group_cols=["stratum_num", "sex", "length_bin"], 
        index_cols=["group"],
        strat_cols=["stratum_num"], 
        value_col="proportion_overall"
    )

    # Compute the within-grouped proportions
    # ---- Execute calculation
    length_proportions_df = calculate_within_group_proportions(
        proportions_dict, 
        group_cols=["stratum_num", "sex"]
    )
    
    # ---- Convert into a table for just the within-grouped proportions
    length_proportions_group = utils.create_pivot_table(
        length_proportions_df, 
        index_cols=["group", "sex", "length_bin"],
        strat_cols=["stratum_num"], 
        value_col="within_group_proportion"
    )
    
    # ---- Create table for the combined within-grouped proportions needed for later calculations
    length_proportions_all = utils.create_pivot_table(
        length_proportions_df, 
        index_cols=["group", "length_bin"], 
        strat_cols=["stratum_num"], 
        value_col="proportion"
    )

    # Convert the binned weights (with `sex="all"` included) into a table
    binned_weights = utils.create_pivot_table(
        binned_weight_table, 
        ["sex", "length_bin"], 
        [], 
        "weight_fitted"
    )

    # Create table for the overall sexed group proportions
    group_proportions = utils.create_grouped_table(
        proportions_dict, 
        group_cols=["stratum_num", "sex"], 
        index_cols=["group", "sex"],
        strat_cols=["stratum_num"], 
        value_col="proportion_overall"
    )

    # Compute the re-weighted proportions from the mixture of the different groups
    adjusted_proportions = calculate_adjusted_proportions(
        group_keys, 
        aggregate_proportions, 
        group_proportions
    )

    # Calculate final weights
    fitted_weight_df = calculate_grouped_weights(
        binned_weights, 
        length_proportions_group, 
        length_proportions_all, 
        aggregate_proportions, 
        adjusted_proportions, 
        group_keys
    )

    return fitted_weight_df


# # TODO: keeping inputs/outputs are dataframes for now,
# #       think about changing to xarray dataarray later
# #       When we are ready to use xarray, output can be a dataarray:
# #       -- dimensions: stratum, sex
# #       -- values in the dataarray are the averaged weights for each sex and stratum combination
# # The current `fit_length_weights` function
# def stratum_averaged_weight(
#     df_length_weight: pd.DataFrame,
#     dict_df_bio: Dict[pd.DataFrame],
# ) -> pd.DataFrame:
#     """
#     Calculate the weight proportion for each length-bin across all and sexed fish

#     Parameters
#     ----------
#     proportions_dict: dict
#         Dictionary containing multiple dataframes with aged and unaged quantized/binned proportions
#     length_weight_dict: dict
#         Dictionary containing length-weight regression terms and fitted values
#     stratum_col: str
#         Name of stratum column
#     """
#     df_fitted_weight: pd.DataFrame
#     return df_fitted_weight


# # TODO: keeping inputs/outputs are dataframes for now,
# #       think about changing to xarray dataarray later
# #       When we are ready to use xarray, consider the following:
# #       -- dimensions: stratum, sex, length, age
# #       -- types of proportions:
# #          -- Aged fish by stratum/sex/age/length
# #          -- Unaged fish by stratum/sex/length
# #          -- Combined aged+unaged proportions by sex
# #          -- Overall aged vs unaged proportions
# #       -- Some of the above can be derived from others
# # The current weight_proportions()
# def weight_proportions(
#     df_catch: pd.DataFrame,
#     dict_df_weight_proportion: dict,  # weight proportions
#     df_length_weight: pd.DataFrame,  # length-weight regression
# ) -> Dict[pd.DataFrame]:
#     pass


# def assemble_proportions(
#     dict_df_number_proportion: Dict[pd.DataFrame],
#     dict_df_weight_proportion: Dict[pd.DataFrame],
# ) -> xr.Dataset:
#     """
#     Assemble xr.Dataset of number and weight proportions from dictionaries of dataframes.

#     Parameters
#     ----------
#     dict_df_number_proportion : dict
#         Dictionary containing multiple dataframes with aged and unaged number proportions
#     dict_df_weight_proportion : dict
#         Dictionary containing multiple dataframes with aged and unaged weight proportions

#     Returns
#     -------
#     xr.Dataset
#         Dataset containing proportions across stratum, sex, age_bin, and length bin.

#         # TODO: sketch of ds_proportions
#         - dimensions: stratum, sex, length_bin, age_bin
#         - variables:
#           - abundance_unaged: (stratum, sex, length_bin)
#           - abundance_aaged: (stratum, sex, length_bin, age_bin)
#           - biomass_aged: (stratum, sex, length_bin, age_bin)
#           - biomass_unaged: (stratum, sex, length_bin)


#     """
#     ds_proportions: xr.Dataset
#     return ds_proportions
