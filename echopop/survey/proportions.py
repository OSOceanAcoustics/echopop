import functools
from itertools import product
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from .. import utils


def compute_binned_counts(
    data: pd.DataFrame,
    groupby_cols: List[str],
    count_col: str,
    agg_func: str = "size",
) -> pd.DataFrame:
    """
    Compute binned counts with grouping.

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
        Column-value pairs to exclude. Format: {column: value_to_exclude}

    Returns
    -------
    pd.DataFrame
        Grouped counts
    """
    df = data.copy()

    # Apply aggregation
    if agg_func == "size":
        return (
            df.groupby(groupby_cols, observed=False)[count_col]
            .size()
            .to_frame("count")
            .reset_index()
        )
    else:
        return (
            df.groupby(groupby_cols, observed=False)[count_col]
            .agg(agg_func)
            .to_frame("count")
            .reset_index()
        )


def number_proportions(
    data: Union[Dict[str, pd.DataFrame], pd.DataFrame],
    group_columns: List[str] = ["stratum_num"],
    exclude_filters: Dict[str, Any] = {},
) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Calculate number proportions from one or more count DataFrames.

    This function computes both within-group proportions and overall proportions
    across all provided dataframes. It can handle different DataFrame structures,
    as long as they all have a 'count' column and share the grouping columns.

    Parameters
    ----------
    data : Union[Dict[str, pd.DataFrame], pd.DataFrame]
        Either a dictionary of DataFrames or a single DataFrame containing count data.
        When providing a dictionary, keys will be used as column aliases if column_aliases is None.
    group_columns : List[str], default ["stratum_num"]
        Columns to group by for calculating totals.
    exclude_filters : Dict[str, Any], default {}}
        Filters to exclude rows from dataframes. This should match the same keys. When supplied an
        empty DataFrame, then no filters are applied. When a filter is expected to be applied to
        only one DataFrame, then all other key-specific dictionaries should be empty.

    Returns
    -------
    Union[pd.DataFrame, Dict[str, pd.DataFrame]]
        If only one DataFrame is provided, returns that DataFrame with added proportion columns.
        If multiple DataFrames are provided, returns a dictionary with keys based on column_aliases
        or default names.

    Notes
    -----
    The function calculates two types of proportions:
    - "proportion": Within-group proportion (count divided by total count within the same DataFrame
    /group)
    - "proportion_overall": Overall proportion (count divided by the sum of counts across all
    DataFrames)

    Examples
    --------
    >>> # Single DataFrame
    >>> result = number_proportions(aged_counts_df)
    >>>
    >>> # Dictionary of DataFrames (keys become aliases)
    >>> data_dict = {"aged": aged_counts_df, "unaged": unaged_counts_df}
    >>> result = number_proportions(data_dict)
    """

    # Handle DataFrame input
    if isinstance(data, pd.DataFrame):
        data = {"data": data}
        exclude_filters = {"data": exclude_filters}

    # Get the column aliases if not supplied
    column_aliases = list(data.keys())

    # Fill out filter dictionary
    data_keys = set(list(data.keys()) + list(exclude_filters.keys()))
    # ---- Complete the dictionary
    exclude_filters = {
        k: exclude_filters[k] if k in data and k in exclude_filters else {} for k in data_keys
    }
    # ---- Apply the filters to the entire dictionary
    data_dict = {
        k: utils.apply_filters(data[k], exclude_filter=exclude_filters[k]) for k in data_keys
    }

    # Transform to a list of dataframes
    df_list = list(data_dict[k] for k in column_aliases)

    # Apply filters if provided
    df_list = [utils.apply_filters(df, exclude_filter=exclude_filters) for df in df_list]

    if not df_list:
        raise ValueError("At least one DataFrame must be provided")

    # Ensure all inputs have "count" column
    for i, df in enumerate(df_list):
        if "count" not in df.columns:
            raise ValueError(f"DataFrame {i} does not have a 'count' column")

    # Initialize the DataFrame with unique group combinations
    count_total = pd.DataFrame(
        product(
            *(
                pd.concat([df[col] for df in df_list if col in df.columns]).unique()
                for col in group_columns
            )
        ),
        columns=group_columns,
    )  # Create dynamic column names based on number of DataFrames

    # Ensure we have enough aliases for all DataFrames
    if len(column_aliases) >= len(df_list):
        total_cols = [f"total_{i}" for i in column_aliases[: len(df_list)]]
    else:
        total_cols = [f"total_df_{i}" for i in range(len(df_list))]

    # Map totals for each DataFrame and generate merge operations
    merge_operations = map(
        lambda i: count_total.merge(
            df_list[i]
            .groupby(group_columns)["count"]
            .sum()
            .reset_index()
            .rename(columns={"count": total_cols[i]}),
            on=group_columns,
            how="left",
        ),
        range(len(df_list)),
    )

    # Reduce merge operations
    count_total = functools.reduce(
        lambda left, right: pd.merge(left, right, on=group_columns, how="left"),
        merge_operations,
        count_total,
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
        count_total_ridx = count_total_idx.reindex(
            df.index.values
        )  # Get the alias for this dataframe (for column naming)
        alias = column_aliases[i]

        # Calculate within-group proportion
        df["proportion"] = (df["count"] / count_total_ridx[total_col]).astype(float).fillna(0.0)

        # Calculate overall proportion
        df["proportion_overall"] = (
            (df["count"] / count_total_ridx["total_overall"]).astype(float).fillna(0.0)
        )

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
    contrast_vars: Optional[Union[str, List[str]]] = None,
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
    if "_global_" in interpolators:
        # Apply global interpolation
        global_interp = interpolators["_global_"]
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
    interpolate_regression: bool = True,
    table_index: List[str] = ["length_bin"],
    table_cols: List[str] = [],
    length_weight_dataset: Optional[pd.DataFrame] = None,
    include_filter: Optional[Dict[str, Any]] = None,
    contrast_vars: Optional[Union[str, List[str]]] = None,
) -> pd.DataFrame:
    """
    Process length-weight data with optional interpolation to create a weight table.

    Creates a pivot table of weights by length bins and optionally other stratification
    variables. Can use either direct weights (when `interpolate_regression=False`) or interpolated
    weights based on length-weight relationships.

    Parameters
    ----------
    length_dataset : pd.DataFrame
        Dataset with length measurements, must contain 'length_bin' and 'length_count'
        columns when using interpolation, or direct 'weight' values when not interpolating
    interpolate_regression : bool, default True
        Whether to use interpolation for weights. If True, interpolates weights based on
        length-weight relationships. Values from weights fitted to binned length values for a
        different dataset are used to for interpolating the lengths in `length_dataset`. If False,
        uses existing weight values in length_dataset.
    table_index : List[str], default ["length_bin"]
        Columns to use as index in the pivot table
    table_cols : List[str], default []
        Variable(s) to stratify the final pivot table by. These will be included as columns in the
        pivot table.
    length_weight_dataset : pd.DataFrame, optional
        Dataset with length-weight relationships. Required when interpolate=True. Must contain
        'length_bin' and 'weight_fitted' columns for interpolation. Not used when
        interpolate_regression=False.
    include_filter : Dict[str, Any], optional
        Filter to apply to both datasets (e.g., to include only certain sexes)
    contrast_vars : str, List[str], or None
        Variable(s) to use for contrast in interpolation. If None or empty list, a global
        interpolator is used.

    Returns
    -------
    pd.DataFrame
        Pivot table with weights stratified by specified variables

    Raises
    ------
    ValueError
        If interpolate_regression=True but length_weight_dataset is None

    Notes
    -----
    The function expects 'length_bin' objects to have a 'mid' property that
    returns the midpoint of each length interval.

    When interpolate_regression=False, rows with missing weights will be dropped.

    Examples
    --------
    >>> # Using interpolation with sex-specific length-weight relationships
    >>> weights_by_sex = binned_weights(
    ...     length_dataset=length_freq_df,
    ...     length_weight_dataset=length_weight_model,
    ...     interpolate_regression=True,
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
    if interpolate_regression and length_weight_dataset is None:
        raise ValueError("length_weight_dataset must be provided when interpolate_regression=True")

    # Apply filters if provided
    if include_filter:
        length_dataset = utils.apply_filters(length_dataset, include_filter)

    # Working copy of the dataset
    result_dataset = length_dataset.copy()

    if interpolate_regression:
        # Apply filters if provided
        if include_filter:
            # ---- This is applied here since `length_weight_dataset` is optional
            length_weight_dataset = utils.apply_filters(length_weight_dataset, include_filter)

        # Pivot to a long format for the interpolator, if needed
        if utils.is_df_wide(length_weight_dataset):
            length_weight_dataset_valid = length_weight_dataset.stack().reset_index(
                name="weight_fitted"
            )
        else:
            length_weight_dataset_valid = length_weight_dataset.copy()

        # Extract length from the interval categories
        length_weight_dataset_valid.loc[:, "length"] = (
            length_weight_dataset_valid.loc[:, "length_bin"].map(lambda x: x.mid).astype(float)
        )

        # Create interpolators
        interpolators = utils.group_interpolator_creator(
            grouped_data=length_weight_dataset_valid,
            independent_var="length",
            dependent_var="weight_fitted",
            contrast_vars=contrast_vars,
        )

        # Apply interpolation
        result_dataset = apply_weight_interpolation(
            target_df=result_dataset,
            interpolators=interpolators,
            dependent_var="weight",
            independent_var="length",
            count_col="length_count",
            contrast_vars=contrast_vars,
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
        weight_table = (
            result_dataset.pivot_table(
                index=table_index,
                columns=pivot_columns,
                values="weight",
                aggfunc="sum",
                observed=False,
            )
            .astype(float)
            .fillna(0.0)
        )
    else:
        # Without columns (single series)
        weight_table = (
            result_dataset.pivot_table(
                index=table_index, values="weight", aggfunc="sum", observed=False
            )
            .astype(float)
            .fillna(0.0)
        )

    return weight_table


def calculate_adjusted_proportions(
    group_keys: List[str],
    aggregate_table: pd.DataFrame,
    group_proportions_table: pd.DataFrame,
    group_by=List[str],
) -> pd.DataFrame:
    """
    Calculate adjusted proportions across multiple groups.

    Parameters
    ----------
    group_keys : List[str]
        List of group keys/names (e.g., ["aged", "unaged"])
    aggregate_table : pd.DataFrame
        Table with aggregate proportions by group
    group_proportions_table : pd.DataFrame
        Table with grouping-specific proportions by group

    Returns
    -------
    pd.DataFrame
        DataFrame with adjusted proportions with multi-index (group, group_by variable)

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
            # aggregate_table.loc[group] + sex_proportions_table.loc[first_group]
            aggregate_table.loc[group]
            + group_proportions_table.loc[first_group]
        )

    # Calculate first group proportion using all other adjusted proportions
    if len(group_keys) > 1:
        # adjusted_props[first_group] = sex_proportions_table.loc[first_group] / (
        #     sex_proportions_table.loc[first_group]
        #     + pd.concat([adjusted_props[g] for g in group_keys[1:]], axis=0)
        # )
        adjusted_props[first_group] = group_proportions_table.loc[first_group] / (
            group_proportions_table.loc[first_group]
            + pd.concat([adjusted_props[g] for g in group_keys[1:]], axis=0)
        )
    else:
        # If only one group, it gets 100% of the proportions
        # adjusted_props[first_group] = sex_proportions_table.loc[first_group]
        adjusted_props[first_group] = group_proportions_table.loc[first_group]

    # Combine into a single DataFrame with proper multi-index
    adjusted_proportions = pd.concat(
        # [adjusted_props[group] for group in group_keys], keys=group_keys, names=["group", "sex"]
        [adjusted_props[group] for group in group_keys],
        keys=group_keys,
        names=["group"] + group_by,
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
        Proportions by group, groupings, and length bin
    length_proportions_pvt_all : pd.DataFrame
        Aggregated proportions by group and length bin
    aggregate_table : pd.DataFrame
        Table with aggregate proportions by group
    adjusted_proportions : pd.DataFrame
        Adjusted proportions for each group and groupings
    group_keys : List[str]
        List of group keys/names (e.g., ["aged", "unaged"])

    Returns
    -------
    pd.DataFrame
        DataFrame with average weights by a grouping category defined by groupings

    Notes
    -----
    The function calculates weighted average weights for different grouping categories:
    - For "all", it uses aggregate proportions across all groups
    - For individual groupings, it uses sex-specific adjusted proportions

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

    # Initialize the dictionary to contain the grouped weights
    weights_dict = {}

    # Gather the index names
    non_length_level = [lvl for lvl in binned_weight_table_pvt.index.names if lvl != "length_bin"][
        0
    ]

    # Define the grouping categories
    grouping_categories = (
        binned_weight_table_pvt.index.get_level_values(non_length_level).unique().to_list()
    )

    # Weight for all sexes combined
    if "all" in grouping_categories and "all" in binned_weight_table_pvt.index:
        weight_all = binned_weight_table_pvt.loc["all"]["weight_fitted"].dot(
            sum(
                length_proportions_pvt_all.loc[group] * aggregate_table.loc[group]
                for group in group_keys
            )
        )
        weights_dict["all"] = weight_all

    # Weight for each group individually
    if non_length_level in adjusted_proportions.index.names:
        for grouping in [g for g in grouping_categories if g != "all"]:
            if grouping in adjusted_proportions.index.get_level_values(non_length_level):
                weight_base = binned_weight_table_pvt.loc[
                    "all" if "all" in binned_weight_table_pvt.index else grouping
                ]["weight_fitted"]
                weight_grouping = weight_base.dot(
                    sum(
                        length_proportions_pvt.loc[(group, grouping)]
                        * adjusted_proportions.loc[(group, grouping)]
                        for group in group_keys
                        if (group, grouping) in length_proportions_pvt.index
                    )
                )
                weights_dict[grouping] = weight_grouping

    # Combine into final DataFrame
    weights_df = pd.DataFrame(weights_dict)
    # ---- Apply name change
    weights_df.columns.name = non_length_level
    return weights_df


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
                within_group_proportion=df["proportion"]
                / df.groupby(group_cols)["proportion"].transform("sum")
            )
            .groupby(group_cols + ["length_bin"], observed=False)[
                ["within_group_proportion", "proportion"]
            ]
            .sum()
            .reset_index()
            .assign(group=key)
        )
        for key, df in proportions_dict.items()
    ]
    return pd.concat(series, axis=0)


def stratum_averaged_weight(
    proportions_dict: Dict[str, pd.DataFrame],
    binned_weight_table: pd.DataFrame,
    stratify_by: List[str] = ["stratum_num"],
    group_by: List[str] = [],
) -> pd.DataFrame:
    """
    Calculate stratum-specific average weights across multiple datasets with different proportions.

    This function combines length-at-age distributions from different groups (e.g., aged and unaged
    fish) to produce weighted average weights by sex and stratum. It accounts for different
    proportions of individuals in each group and correctly weights the contribution of each group
    to the final average.

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
        DataFrame with average weights defined for each stratum with optional groupings.

    Notes
    -----
    - This function assumes that the length bins in proportions_dict and binned_weight_table match
    - The function requires an 'all' group category in binned_weight_table for calculating combined
      weights
    - Missing strata or group categories will be excluded from the final results

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
        # group_cols=[stratum_col, "sex", "length_bin"],
        group_cols=stratify_by + group_by + ["length_bin"],
        index_cols=["group"],
        # strat_cols=[stratum_col],
        strat_cols=stratify_by,
        value_col="proportion_overall",
    )

    # Compute the within-grouped proportions
    # ---- Execute calculation
    length_proportions_df = calculate_within_group_proportions(
        # proportions_dict, group_cols=[stratum_col, "sex"]
        proportions_dict,
        group_cols=stratify_by + group_by,
    )

    # ---- Convert into a table for just the within-grouped proportions
    length_proportions_group = utils.create_pivot_table(
        length_proportions_df,
        # index_cols=["group", "sex", "length_bin"],
        # strat_cols=[stratum_col],
        index_cols=["group"] + group_by + ["length_bin"],
        strat_cols=stratify_by,
        value_col="within_group_proportion",
    )

    # ---- Create table for the combined within-grouped proportions needed for later calculations
    length_proportions_all = utils.create_pivot_table(
        length_proportions_df,
        index_cols=["group", "length_bin"],
        # strat_cols=[stratum_col],
        strat_cols=stratify_by,
        value_col="proportion",
    )

    # Convert the binned weights (with `sex="all"` included) into the correctly formatted table
    binned_weights = binned_weight_table.unstack("length_bin").to_frame("weight_fitted")

    # Create table for the overall sexed group proportions
    group_proportions = utils.create_grouped_table(
        proportions_dict,
        # group_cols=[stratum_col, "sex"],
        group_cols=stratify_by + group_by,
        index_cols=["group"] + group_by,
        # strat_cols=[stratum_col],
        strat_cols=stratify_by,
        value_col="proportion_overall",
    )

    # Compute the re-weighted proportions from the mixture of the different groups
    adjusted_proportions = calculate_adjusted_proportions(
        group_keys, aggregate_proportions, group_proportions, group_by
    )

    # Calculate final weights
    fitted_weight_df = calculate_grouped_weights(
        binned_weight_table_pvt=binned_weights,
        length_proportions_pvt=length_proportions_group,
        length_proportions_pvt_all=length_proportions_all,
        aggregate_table=aggregate_proportions,
        adjusted_proportions=adjusted_proportions,
        group_keys=group_keys,
    )

    return fitted_weight_df


def aggregate_stratum_weights(input_data, stratum_col="stratum_num"):
    """
    Aggregate weights by stratum across all groups in the input.

    This function sums weights by stratum for each group in the input data,
    handling both single DataFrame and dictionary of DataFrames cases.

    Parameters
    ----------
    input_data : Union[Dict[str, pd.DataFrame], pd.DataFrame]
        Either a DataFrame with multi-level columns including stratum_num
        or a dictionary of such DataFrames
    stratum_col : str, default "stratum_num"
        Column name for stratum identifier

    Returns
    -------
    pd.DataFrame
        DataFrame with stratum weights for each group, with stratum_num as index
        and group names as columns

    Examples
    --------
    >>> # With a dictionary of DataFrames
    >>> stratum_summary = aggregate_stratum_weights(dict_df_weight_distr)
    >>>
    >>> # With a single DataFrame
    >>> aged_summary = aggregate_stratum_weights(dict_df_weight_distr["aged"])
    """
    results = {}

    # Handle case where input is a single DataFrame
    if isinstance(input_data, pd.DataFrame):
        dict_df_weight_distr = {"data": input_data}
    else:
        dict_df_weight_distr = input_data

    for group_name, df in dict_df_weight_distr.items():
        # Find which level contains stratum_col
        stratum_level = None
        for i, name in enumerate(df.columns.names):
            if name == stratum_col:
                stratum_level = i
                break

        if stratum_level is None:
            print(f"Warning: No {stratum_col} level found in {group_name}")
            continue

        # Aggregate weights by stratum - this returns a Series indexed by stratum_col
        stratum_weights = df.T.groupby(level=stratum_level).sum().T.sum()

        # Store the series with appropriate name
        results[f"{group_name}"] = stratum_weights

    # Combine all results into a DataFrame
    if not results:
        return pd.DataFrame(index=pd.Index([], name=stratum_col))

    final_df = pd.DataFrame(results)

    # Fill NaN values with 0
    final_df.fillna(0, inplace=True)

    return final_df


def scale_weights_by_stratum(
    weights_df: Union[pd.Series, pd.DataFrame],
    reference_weights_df: pd.DataFrame,
    stratum_col: str = "stratum_num",
):
    """
    Scale weights in a DataFrame using reference weights by stratum.

    This function adjusts the weights in the input DataFrame to match the
    reference weight distribution by stratum while maintaining the relative
    proportions within each stratum.

    Parameters
    ----------
    weights_df : Union[pd.Series, pd.DataFrame]
        Either a Series or DataFrame with one- or multi-level columns including the stratum name
    reference_weights_df : pd.DataFrame
        DataFrame with stratum as index or column and weight column
    stratum_col : str, default "stratum_num"
        Column name for stratum identifier

    Returns
    -------
    pd.DataFrame
        DataFrame with weights scaled to match reference weight distribution

    Examples
    --------
    >>> standardized_weights = scale_weights_by_stratum(
    ...     weights_df=dict_df_weight_distr["unaged"],
    ...     reference_weights_df=stratum_weights,
    ...     stratum_col="stratum_ks"
    ... )

    Notes
    -----
    The reference_weights_df must have either the specified stratum_col as index or as a
    column along with a 'weight' column.
    """
    # Create copy
    reference_copy = reference_weights_df.copy()

    # Convert to a DataFrame, if needed
    if isinstance(reference_copy, pd.Series):
        reference_copy = reference_copy.to_frame()

    # Set index of reference_weights to stratum_num if it's not already
    if reference_copy.index.name is None:
        if stratum_col in reference_copy.columns:
            reference_copy.set_index(stratum_col, inplace=True)
        else:
            raise ValueError(
                f"reference_weights_df must have the defined `stratum_col` ({stratum_col}) "
                f"column or index!"
            )

    # Check for appropriate indexing
    if stratum_col not in (set([reference_copy.index.name]) | set(reference_copy.index.names)):
        raise ValueError(
            f"reference_weights_df must have the defined `stratum_col` ({stratum_col}) "
            f"column or index!"
        )

    # Check that stratum indices match
    if stratum_col not in (set([weights_df.columns.name]) | set(weights_df.columns.names)):
        raise ValueError(
            f"weights_df must have the defined `stratum_col` ({stratum_col}) as a column index!"
        )

    # Sum weights by stratum and sex
    summed_weights = weights_df.sum()

    # Calculate within-stratum proportions for each sex
    strata_totals = summed_weights.unstack(stratum_col).sum(axis=0)

    # Simple standardization: divide by strata totals and multiply by reference weights
    scaled = (
        (summed_weights / strata_totals).unstack(stratum_col) * reference_copy["weight"]
    ).fillna(0.0)

    # Fill any NaN values with 0
    return scaled


def weight_proportions(
    weight_data: Dict[str, pd.DataFrame],
    catch_data: pd.DataFrame,
    group: str,
    stratum_col: str = "stratum_num",
) -> pd.DataFrame:
    """
    Calculate stratified weight proportions of a dataset relative to total stratified weights.

    This function computes the proportional weight distribution of a specific group
    relative to the total weights across all strata, using haul weights from catch data.

    Parameters
    ----------
    weight_data : Dict[str, pd.DataFrame]
        Dictionary of DataFrames containing weight distributions for all groups
    catch_data : pd.DataFrame
        DataFrame with catch data including stratum and haul_weight columns
    group : str
        Identifier for the group being analyzed (e.g., 'aged', 'unaged')
    stratum_col : str, default "stratum_num"
        Column name for stratum identifier

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated weight proportions

    Examples
    --------
    >>> weight_props = weight_proportions(
    ...     weight_data=dict_df_weight_distr,
    ...     catch_data=dict_df_bio_binned_ks["catch"],
    ...     group="aged"
    ... )
    """
    # Compute the total weights per stratum from the biological data
    stratum_weights = catch_data.groupby([stratum_col])["weight"].sum().reset_index(name="weight")

    # Compute the total weights among the different groups
    stratum_summary = aggregate_stratum_weights(weight_data, stratum_col)

    # Get the total stratified weights across groups
    total_stratum_weights = (
        stratum_weights.set_index([stratum_col])["weight"] + stratum_summary[group]
    )

    # Prepare the data: reindex the DataFrame
    data_pvt = (
        weight_data[group]
        .stack(list(range(weight_data[group].columns.nlevels)), future_stack=True)
        .unstack(stratum_col)
    )

    # Compute the weight proportions relative to the global stratified total weights
    return data_pvt / total_stratum_weights


def scale_weight_proportions(
    weight_data: pd.DataFrame,
    reference_weight_proportions: pd.DataFrame,
    catch_data: pd.DataFrame,
    number_proportions: Dict[str, pd.DataFrame],
    binned_weights: pd.DataFrame,
    group: str,
    group_columns: List[str],
    stratum_col: str = "stratum_num",
) -> pd.DataFrame:
    """
    Calculate detailed weight proportions with adjustments for fitted weights and reference data.

    This function computes weight proportions considering reference data, length-binned
    proportions, and fitted weight tables to generate comprehensive weight distribution.

    Parameters
    ----------
    weight_data : pd.DataFrame
        DataFrame containing weight distribution data for the specific group.
    reference_weight_proportions : pd.DataFrame
        DataFrame with reference data proportions to compare against.
    catch_data : pd.DataFrame
        DataFrame containing catch data with stratum_num and haul_weight columns.
    number_proportions : Dict[str, pd.DataFrame]
        Dictionary containing proportion data by various grouping factors.
    binned_weights : pd.DataFrame
        DataFrame containing fitted weights by length bins.
    group : str
        Identifier for the group being analyzed (e.g., 'unaged').
    group_columns : List[str]
        List of column names used for grouping (e.g., ['sex']).
    stratum_col : str, default "stratum_num"
        Column name for stratum identifier.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated detailed weight proportions.

    Examples
    --------
    >>> props = scale_weight_proportions(
    ...     weight_data=standardized_unaged_sex_weights,
    ...     reference_data=weight_proportions,
    ...     catch_data=catch_data,
    ...     proportion_dict=proportion_dict,
    ...     binned_weight_table=binned_weight_table,
    ...     group="unaged",
    ...     group_columns=["sex"]
    ... )
    """

    # Compute the total weights per stratum from the biological data
    stratum_weights = catch_data.groupby([stratum_col])["weight"].sum().to_frame()

    # Compute the total weights among the different groups
    stratum_summary = aggregate_stratum_weights(weight_data, stratum_col)

    # Get the total stratified weights across groups
    total_stratum_weights = stratum_weights["weight"] + stratum_summary["data"]

    # Calculate the overall weight proportions relative to the total stratified weights
    weight_proportions_grouped_overall = weight_data / total_stratum_weights

    # Compute the within-group proportions
    weight_proportions_grouped = (
        (weight_proportions_grouped_overall / weight_proportions_grouped_overall.sum())
        .astype(float)
        .fillna(0.0)
    )

    # Get the length-binned number proportions
    length_proportions_df = calculate_within_group_proportions(
        number_proportions, group_cols=[stratum_col] + group_columns
    )

    # Reindex the number proportions to generate pivot table
    length_proportions_all = utils.create_pivot_table(
        length_proportions_df,
        index_cols=["group", "length_bin"],
        strat_cols=[stratum_col],
        value_col="proportion",
    )

    # Calculate the average weights per length bin within each stratum
    fitted_weights = length_proportions_all.loc[group].T * binned_weights.to_numpy()

    # Compute the average within-group weight proportions
    fitted_weight_proportions = (
        (fitted_weights.T / fitted_weights.sum(axis=1)).astype(float).fillna(0.0)
    )

    # Get the overall weight proportions per stratum for the reference data
    stratum_reference_proportions = reference_weight_proportions.sum()

    # Compute the complementary proportions
    stratum_data_proportions = 1 - stratum_reference_proportions

    # Distribute the relative proportions of the defined group relative to total stratified weights
    weight_length_proportions_overall = stratum_data_proportions * fitted_weight_proportions

    # Get all group keys
    group_keys_list = list(
        weight_data.index.to_frame(index=False)[group_columns]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )

    # Distribution the overall sex proportions over the within-group proportions to get the overall
    # proportions
    weight_proportions_overall = (
        pd.concat(
            [weight_length_proportions_overall]
            * len(group_keys_list),  # Create exactly the right number of copies
            keys=group_keys_list,
            names=group_columns,
        )
        * weight_proportions_grouped
    )

    return weight_proportions_overall


def get_nasc_proportions_slice(
    number_proportions: pd.DataFrame,
    ts_length_regression_parameters: Dict[str, float],
    stratify_by: List[str],
    include_filter: Dict[str, Any] = {},
    exclude_filter: Dict[str, Any] = {},
) -> pd.Series:
    """
    Calculate NASC proportions for a specific slice of the population

    This function computes weighted acoustic proportions by applying target strength-length
    regression to convert length-based number proportions into acoustic backscatter
    proportions. The target strength weighting accounts for the fact that larger fish
    contribute disproportionately to acoustic backscatter.

    Parameters
    ----------
    number_proportions : pd.DataFrame
        DataFrame containing number proportions with columns:
        - length_bin: length bins (with .mid attribute for midpoint values)
        - proportion: number proportions
        - Additional grouping columns (age_bin, sex, etc.)
    stratify_by : List[str]
        Column names for stratification (e.g., ["stratum_ks"])
    include_filter : Dict[str, Any], default {}
        Filter criteria to include specific grouping, e.g.:
        {"age_bin": [1], "sex": ["female"]}
    exclude_filter : Dict[str, Any], default {}
        Groups to exclude, e.g., {"length_bin": [small_lengths]}
    ts_length_regression_parameters : Dict[str, float]
        Target strength-length regression parameters:
        - slope: regression slope
        - intercept: regression intercept

    Returns
    -------
    pd.Series
        NASC proportions by strata, representing the fraction of total acoustic
        backscatter contributed by the filtered group

    Notes
    -----
    The target strength calculation follows: TS = slope * log10(length) + intercept
    Backscatter cross-section (sigma_bs) = 10^(TS/10)

    Examples
    --------
    >>> # Get NASC proportions for age-1 fish
    >>> nasc_props = get_nasc_proportions_slice(
    ...     number_proportions=aged_data,
    ...     include_filter={"age_bin": [1]}
    ... )
    """

    # Get length values from length bins
    length_vals = number_proportions["length_bin"].apply(lambda x: x.mid).astype(float).unique()

    # Compute equivalent sigma_bs using target strength-length regression
    sigma_bs_equiv = 10 ** (
        (
            ts_length_regression_parameters["slope"] * np.log10(length_vals)
            + ts_length_regression_parameters["intercept"]
        )
        / 10.0
    )

    # Create pivot table for total population (all groups)
    aggregate_table = utils.create_pivot_table(
        number_proportions,
        index_cols=["length_bin"],
        strat_cols=stratify_by,
        value_col="proportion",
    )

    # Calculate total weighted sigma_bs for all strata
    aggregate_weighted_sigma_bs = (aggregate_table.T * sigma_bs_equiv).T.sum(axis=0)

    # Create pivot table including filter dimensions
    filtered_population_table = utils.create_pivot_table(
        number_proportions,
        index_cols=list(include_filter.keys()) + ["length_bin"],
        strat_cols=stratify_by,
        value_col="proportion",
    )

    # Apply filter to extract target group
    target_group_table = utils.apply_filters(
        filtered_population_table, include_filter=include_filter, exclude_filter=exclude_filter
    )

    # Aggregate target group over length and strata dimensions
    target_group_aggregated = target_group_table.unstack("length_bin").sum().unstack(stratify_by)

    # Calculate weighted sigma_bs for target group
    target_weighted_sigma_bs = (target_group_aggregated.T * sigma_bs_equiv).T.sum(axis=0)

    # Calculate and return proportions
    return target_weighted_sigma_bs / aggregate_weighted_sigma_bs


def get_number_proportions_slice(
    number_proportions: pd.DataFrame,
    stratify_by: List[str],
    exclude_filter: Dict[str, Any] = {},
    include_filter: Dict[str, Any] = {},
) -> pd.Series:
    """
    Extract number proportions for a specific slice of the population

    This function creates pivot tables from number proportion data and applies
    inclusion/exclusion filters to extract proportions for specific population
    segments. Handles multiple stratification dimensions and returns appropriately
    aggregated results.

    Parameters
    ----------
    number_proportions : pd.DataFrame
        DataFrame with number proportions containing:
        - proportion: number proportions
        - length_bin: length bins
        - Additional grouping columns (age_bin, sex, etc.)
    stratify_by : List[str]
        Columns for stratification
    exclude_filter : Dict[str, Any], default {}
        Groups to exclude, e.g., {"length_bin": [small_lengths]}
    include_filter : Dict[str, Any], default {}
        Groups to include, e.g., {"age_bin": [1], "sex": ["female"]}

    Returns
    -------
    pd.Series
        Proportions by strata.

    Notes
    -----
    - Handles missing columns gracefully by intersecting with available data
    - For single stratification variable, sums over all other dimensions
    - For multiple variables, preserves structure based on data complexity

    Examples
    --------
    >>> # Get age-1 proportions by stratum
    >>> age1_props = get_number_proportions_slice(
    ...     aged_data,
    ...     include_filter={"age_bin": [1]}
    ... )

    >>> # Exclude small fish, include only females
    >>> female_large_props = get_number_proportions_slice(
    ...     aged_data,
    ...     exclude_filter={"length_bin": small_length_bins},
    ...     include_filter={"sex": ["female"]}
    ... )
    """

    # Determine index columns from filter keys and length_bin
    index_set = set(list(exclude_filter.keys()) + list(include_filter.keys()) + ["length_bin"])

    # Intersect with actually available columns/indices
    available_columns = set(list(number_proportions.index.names) + list(number_proportions.columns))
    index_set = index_set.intersection(available_columns)

    # Separate stratification columns from index columns
    column_set = set(stratify_by).difference(index_set)

    # Create pivot table
    group_table = utils.create_pivot_table(
        number_proportions,
        index_cols=list(index_set),
        strat_cols=list(column_set),
        value_col="proportion",
    )

    # Apply filters
    filtered_table = utils.apply_filters(
        group_table, include_filter=include_filter, exclude_filter=exclude_filter
    )

    # Handle aggregation if no filter is applied
    if len(include_filter) == 0 and len(exclude_filter) == 0:
        return filtered_table

    # Handle aggregation based on stratification complexity
    if len(stratify_by) == 1:
        # Single stratification variable - sum over all other dimensions
        aggregated_result = filtered_table.unstack("length_bin").sum().unstack(stratify_by)
        return aggregated_result.sum(axis=0)
    else:
        # Multiple stratification variables - preserve structure
        current_dimensions = set(
            list(filtered_table.index.names) + list(filtered_table.columns.names)
        )

        # If already at target dimensions, return as-is
        if len(current_dimensions.difference(set(stratify_by))) == 0:
            return filtered_table
        else:
            # Aggregate step by step
            intermediate_result = filtered_table.unstack("length_bin").sum()

            # Find additional dimensions to unstack
            remaining_dimensions = list(
                set(stratify_by).difference(intermediate_result.index.names)
            )

            # Fill in missing dimensions if needed
            if len(remaining_dimensions) == 0:
                remaining_dimensions = list(set(stratify_by).difference(["length_bin"]))

            # Final unstack
            return intermediate_result.unstack(remaining_dimensions)


def get_weight_proportions_slice(
    weight_proportions: pd.DataFrame,
    stratify_by: List[str],
    include_filter: Dict[str, Any] = {},
    exclude_filter: Dict[str, Any] = {},
    number_proportions: Dict[str, pd.DataFrame] = {},
    length_threshold_min: float = 0.0,
    weight_proportion_threshold: float = 1e-10,
) -> pd.Series:
    """
    Calculate weight proportions for a specific slice of the population with optional thresholding

    This function computes weight proportions for a target population group, with
    optional thresholding based on number proportions. The thresholding helps
    identify cases where estimates may be unreliable due to low sample sizes.

    Parameters
    ----------
    weight_proportions : pd.DataFrame
        DataFrame with weight proportions and hierarchical index structure
    stratify_by : List[str]
        Stratification columns (e.g., ["stratum_ks"])
    include_filter : Dict[str, Any], default {}
        Filter criteria for target group, e.g., {"age_bin": [1]}
    exclude_filter : Dict[str, Any], default {}
        Groups to exclude, e.g., {"length_bin": [small_lengths]}
    number_proportions : Dict[str, pd.DataFrame], default {}
        Number proportions required for thresholding particular length values. This argument is
        expected to be a dictionary with a key paired with at least one DataFrame
        (e.g. {"grouped": df}). When no value is supplied, `number_proportions` defaults to an
        empty dictionary with no thresholding applied to the subsequent number proportions.
    length_threshold_min : float, default 0.0
        Minimum length for threshold calculations (e.g., 10.0 cm). This is only used when number
        proportions are provided.
    weight_proportion_threshold : float, default 1e-10
        Threshold value for proportion comparisons

    Returns
    -------
    pd.Series
        Weight proportions by strata for the target group

    Notes
    -----
    Thresholding Logic:
    1. Calculate target group weight proportions normally
    2. If number_proportions provided, calculate corresponding number thresholds
    3. Apply threshold mask where both weight and number proportions are very small
    4. Set masked values to 0.0 to indicate unreliable estimates

    Examples
    --------
    >>> # Simple weight proportions without thresholding
    >>> weight_props = calculate_weight_proportions_slice(
    ...     weight_data,
    ...     stratify_by=["stratum_ks"],
    ...     include_filter={"age_bin": [1]}
    ... )

    >>> # With thresholding using number proportions
    >>> weight_props = calculate_weight_proportions_slice(
    ...     weight_data,
    ...     stratify_by=["stratum_ks"],
    ...     include_filter={"age_bin": [1]},
    ...     number_proportions={"aged": aged_nums, "unaged": unaged_nums},
    ...     length_threshold_min=10.0
    ... )
    """

    # Check if filter keys are present in weight proportions
    filter_keys_present = all(
        key in weight_proportions.index.names for key in include_filter.keys()
    )

    if not filter_keys_present:
        raise ValueError(
            f"Filter keys {list(include_filter.keys())} not found in weight_proportions index"
        )

    # Calculate basic weight proportions for target group
    target_weight_table = utils.apply_filters(
        weight_proportions, include_filter=include_filter, exclude_filter=exclude_filter
    )

    # Aggregate target group proportions
    target_group_weight_proportions = (
        target_weight_table.unstack("length_bin").sum().unstack(stratify_by).sum(axis=0)
    )

    # Normalize by total weight proportions
    total_weight_proportions = weight_proportions.sum(axis=0)
    proportions_weight = target_group_weight_proportions / total_weight_proportions

    # Apply thresholding if number proportions provided
    if len(number_proportions) > 0:
        # ---- Get all unique length values across datasets
        all_length_vals = np.concatenate(
            [
                df["length_bin"].apply(lambda x: x.mid).astype(float).unique()
                for df in number_proportions.values()
            ]
        )
        length_vals = np.unique(all_length_vals)
        # ---- Create length exclusion filter
        length_exclusion_filter = {
            "length_bin": length_vals[length_vals < length_threshold_min],
            **exclude_filter,
        }
        # ---- Get filtered number proportions for each dataset
        filtered_number_proportions_dict = {
            key: get_number_proportions_slice(
                df,
                stratify_by=stratify_by + ["length_bin"],
                exclude_filter=length_exclusion_filter,
                include_filter=include_filter,
            )
            for key, df in number_proportions.items()
        }
        # ---- Calculate element-wise product across all datasets
        filtered_number_proportions = functools.reduce(
            lambda df1, df2: df1.mul(df2, fill_value=0),
            filtered_number_proportions_dict.values(),
        ).sum()
        # ---- Apply threshold mask
        threshold_mask = (target_group_weight_proportions <= weight_proportion_threshold) & (
            filtered_number_proportions <= weight_proportion_threshold
        )
        # ---- Set masked values to 0
        proportions_weight[threshold_mask] = 0.0

    # Return the sliced weight proportions
    return proportions_weight
