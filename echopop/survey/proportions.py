import functools
import operator
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import xarray as xr

from .. import utils


def compute_binned_counts(
    data: pd.DataFrame,
    groupby_cols: List[str],
    count_col: str,
    agg_func: str = "size",
) -> xr.DataArray:
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
    xr.DataArray
        Grouped counts with coordinates based on the column names provided in 'groupby_cols'.
    """

    # Copy the dataset
    df = data.copy()

    # Apply aggregation
    aggregate = df.groupby(groupby_cols, observed=False)[count_col].agg(agg_func).to_frame("count")

    # Convert to xarray.DataArray
    return aggregate.to_xarray()["count"]


def number_proportions(
    data: Union[xr.DataArray, xr.Dataset],
    group_columns: list = [],
    exclude_filters: dict = {},
) -> Union[xr.DataArray, xr.Dataset]:
    """
    Calculate number proportions from one or more xarray DataArrays or Datasets.

    This function computes both within-group proportions and overall proportions
    across all provided xarray objects. It can handle different xarray structures,
    as long as they all have a 'count' variable and share the grouping columns.

    Parameters
    ----------
    data : Union[xr.DataArray, xr.Dataset]
        Either a single xarray.DataArray or an xarray.Dataset containing multiple named
        xarray.DataArray objects.
    group_columns : list, default ["stratum_num"]
        Dimensions to group by for calculating totals.
    exclude_filters : dict, default {}
        Filters to exclude coordinates from xarray objects. For a single DataArray, this should be
        a dict like {"sex": "unsexed"}; for a Dataset, a dict of dicts keyed by variable name.

    Returns
    -------
    Union[xr.DataArray, Dict[xr.Dataset]]
        If only one DataArray is provided, returns that xarray.DataArray with added proportion
        coordinates.
        If a Dataset is provided, returns a dictionary of xarray.Datasets for each variable
        within 'data'.

    Notes
    -----
    The function calculates two types of proportions:
    - "proportion": Within-group proportion (count divided by total count within the same DataArray
    /group)
    - "proportion_overall": Overall proportion (count divided by the sum of counts across all
    DataArrays within the original Dataset)

    Examples
    --------
    >>> # Single DataArray
    >>> result = number_proportions(aged_counts)
    >>>
    >>> # Dataset
    >>> dataset_input = Dataset({"aged": aged_counts, "unaged": unaged_counts})
    >>> result = number_proportions(dataset_input)
    """

    # Handle xr.DataArray
    if isinstance(data, xr.DataArray):
        data = xr.Dataset({"data": data})
        # ---- Ensure exclude_filters matches the key
        exclude_filters = {"data": exclude_filters}
    elif isinstance(data, xr.Dataset):
        # ---- Check that all keys in exclude_filters exist in the dataset
        missing_keys = [k for k in exclude_filters if k not in data]
        if missing_keys:
            missing_str = "', '".join(missing_keys)
            raise KeyError(f"Variables in 'exclude_filters' not found in 'data': '{missing_str}'.")

    # Apply exclude_filters to each variable, transform to DataFrame, filter, then back to DataArray
    filtered_vars = {}
    for varname, da in data.data_vars.items():
        filt = exclude_filters.get(varname, {})
        # ---- Subset filter with overlapping coordinates
        filt_pruned = {k: v for k, v in filt.items() if k in data.coords}
        if filt:
            da = da.drop_sel(filt_pruned, errors="ignore")
        # ---- Convert back to DataArray
        if da.size == 0:
            # ---- If all filtered out, create empty DataArray with same dims
            filtered_vars[varname] = da.isel({d: slice(0) for d in da.dims})
        else:
            filtered_vars[varname] = da

    # Compute group totals and overall total for each variable
    group_totals = {}
    for varname, da in filtered_vars.items():
        if not group_columns:
            group_total = da.copy()
        else:
            sum_dims = [d for d in da.dims if d not in group_columns]
            if sum_dims:
                group_total = da.sum(dim=sum_dims, skipna=True)
            elif group_columns:
                group_total = da.groupby(group_columns).sum()
            else:
                group_total = da.sum(skipna=True)
        group_totals[varname] = group_total

    # Calculate grand totals
    overall_total = xr.concat(list(group_totals.values()), dim="source").sum(dim="source")

    # Compute proportions for each variable
    output = {}
    for varname, da in filtered_vars.items():
        sum_dims = [d for d in da.dims if d not in group_columns]
        if sum_dims:
            group_total_full = group_totals[varname]
            # ---- Expand group_total_full to match da's shape
            missing_dims = set(sum_dims) - set(group_total_full.dims)
            if missing_dims:
                expand_dict = {d: da.coords[d] for d in missing_dims}
                group_total_full = group_total_full.expand_dims(expand_dict)
            proportion = da / group_total_full
            proportion_overall = da / overall_total
        elif group_columns:
            # ---- Broadcast
            sel_dict = {g: da.coords[g] for g in group_columns}
            proportion = da / group_totals[varname].sel(**sel_dict)
            proportion_overall = da / overall_total.sel(**sel_dict)
        else:
            proportion = da / group_totals[varname]
            proportion_overall = da / overall_total
        # ---- Attach as new DataArrays with appropriate names
        output[varname] = xr.Dataset(
            {
                "count": da,
                "proportion": proportion,
                "proportion_overall": proportion_overall,
            }
        )

    # Return either the dictionary (if multiple Datasets)
    # ---- or just the single Dataset (if only one)
    return list(output.values())[0] if len(output) == 1 else output


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
    length_data: pd.DataFrame,
    interpolate_regression: bool = True,
    group_columns: list = ["length_bin"],
    length_weight_data: Optional[xr.DataArray] = None,
    include_filter: Optional[Dict[str, Any]] = None,
) -> xr.DataArray:
    """
    Process length-weight data and return an xarray.DataArray of weights binned by group_columns.

    This function creates a DataArray of weights by length bins and optionally other grouping
    variables. If interpolate_regression is True, weights are interpolated using a regression
    model provided as a DataArray. If False, weights are taken directly from the length_dataset.
    Filtering is applied using include_filter, which keeps only rows matching the filter.

    Parameters
    ----------
    length_dataset : pd.DataFrame
        Dataset with length measurements. Must contain 'length_bin' and 'length_count'
        columns when using interpolation, or direct 'weight' values when not interpolating.
    interpolate_regression : bool, default True
        Whether to use interpolation for weights. If True, interpolates weights based on
        length-weight relationships from length_weight_dataset. If False, uses existing
        weight values in length_dataset.
    group_columns : list, default ["length_bin"]
        Columns to use as coordinates in the output DataArray.
    length_weight_dataset : xr.DataArray, optional
        DataArray with length-weight relationships. Required when interpolate_regression=True.
        Must contain 'length_bin' as a coordinate.
    include_filter : Dict[str, Any], optional
        Filter to apply to both datasets (e.g., to include only certain sexes).
        Rows are kept if they match the filter.

    Returns
    -------
    xr.DataArray
        DataArray of weights binned by group_columns.

    Raises
    ------
    ValueError
        If interpolate_regression=True but length_weight_dataset is None or missing 'length_bin'.

    Notes
    -----
    - The output DataArray will have coordinates as specified by group_columns.
    - Filtering is performed before binning/interpolation.
    - If interpolate_regression is True, interpolation is performed using the regression model
      in length_weight_dataset, which is first converted to a DataFrame for compatibility.
    - If interpolate_regression is False, weights are taken directly from length_dataset.

    Examples
    --------
    >>> # Using interpolation with regression model
    >>> weights = binned_weights(
    ...     length_dataset=length_freq_df,
    ...     length_weight_dataset=length_weight_da,
    ...     interpolate_regression=True,
    ...     group_columns=["length_bin", "sex"],
    ...     include_filter={"sex": ["female", "male"]}
    ... )
    >>> # Using direct weights (no interpolation)
    >>> weights = binned_weights(
    ...     length_dataset=specimen_df,
    ...     interpolate_regression=False,
    ...     group_columns=["length_bin", "sex"],
    ...     include_filter={"sex": ["female", "male"]}
    ... )
    """

    # Apply include_filter if provided
    if include_filter:
        length_data = utils.apply_filters(length_data, include_filter=include_filter)

    # Check that all group_columns are present in length_dataset
    missing_cols = [col for col in group_columns if col not in length_data.columns]
    if missing_cols:
        missing_str = "', '".join(missing_cols)
        raise KeyError(f"Columns in 'group_columns' not found in 'length_data': '{missing_str}'.")

    # Ensure 'length_bin' is always included in group_columns
    if "length_bin" not in group_columns:
        group_columns = ["length_bin"] + list(group_columns)

    # Straightforward case without interpolation
    if not interpolate_regression:
        return (
            length_data.groupby(group_columns, observed=False)["weight"]
            .sum()
            .to_xarray()
            .fillna(0.0)
        )

    # Verify compatible arguments
    if length_weight_data is None:
        raise ValueError(
            "'length_weight_data' must be provided when 'interpolate_regression'=True."
        )

    # Validate coordinates
    if "length_bin" not in length_weight_data.coords:
        raise ValueError("'length_weight_data' must have 'length_bin' as a coordinate.")

    # Match the inclusion filter, if needed
    if include_filter:
        length_weight_data = length_weight_data.sel(include_filter)

    # Convert to a DataFrame for interpolation
    length_weight_data_cnv = length_weight_data.to_dataframe().reset_index()

    # Extract the defining lengths for each interval category of length_bin
    length_weight_data_cnv["length"] = (
        length_weight_data_cnv["length_bin"].map(lambda x: x.mid).astype(float)
    )

    # Define contrast variables
    contrast_vars = list(set(list(length_weight_data.coords.keys())) - set(["length_bin"]))

    # Create interpolation generators
    interpolators = utils.group_interpolator_creator(
        grouped_data=length_weight_data_cnv,
        independent_var="length",
        dependent_var="weight_fitted",
        contrast_vars=contrast_vars,
    )

    # Apply the interpolation
    result = apply_weight_interpolation(
        target_df=length_data.copy(),
        interpolators=interpolators,
        dependent_var="weight",
        independent_var="length",
        count_col="length_count" if "length_count" in length_data.columns else None,
        contrast_vars=contrast_vars,
    )

    # Return the grouped sums
    return (
        result.groupby(group_columns, observed=False)["weight"]
        .sum()
        .to_xarray()
        .fillna(0.0)
        .astype(float)
    )


def calculate_adjusted_proportions(
    group_keys: List[str],
    aggregate_proportions: xr.DataArray,
    group_proportions: xr.DataArray,
) -> xr.DataArray:
    """
    Calculate adjusted proportions across multiple groups.

    Parameters
    ----------
    group_keys : List[str]
        List of group keys/names (e.g., ["aged", "unaged"])
    aggregate_table : xr.DataArray
        Table with aggregate proportions by group
    group_proportions_table : xr.DataArray
        Table with grouping-specific proportions by group

    Returns
    -------
    xr.DataArray
        DataArray with adjusted proportions

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

    # Initialize the dictionary for the group-specific DataArrays
    adjusted_props = {}
    # ---- Use the first group as the reference
    first_group = group_keys[0]

    # For all other groups, calculate relative to the first group
    for group in group_keys[1:]:
        adjusted_props[group] = aggregate_proportions.sel(group=group) / (
            aggregate_proportions.sel(group=group) + group_proportions.sel(group=first_group)
        )

    # Calculate first group proportion using all other adjusted proportions
    if len(group_keys) > 1:
        other_adjusted = sum(adjusted_props[g] for g in group_keys[1:])
        adjusted_props[first_group] = group_proportions.sel(group=first_group) / (
            group_proportions.sel(group=first_group) + other_adjusted
        )
    else:
        # ---- If only one group, it gets 100% of the proportions
        adjusted_props[first_group] = group_proportions.sel(group=first_group)

    # Combine into a single DataArray with 'group' dimension
    return xr.concat(
        [adjusted_props[group].expand_dims({"group": [group]}) for group in group_keys], dim="group"
    )


def calculate_grouped_weights(
    length_weight_data: xr.DataArray,
    within_group_proportions: xr.DataArray,
    within_group_proportions_all: xr.DataArray,
    aggregate_proportions: xr.DataArray,
    adjusted_proportions: xr.DataArray,
    group_keys: List[str],
) -> xr.DataArray:
    """
    Calculate average weights by sex category across groups.

    Parameters
    ----------
    length_weight_data : xr.DataArray
        DataArray with weight values distributed across length bins and optionally other coordinates
    within_group_proportions : xr.DataArray
        Within-group proportions distributed over length and optionally other coordinates
    within_group_proportions_all : pd.DataFrame
        Across-group/overall proportions distributed over length and optionally other coordinates
    aggregate_proportions : xr.DataArray
        Table with aggregate proportions by group
    adjusted_proportions : xr.DataArray
        Adjusted proportions for each group and groupings
    group_keys : List[str]
        List of group keys/names (e.g., ["aged", "unaged"])

    Returns
    -------
    xr.DataArray
        DataArray with average weights by a grouping category defined by groupings

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
    weight_das = {}

    # Identify non-length dimensions (grouping dimensions)
    non_length_dims = [d for d in length_weight_data.dims if d != "length_bin"]

    # Get grouping categories from the first non-length dimension
    if not non_length_dims:
        raise ValueError("'length_weight_data' must have at least one dimension not 'length_bin'.")
    elif len(non_length_dims) > 1:
        # ---- Current dimensions
        all_dims = ", ".join(f"'{c}'" for c in non_length_dims)
        raise ValueError(
            f"'length_weight_data' must only have one additional dimension besides 'length_bin'. "
            f"Got: {all_dims}."
        )

    # Weight for 'all' category combined
    grouping_dim = non_length_dims[0]
    grouping_vals = length_weight_data[grouping_dim].values
    if "all" in grouping_vals:
        weight_all_values = length_weight_data.sel({grouping_dim: "all"})
        weighted_sum = sum(
            (within_group_proportions_all.sel(group=g) * weight_all_values).sum(dim="length_bin")
            * aggregate_proportions.sel(group=g)
            for g in group_keys
        )
        weight_das["all"] = weighted_sum

    # Weight for each grouping category individually
    if grouping_dim in adjusted_proportions.dims:
        for grouping in grouping_vals:
            # ---- Skip if "all" to avoid double-processing
            if grouping == "all":
                continue
            if grouping in adjusted_proportions.coords[grouping_dim].values:
                # ---- Use 'all' weights if available, otherwise use category-specific weights
                if "all" in grouping_vals:
                    weight_base = length_weight_data.sel({grouping_dim: "all"})
                else:
                    weight_base = length_weight_data.sel({grouping_dim: grouping})
                # ---- Weighted sums
                weighted_sum = sum(
                    (
                        within_group_proportions.sel(group=g, **{grouping_dim: grouping})
                        * weight_base
                    ).sum(dim="length_bin")
                    * adjusted_proportions.sel(group=g, **{grouping_dim: grouping})
                    for g in group_keys
                    if grouping in within_group_proportions.sel(group=g).coords[grouping_dim].values
                )
                weight_das[grouping] = weighted_sum

    # Combine into final DataArray
    if "all" in grouping_vals:
        ordered_cats = ["all"] + [c for c in grouping_vals if c != "all"]
    else:
        ordered_cats = list(grouping_vals)

    # Format and return
    result = xr.concat(
        [weight_das[cat] for cat in ordered_cats if cat in weight_das],
        dim=xr.DataArray(
            [cat for cat in ordered_cats if cat in weight_das],
            dims=[grouping_dim],
            name=grouping_dim,
        ),
    ).fillna(0.0)
    return result


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
    number_proportions: Dict[str, xr.Dataset],
    length_weight_data: xr.DataArray,
    group_columns: List[str] = [],
) -> xr.DataArray:
    """
    Calculate stratum- and group-specific average weights using xarray inputs.

    This function combines length-at-age or length-only distributions from multiple groups
    (e.g., aged and unaged fish) to produce weighted average weights by group and stratum,
    using xarray objects for all calculations. It accounts for the proportions of each group
    and correctly weights their contributions to the final average.

    Parameters
    ----------
    number_proportions : Dict[str, xr.Dataset]
        Dictionary mapping group names to xarray Datasets containing proportion data.
        Each Dataset must have at least one variable representing proportions and
        coordinates for all relevant grouping and binning dimensions.

    length_weight_data : xr.DataArray
        xarray DataArray containing fitted weights, with coordinates matching the
        grouping and binning dimensions used in the proportion Datasets.

    group_columns : List[str], optional (default is [])
        List of dimension names to use for grouping (e.g., strata, sex, etc.).
        If None, all non-binning dimensions present in the data will be used.

    Returns
    -------
    xr.DataArray
        DataArray with average weights for each stratum and group, indexed by group_columns and
        groupings.

    Notes
    -----

    - The function assumes that the length bins in number_proportions and length_weight_data match.

    - The function requires an 'all' group category in length_weight_data for calculating combined
      weights.

    - Missing strata or group categories will be excluded from the final results.

    - All calculations are performed using xarray and pandas interconversion for flexibility.

    Examples
    --------
    >>> # Prepare xarray Datasets for proportions (see proportions.number_proportions)
    >>> number_proportions = {'aged': ds_aged, 'unaged': ds_unaged}
    >>> # Prepare xarray DataArray for weights (see biology.length_binned_weights)
    >>> length_weight_data = da_binned_weights
    >>> # Calculate average weights
    >>> weights = stratum_averaged_weight(number_proportions, length_weight_data,
    group_columns=['stratum_num', 'sex'])
    >>> print(weights)
    <xarray.DataArray ...>
    """

    # Get variable names for number proportions
    proportion_vars = list(number_proportions.keys())

    # Identify non-length dimensions from length_weight_data
    grps = [d for d in length_weight_data.dims if d != "length_bin"]

    # Reorganize the individual DataArrays
    das_list = [ds["proportion"] for ds in number_proportions.values()]
    das_aligned = xr.align(*das_list, join="inner")

    # Find shared dimensions
    shared_dims = set().union(*map(lambda da: da.dims, das_list))

    # Extract proportion_overall and add group dimension
    proportion_arrays = []
    for group_name, ds in number_proportions.items():
        # ---- Get the grouped array
        grouped_array = ds["proportion_overall"].sum(
            dim=[d for d in ds.dims if d not in group_columns]
        )
        # ---- Expand the dimensions
        proportion_arrays.append(grouped_array.expand_dims({"group": [group_name]}))
    aggregate_proportions = xr.concat(proportion_arrays, dim="group")

    # Calculate the within-group proportions
    within_grp_props_list = []
    for da in das_aligned:
        sum_dims = [d for d in da.dims if d not in shared_dims]
        within_grp_props_list.append(da.sum(dim=sum_dims))
    within_grp_props_orig = xr.concat(
        within_grp_props_list,
        dim=xr.IndexVariable("group", list(number_proportions.keys())),
        join="outer",
    ).sum(dim=shared_dims - set([*["length_bin"], *group_columns, *grps]))
    # ---- Get the grouped aggregates for normalizing
    within_grp_props_norm = xr.concat(
        within_grp_props_list,
        dim=xr.IndexVariable("group", list(number_proportions.keys())),
    ).sum(dim=shared_dims - set([*group_columns, *grps]))
    within_grp_props = within_grp_props_orig / within_grp_props_norm

    # Generalize the overall groups
    within_grp_props_all = within_grp_props.sum(dim=grps)

    # Create list of arrays for overall proportions
    das_list_overall = [ds["proportion_overall"] for ds in number_proportions.values()]
    das_aligned_overall = xr.align(*das_list_overall, join="inner")

    # Compute the grouped overall proportions
    within_grp_props_overall_list = []
    for da in das_aligned_overall:
        sum_dims = [d for d in da.dims if d not in shared_dims]
        within_grp_props_overall_list.append(da.sum(dim=sum_dims))
    within_grp_props_overall = xr.concat(
        within_grp_props_overall_list,
        dim=xr.IndexVariable("group", list(number_proportions.keys())),
        join="outer",
    ).sum(dim=shared_dims - set([*group_columns, *grps]))

    # Compute the re-weighted proportions from the mixture of the different groups
    adjusted_proportions = calculate_adjusted_proportions(
        proportion_vars, aggregate_proportions, within_grp_props_overall
    )

    # Calculate final weights
    fitted_weights = calculate_grouped_weights(
        length_weight_data,
        within_grp_props,
        within_grp_props_all,
        aggregate_proportions,
        adjusted_proportions,
        proportion_vars,
    )

    return fitted_weights


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


def weight_proportions(
    weight_data: xr.DataArray,
    catch_data: pd.DataFrame,
    group_columns: str = [],
) -> xr.DataArray:
    """
    Calculate stratified weight proportions using xarray and pandas inputs.

    This function computes the proportional weight distribution for each group
    relative to the total weights across all strata or groupings, using an xarray
    DataArray for biological weights and a pandas DataFrame for catch weights.
    All grouping is handled dynamically based on the provided group_columns.

    Parameters
    ----------
    weight_data : xr.DataArray
        xarray DataArray containing weight distributions for all groups, with
        dimensions including those in group_columns.
    catch_data : pd.DataFrame
        DataFrame with catch data, including columns for group_columns and "weight".
    group_columns : list of str
        List of dimension/column names to group by (e.g., strata, sex, etc.).

    Returns
    -------
    xr.DataArray
        DataArray containing the calculated weight proportions, indexed by all
        grouping dimensions present in the input data.

    Notes
    -----
    - No column or dimension names are hard-coded; all logic is dynamic.
    - The function assumes that grouping columns are consistent between weight_data
      and catch_data.
    - Missing or extra group/category combinations are handled automatically by xarray/pandas.

    Examples
    --------
    >>> weight_props = weight_proportions(
    ...     weight_data=ds_da_weight_dist["aged"],
    ...     catch_data=dict_df_bio["catch"],
    ...     group_columns=["stratum_ks"]
    ... )
    >>> print(weight_props)
    <xarray.DataArray ...>
    """

    # Compute the grouped catch weights
    catch_weights = xr.DataArray(catch_data.groupby(group_columns)["weight"].sum())

    # Compute the total weights per stratum from the biological data
    group_weights = weight_data.sum(dim=[d for d in weight_data.dims if d not in group_columns])

    # Get the overall sums
    total_weights = catch_weights + group_weights

    # Compute the weight proportions for the array
    arr = weight_data / total_weights
    arr.name = "proportion_overall"
    return arr.to_dataset()


def fitted_weight_proportions(
    weight_data: xr.DataArray,
    reference_weight_proportions: xr.Dataset,
    catch_data: pd.DataFrame,
    number_proportions: xr.DataArray,
    binned_weights: xr.DataArray,
    stratum_dim: List[str] = [],
) -> xr.Dataset:
    """
    Calculate weight proportions based on fitted weights with adjustments based on reference
    values.

    This function computes comprehensive weight proportions for a group, accounting for reference
    proportions, length-binned proportions, and fitted weight tables. All grouping and dimension
    logic is handled dynamically based on the provided xarray objects and stratum_dim.

    Parameters
    ----------
    scaled_weight_data : xr.DataArray
        xarray DataArray of scaled weights for the group, with dimensions including stratum_dim.
    reference_weight_proportions : xr.Dataset
        xarray Dataset of reference weight proportions for comparison, with compatible dimensions.
    catch_data : pd.DataFrame
        DataFrame containing catch data, including columns for stratum_dim and "weight".
    number_proportions : xr.DataArray
        xarray DataArray of number proportions by relevant grouping factors.
    binned_weights : xr.DataArray
        xarray DataArray of fitted weights by length bins (and possibly other groupings).
    stratum_dim : list of str
        Stratification dimension name (e.g., ['stratum_ks']).

    Returns
    -------
    xr.DataArray
        DataArray containing the calculated detailed weight proportions, indexed by all grouping
        and binning dimensions present in the input data.

    Examples
    --------
    >>> props = fitted_weight_proportions(
    ...     scaled_weight_data=da_scaled_unaged_weights,
    ...     reference_weight_proportions=ds_da_weight_proportion["aged"],
    ...     catch_data=dict_df_bio["catch"],
    ...     number_proportions=dict_ds_number_proportion["unaged"],
    ...     binned_weights=da_binned_weights_all,
    ...     stratum_dim=["stratum_ks"]
    ... )
    >>> print(props)
    <xarray.DataArray ...>
    """

    # Compute the total strata weights
    subgroup_weights = weight_data.groupby(stratum_dim).sum(dim="length_bin").astype(float)

    # Calculate the total grouped weights
    group_weights = xr.DataArray(
        subgroup_weights.to_numpy().sum(axis=1),
        dims=stratum_dim,
        coords={col: weight_data.coords[col] for col in stratum_dim},
    )

    # Rescale the catch weights
    scaled_weights = (subgroup_weights / group_weights) * catch_data.groupby(stratum_dim)[
        "weight"
    ].sum().to_xarray().fillna(0.0)

    # Calculate the grouped weight proportions
    weight_grouped_props = (
        (
            scaled_weights
            / scaled_weights.sum(dim=[d for d in scaled_weights.coords if d not in stratum_dim])
        )
        .fillna(0.0)
        .squeeze()
    )

    # Calculate the number proportions based on strata
    length_bin_props = number_proportions["proportion"].sum(
        dim=[d for d in number_proportions.coords if d not in stratum_dim + ["length_bin"]]
    )

    # Calculate the average weights per length bin within each stratum
    fitted_length_weights = (length_bin_props * binned_weights).fillna(0.0)

    # Calculate the weight proportions
    # ---- Get the grouped total fitted weights
    total_fitted_weights = fitted_length_weights.sum(
        dim=[d for d in fitted_length_weights.coords if d not in stratum_dim]
    )
    # ---- Compute proportions
    fitted_weight_props = (fitted_length_weights / total_fitted_weights).fillna(0.0).squeeze()

    # Get the overall weight proportions per group for the reference data
    grouped_reference_proportions = reference_weight_proportions.sum(
        dim=[d for d in reference_weight_proportions.dims if d not in stratum_dim]
    )

    # Calculate the complementary proportions
    compl_data_props = 1 - grouped_reference_proportions

    # Distribute the within-group proportions to get the overall proportions
    return compl_data_props * fitted_weight_props * weight_grouped_props


def get_nasc_proportions_slice(
    number_proportions: xr.DataArray,
    ts_length_regression_parameters: Dict[str, float],
    include_filter: Dict[str, Any] = {},
    exclude_filter: Dict[str, Any] = {},
    group_columns: List[str] = [],
) -> xr.DataArray:
    """
    Extract and aggregate NASC proportions for specified groups from an xarray.DataArray, with
    dynamic filtering.

    This function selects and aggregates NASC proportions from an xarray.DataArray according to
    user-specified grouping columns and filters, and returns the result as an xarray.DataArray.
    All grouping and filtering logic is handled dynamically based on the provided arguments.

    Parameters
    ----------
    number_proportions : xr.DataArray
        xarray.DataArray of number proportions by relevant grouping factors, with dimensions
        including those in group_columns and binning columns (e.g., 'sex').
    exclude_filter : Dict[str, Any], optional
        Dictionary specifying values to exclude for filtering.
    include_filter : Dict[str, Any], optional
        Dictionary specifying values to include for filtering.
    ts_length_regression_parameters : Dict[str, float]
        Target strength-length regression parameters:
        - slope: regression slope
        - intercept: regression intercept
    group_columns : list of str, optional
        List of dimension names to use for grouping (e.g., ['stratum_num']).

    Returns
    -------
    xr.DataArray
        DataArray containing the aggregated NASC proportions for each group/stratum,
        indexed by group_columns.

    Notes
    -----
    - No dimension names are hard-coded; all grouping and filtering is dynamic.

    - The function assumes that grouping and binning dimensions are consistent with the DataArray
      structure.

    - Filtering is applied dynamically using include_filter and exclude_filter.

    - Output is an xarray.DataArray indexed by group_columns.

    Examples
    --------
    >>> props = get_nasc_proportions_slice(
    ...     number_proportions=aged_data,
    ...     group_columns=["stratum_ks"],
    ...     include_filter={"age_bin": [1]}
    ... )
    >>> print(props)
    <xarray.DataArray ...>
    """

    # Get length values from length bins
    length_vals = number_proportions["length_bin"].to_series().map(lambda x: x.mid).astype(float)

    # Compute equivalent sigma_bs using target strength-length regression
    sigma_bs_equiv = xr.DataArray(
        10
        ** (
            (
                ts_length_regression_parameters["slope"] * np.log10(length_vals)
                + ts_length_regression_parameters["intercept"]
            )
            / 10.0
        ),
        dims=["length_bin"],
        coords={"length_bin": length_vals.index},
    )

    # Get length-binned proportions aggregated based on group_columns
    aggregate_array = number_proportions["proportion"].sum(
        dim=[d for d in number_proportions.coords if d not in [*group_columns, "length_bin"]]
    )

    # Calculate total weighted sigma_bs for all strata
    aggregate_weighted_sigma_bs = (aggregate_array * sigma_bs_equiv.T).sum(
        dim=[d for d in aggregate_array.coords if d not in group_columns]
    )

    # Apply filters, if any
    filtered_proportions = number_proportions["proportion"]
    if include_filter:
        filtered_proportions = filtered_proportions.sel(include_filter)
    if exclude_filter:
        filtered_proportions = filtered_proportions.drop_sel(exclude_filter)
    # ---- Drop any singleton coordinates
    grouped_props = filtered_proportions.squeeze(drop=True)

    # Weight sigma_bs
    weighted_sigma_bs = (grouped_props * sigma_bs_equiv.T).sum(
        dim=[d for d in grouped_props.coords if d not in group_columns]
    )

    # Calculate the proportions
    return weighted_sigma_bs / aggregate_weighted_sigma_bs


def get_number_proportions_slice(
    number_proportions: xr.DataArray,
    stratum_dim: List[str] = [],
    exclude_filter: Dict[str, Any] = {},
    include_filter: Dict[str, Any] = {},
) -> xr.DataArray:
    """
    Extract and aggregate number proportions for specified groups from an xarray.DataArray, with
    dynamic filtering.

    This function selects and aggregates number proportions from an xarray.DataArray according to
    user-specified grouping columns and filters, and returns the result as an xarray.DataArray.
    All grouping and filtering logic is handled dynamically based on the provided arguments.

    Parameters
    ----------
    number_proportions : xr.DataArray
        xarray.DataArray of number proportions by relevant grouping factors, with dimensions
        including those in group_columns and binning columns (e.g., 'length_bin').
    stratum_dim : list of str
        Stratification dimension name (e.g., ['stratum_ks']).
    exclude_filter : Dict[str, Any], optional
        Dictionary specifying values to exclude for filtering.
    include_filter : Dict[str, Any], optional
        Dictionary specifying values to include for filtering.

    Returns
    -------
    xr.DataArray
        DataArray containing the aggregated number proportions for each group/stratum,
        indexed by group_columns.

    Notes
    -----

    - No dimension names are hard-coded; all grouping and filtering is dynamic.

    - The function assumes that grouping and binning dimensions are consistent with the DataArray
      structure.

    - Filtering is applied dynamically using include_filter and exclude_filter.

    - Output is an xarray.DataArray indexed by group_columns.

    Examples
    --------
    >>> props = get_number_proportions_slice(
    ...     number_proportions=dict_ds_number_proportion["aged"],
    ...     group_columns=["stratum_ks"],
    ...     include_filter={"age_bin": [1]}
    ... )
    >>> print(props)
    <xarray.DataArray ...>
    """

    # Determine index columns from filter keys and length_bin
    filter_indices = set(list(exclude_filter.keys()) + list(include_filter.keys()) + ["length_bin"])

    # Intersect with actually available columns/indices
    index_set = set(number_proportions.coords).intersection(filter_indices)

    # Get length-binned proportions aggregated based on strata
    aggregate_array = number_proportions["proportion"].sum(
        dim=[d for d in number_proportions.coords if d not in [*stratum_dim, *index_set]]
    )

    # Apply filters, if any
    if include_filter:
        # ---- Parse existing labels
        to_keep = {k: v for k, v in include_filter.items() if k in aggregate_array.coords}
        if to_keep:
            aggregate_array = aggregate_array.sel(to_keep)
    if exclude_filter:
        # ---- Parse existing labels
        to_drop = {k: v for k, v in exclude_filter.items() if k in aggregate_array.coords}
        if to_drop:
            aggregate_array = aggregate_array.drop_sel(to_drop)
    # ---- Drop any singleton coordinates
    grouped_props = aggregate_array.squeeze(drop=True)

    # Aggregate further over the defined strata
    return grouped_props.sum(dim=[d for d in grouped_props.coords if d not in stratum_dim])


def get_weight_proportions_slice(
    weight_proportions: xr.Dataset,
    stratum_dim: List[str] = [],
    include_filter: Dict[str, Any] = {},
    exclude_filter: Dict[str, Any] = {},
    number_proportions: Dict[str, xr.Dataset] = {},
    length_threshold_min: Optional[float] = None,
    weight_proportion_threshold: float = 1e-10,
) -> xr.DataArray:
    """
    Calculate weight proportions for a specific slice of the population using xarray, with optional
    thresholding.

    This function computes weight proportions for a target population group from an
    xarray.Dataset, with optional thresholding based on number proportions. Thresholding helps
    identify and mask unreliable estimates due to low sample sizes or small proportions. Number
    proportions must be provided as a dictionary of xarray.Dataset objects, each containing a
    "proportion" variable.

    Parameters
    ----------
    weight_proportions : xr.Dataset
        Dataset with weight proportions and all relevant grouping dimensions.
    stratum_dim : list of str
        Stratification dimension name (e.g., ['stratum_ks']).
    include_filter : Dict[str, Any], default {}
        Dictionary of filter criteria for the target group (e.g., {"age_bin": [1]}).
    exclude_filter : Dict[str, Any], default {}
        Dictionary of groups to exclude (e.g., {"length_bin": [small_lengths]}).
    number_proportions : Dict[str, xr.Dataset], default {}
        Dictionary of xarray.Dataset objects for thresholding. Each Dataset must contain a
        "proportion" variable and relevant grouping/binning dimensions.
    length_threshold_min : float or None, default None
        Minimum length for threshold calculations (e.g., 10.0 cm). Used only if number_proportions
        are provided. If None, no length-based thresholding is applied.
    weight_proportion_threshold : float, default 1e-10
        Threshold value for masking small proportions.

    Returns
    -------
    xr.DataArray
        Weight proportions by group_columns for the target group, with unreliable estimates masked
        to 0.

    Notes
    -----
    Thresholding Logic:

    1. Calculate target group weight proportions.
    2. If number_proportions are provided, calculate corresponding number thresholds.
    3. If length_threshold_min is set, exclude length bins below this threshold.
    4. Apply a mask where both weight and number proportions are very small.
    5. Set masked values to 0.0 to indicate unreliable estimates.

    Examples
    --------
    >>> weight_props = get_weight_proportions_slice(
    ...     weight_proportions,
    ...     group_columns=["stratum_ks"],
    >>> weight_props = get_weight_proportions_slice(
    ...     weight_proportions,
    ...     group_columns=["stratum_ks"],
    ...     include_filter={"age_bin": [1]},
    ...     number_proportions={"aged": ds_aged, "unaged": ds_unaged},
    ...     length_threshold_min=10.0
    ... )
    >>> print(weight_props)
    <xarray.DataArray ...>
    """

    # Determine index columns from filter keys and length_bin
    index_set = set(list(exclude_filter.keys()) + list(include_filter.keys()) + ["length_bin"])

    # Get length-binned proportions aggregated based on strata
    aggregate_array = weight_proportions["proportion_overall"].sum(
        dim=[d for d in weight_proportions.coords if d not in [*stratum_dim, *index_set]]
    )

    # Apply filters, if any
    if include_filter:
        aggregate_array = aggregate_array.sel(include_filter)
    if exclude_filter:
        aggregate_array = aggregate_array.drop_sel(exclude_filter)
    # ---- Drop any singleton coordinates
    grouped_props = aggregate_array.squeeze(drop=True)

    # Aggregate target stratum proportions
    target_group_weight_proportions = grouped_props.sum(
        dim=[d for d in grouped_props.dims if d not in stratum_dim], skipna=True
    )

    # Normalize by the unfiltered proportions
    total_weight_proportions = weight_proportions["proportion_overall"].sum(
        dim=[d for d in weight_proportions.dims if d not in stratum_dim], skipna=True
    )
    proportions_weight = target_group_weight_proportions / total_weight_proportions

    # Apply thresholding if number proportions provided
    if len(number_proportions) > 0:
        # ---- Optional length-based threshold filter
        if length_threshold_min:
            # ---- Get all unique length values across datasets
            all_length_vals = np.unique(
                np.concatenate(
                    [
                        [bin.mid for bin in ds["length_bin"].values]
                        for ds in number_proportions.values()
                    ]
                )
            )
            # ---- Create length exclusion filter
            length_exclusion_filter = {
                "length_bin": all_length_vals[all_length_vals < length_threshold_min],
                **exclude_filter,
            }
            # ---- Get filtered number proportions for each dataset
            filtered_number_proportions_grp = {
                key: get_number_proportions_slice(
                    ds,
                    stratum_dim=stratum_dim + ["length_bin"],
                    exclude_filter=length_exclusion_filter,
                    include_filter=include_filter,
                )
                for key, ds in number_proportions.items()
            }
            # ---- Calculate element-wise product across all arrays
            filtered_number_proportions = functools.reduce(
                operator.mul, filtered_number_proportions_grp.values()
            ).sum(
                dim=[
                    d
                    for d in next(iter(filtered_number_proportions_grp.values())).dims
                    if d not in stratum_dim
                ]
            )
            # ---- Align/reindex
            proportions_weight = proportions_weight.reindex_like(
                filtered_number_proportions
            ).fillna(0.0)
            target_group_weight_proportions = target_group_weight_proportions.reindex_like(
                filtered_number_proportions
            ).fillna(0.0)
            # ---- Apply threshold mask
            threshold_mask_num = filtered_number_proportions <= weight_proportion_threshold
            # ---- Set masked values to 0
            proportions_weight = xr.where(threshold_mask_num.values, 0, proportions_weight)

    # Apply weight thresholding
    threshold_mask_wgts = target_group_weight_proportions <= weight_proportion_threshold
    proportions_weight = xr.where((threshold_mask_wgts.values), 0, proportions_weight)

    # Output the masked proportions array
    return proportions_weight
