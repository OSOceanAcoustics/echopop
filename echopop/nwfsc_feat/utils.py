from functools import reduce
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import interpolate as interp


def binned_distribution(
    bins: npt.NDArray[np.number], return_dataframe: bool = True
) -> Union[pd.DataFrame, Tuple[npt.NDArray[np.number], npt.NDArray[np.number]]]:
    """
    Create centered bins for data binning operations.

    This function takes an array of bin edges and creates centered bins by calculating
    the mean bin width and extending the bins to create proper intervals for binning.
    The centered bins can be used with pandas.cut() for data discretization.

    Parameters
    ----------
    bins : npt.NDArray[np.number]
        Array of bin edge values. Must be 1-dimensional and contain at least 2 elements.
        Values should be in ascending order for proper binning behavior.
    return_dataframe : bool, default True
        If True, returns a DataFrame with bins and intervals.
        If False, returns a tuple of (bins, centered_bins) arrays.
        [!!! NOTE: THIS HAS BEEN ADDED AS AN ARGUMENT FOR PRIMARILY TESTING PURPOSES BEFORE
        DETERMINING THE MOST APPROPRIATE OUTPUT]

    Returns
    -------
    pd.DataFrame or tuple
        If return_dataframe is True:
            DataFrame with columns:
            - 'bin': Original bin values
            - 'interval': pd.Interval objects representing the binning intervals
        If return_dataframe is False:
            Tuple containing:
            - bins: Original bin array
            - centered_bins: Array of centered bin edges for use with pd.cut()

    Raises
    ------
    ValueError
        If bins array has fewer than 2 elements or is not 1-dimensional.
    TypeError
        If bins is not a numpy array or cannot be converted to one.

    Examples
    --------
    >>> import numpy as np
    >>> bins = np.array([1, 2, 3, 4, 5])
    >>> result = binned_distribution(bins)
    >>> print(result.columns)
    Index(['bin', 'interval'], dtype='object')

    >>> bins = np.linspace(0, 10, 11)
    >>> bins_array, centered_array = binned_distribution(bins, return_dataframe=False)
    >>> len(centered_array) == len(bins_array) + 1
    True

    Notes
    -----
    The function calculates the bin width as the mean of half the differences between
    consecutive bin values. This approach works well for both evenly and unevenly
    spaced bins.

    The centered bins extend beyond the original range by one bin width on each side,
    ensuring that all original bin values fall within the created intervals.
    """

    # Compute binwidth as mean of half the differences
    binwidth = np.mean(np.diff(bins) / 2.0)

    # Create centered bins by extending the range
    centered_bins = np.concatenate([[bins[0] - binwidth], bins + binwidth])

    if return_dataframe:
        # Generate DataFrame with bins and intervals
        intervals = pd.cut(bins, centered_bins)
        return pd.DataFrame({"bin": bins, "interval": intervals})
    else:
        return bins, centered_bins


def binify(
    data: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    bins: npt.NDArray[np.number],
    bin_column: str,
) -> None:
    """
    Apply binning to biological data using predefined bin distributions.

    This function bins continuous variables (like length or age) in biological datasets
    using bin edge arrays. It creates interval distributions internally and can handle single
    DataFrames or dictionaries of DataFrames, automatically skipping DataFrames that
    don't contain the target column. The data is modified in place.

    Parameters
    ----------
    data : pd.DataFrame or dict of pd.DataFrame
        Target data to bin. Can be a single DataFrame or dictionary of DataFrames.
        Data will be modified in place.
    bins : npt.NDArray[np.number]
        Array of bin edge values. Must be 1-dimensional and contain at least 2 elements.
        Values should be in ascending order for proper binning behavior.
    bin_column : str
        Name of the column in data to apply binning to (e.g., 'length', 'age').

    Returns
    -------
    None
        Data is modified in place, nothing is returned.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from echopop.nwfsc_feat.utils import binify
    >>> 
    >>> # Create sample data
    >>> bio_data = pd.DataFrame({
    ...     'length': [25.5, 30.2, 35.8, 40.1, 45.3],
    ...     'weight': [150, 220, 310, 420, 580]
    ... })
    >>> 
    >>> # Create bin edges
    >>> length_bins = np.array([20, 30, 40, 50])
    >>> 
    >>> # Apply binning (modifies bio_data in place)
    >>> binify(bio_data, length_bins, 'length')
    >>> print('length_bin' in bio_data.columns)
    True
    >>> 
    >>> # Works with numpy linspace too
    >>> age_bins = np.linspace(start=1., stop=22., num=22)
    >>> bio_data['age'] = [5, 8, 12, 15, 18]
    >>> binify(bio_data, age_bins, 'age')
    >>> print('age_bin' in bio_data.columns)
    True
    >>> 
    >>> # Works with dictionary of DataFrames too
    >>> data_dict = {'catch': bio_data.copy(), 'length': bio_data.copy()}
    >>> binify(data_dict, length_bins, 'length')
    >>> print('length_bin' in data_dict['catch'].columns)
    True
    """
    # Create bin distribution internally using binned_distribution function
    bin_distribution = binned_distribution(bins, return_dataframe=True)
    
    # Extract bin categories
    try:
        bin_intervals = bin_distribution["interval"].cat.categories
    except AttributeError:
        raise ValueError("Failed to create proper interval categories from bins")

    # Format new bin column name
    bin_column_name = f"{bin_column}_bin"

    def _apply_binning(df: pd.DataFrame) -> None:
        """Apply binning to a single DataFrame in place, skip if column missing."""
        if bin_column not in df.columns:
            return  # Skip DataFrames without target column

        df[bin_column_name] = pd.cut(df[bin_column], bins=bin_intervals)

    # Handle different input types
    if isinstance(data, pd.DataFrame):
        # Single DataFrame
        _apply_binning(data)

    elif isinstance(data, dict):
        # Dictionary of DataFrames - modify each DataFrame in place
        for key, df in data.items():
            if isinstance(df, pd.DataFrame):
                _apply_binning(df)  # Automatically skips if column missing

    else:
        raise TypeError(f"data must be DataFrame or dict of DataFrames, got {type(data)}")


def apply_filters(
    df: pd.DataFrame,
    include_filter: Optional[Dict[str, Any]] = None,
    exclude_filter: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    Apply inclusion and exclusion filters to a DataFrame.

    Filters rows based on column values, supporting both inclusion and
    exclusion criteria with single values or lists of values.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to filter
    include_filter : Dict[str, Any], optional
        Dictionary of column:value(s) pairs. Rows will be kept if they match.
        If value is a list, rows matching any value in the list will be kept.
    exclude_filter : Dict[str, Any], optional
        Dictionary of column:value(s) pairs. Rows will be excluded if they match.
        If value is a list, rows matching any value in the list will be excluded.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame

    Examples
    --------
    >>> # Keep only females and males
    >>> apply_filters(df, include_filter={"sex": ["female", "male"]})

    >>> # Exclude unsexed specimens and small fish
    >>> apply_filters(df, exclude_filter={"sex": "unsexed", "length": 10})
    """

    # Initialize boolean mask
    mask = pd.Series(True, index=df.index)

    # Inclusion filter
    if include_filter:
        include_conditions = [
            df[k].isin(v) if isinstance(v, list) else df[k] == v for k, v in include_filter.items()
        ]
        mask &= reduce(lambda a, b: a & b, include_conditions)

    # Exclusion filter
    if exclude_filter:
        exclude_conditions = [
            ~df[k].isin(v) if isinstance(v, list) else df[k] != v for k, v in exclude_filter.items()
        ]
        mask &= reduce(lambda a, b: a & b, exclude_conditions)

    # Return masked DataFrame
    return df[mask]


def group_interpolator_creator(
    grouped_data: pd.DataFrame,
    independent_var: str,
    dependent_var: str,
    contrast_vars: Optional[Union[str, List[str]]] = None,
) -> Dict:
    """
    Create interpolator functions grouped by one or more contrast variables.

    Generates scipy interpolation functions for length-weight relationships,
    either as a single global interpolator or grouped by specified contrast
    variables (e.g., sex, age class).

    Parameters
    ----------
    grouped_data : pd.DataFrame
        Data containing independent and dependent variables
    independent_var : str
        Column name of the independent variable (e.g., 'length')
    dependent_var : str
        Column name of the dependent variable (e.g., 'weight_fitted')
    contrast_vars : str, List[str], or None, optional
        Column name(s) to group by (e.g., 'sex' or ['sex', 'age_bin']).
        If None or empty list, a single global interpolator is created.

    Returns
    -------
    Dict
        Dictionary of interpolator functions.
        When contrast_vars is provided, keys are contrast variable values.
        When contrast_vars is None or empty, contains a single entry with key '_global_'.

    Notes
    -----
    The interpolation is linear and extrapolates using the endpoints
    when values outside the range are requested.

    Requires at least 2 points per group for valid interpolation.
    """
    # Check if we have contrast variables to group by
    if contrast_vars is None or (isinstance(contrast_vars, list) and len(contrast_vars) == 0):
        # Create a single global interpolator
        if len(grouped_data) < 2:
            # Not enough points for interpolation
            return {"_global_": None}

        # Sort the data
        sorted_data = grouped_data.sort_values(by=independent_var)

        # Create the interpolator
        global_interp = interp.interp1d(
            sorted_data[independent_var],
            sorted_data[dependent_var],
            kind="linear",
            bounds_error=False,
            fill_value=(sorted_data[dependent_var].iloc[0], sorted_data[dependent_var].iloc[-1]),
        )

        # Return a dictionary with a special key for the global interpolator
        return {"_global_": global_interp}

    # Convert single string to list for consistent handling
    if isinstance(contrast_vars, str):
        contrast_vars = [contrast_vars]

    # Interpolator generation helper function
    def interpolator_factory(sub_group):
        if len(sub_group) < 2:
            # Need at least 2 points for interpolation
            return None

        # Sort the grouped values
        sub_group_sort = sub_group.sort_values(by=independent_var)

        # Return the interpolation object for the specific sub-group
        return interp.interp1d(
            sub_group_sort[independent_var],
            sub_group_sort[dependent_var],
            kind="linear",
            bounds_error=False,
            fill_value=(
                sub_group_sort[dependent_var].iloc[0],
                sub_group_sort[dependent_var].iloc[-1],
            ),
        )

    # Produce a dictionary of interpolator functions
    interpolators = (
        grouped_data.groupby(contrast_vars)
        .apply(interpolator_factory, include_groups=False)
        .to_dict()
    )

    return interpolators


def create_grouped_series(
    proportions_dict: Dict[str, pd.DataFrame], group_cols: List[str], value_col: str
) -> pd.DataFrame:
    """
    Create and combine grouped series from a dictionary of DataFrames.

    Parameters
    ----------
    proportions_dict : Dict[str, pd.DataFrame]
        Dictionary of DataFrames with proportion data, keyed by group name
    group_cols : List[str]
        Columns to group by (e.g., ["stratum_num", "sex"])
    value_col : str
        Column name containing values to be aggregated (e.g., "proportion_overall")

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with aggregated values and original dictionary keys as 'group' column

    Examples
    --------
    >>> props_dict = {'aged': aged_df, 'unaged': unaged_df}
    >>> grouped = create_grouped_series(
    ...    props_dict,
    ...    ["stratum_num", "sex"],
    ...    "proportion_overall"
    ... )
    """
    series = [
        df.groupby(group_cols, observed=False)[value_col].sum().reset_index().assign(group=key)
        for key, df in proportions_dict.items()
    ]
    return pd.concat(series, axis=0)


def create_pivot_table(
    df: pd.DataFrame,
    index_cols: List[str],
    strat_cols: List[str],
    value_col: str,
) -> pd.DataFrame:
    """
    Create a pivot table with standard settings.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to pivot
    index_cols : List[str]
        Column names to use as pivot table index
    strat_cols : List[str]
        Column names to use as pivot table columns
    value_col : str
        Column name containing values to aggregate

    Returns
    -------
    pd.DataFrame
        Pivot table with values aggregated by sum, converted to float, and NaN filled with 0

    Examples
    --------
    >>> pivot = create_pivot_table(
    ...     df=grouped_data,
    ...     index_cols=["group", "sex"],
    ...     strat_cols=["stratum_num"],
    ...     value_col="proportion_overall"
    ... )
    """
    return (
        df.pivot_table(
            index=index_cols, columns=strat_cols, values=value_col, aggfunc="sum", observed=False
        )
        .astype(float)
        .fillna(0.0)
    )


def create_grouped_table(
    proportions_dict: Dict[str, pd.DataFrame],
    group_cols: List[str],
    index_cols: List[str],
    strat_cols: List[str],
    value_col: str,
) -> pd.DataFrame:
    """
    Create grouped series and then convert to pivot table.

    Parameters
    ----------
    proportions_dict : Dict[str, pd.DataFrame]
        Dictionary of DataFrames with proportion data, keyed by group name
    group_cols : List[str]
        Columns to group by (e.g., ["stratum_num", "sex"])
    index_cols : List[str]
        Column names to use as pivot table index
    strat_cols : List[str]
        Column names to use as pivot table columns
    value_col : str
        Column name containing values to aggregate

    Returns
    -------
    pd.DataFrame
        Pivot table with aggregated values by group

    Examples
    --------
    >>> props_dict = {'aged': aged_df, 'unaged': unaged_df}
    >>> grouped_table = create_grouped_table(
    ...     props_dict,
    ...     group_cols=["stratum_num", "sex"],
    ...     index_cols=["group"],
    ...     strat_cols=["stratum_num"],
    ...     value_col="proportion_overall"
    ... )
    """
    # First create grouped series
    grouped_df = create_grouped_series(proportions_dict, group_cols, value_col)

    # Then convert to pivot table
    return create_pivot_table(grouped_df, index_cols, strat_cols, value_col)
