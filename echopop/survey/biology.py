from typing import List

import numpy as np
import pandas as pd


def fit_length_weight_regression(data: pd.DataFrame) -> pd.Series:
    """
    Fit a log-linear length-weight regression to biological data.

    This function fits a linear regression to log10-transformed length and weight
    data, following the standard allometric relationship: log(weight) = a + b*log(length).
    The function can be used standalone or as part of a groupby operation.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing 'length' and 'weight' columns.

    Returns
    -------
    pd.Series
        Series with regression coefficients:
        - 'slope': The slope coefficient (b) from log(weight) = a + b*log(length)
        - 'intercept': The intercept coefficient (a) from log(weight) = a + b*log(length)

    Examples
    --------
    >>> # Standalone usage
    >>> coeffs = fit_length_weight_regression(length_weight_dataset)

    >>> # With groupby
    >>> grouped_coeffs = data.groupby('sex').apply(fit_length_weight_regression)
    """
    # Remove missing values
    clean_data = data.dropna(subset=["length", "weight"])

    # Fit the length-weight log-linear regression
    regression_coefficients = pd.Series(
        np.polyfit(np.log10(clean_data["length"]), np.log10(clean_data["weight"]), 1),
        index=["slope", "intercept"],
    )

    return regression_coefficients


def quantize_length_data(df, group_columns: List[str]):
    """
    Process DataFrame to ensure it has 'length' and 'length_count' columns.

    Aggregates fish length data by grouping variables and length, either counting occurrences (if
    no length_count exists) or summing existing counts.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing length data. Must have a 'length' column.
        Optionally can have 'length_count' column with existing counts.
    group_columns : List[str]
        List of column names to group by before aggregating length counts.
        Common examples: ['stratum', 'haul_num', 'sex']

    Returns
    -------
    pd.DataFrame
        DataFrame with grouping columns, 'length', and 'length_count' columns.
        Index will be a MultiIndex of group_columns + ['length'].

    Examples
    --------
    >>> # Data without existing counts - will count occurrences
    >>> specimen_df = pd.DataFrame({
    ...     'stratum': [1, 1, 1, 2, 2],
    ...     'sex': ['M', 'M', 'F', 'F', 'M'],
    ...     'length': [20.5, 20.5, 22.1, 18.3, 20.5]
    ... })
    >>> result = quantize_length_data(specimen_df, ['stratum', 'sex'])
    >>> print(result)
                           length_count
    stratum sex length
    1       F   22.1               1
            M   20.5               2
    2       F   18.3               1
            M   20.5               1

    >>> # Data with existing counts - will sum them
    >>> length_df = pd.DataFrame({
    ...     'stratum': [1, 1, 2],
    ...     'length': [20.5, 22.1, 20.5],
    ...     'length_count': [5, 3, 2]
    ... })
    >>> result = quantize_length_data(length_df, ['stratum'])
    >>> print(result)
                   length_count
    stratum length
    1       20.5              5
            22.1              3
    2       20.5              2

    Notes
    -----
    This function automatically detects whether to count fish (size operation) or sum existing
    counts based on the presence of a 'length_count' column.

    The resulting DataFrame will have a MultiIndex with group_columns + ['length'] and a single
    'length_count' column containing the aggregated counts.
    """

    # Create copy
    df = df.copy()

    # Define which column should be quantized
    if "length_count" not in df.columns:
        # ---- Column name
        sum_var_column = "length"
        # ---- Operation
        var_operation = "size"
    else:
        # ---- Column name
        sum_var_column = "length_count"
        # ---- Operation
        var_operation = "sum"

    # Aggregate and return
    return df.groupby(group_columns + ["length"]).agg(length_count=(sum_var_column, var_operation))
