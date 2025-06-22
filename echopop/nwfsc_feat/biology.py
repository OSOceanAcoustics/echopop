from typing import Union

import numpy as np
import pandas as pd

from .utils import binned_distribution, create_pivot_table


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


def length_binned_weights(
    data: pd.DataFrame,
    length_bins: np.ndarray,
    regression_coefficients: Union[pd.Series, pd.DataFrame],
    impute_bins: bool = True,
    minimum_count_threshold: int = 0,
) -> pd.DataFrame:
    """
    Compute length-binned average weights as a pivot table.

    This function calculates fitted weights for length bins by combining modeled weights
    (from length-weight regression) with observed mean weights. For bins with sufficient
    sample sizes, observed means are used; for bins with low sample sizes, modeled weights
    are used if imputation is enabled. The result is always returned as a wide-format
    pivot table.

    Parameters
    ----------
    data : pd.DataFrame
        Specimen data containing 'length', 'weight', and 'length_bin' columns.
        The 'length_bin' column must already exist in the data.
    length_bins : np.ndarray
        Array of bin edge values for length binning. Used to create a complete
        bin distribution that ensures all bins are represented.
    regression_coefficients : pd.Series or pd.DataFrame
        Length-weight regression coefficients from fit_length_weight_regression().
        If DataFrame, represents grouped coefficients (e.g., by sex).
    impute_bins : bool, default True
        Whether to use modeled weights for bins with insufficient data.
        If False, only observed means are used regardless of sample size.
    minimum_count_threshold : int, default 0
        Minimum number of specimens required to use observed mean instead of modeled weight.
        Only relevant when impute_bins=True.

    Returns
    -------
    pd.DataFrame
        Pivot table with length bins as index and grouping variables as columns.
        If no grouping variable is present, column will be named "all".
        Values are fitted weights for each bin-group combination.

    Examples
    --------
    >>> # Single coefficient set - returns pivot table with "all" column
    >>> coeffs = fit_length_weight_regression(specimen_data)
    >>> length_bins = np.array([10, 15, 20, 25, 30])
    >>> fitted = length_binned_weights(specimen_data, length_bins, coeffs)
    >>> # fitted.columns == ["all"]

    >>> # Grouped coefficients (e.g., by sex) - returns pivot table with sex columns
    >>> sex_coeffs = specimen_data.groupby('sex').apply(fit_length_weight_regression)
    >>> fitted = length_binned_weights(specimen_data, length_bins, sex_coeffs,
    ...                               minimum_count_threshold=5)
    >>> # fitted.columns might be ["female", "male"]
    
    >>> # No imputation - use only observed means
    >>> fitted = length_binned_weights(specimen_data, length_bins, coeffs,
    ...                               impute_bins=False)
    """

    # Make a copy to avoid modifying original data
    data = data.copy()

    # Create length distribution from bins
    length_distribution = binned_distribution(length_bins)

    # Handle different coefficient input types
    if isinstance(regression_coefficients, pd.Series):
        # Single set of coefficients - convert to DataFrame for consistent processing
        regression_df = pd.DataFrame([regression_coefficients])
        # Reset index to avoid grouping complications
        regression_df.reset_index(drop=True, inplace=True)
    else:
        # Already a DataFrame from groupby operation
        regression_df = regression_coefficients.reset_index()

    # Initialize fitted weights dataframe
    weight_fitted_df = length_distribution.copy()

    # Cross merge with regression coefficients
    weight_fitted_df = weight_fitted_df.merge(regression_df, how="cross").rename(
        columns={"interval": "length_bin"}
    )

    # Predict weight per bin using allometric relationship: weight = 10^intercept * length^slope
    weight_fitted_df["weight_modeled"] = (
        10.0 ** weight_fitted_df["intercept"] * weight_fitted_df["bin"] ** weight_fitted_df["slope"]
    )  # Get the column names if any grouping is required
    cols = [name for name in regression_coefficients.index.names if name is not None] + [
        "length_bin"
    ]

    # Quantize weight counts per length bin
    binned_weight_distribution = (
        data.groupby(cols, observed=False)["length"].size().fillna(0).to_frame("count")
    )

    # Quantize weight means per length bin
    binned_weight_distribution["weight_mean"] = (
        data.groupby(cols, observed=False)["weight"].mean().fillna(0)
    )

    # Merge with the fitted weights
    binned_weight_distribution["weight_modeled"] = weight_fitted_df.set_index(cols)[
        "weight_modeled"
    ]

    # Create distribution mask based on imputation settings
    if impute_bins:
        # Use modeled weights when count is below threshold
        if minimum_count_threshold > 0:
            distribution_mask = binned_weight_distribution["count"] < minimum_count_threshold
        else:
            # When threshold is 0, use modeled weights for empty bins only
            distribution_mask = binned_weight_distribution["count"] == 0
    else:
        # No imputation - always use observed means (mask is all False)
        distribution_mask = pd.Series(False, index=binned_weight_distribution.index)    # Apply mask to determine final fitted weights
    binned_weight_distribution["weight_fitted"] = np.where(
        distribution_mask,
        binned_weight_distribution["weight_modeled"],
        binned_weight_distribution["weight_mean"],
    )

    # Reset the index to get grouping columns as regular columns
    long_format_df = binned_weight_distribution.reset_index()
    
    # Determine grouping columns and create pivot table
    grouping_cols = [name for name in regression_coefficients.index.names if name is not None]
    
    if grouping_cols:
        # Create pivot table with grouping columns as pivot columns
        pivot_df = create_pivot_table(
            df=long_format_df,
            index_cols=["length_bin"],
            strat_cols=grouping_cols,
            value_col="weight_fitted"
        )
    else:
        # No grouping variables - create pivot table with "all" column
        long_format_df["all"] = "all"
        pivot_df = create_pivot_table(
            df=long_format_df,
            index_cols=["length_bin"],
            strat_cols=["all"],
            value_col="weight_fitted"
        )
    
    return pivot_df
