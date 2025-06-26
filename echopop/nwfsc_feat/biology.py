from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .utils import binned_distribution


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
    Compute length-binned average weights using regression coefficients and observed data.

    This function calculates fitted weights for length bins by combining modeled weights
    (from length-weight regression) with observed mean weights. For bins with sufficient
    sample sizes, observed means are used; for bins with low sample sizes, modeled weights
    are used if imputation is enabled.

    Parameters
    ----------
    data : pd.DataFrame
        Specimen data containing 'length', 'weight', and 'length_bin' columns.
        The 'length_bin' column must already exist in the data.
    length_distribution : pd.DataFrame
        DataFrame with length bin information, containing 'bin' and 'interval' columns.
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
        DataFrame with fitted weights for each length bin and group combination.
        Contains grouping columns (if any), 'length_bin', and 'weight_fitted'.

    Examples
    --------
    >>> # Single coefficient set
    >>> coeffs = fit_length_weight_regression(specimen_data)
    >>> fitted = compute_binned_weights(specimen_data, length_dist, coeffs)

    >>> # Grouped coefficients (e.g., by sex)
    >>> sex_coeffs = specimen_data.groupby('sex').apply(fit_length_weight_regression)
    >>> fitted = compute_binned_weights(specimen_data, length_dist, sex_coeffs,
    ...                                minimum_count_threshold=5)

    >>> # No imputation - use only observed means
    >>> fitted = compute_binned_weights(specimen_data, length_dist, coeffs,
    ...                                impute_bins=False)
    """
    # Make a copy to avoid modifying original data
    data = data.copy()  # Create length distribution from bins
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
        distribution_mask = pd.Series(False, index=binned_weight_distribution.index)

    # Apply mask to determine final fitted weights
    binned_weight_distribution["weight_fitted"] = np.where(
        distribution_mask,
        binned_weight_distribution["weight_modeled"],
        binned_weight_distribution["weight_mean"],
    )

    # Reset the index and pare down the output columns
    return binned_weight_distribution.reset_index().filter(cols + ["weight_fitted"])


def set_population_metrics(
    df_nasc: pd.DataFrame,
    metrics: List[str] = ["abundance", "biomass", "biomass_density"],
    stratify_by: str = "stratum_num",
    df_average_weight: Optional[Union[pd.DataFrame, float]] = None,
) -> None:
    """
    Convert acoustically derived number densities into population metrics (in-place)

    This function calculates abundance, biomass, and biomass density from acoustic
    number density estimates. The calculations are performed in-place, adding new
    columns to the input DataFrame. For biomass calculations, average weights are
    required and can be provided as a scalar value or DataFrame with group-specific weights.

    Parameters
    ----------
    df_nasc : pd.DataFrame
        DataFrame containing acoustic data. Must include:
        - 'number_density': Number of animals per unit area (animals/nmi²)
        - 'area_interval': Area of each sampling interval (nmi²)
        Additional columns depend on stratify_by specification.
    metrics : List[str], default ["abundance", "biomass", "biomass_density"]
        Population metrics to calculate. Options:
        - 'abundance': Total number of animals in each interval
        - 'biomass': Total biomass in each interval (requires df_average_weight)
        - 'biomass_density': Biomass per unit area (requires df_average_weight)
    stratify_by : str, default "stratum_num"
        Column name to use for stratification when aligning with df_average_weight.
        If provided and df_average_weight is a DataFrame, df_nasc will be indexed
        by this column for weight alignment.
    df_average_weight : Optional[Union[pd.DataFrame, float]], default None
        Average weight data for biomass calculations:
        - float: Single weight value applied to all observations
        - pd.DataFrame: Group-specific weights. If stratify_by is specified,
          this DataFrame should be indexed or contain the stratification column.
        Required when 'biomass' or 'biomass_density' are in metrics.

    Returns
    -------
    None
        Function modifies df_nasc in-place by adding the requested metric columns.

    Raises
    ------
    KeyError
        If required columns ('number_density', 'area_interval') are missing from df_nasc.
    ValueError
        If biomass metrics are requested but df_average_weight is not provided.

    Examples
    --------
    >>> # Basic usage with scalar average weight
    >>> nasc_df = pd.DataFrame({
    ...     'number_density': [100, 200, 150],
    ...     'area_interval': [2.5, 3.0, 2.8],
    ...     'stratum_num': [1, 2, 3]
    ... })
    >>> set_population_metrics(nasc_df, df_average_weight=0.5)
    >>> print(nasc_df.columns)
    Index(['number_density', 'area_interval', 'stratum_num', 'abundance', 'biomass',
    'biomass_density'])

    >>> # With grouped average weights
    >>> weight_df = pd.DataFrame({
    ...     'stratum_num': [1, 2, 3],
    ...     'avg_weight': [0.4, 0.6, 0.5]
    ... }).set_index('stratum_num')
    >>> nasc_df = pd.DataFrame({
    ...     'number_density': [100, 200, 150],
    ...     'area_interval': [2.5, 3.0, 2.8],
    ...     'stratum_num': [1, 2, 3]
    ... })
    >>> set_population_metrics(nasc_df,
    ...                        stratify_by='stratum_num',
    ...                        df_average_weight=weight_df['avg_weight'])

    >>> # Calculate only specific metrics
    >>> set_population_metrics(nasc_df, metrics=['abundance'], df_average_weight=0.5)
    """

    # Handle stratification and weight alignment
    reset_index_needed = False
    if stratify_by in df_nasc.columns and isinstance(df_average_weight, (pd.DataFrame, pd.Series)):
        # Set the stratification column as index for alignment
        df_nasc.set_index(stratify_by, inplace=True)
        reset_index_needed = True

        # Ensure weights are properly aligned with the NASC data
        if isinstance(df_average_weight, (pd.DataFrame, pd.Series)):
            df_average_weight = df_average_weight.copy()

            # If weights have the stratify_by column but not as index, set it
            if hasattr(df_average_weight, "columns") and stratify_by in df_average_weight.columns:
                df_average_weight.set_index(stratify_by, inplace=True)

            # Reindex to match NASC data
            df_average_weight = df_average_weight.reindex(df_nasc.index)

    # Abundance
    if "abundance" in metrics:
        df_nasc["abundance"] = np.round(df_nasc["area_interval"] * df_nasc["number_density"])

    # Biomass
    if "biomass" in metrics:
        # ---- Temporary abundance, if not already present
        if "abundance" not in df_nasc.columns:
            abundance_tmp = np.round(df_nasc["area_interval"] * df_nasc["number_density"])
        else:
            abundance_tmp = df_nasc["abundance"]
        # ---- Complete calculation
        df_nasc["biomass"] = abundance_tmp * df_average_weight

    # Biomass density
    if "biomass_density" in metrics:
        df_nasc["biomass_density"] = df_nasc["number_density"] * df_average_weight

    # Reset the index if it was modified
    if reset_index_needed:
        df_nasc.reset_index(inplace=True)
