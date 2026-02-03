from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
import xarray as xr

from ...utils import binned_distribution


def length_binned_weights(
    data: pd.DataFrame,
    length_bins: np.ndarray,
    regression_coefficients: Union[pd.Series, pd.DataFrame],
    impute_bins: bool = True,
    minimum_count_threshold: int = 0,
) -> xr.DataArray:
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
    xr.DataArray
        Fitted weights indexed by length bin and grouping dimensions.

        Dimensions
        ----------
        length_bin
            Length-bin intervals.
        *group_dims
            One dimension per grouping variable implied by
            ``regression_coefficients`` (e.g., ``sex``). Present only when
            grouped coefficients are supplied.

        Data
        ----
        weight_fitted
            Fitted weight values for each coordinate combination.

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
        # ---- Single set of coefficients - convert to DataFrame for consistent processing
        regression_df = pd.DataFrame([regression_coefficients])
        # ---- Reset index to avoid grouping complications
        regression_df.reset_index(drop=True, inplace=True)
        # ---- Assign metadata variables required for the output
        group_cols = []
    else:
        # ---- Already a DataFrame from groupby operation
        regression_df = regression_coefficients.reset_index()
        group_cols = [name for name in regression_coefficients.index.names if name is not None]

    # Initialize fitted weights dataframe
    weight_fitted_df = length_distribution.copy()

    # Cross merge with regression coefficients
    weight_fitted_df = weight_fitted_df.merge(regression_df, how="cross").rename(
        columns={"interval": "length_bin"}
    )

    # Predict weight per bin using allometric relationship: weight = 10^intercept * length^slope
    weight_fitted_df["weight_modeled"] = (
        10.0 ** weight_fitted_df["intercept"] * weight_fitted_df["bin"] ** weight_fitted_df["slope"]
    )

    # Get the column names if any grouping is required
    index_cols = group_cols + ["length_bin"]

    # Quantize weight counts per length bin
    binned_weight_distribution = (
        data.groupby(index_cols, observed=False)["length"].size().fillna(0).to_frame("count")
    )

    # Quantize weight means per length bin
    binned_weight_distribution["weight_mean"] = (
        data.groupby(index_cols, observed=False)["weight"].mean().fillna(0)
    )

    # Merge with the fitted weights
    binned_weight_distribution["weight_modeled"] = weight_fitted_df.set_index(index_cols)[
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

    # Reset index and prepare output
    result = binned_weight_distribution.reset_index()

    # Mutate, if needed, or otherwise return the pivoted DataFrame
    pivot_result = result.pivot(index="length_bin", columns=group_cols, values="weight_fitted")
    if isinstance(pivot_result, pd.Series):
        pivot_result = pivot_result.to_frame("all")

    # Get dimensions and coordinates
    if group_cols:
        dims = ["length_bin"] + group_cols
        coords = {
            "length_bin": pivot_result.index,
            **{
                col: (
                    pivot_result.columns.get_level_values(col)
                    if isinstance(pivot_result.columns, pd.MultiIndex)
                    else pivot_result.columns
                )
                for col in group_cols
            },
        }
        output = pivot_result.values
    else:
        dims = ["length_bin"]
        # ---- Convert MultiIndex to Index with no name to avoid conflict
        if isinstance(pivot_result.index, pd.MultiIndex) and len(pivot_result.index.names) == 1:
            coords = {"length_bin": pivot_result.index.get_level_values(0)}
        else:
            coords = {"length_bin": pivot_result.index}
        output = pivot_result.values.squeeze()

    # Convert to an xarray.DataArray
    return xr.DataArray(
        output,
        dims=dims,
        coords=coords,
        name="weight_fitted",
    )


def compute_abundance(
    dataset: pd.DataFrame,
    exclude_filter: Dict[str, str] = {},
    number_proportions: Optional[Dict[str, xr.Dataset]] = None,
) -> None:
    """
    Convert number density estimates in abundances.

    Parameters
    ----------
    dataset : pd.DataFrame
        DataFrame containing transect data with number densities already computed. Must include:
        - 'number_density': Number of animals per unit area (animals/nmi²)
        - 'area_interval': Area of each sampling interval (nmi²)
    stratify_by : List[str], default []
        Column name to use for stratification when aligning with the number proportions dictionary,
        if specified.
    group_by : List[str], default []
        Grouping columns to apply to number density and abundances (e.g. sex). This will produce
        additional columns formatted as `number_density_{group}` and `abundance_{group}`.
    exclude_filter : Dict[str, str], default {}
        Dictionary specifying which population segments to exclude and redistribute. Keys are
        column/index names, values are the categories to exclude.
    number_proportions : Optional[Dict[str, pd.DataFrame]], default None
        When provided, the number densities and abundances will be reapportioned according to the
        defined number proportions. If `group_by` is not empty, then the resulting DataFrame will
        include additional columns for each group as well as the overall number density and
        abundance.

    Returns
    -------
    None
        Function modifies the transect DataFrame in-place.

    Note
    ----
    When `number_proportions` is provided, the function first calculates base abundance as
    area_interval * number_density, then applies proportions to create grouped abundance columns.
    The original 'abundance' column represents the total.
    """

    # If no grouping, run the simple abundance calculation
    dataset["abundance"] = dataset["area_interval"] * dataset["number_density"]

    # Compute grouped values, if needed
    if number_proportions is not None:
        # ---- Select overall proportions from each dataset
        overall_props = {k: ds["proportion_overall"] for k, ds in number_proportions.items()}
        # ---- Get the set of coordinate names for each DataArray
        coord_sets = [set(da.coords) for da in overall_props.values()]
        # ---- Find the intersection (shared coordinates)
        shared_coords = set.intersection(*coord_sets) - {"length_bin"}
        # ---- Find overlapping indices
        idx_names = list(shared_coords.intersection(set(dataset.columns)))
        nonidx_names = [id for id in list(shared_coords) if id not in idx_names]
        # ---- Create grouped table from number proportions
        grouped_proportions = sum(
            da.sum(dim=[d for d in da.dims if d not in shared_coords])
            for da in overall_props.values()
        )
        # ---- Apply exclusion filter, if required
        if exclude_filter:
            # ---- Parse existing labels
            to_drop = {
                k: v
                for k, v in exclude_filter.items()
                if k in grouped_proportions.coords and v in grouped_proportions.coords[k].values
            }
            if to_drop:
                grouped_proportions_excl = grouped_proportions.drop_sel(to_drop)
            else:
                grouped_proportions_excl = grouped_proportions
        else:
            grouped_proportions_excl = grouped_proportions
        # ---- Refine if no grouping
        if len(nonidx_names) == 0:
            number_density_vals = dataset["number_density"].values
            abundance_vals = dataset["abundance"].values
        else:
            number_density_vals = dataset["number_density"].values[:, None]
            abundance_vals = dataset["abundance"].values[:, None]
        # ---- Reindex the DataArray
        ridx = {k: dataset[k] for k in idx_names}
        grouped_proportions_ridx = grouped_proportions_excl.reindex(ridx).fillna(0.0)
        # ---- Set the index of the tabular data
        dataset.set_index(idx_names, inplace=True)
        # ---- Compute number density
        grouped_number_density = number_density_vals * grouped_proportions_ridx
        # ---- Compute abundance
        grouped_abundance = abundance_vals * grouped_proportions_ridx
        # ---- Add the number densities to the dataset
        dataset[
            [f"number_density_{c}" for c in grouped_number_density.coords[nonidx_names[0]].values]
        ] = grouped_number_density
        # ---- Add abundances to the dataset
        dataset[[f"abundance_{c}" for c in grouped_abundance.coords[nonidx_names[0]].values]] = (
            grouped_abundance
        )
        # ---- Reset the index
        dataset.reset_index(inplace=True)


# !!! ----> Rename `dataset` to not be `dataset` since it is a DataFrame, not an xarray object
def matrix_multiply_grouped_table(
    dataset: pd.DataFrame,
    table: xr.DataArray,
    variable: str,
    output_variable: str,
    group: Optional[str] = None,
) -> None:
    """
    Multiply multiple data columns by identically indexed grouped table columns from a separate
    DataFrame.

    Parameters
    ----------
    dataset : pd.DataFrame
        DataFrame containing the data to be multiplied.
    table : xr.DataArray
        DataArray containing the grouped table with columns to multiply against.
    variable : str
        The name of the variable in `dataset` that will be multiplied by the grouped table.
        This variable should have columns named as `{variable}_{suffix}` where `suffix` is a
        group identifier.
    output_variable : str
        The name of the output variable/column that will be created in `dataset` to store the
        results.
    group : Optional[str], default None
        The specific group to filter the grouped table by. If None, all groups will be used.
    """

    # Create pattern for column filtering
    prefix_pattern = variable + "_"

    # Get the overlapping columns
    variable_columns = dataset.filter(like=prefix_pattern).columns

    # Gather the suffixes corresponding to the target group
    suffixes = variable_columns.str.replace(prefix_pattern, "", regex=False).to_list()

    # Reindex table
    ridx = {k: dataset.reset_index()[k] for k in dataset.index.names}
    table_ridx = table.reindex(ridx).fillna(0.0)

    # Apply an inclusion filter
    table_groups = table_ridx.sel({group: suffixes})

    # Get columns that also exist for dataset
    variable_overlap = [prefix_pattern + col for col in table_groups[group].values]

    # Run the multiplication
    table_matrix = dataset.filter(variable_overlap).to_numpy() * table_groups.T

    # Set up new columns
    new_columns = [f"{output_variable}_{c}" for c in table_groups[group].values]

    # Add the output variables
    dataset[new_columns] = table_matrix.fillna(0.0).values

    # Calculate the remainder comprising the ungrouped values
    remainder = dataset[variable] - dataset[variable_overlap].fillna(0.0).sum(axis=1)

    # Calculate the output variable for the ungrouped/excluded values
    remainder_matrix = remainder * table_ridx.sel({group: "all"}).fillna(0.0)

    # Compute the overall output variable
    dataset[output_variable] = dataset[new_columns].sum(axis=1) + remainder_matrix


def compute_biomass(
    dataset: pd.DataFrame,
    stratum_weights: Optional[Union[xr.DataArray, float]] = None,
) -> None:
    """
    Convert number density and abundance estimates to biomass density and total biomass,
    respectively.

    Parameters
    ----------
    dataset : pd.DataFrame
        DataFrame containing transect data with number densities already computed. Must include:
        - 'number_density': Number of animals per unit area (animals/nmi²)
        - 'abundance': Total number of animals in each interval (animals)
    stratify_by : List[str], default []
        Column name to use for stratification when aligning with the number proportions dictionary,
        if specified.
    group_by : List[str], default []
        Grouping columns to apply to number density and abundances (e.g. sex). This will produce
        additional columns formatted as `number_density_{group}` and `abundance_{group}`.
    df_average_weight : Optional[Dict[str, pd.DataFrame]], default None
        DataFrame or Series containing the average weight data for defined strata. This DataFrame
        must be:
        - Indexed by the same column as `stratify_by` if specified
        - If a DataFrame, the columns should correspond to the grouping defined in `group_by`. For
        example, if `group_by=['sex']`, then the DataFrame should have a column for each. If there
        are additional groups that are present in the DataFrame but will not be used, then a
        global set of values named 'all' should be present.

    Returns
    -------
    None
        Function modifies the transect DataFrame in-place.

    Examples
    --------
    >>> # Calculate biomass from abundance using weight table
    >>> weight_table = pd.DataFrame({'female': [0.5, 0.6], 'male': [0.4, 0.5], 'all': [0.45, 0.55]})
    >>> dataset = pd.DataFrame({'abundance_female': [100, 200], 'abundance_male': [150, 250], \
        'abundance': [250, 450]})
    >>> compute_biomass(
    ...     dataset=dataset,
    ...     stratify_by=['stratum'],
    ...     group_by=['sex'],
    ...     df_average_weight=df_stratum_weights,
    ... )
    >>> # Creates 'biomass_female', 'biomass_male' in-place
    """

    # Find overlapping indices
    idx_names = list(set(stratum_weights.coords).intersection(set(dataset.columns)))
    nonidx_names = [id for id in set(stratum_weights.coords) if id not in idx_names]

    # Set index
    dataset.set_index(idx_names, inplace=True)

    # !!! ----> Rename `dataset` to not be `dataset` since it is a DataFrame, not an xarray object
    # If grouped beyond just index
    if len(nonidx_names) > 0:
        # ---- Compute the biomass densities across groups
        matrix_multiply_grouped_table(
            dataset,
            table=stratum_weights,
            group=nonidx_names[0],
            variable="number_density",
            output_variable="biomass_density",
        )
        # ---- Compute the biomass densities across groups
        matrix_multiply_grouped_table(
            dataset,
            table=stratum_weights,
            group=nonidx_names[0],
            variable="abundance",
            output_variable="biomass",
        )
    # ---- Ungrouped
    else:
        # ---- Reindex table
        ridx = {k: dataset.reset_index()[k] for k in dataset.index.names}
        stratum_weights_ridx = stratum_weights.reindex(ridx).fillna(0.0)
        # ---- Compute biomass densities
        dataset["biomass_density"] = dataset["number_density"] * stratum_weights_ridx.sel(
            {nonidx_names[0]: "all"}
        )
        # ---- Compute biomass
        dataset["biomass"] = dataset["abundance"] * stratum_weights.sel({nonidx_names[0]: "all"})

    # Reset the index
    dataset.reset_index(inplace=True)


def remove_specimen_hauls(
    biodata_dict: Dict[str, pd.DataFrame],
) -> None:
    """
    Remove hauls from the catch data where all samples were individually processed.

    This function filters the catch data to exclude hauls that don't have corresponding length
    frequency data, ensuring consistency between catch weights and length samples. This prevents
    double-counting when specimen data represents the entire catch.

    Parameters
    ----------
    biodata_dict : Dict[str, pd.DataFrame]
        Dictionary containing biological data with keys typically including 'catch' and 'length'.
        The 'length' DataFrame should contain a 'haul_num' column identifying which hauls have
        length frequency data. The 'catch' DataFrame will be filtered to match.

    Returns
    -------
    None
        Function modifies biodata_dict in place by filtering the 'catch' DataFrame.

    Notes
    -----
    This function addresses a common issue in fisheries data where:
    - Some hauls have bulk catch weights but no individual fish measurements
    - Other hauls have detailed individual fish data that represents the entire catch
    - Including both would lead to double-counting of biomass

    By filtering catch data to only include hauls with length frequency data,
    the function ensures that biomass estimates are based on consistent sampling methods.

    The function is typically called after loading biological data but before
    computing length-weight relationships or abundance estimates.
    """

    # Get unique haul numbers
    haul_numbers = biodata_dict["length"]["haul_num"].unique()

    # Find incompatible hauls
    biodata_dict["catch"] = biodata_dict["catch"].loc[
        biodata_dict["catch"]["haul_num"].isin(haul_numbers)
    ]
