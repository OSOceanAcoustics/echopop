import warnings
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from . import utils

warnings.simplefilter("always")


def partition_transect_data(
    dataset: pd.DataFrame,
    partition_dict: Dict[str, Union[pd.DataFrame, pd.Series]],
) -> pd.DataFrame:
    """
    Partition NASC, abundance (and number density), and biomass (and biomass density) transect
    values across indexed groups.

    Parameters
    ----------
    dataset : pd.DataFrame
        DataFrame containing transect data with number densities already computed. Must include:
        - A column that matches the index of the `partition_dict` DataFrames/Series.
    partition_dict : Dict[str, Union[pd.DataFrame, pd.Series]]
        Dictionary containing partitioning data for NASC, abundance, and biomass. Valid key names
        are limited to 'abundance', 'biomass', and 'nasc'. Each key maps to different variables:
        - 'abundance': pd.Series or DataFrame with abundance proportions. This partitions the
        number density and abundance estimates.
        - 'biomass': pd.Series or DataFrame with biomass proportions. This partitions the biomass
        density and biomass estimates.
        - 'nasc': pd.Series or DataFrame with NASC proportions.

    Returns
    -------
    pd.DataFrame
        DataFrame with the same structure as the input dataset, but with defined partitions applied
        to NASC and the other abundance/biomass columns.

    Notes
    -----
    The DataFrames/Series in `partition_dict` must have indices that correspond to values in the
    `dataset`. Missing indices will result in NaN values for those rows. The function automatically
    identifies and uses the appropriate index columns for merging.
    """

    # Create copy
    dataset = dataset.copy()

    # Get the index names
    index_names = list(set().union(*[set(df.index.names) for df in partition_dict.values()]))

    # Set the index of the input dataset
    dataset.set_index(index_names, inplace=True)

    # NASC, if present
    if "nasc" in partition_dict:
        dataset["nasc"] = dataset["nasc"] * (1 - partition_dict["nasc"].reindex(dataset.index))

    # Abundance and number density, if present
    if "abundance" in partition_dict:
        # ---- Get the inverse proportions
        abundance_proportions = 1 - partition_dict["abundance"].reindex(dataset.index)
        # ---- Map the appropriate columns for abundance
        abundance_names = dataset.filter(like="abundance").columns
        # ---- Adjust abundances
        dataset[abundance_names] = dataset[abundance_names].mul(abundance_proportions, axis=0)
        # ---- Map the appropriate columns for number density
        number_density_names = dataset.filter(like="number_density").columns
        # ---- Adjust number densities
        dataset[number_density_names] = dataset[number_density_names].mul(
            abundance_proportions, axis=0
        )

    # Biomass and biomass density, if present
    if "biomass" in partition_dict:
        # ---- Get the inverse proportions
        biomass_proportions = 1 - partition_dict["biomass"].reindex(dataset.index)
        # ---- Map the appropriate columns for biomass and biomass density
        biomass_names = dataset.filter(like="biomass").columns
        # ---- Adjust biomass
        dataset[biomass_names] = (biomass_proportions * dataset[biomass_names].T).T

    # Reset the index
    dataset.reset_index(inplace=True)

    # Return the partitioned dataset
    return dataset


def mesh_biomass_to_nasc(
    mesh_data_df: pd.DataFrame,
    biodata: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    mesh_biodata_link: Dict[str, str],
    stratum_weights_df: pd.DataFrame,
    stratum_sigma_bs_df: pd.DataFrame,
    group_by: List[str] = [],
) -> None:
    """
    Convert biomass estimates distributed across a grid or mesh into grouped (e.g. sex-specific)
    abundance and NASC values.

    This function takes kriged biomass density estimates from a mesh/grid and converts them to
    abundance and NASC (Nautical Area Scattering Coefficient) values by applying biological
    proportion data and acoustic backscattering coefficients. The function modifies the input
    mesh DataFrame in place.

    Parameters
    ----------
    mesh_data_df : pd.DataFrame
        DataFrame containing kriged biomass density estimates with spatial mesh information.
        Must contain columns for linking to biodata and either 'biomass', or 'biomass_density'
        and 'area', columns.
    biodata : pd.DataFrame or Dict[str, pd.DataFrame]
        Biological proportion data. Can be a single DataFrame or dictionary of DataFrames
        containing proportion data for different biological groups (e.g., aged/unaged).
    mesh_biodata_link : Dict[str, str]
        Dictionary mapping column names from mesh_data_df to biodata column names for joining.
        Example: {"geostratum_ks": "stratum_ks"}
    stratum_weights_df : pd.DataFrame
        DataFrame containing average weights per stratum for converting biomass to abundance.
        Index should match the linked columns from mesh_biodata_link.
    stratum_sigma_bs_df : pd.DataFrame
        DataFrame containing sigma_bs (backscattering coefficient) values per stratum.
        Must have 'sigma_bs' column with index matching linked columns.
    group_by : List[str], default []
        List of column names to group biological data by (e.g., ["sex", "age_bin"]).

    Returns
    -------
    None
        Function modifies mesh_data_df in place, adding new columns for biomass, abundance,
        and NASC values.

    Raises
    ------
    KeyError
        If required link columns are missing from biodata or mesh DataFrames.

    Notes
    -----
    The function performs the following steps:
    1. Computes biomass from biomass_density x area if not present
    2. Links mesh data to biological proportions using mesh_biodata_link
    3. Applies proportions to distribute total biomass across biological groups
    4. Converts biomass to abundance using stratum weights
    5. Calculates NASC using abundance and sigma_bs values

    The NASC calculation uses the formula: NASC = abundance x sigma_bs x 4π
    """

    # Convert to a Dictionary if needed
    if isinstance(biodata, pd.DataFrame):
        biodata = {"": biodata}

    # Get the column indices that will be used for setting shared indices
    column_names = list(
        set(
            col
            for df in biodata.values()
            if isinstance(df.columns, pd.Index)
            for col in df.columns.names
        )
    )

    # Check
    # ---- Setdiff: biodata
    missing_biodata_columns = set(list(mesh_biodata_link.values())).difference(column_names)
    if len(missing_biodata_columns) > 0:
        raise KeyError(
            f"The following link columns were missing from the biological data: "
            f"{', '.join(col for col in missing_biodata_columns)}."
        )
    # ---- Setdiff: mesh
    missing_mesh_columns = set(list(mesh_biodata_link.keys())).difference(mesh_data_df.columns)
    if len(missing_mesh_columns) > 0:
        raise KeyError(
            f"The following link columns were missing from the mesh DataFrame: "
            f"{', '.join(col for col in missing_mesh_columns)}."
        )

    # Compute biomass if not already in DataFrame
    if "biomass" not in mesh_data_df.columns:
        mesh_data_df["biomass"] = mesh_data_df["biomass_density"] * mesh_data_df["area"]

    # Set the link column names
    mesh_data_df.rename(columns=mesh_biodata_link, inplace=True)

    # Map the link columns
    mesh_index = [v for v in mesh_biodata_link.values() if v in mesh_data_df.columns]

    # Index
    mesh_data_df.set_index(mesh_index, inplace=True)

    # Stack the proportions to be as a function of stratum
    biodata_stk = {
        key: df.unstack(group_by).sum().unstack(group_by).reindex(mesh_data_df.index)
        for key, df in biodata.items()
    }

    # Get the column names of the biodata
    biodata_columns = sorted(
        set(
            col
            for df in biodata_stk.values()
            if isinstance(df, pd.DataFrame)
            for col in df.columns.tolist()
        )
    )
    # ---- Fill in if empty
    if len(biodata_columns) == 0:
        biodata_columns = set(["apportioned"])

    # Stack the proportions
    stacked_proportions = pd.concat(list(biodata_stk.values()), axis=1)

    # If only a single expected output, summarize
    if len(biodata_columns) == 1:
        stacked_proportions = stacked_proportions.sum(axis=1).to_frame(list(biodata_columns)[0])

    # Multiply by biomass column, broadcasting over all columns
    stacked_proportions = stacked_proportions.mul(mesh_data_df["biomass"], axis=0)

    # Group columns by name and sum (to handle duplicate columns from stacking)
    stacked_proportions = stacked_proportions.T.groupby(stacked_proportions.columns).sum().T

    # Prefix columns
    stacked_proportions = stacked_proportions.add_prefix("biomass_")

    # Add to the mesh DataFrame
    mesh_data_df[stacked_proportions.columns.tolist()] = stacked_proportions

    # Get the column names of the resulting biomass values
    biomass_columns = stacked_proportions.columns.tolist() + ["biomass"]

    # Reindex the stratum weights
    stratum_weights_df_idx = stratum_weights_df.reindex(mesh_data_df.index)

    # Calculate abundance
    mesh_data_df[["abundance"] + [f"abundance_{name}" for name in biodata_columns]] = mesh_data_df[
        biomass_columns
    ].div(stratum_weights_df_idx, axis=0)

    # Index the stratum sigma_bs
    stratum_sigma_bs_df_idx = stratum_sigma_bs_df.reindex(mesh_data_df.index)

    # Calculate NASC
    mesh_data_df["nasc"] = (
        mesh_data_df["abundance"] * stratum_sigma_bs_df_idx["sigma_bs"] * 4.0 * np.pi
    )

    # Rename the aligned columns
    # ---- Create inverted link dictionary
    inverted_link = {v: k for k, v in mesh_biodata_link.items()}

    # Reset the index
    mesh_data_df.reset_index(inplace=True)

    # Rename
    mesh_data_df.rename(columns=inverted_link, inplace=True)


def distribute_kriged_estimates(
    mesh_data_df: pd.DataFrame,
    proportions: Union[Dict[str, Any], pd.DataFrame],
    variable: str,
    group_by: List[str],
    stratify_by: List[str],
    mesh_proportions_link=Dict[str, str],
) -> Dict[str, pd.DataFrame]:
    """
    Distribute kriged estimates (e.g. abundance, biomass) across biological groups and other metrics
    such as age and length bins.

    This function takes kriged estimates from a mesh and distributes them across different
    biological categories (e.g., sex, age, length) using proportion data. The distribution
    is performed by multiplying the kriged estimates with the corresponding proportions.

    Parameters
    ----------
    mesh_data_df : pd.DataFrame
        DataFrame containing kriged estimates with spatial mesh information.
        Must contain the variable column specified in the 'variable' parameter.
    proportions : Dict[str, Any] or pd.DataFrame
        Proportion data for distributing the kriged estimates. Can be a single DataFrame
        or dictionary of DataFrames (e.g., {"aged": df_aged, "unaged": df_unaged}).
    variable : str
        Name of the column in mesh_data_df containing the values to distribute
        (e.g., "biomass", "abundance").
    group_by : List[str]
        List of column names that define the biological groups for distribution
        (e.g., ["sex", "age_bin", "length_bin"]).
    stratify_by : List[str]
        List of column names used for stratification (e.g., ["stratum_ks"]).
    mesh_proportions_link : Dict[str, str]
        Dictionary mapping column names from mesh_data_df to proportions column names
        for joining. Example: {"geostratum_ks": "stratum_ks"}

    Returns
    -------
    Dict[str, pd.DataFrame]
        Distributed estimates. Returns a dictionary. Each DataFrame contains the original
        biological group structure with values distributed according to proportions.

    Notes
    -----
    The function performs the following steps:
    1. Converts single DataFrame to dictionary format if needed
    2. Aggregates mesh data by stratification variables
    3. Creates pivot tables from proportions data
    4. Multiplies aggregated mesh values by proportions
    5. Returns results in the same format as input proportions

    The distribution preserves the biological group structure while scaling values
    according to the proportion data provided.
    """

    # Convert to compatible format
    if isinstance(proportions, pd.DataFrame):
        proportions = {"distributed_mesh": proportions}

    # Verify that `variable` exists
    if variable not in mesh_data_df.columns:
        raise ValueError(f"Variable '{variable}' not found in `mesh_data_df`.")

    # Sum variables over indices
    mesh_data_pvt = mesh_data_df.rename(columns=mesh_proportions_link).pivot_table(
        index=stratify_by, values=variable, aggfunc="sum", observed=False
    )

    # Parse the additional columns that are required for grouping
    proportions_group_columns = {
        k: [c for c in (list(v.index.names) + list(v.columns)) if c in group_by]
        for k, v in proportions.items()
    }

    # Convert to DataFrame(s) to pivot table(s)
    proportions_grouped_pvt = {
        k: (
            df
            if utils.is_pivot_table(df)
            else utils.create_pivot_table(
                df,
                index_cols=proportions_group_columns[k],
                strat_cols=stratify_by,
                value_col="proportion_overall",
            )
        )
        for k, df in proportions.items()
    }

    # Distribute the variable over each table
    apportioned_grouped_pvt = {
        k: df.mul(mesh_data_pvt[variable]).fillna(0.0) for k, df in proportions_grouped_pvt.items()
    }

    # Return output in the same format as the proportions input
    return apportioned_grouped_pvt


def impute_kriged_table(
    reference_table_df: pd.DataFrame,
    initial_table_df: pd.DataFrame,
    standardized_table_df: pd.DataFrame,
    group_by: List[str],
    impute_variable: List[str],
    table_reference_indices: List[str],
) -> pd.DataFrame:
    """
    Perform nearest-neighbor imputation of missing values in population apportionment tables.

    This function fills missing values in standardized population tables by finding the nearest
    non-zero reference values and using them to impute proportional values based on the
    reference table structure.

    Parameters
    ----------
    reference_table_df : pd.DataFrame
        Reference table used to determine which values should be imputed and to calculate
        imputation ratios. Typically contains aged fish data with complete age structure.
    initial_table_df : pd.DataFrame
        Original population table before standardization. Used to determine which cells
        have actual data vs. missing values.
    standardized_table_df : pd.DataFrame
        Standardized population table that needs imputation. This will be modified with
        imputed values.
    group_by : List[str]
        List of column names that define the grouping variables (e.g., ["sex"]).
    impute_variable : List[str]
        List of variables that define the dimension being imputed (e.g., ["age_bin"]).
    table_reference_indices : List[str]
        List of index names that are in the reference table but not in the population table.

    Returns
    -------
    pd.DataFrame
        Copy of standardized_table_df with missing values imputed using nearest-neighbor
        approach based on reference table proportions.

    Notes
    -----
    The imputation algorithm:
    1. Identifies cells in the reference table that are zero (missing in reference)
    2. Finds corresponding non-zero cells in the initial table
    3. For each missing cell, finds the nearest non-zero reference cell
    4. Imputes values using proportional scaling based on reference ratios

    This approach preserves the biological structure while filling gaps in the data
    using the most similar available information.
    """

    # Sum the reference proportions across the impute variable
    reference_stk = (
        reference_table_df.sum(axis=1).unstack(impute_variable).sum(axis=1).unstack(group_by)
    )

    # Unstack the grouping and imputate variable
    reference_group_stk = reference_table_df.sum(axis=1).unstack(group_by + impute_variable)

    # Get the mask for all 0.0's and non-zeros
    ref_zero_mask = reference_stk == 0.0
    ref_nonzero_mask = reference_stk != 0.0

    # Gather the indices for each column
    ref_zero_indices = {
        col: reference_stk.index[ref_zero_mask[col]].tolist() for col in reference_stk.columns
    }
    ref_nonzero_indices = {
        col: reference_stk.index[ref_nonzero_mask[col]].tolist() for col in reference_stk.columns
    }

    # Create translation for row numbers
    interval_to_numeric = {interval: i for i, interval in enumerate(reference_stk.index)}

    # Convert to row numbers
    ref_zero_rows = {
        col: np.array([interval_to_numeric[ival] for ival in intervals])
        for col, intervals in ref_zero_indices.items()
    }
    ref_nonzero_rows = {
        col: np.array([interval_to_numeric[ival] for ival in intervals])
        for col, intervals in ref_nonzero_indices.items()
    }

    # Check keys
    if not all(col in reference_stk.columns for col in initial_table_df.columns):
        # ---- Get difference
        missing_columns = set(initial_table_df.columns).difference(reference_stk.columns)
        # ---- Print warning
        warnings.warn(
            f"The following columns are missing from the reference table and will therefore be "
            f"skipped during imputation: {', '.join(missing_columns)}. ",
            stacklevel=2,
        )

    # Apply the non-zero reference indices to each column
    table_nonzeros_mask = {
        col: initial_table_df[col].loc[ref_zero_indices[col]] != 0.0
        for col in reference_stk.columns
    }

    # Get the actual indices of the masked values
    nonzero_reference_to_table_indices = {
        col: np.array(ref_zero_indices[col])[table_nonzeros_mask[col]]
        for col in table_nonzeros_mask
    }

    # Convert to numeric indices for compatibility with `iloc`
    nonzero_reference_to_table_rows = {
        col: np.array([interval_to_numeric[ival] for ival in intervals])
        for col, intervals in nonzero_reference_to_table_indices.items()
    }

    # Create copy of table
    standardized_table_copy = standardized_table_df.copy()

    # Define the column indices to ensure consistent ordering across the table
    # ---- Original table indices
    column_indices = table_reference_indices + group_by
    # ---- Reverse table indices
    column_indices_rev = list(reversed(column_indices))

    # Check if the columns are in the expected order
    if list(standardized_table_copy.columns.names) == column_indices:
        # ---- Swap the indices
        standardized_table_copy.columns = standardized_table_copy.columns.reorder_levels(
            column_indices_rev
        )

    # Get the nearest-neighbor rows and recompute the indices
    imputed_rows = {
        col: arr[
            np.argmin(
                np.abs(
                    ref_zero_rows[col][table_nonzeros_mask[col]][:, np.newaxis]
                    - ref_nonzero_rows[col]
                ),
                axis=1,
            )
        ]
        for col, arr in ref_nonzero_rows.items()
    }

    # Impute to replace these values
    imputed_values = {
        col: (
            initial_table_df[col].loc[nonzero_reference_to_table_indices[col]].to_numpy()
            * reference_group_stk.iloc[imputed_rows[col]][col].T
            / reference_stk.iloc[imputed_rows[col]][col]
        ).T
        for col in reference_stk.columns
    }

    # Update the standardized values
    for col in initial_table_df.columns:
        if col in reference_stk.columns:
            standardized_table_copy.iloc[
                nonzero_reference_to_table_rows[col], standardized_table_copy.columns.get_loc(col)
            ] = imputed_values[col]

    # Return to the original column index order
    standardized_table_copy = standardized_table_copy.reorder_levels(column_indices, axis=1)

    # Return the imputed standardized table
    return standardized_table_copy


def standardize_kriged_estimates(
    population_table: pd.DataFrame,
    reference_table: pd.DataFrame,
    group_by: List[str],
    impute: bool = True,
    impute_variable: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Standardize kriged population estimates using reference proportions and optionally impute
    missing values.

    This function standardizes population estimates (e.g., unaged fish distributions) by applying
    reference proportions (e.g., from aged fish data) to ensure consistent age structure across
    datasets. It can also perform nearest-neighbor imputation for missing values.

    Parameters
    ----------
    population_table : pd.DataFrame
        Population estimates to be standardized. Should be a pivot table or similar structure
        with biological groups as index/columns.
    reference_table : pd.DataFrame
        Reference population data used for standardization. Should have the same or compatible
        structure as population_table with additional detail (e.g., age information).
    group_by : List[str]
        List of column names that define the grouping variables for standardization
        (e.g., ["sex"] to standardize within each sex).
    impute : bool, default True
        Whether to perform nearest-neighbor imputation for missing values after standardization.
    impute_variable : List[str], optional
        List of variables to use for imputation. Required if impute=True.
        Typically refers to the dimension being imputed (e.g., ["age_bin"]).

    Returns
    -------
    pd.DataFrame
        Standardized population estimates with the same structure as population_table
        but adjusted according to reference proportions. If impute=True, missing values
        are filled using nearest-neighbor imputation.

    Notes
    -----
    The standardization process:
    1. Computes reference proportions from the reference table
    2. Applies these proportions to the population table
    3. Redistributes values to match reference age/size structure
    4. Optionally imputes missing values using nearest neighbors

    This is commonly used to distribute unaged fish data according to the age
    structure observed in aged fish samples from the same area/stratum.
    """

    # Get the table indices
    table_indices = population_table.index.names

    # Get the shared index names across the reference tables
    reference_index_names = list(reference_table.index.names)

    # Get the indices for each reference table required for summation that must be unstacked
    reference_stack_indices = list(set(reference_index_names).difference(table_indices))

    # Reorient the columns accordingly
    population_table_reshape = population_table.sum(axis=1).unstack(group_by)

    # Standardize the apportioned table values
    standardized_table = (
        population_table_reshape
        * reference_table.sum(axis=1).unstack(reference_stack_indices + group_by)
        / reference_table.unstack(reference_stack_indices).sum(axis=1).unstack(group_by)
    ).fillna(0.0)

    # If imputation is requested, perform it
    if not impute:
        return standardized_table
    else:
        # ---- Impute
        standardized_table_imputed = impute_kriged_table(
            reference_table_df=reference_table,
            initial_table_df=population_table_reshape,
            standardized_table_df=standardized_table,
            table_reference_indices=reference_stack_indices,
            group_by=group_by,
            impute_variable=impute_variable,
        )
        return standardized_table_imputed


def combine_population_tables(
    population_table: Dict[str, pd.DataFrame],
    table_names: List[str],
    table_index: List[str],
    table_columns: List[str],
) -> pd.DataFrame:
    """
    Combine and sum population estimates across defined tables to yield a single population table.

    This function takes multiple population estimate tables and combines them into a single
    consolidated table by stacking, pivoting, and summing the values across common dimensions.

    Parameters
    ----------
    population_table : Dict[str, pd.DataFrame]
        Dictionary of population estimate tables to combine. Keys are table names (e.g., "aged",
        "unaged", "speciesA", "speciesB) and values are DataFrames with population data.
    table_names : List[str]
        List of table names (keys) from population_table to include in the combination. Tables not
        in this list will be ignored.
    table_index : List[str]
        List of column names to use as the index in the final combined table (e.g., ["length_bin"]).
    table_columns : List[str]
        List of column names to use as columns in the final combined table (e.g., ["age_bin",
        "sex"]).

    Returns
    -------
    pd.DataFrame
        Combined population table with table_index as index and table_columns as columns.
        Values are summed across all input tables for each index/column combination.

    Notes
    -----
    The combination process:
    1. Filters population_table to include only tables listed in table_names
    2. Stacks each table to convert to long format with 'value' column
    3. Creates pivot tables with consistent index and column structure
    4. Sums all pivot tables element-wise to get final combined table
    """

    # Subset the table to only include the tables-of-interest
    tables_to_combine = {
        k: v.stack(v.columns.names, future_stack=True).to_frame("value")
        for k, v in population_table.items()
        if k in table_names
    }

    # Create a pivot table for each table with identical indices
    compatible_tables = {
        k: utils.create_pivot_table(
            v, index_cols=table_index, strat_cols=table_columns, value_col="value"
        )
        for k, v in tables_to_combine.items()
    }

    # Sum the compatible tables
    return sum(compatible_tables.values())


def redistribute_population_table(
    population_table: pd.DataFrame,
    exclusion_filter: Dict[str, Any],
    group_by: List[str],
) -> pd.DataFrame:
    """
    Redistribute population estimates across groups after excluding a specific subset.

    This function removes specified population segments (e.g., age-1 fish) and redistributes
    their population estimates proportionally across the remaining groups. This is commonly used
    when certain age classes need to be excluded from final estimates but their contributions
    should be reallocated to other age groups.

    Parameters
    ----------
    population_table : pd.DataFrame
        Population estimates table with biological groups as index/columns. Typically a pivot table
        with population values.
    exclusion_filter : Dict[str, Any]
        Dictionary specifying which population segments to exclude and redistribute. Keys are
        column/index names, values are the categories to exclude. Example: {"age_bin": [1]} to
        exclude age-1 fish.
    group_by : List[str]
        List of column names that define the grouping variables for redistribution (e.g., ["sex"]
        to redistribute within each sex separately).

    Returns
    -------
    pd.DataFrame
        Population table with excluded segments removed and their values redistributed
        proportionally across remaining groups. Has the same structure as input table.

    Raises
    ------
    UserWarning
        If the redistributed table sums don't match the original table sums within tolerance.

    Notes
    -----
    The redistribution algorithm:
    1. Identifies population segments matching the exclusion filter
    2. Calculates the total excluded biomass/abundance for each group
    3. Removes excluded segments (sets them to zero)
    4. Redistributes excluded totals proportionally across remaining segments
    5. Ensures total population remains constant after redistribution

    This maintains the biological realism of population estimates while allowing for policy-based
    exclusions (e.g., removing age-1 fish from assessments).
    """

    # Find any columns that are not in the group_by list
    # ---- Get column names
    column_names = population_table.columns.names
    # ---- Identify extra columns that are not in the group_by list
    extra_columns = [col for col in column_names if col not in group_by]
    # ---- Stack the population table
    stacked_table = population_table.stack(extra_columns, future_stack=True)

    # If no appropriate filter is defined, then nothing is redistributed
    if len(exclusion_filter) == 0:
        return population_table

    # Apply inverse of exclusion filter to get the values being excluded
    excluded_grouped_table = utils.apply_filters(stacked_table, include_filter=exclusion_filter)

    # Replace the excluded values in the full table with 0.
    filtered_grouped_table = utils.apply_filters(
        stacked_table, exclude_filter=exclusion_filter, replace_value=0.0
    )

    # Get the sums for each group across the excluded and filtered tables
    # ---- Filtered/included
    filtered_grouped_sum = filtered_grouped_table.sum()
    # ---- Excluded
    excluded_grouped_sum = excluded_grouped_table.sum()
    # ---- Handling if `pandas.Series`
    if isinstance(excluded_grouped_sum, pd.Series):
        excluded_grouped_sum = excluded_grouped_sum.reindex(
            index=filtered_grouped_sum.index, method="ffill"
        )

    # Get the redistributed values that will be added to the filtered table values
    adjustment_table = (
        filtered_grouped_table * excluded_grouped_sum / filtered_grouped_sum
    ).fillna(0.0)

    # Add the adjustments to the filtered table
    filtered_grouped_table += adjustment_table

    # Check sums
    if len(set(exclusion_filter) & set(group_by)) > 0:
        check_sum = (
            filtered_grouped_table.stack(group_by).sum() - stacked_table.stack(group_by).sum()
        )
    else:
        check_sum = filtered_grouped_table.sum() - stacked_table.sum()
    if np.any(check_sum > 1e-6):
        # ---- Raise a warning with the indices where the sums do not match
        warnings.warn(
            "The sums of the table with the redistributed estimates do not match the original "
            "table filtered table do not match the original table."
        )

    # Restore the original column structure
    redistributed_table = filtered_grouped_table.unstack(extra_columns).reorder_levels(
        column_names, axis=1
    )

    # Return the the redistributed table
    return redistributed_table
