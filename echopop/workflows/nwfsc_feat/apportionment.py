import warnings
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import xarray as xr

from ... import utils

warnings.simplefilter("always")


def remove_group_from_estimates(
    transect_data: pd.DataFrame,
    group_proportions: xr.Dataset,
) -> pd.DataFrame:
    """
    Partition NASC, abundance (and number density), and biomass (and biomass density) transect
    values across indexed groups.

    Parameters
    ----------
    transect_data : pd.DataFrame
        DataFrame containing transect data with number densities already computed. Must include:
        - A column that matches the index of the `group_proportions` DataArrays.
    group_proportions : xr.Dataset
        xarray.Dataset containing partitioning data for NASC, abundance, and biomass.
        Valid variable names are limited to 'abundance', 'biomass', and 'nasc'.
        Each variable must be a 1D DataArray with dimensions matching grouping columns in
        `transect_data`.

    Returns
    -------
    pd.DataFrame
        DataFrame with the same structure as the input dataset, but with defined partitions applied
        to NASC and the other abundance/biomass columns.

    Notes
    -----
    The DataArrays in `group_proportions` must have dimensions that correspond to columns in the
    `transect_data`. Missing indices will result in NaN values for those rows. The function
    automatically identifies and uses the appropriate columns for merging.
    """

    # Create copy
    transect_data = transect_data.copy()

    # Empty Dataset not allowed
    if len(group_proportions.data_vars) == 0:
        raise ValueError("Input 'group_proportions' Dataset cannot be empty.")

    # Validate that there is a matching column with the Dataset coordinates
    # ---- Get coordinates
    grp_coords = list(group_proportions.coords.keys())
    if not set(grp_coords) <= set(transect_data.columns):
        # ---- Format error
        missing_coords = ", ".join(f"'{c}'" for c in grp_coords)
        raise KeyError(
            f"The following coordinates from 'group_proportions' could not be found in "
            f"'transect_data': {missing_coords}."
        )

    # Index the transect data based on the DataArray dimensions
    transect_data.set_index(grp_coords, inplace=True)
    transect_data_cnv = transect_data.to_xarray()

    # Align coordinates
    group_proportions_aligned = group_proportions.sel(
        {coord: transect_data_cnv[coord] for coord in grp_coords}
    )

    # Adjust NASC, if present; otherwise, drop to avoid partial evaluation
    if "nasc" in group_proportions_aligned:
        transect_data["nasc"] = transect_data_cnv["nasc"] * (1 - group_proportions_aligned["nasc"])
    else:
        transect_data.drop(columns=["nasc"], inplace=True)

    # Adjust number density and abundance, if present; otherwise, drop to avoid partial evaluation
    # ---- Get variables
    abundance_columns = [
        v for v in transect_data_cnv.data_vars if ("abundance" in v or "number_density" in v)
    ]
    if "abundance" in group_proportions_aligned:
        transect_data[abundance_columns] = (
            transect_data_cnv[abundance_columns] * (1 - group_proportions_aligned["abundance"])
        ).to_pandas()
    else:
        transect_data.drop(columns=abundance_columns, inplace=True)

    # Adjust biomass and biomass density, if present; otherwise, drop to avoid partial evaluation
    # ---- Get variables
    biomass_columns = [
        v for v in transect_data_cnv.data_vars if ("biomass" in v or "biomass_density" in v)
    ]
    if "biomass" in group_proportions_aligned:
        transect_data[biomass_columns] = (
            transect_data_cnv[biomass_columns] * (1 - group_proportions_aligned["biomass"])
        ).to_pandas()
    else:
        transect_data.drop(columns=biomass_columns, inplace=True)

    # Restore original indexing
    transect_data.reset_index(inplace=True)

    # Return the partitioned transect data
    return transect_data


def mesh_biomass_to_nasc(
    mesh_data: pd.DataFrame,
    biodata: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    mesh_biodata_link: Dict[str, str],
    stratum_weights: pd.DataFrame,
    stratum_sigma_bs: pd.DataFrame,
    group_columns: List[str] = [],
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
    if isinstance(biodata, xr.Dataset):
        biodata = {"": biodata}

    # Get the coordinates that will be used for setting shared indices
    coord_names = set()
    for ds in biodata.values():
        coord_names.update(ds.coords.keys())

    # Check
    # ---- Setdiff: biodata
    missing_biodata_columns = set(list(mesh_biodata_link.values())).difference(coord_names)
    if len(missing_biodata_columns) > 0:
        raise KeyError(
            f"The following link columns were missing from the biological data: "
            f"{', '.join(col for col in missing_biodata_columns)}."
        )
    # ---- Setdiff: mesh
    missing_mesh_columns = set(list(mesh_biodata_link.keys())).difference(mesh_data.columns)
    if len(missing_mesh_columns) > 0:
        raise KeyError(
            f"The following link columns were missing from the mesh DataFrame: "
            f"{', '.join(col for col in missing_mesh_columns)}."
        )
    # ---- Biomass column
    if "biomass" not in mesh_data.columns:
        raise KeyError(
            "The input kriging mesh DataFrame does not have a column 'biomass'. Please compute the "
            "kriged biomass."
        )

    # Set the link column names
    mesh_data.rename(columns=mesh_biodata_link, inplace=True)

    # Map the link columns
    mesh_index = list(mesh_biodata_link.values())

    # Index
    mesh_data.set_index(mesh_index, inplace=True)
    mesh_data.sort_index(inplace=True)

    # Stack the proportions to be as a function of group
    biodata_cnv = {
        k: ds["proportion_overall"]
        .sum(dim=[v for v in ds.coords.keys() if v not in group_columns])
        .to_pandas()
        for k, ds in biodata.items()
    }

    # Get the column names of the biodata
    biodata_columns = sorted(
        set(
            col
            for df in biodata_cnv.values()
            if isinstance(df, pd.DataFrame)
            for col in df.columns.tolist()
        )
    )
    # ---- Fill in if empty
    if len(biodata_columns) == 0:
        biodata_columns = set(["apportioned"])

    # Stack the proportions
    stacked_proportions = pd.concat(list(biodata_cnv.values()), axis=1)

    # If only a single expected output, summarize
    if len(biodata_columns) == 1:
        stacked_proportions = stacked_proportions.sum(axis=1).to_frame(list(biodata_columns)[0])

    # Multiply by biomass column, broadcasting over all columns
    stacked_proportions = stacked_proportions.mul(mesh_data["biomass"], axis=0).dropna()

    # Group columns by name and sum (to handle duplicate columns from stacking)
    stacked_proportions = stacked_proportions.T.groupby(stacked_proportions.columns).sum().T

    # Prefix columns
    stacked_proportions = stacked_proportions.add_prefix("biomass_")

    # Add to the mesh
    mesh_data[stacked_proportions.columns.tolist()] = stacked_proportions

    # Get the column names of the resulting biomass values
    biomass_columns = stacked_proportions.columns.tolist() + ["biomass"]

    # Reindex the stratum weights
    stratum_weights_idx = stratum_weights.to_series().reindex(mesh_data.index)

    # Calculate abundance
    mesh_data[[f"abundance_{name}" for name in biodata_columns] + ["abundance"]] = mesh_data[
        biomass_columns
    ].div(stratum_weights_idx, axis=0)

    # Index the stratum sigma_bs
    stratum_sigma_bs_idx = stratum_sigma_bs.reindex(mesh_data.index)

    # Calculate NASC
    mesh_data["nasc"] = mesh_data["abundance"] * stratum_sigma_bs_idx["sigma_bs"] * 4.0 * np.pi

    # Rename the aligned columns
    # ---- Create inverted link dictionary
    inverted_link = {v: k for k, v in mesh_biodata_link.items()}

    # Reset the index
    mesh_data.reset_index(inplace=True)

    # Rename
    mesh_data.rename(columns=inverted_link, inplace=True)


def impute_kriged_table(
    reference_table_df: pd.DataFrame,
    initial_table_df: pd.DataFrame,
    standardized_table_df: pd.DataFrame,
    group_by: List[str],
    impute_variable: List[str],
    table_reference_indices: List[str],
) -> pd.DataFrame:
    """
    Convert biomass estimates distributed across a grid or mesh into grouped (e.g. sex-specific)
    abundance and NASC values.

    This function takes kriged biomass density estimates from a mesh/grid and converts them to
    abundance and NASC (Nautical Area Scattering Coefficient) values by applying biological
    proportion data and acoustic backscattering coefficients. The function modifies the input
    mesh DataFrame in place.

    Parameters
    ----------
    mesh_data : pd.DataFrame
        DataFrame containing kriged biomass density estimates with spatial mesh information.
        Must contain columns for linking to biodata and either 'biomass', or 'biomass_density'
        and 'area', columns.
    biodata : Union[xr.Dataset, Dict[str, xr.Dataset]]
        Biological proportion data. Can be a single xarray.Dataset or a dictionary of Datasets
        containing proportion data for different biological groups (e.g., aged/unaged).
        Each Dataset must contain a variable named "proportion_overall".
    mesh_biodata_link : Dict[str, str]
        Dictionary mapping column names from mesh_data to biodata column names for joining.
        Example: {"geostratum_ks": "stratum_ks"}
    stratum_weights : pd.DataFrame
        DataFrame containing average weights per stratum for converting biomass to abundance.
        Index should match the linked columns from mesh_biodata_link.
    stratum_sigma_bs : pd.DataFrame
        DataFrame containing sigma_bs (backscattering coefficient) values per stratum.
        Must have 'sigma_bs' column with index matching linked columns.
    group_columns : List[str], default []
        List of column names to group biological data by (e.g., ["sex", "age_bin"]).

    Returns
    -------
    None
        Function modifies mesh_data in place, adding new columns for biomass, abundance,
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


def distribute_unaged_from_aged(
    population_table: xr.DataArray,
    reference_table: xr.DataArray,
    group_columns: List[str] = [],
    impute: bool = True,
    impute_variable: Optional[List[str]] = None,
) -> xr.DataArray:
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

    # Get original table coordinates
    table_coords = list(population_table.coords.keys())
    table_noncoords = list(set(table_coords) - set(group_columns + ["length_bin"]))

    # Get reference coordinates
    reference_coords = list(reference_table.coords.keys())

    # Aggregate the grouped dimensions
    grouped_table = (
        population_table.sum(
            dim=[v for v in population_table.coords.keys() if v not in set(table_coords)]
        )
        .to_dataframe(name="value")["value"]
        .unstack(table_noncoords)
    )

    # Aggregate the reference
    grouped_reference = (
        reference_table.sum(
            dim=[v for v in population_table.coords.keys() if v not in set(reference_coords)]
        )
        .to_dataframe(name="value")["value"]
        .unstack(table_noncoords)
    )

    # Get the shared index names across the reference tables
    reference_index_names = list(grouped_reference.index.names)

    # Get the indices for each reference table required for summation that must be unstacked
    reference_stack_indices = list(set(reference_index_names).difference(table_coords))

    # Reorient the columns accordingly
    population_table_reshape = grouped_table.sum(axis=1).unstack(group_columns)

    # Standardize the apportioned table values
    standardized_table = (
        population_table_reshape
        * grouped_reference.sum(axis=1).unstack(reference_stack_indices + group_columns)
        / grouped_reference.unstack(reference_stack_indices).sum(axis=1).unstack(group_columns)
    ).fillna(0.0)

    # If imputation is requested, perform it
    if not impute:
        result = standardized_table.unstack().to_xarray()
        if population_table.name:
            result.name = population_table.name
        return result
    else:
        # ---- Impute
        standardized_table_imputed = impute_kriged_table(
            reference_table_df=grouped_reference,
            initial_table_df=population_table_reshape,
            standardized_table_df=standardized_table,
            table_reference_indices=reference_stack_indices,
            group_by=group_columns,
            impute_variable=impute_variable,
        )
        result = standardized_table_imputed.unstack().to_xarray()
        if population_table.name:
            result.name = population_table.name
        return result


def sum_population_tables(
    population_table: Dict[str, xr.Dataset],
) -> xr.DataArray:
    """
    Combine and sum population estimates across defined tables to yield a single population table.

    This function takes a dictionary of xarray Datasets, each representing a population estimate
    table (e.g., "aged", "unaged"), and combines them into a single consolidated xarray DataArray by
    aggregating over non-grouped dimensions and summing values across all input tables.
    Parameters
    ----------
    population_table : Dict[str, xr.Dataset]
        Dictionary of population estimate tables to combine. Keys are table names (e.g., "aged",
        "unaged", "speciesA", "speciesB") and values are xarray Datasets with population data. Each
        Dataset should have compatible coordinate/dimension names.
    group_columns : List[str], optional
        List of coordinate names to retain as dimensions in the final combined table (e.g.,
        ["sex", "age_bin"]). All other coordinates will be summed over.

    Returns
    -------
    xr.DataArray
        Combined population table as an xarray DataArray, with group_columns as dimensions.
        Values are summed across all input tables for each group_columns combination.

    Notes
    -----
    The combination process:
    1. Identifies all unique coordinate names across the input Datasets.
    2. Determines which coordinates are not in group_columns and sums over those dimensions.
    3. Converts each aggregated Dataset to a DataFrame, aligns, and sums them element-wise.
    4. Returns the result as an xarray DataArray with group_columns as dimensions.

    Examples
    --------
    >>> combined = sum_population_table(
    ...     population_table={"aged": ds_aged, "unaged": ds_unaged},
    ...     group_columns=["sex", "age_bin", "length_bin"]
    ... )
    >>> print(combined)
    <xarray.DataArray ...>

    """

    # Get all unique coordinate names
    table_coords = {coord for ds in population_table.values() for coord in ds.coords.keys()}

    # Aggregate the grouped dimensions
    grouped_tables = {
        k: da.sum(dim=[v for v in da.coords.keys() if v not in set(table_coords)])
        .to_dataframe(name="value")["value"]
        .unstack()
        for k, da in population_table.items()
    }

    # Stack the tables
    stacked_tables = {
        k: v.stack(v.columns.names, future_stack=True).to_frame("value")
        for k, v in grouped_tables.items()
    }

    # Find common indices across all DataFrames
    common_indices = set.intersection(*(set(df.index.names) for df in stacked_tables.values()))

    # Collapse excess indices for all DataFrames (no explicit loop)
    def collapse_to_common(df):
        return df.groupby(level=list(common_indices), observed=True).sum()

    reduced_data = {k: collapse_to_common(df) for k, df in stacked_tables.items()}

    # Sum the compatible tables
    table_summed = sum(reduced_data.values())

    # Return the full DataArray
    result = table_summed["value"].to_xarray()
    result.name = "estimate"
    return result


def reallocate_excluded_estimates(
    population_table: xr.DataArray,
    exclusion_filter: Dict[str, Any],
    group_columns: List[str] = [],
) -> Union[xr.Dataset, xr.DataArray]:
    """
    Redistribute population estimates across groups after excluding a specific subset.

    This function removes specified population segments (e.g., age-1 fish) and redistributes
    their population estimates proportionally across the remaining groups. This is commonly used
    when certain age classes need to be excluded from final estimates but their contributions
    should be reallocated to other age groups.

    Parameters
    ----------
    population_table : xr.DataArray
        Population estimates as an xarray DataArray, typically with biological groups as dimensions.
    exclusion_filter : Dict[str, Any]
        Dictionary specifying which population segments to exclude and redistribute. Keys are
        dimension names, values are the categories to exclude. Example: {"age_bin": [1]} to
        exclude age-1 fish.
    group_columns : List[str], optional
        List of dimension names that define the grouping variables for redistribution (e.g., ["sex"]
        to redistribute within each sex separately). If empty, redistribution is performed over all
        dimensions.

    Returns
    -------
    Union[xr.Dataset, xr.DataArray]
        If `group_columns` is not empty, returns an xarray.Dataset with excluded segments removed
        and their values redistributed proportionally across remaining groups, preserving the
        'group_columns' as dimensions. If `group_columns` is empty, returns an xarray.DataArray
        with the same structure as the input, but with excluded segments removed and their values
        redistributed.

    Raises
    ------
    UserWarning
        If the redistributed table sums don't match the original table sums within tolerance.

    Notes
    -----
    The redistribution algorithm:
    1. Identifies population segments matching the exclusion filter.
    2. Calculates the total excluded biomass/abundance for each group.
    3. Removes excluded segments (sets them to zero).
    4. Redistributes excluded totals proportionally across remaining segments.
    5. Ensures total population remains constant after redistribution.

    This maintains the biological realism of population estimates while allowing for policy-based
    exclusions (e.g., removing age-1 fish from assessments).

    Examples
    --------
    >>> result = reallocate_excluded_estimates(
    ...     population_table=da_population,
    ...     exclusion_filter={"age_bin": [1]},
    ...     group_columns=["sex", "age_bin"]
    ... )
    >>> print(result)
    <xarray.Dataset ...>

    >>> result = reallocate_excluded_estimates(
    ...     population_table=da_population,
    ...     exclusion_filter={"age_bin": [1]},
    ...     group_columns=[]
    ... )
    >>> print(result)
    <xarray.DataArray ...>
    """

    # If no appropriate filter is defined, then nothing is redistributed
    if len(exclusion_filter) == 0:
        return population_table

    # Convert to DataFrame
    population_est = population_table.to_series().unstack(group_columns)

    # Apply inverse of exclusion filter to get the values being excluded
    population_excluded = utils.apply_filters(population_est, include_filter=exclusion_filter)

    # Replace the excluded values in the full table with 0.
    population_masked = utils.apply_filters(
        population_est, exclude_filter=exclusion_filter, replace_value=0.0
    )

    # Get the sums for each group across the excluded and filtered tables
    # ---- Filtered/included
    population_masked_sum = population_masked.sum()
    # ---- Excluded
    population_excluded_sum = population_excluded.sum()
    # ---- Handling if `pandas.Series`
    if isinstance(population_excluded_sum, pd.Series):
        population_excluded_sum = population_excluded_sum.reindex(
            index=population_masked_sum.index, method="ffill"
        )

    # Get the redistributed values that will be added to the filtered table values
    population_adjusted = (
        population_masked * population_excluded_sum / population_masked_sum
    ).fillna(0.0)

    # Add the adjustments to the masked table
    population_masked += population_adjusted

    # Unstack and return to DataArray
    result = population_masked.to_xarray().squeeze()
    if isinstance(result, xr.DataArray) and not result.name:
        result.name = population_table.name
    return result


def distribute_population_estimates(
    data: pd.DataFrame,
    proportions: Union[Dict[str, xr.Dataset], xr.Dataset],
    variable: str,
    group_columns: List[str] = [],
    data_proportions_link: Optional[Dict[str, str]] = None,
) -> Union[xr.DataArray, Dict[str, xr.DataArray]]:
    """
    Distribute population estimates (e.g. abundance, biomass) using proportions grouped by metrics
    such as age and length bins.

    This function takes population estimates and distributes them across different biological
    categories (e.g. sex, age, length) using calculated proportions. The distribution is performed
    by multiplying the population estimates with the corresponding proportions.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing population estimates. Must contain the columns specified by
        'variable' and 'group_columns'. If 'data_proportions_link' is defined, then the
        key value must also be present in the DataFrame.
    proportions : Union[Dict[str, xr.Dataset], xr.Dataset]
        Proportion data for distributing the population estimates. Can be a single xarray.Dataset
        or a dictionary of Datasets (e.g., {"aged": ds_aged, "unaged": ds_unaged}).
        Each Dataset must contain a variable named "proportion_overall".
    variable : str
        Name of the column in data containing the values to distribute (e.g. "biomass",
        "abundance").
    group_columns : List[str]
        List of column names that define the biological groups for distribution and stratification
        (e.g. ["sex", "age_bin", "length_bin", "stratum_ks"]).
    data_proportions_link : Optional[Dict[str, str]], default None
        Dictionary mapping column names from 'data' to those in 'proportions'. For
        instance, the dictionary `{'stratum_A': 'stratum_B'}` links 'stratum_A' within
        'data' with 'stratum_B' in the 'proportions' Dataset(s).

    Returns
    -------
    Union[xr.DataArray, Dict[str, xr.DataArray]]
        Distributed estimates. If a single group is present, returns an xarray.DataArray.
        Otherwise, returns a dictionary of xarray.DataArray objects, each containing the original
        biological group structure with values distributed according to proportions.

    Notes
    -----
    The DataArrays in `proportions` must have dimensions that correspond to columns in the
    `data`. Missing indices will result in NaN values for those rows. The function
    automatically identifies and uses the appropriate columns for merging and grouping.
    """

    # Create copy
    data = data.copy()

    # Type-check and normalize proportions input
    if isinstance(proportions, xr.Dataset):
        proportions = {"data": proportions}

    # Verify that `variable` exists
    if variable not in data.columns:
        raise ValueError(f"Variable '{variable}' missing from `data`.")

    # Rename stratum naming, if defined
    if data_proportions_link is not None:
        data_stratum = next(iter(data_proportions_link))
        if data_stratum not in data.columns:
            raise KeyError(f"Stratum '{data_stratum}' missing from `data`.")
        data.rename(columns=data_proportions_link, inplace=True)

    # Validate group_columns overlapping with data and proportions
    # ---- proportions
    proportions_group_columns = set()
    for ds in proportions.values():
        proportions_group_columns.update(ds.coords.keys())
    # ---- Only require group_columns to be present in at least one of the sources
    missing = [
        col
        for col in group_columns
        if col not in data.columns and col not in proportions_group_columns
    ]
    if missing:
        missing_str = ", ".join(f"'{c}'" for c in missing)
        raise KeyError(
            f"The following columns from 'group_columns' could not be found in either "
            f"'data' or 'proportions': {missing_str}."
        )

    # Sum data variable over indices
    data_group_columns = [col for col in group_columns if col in data.columns]
    data_array = data.pivot_table(
        index=data_group_columns, values=variable, aggfunc="sum", observed=False
    ).to_xarray()

    # Parse the additional columns that are required for grouping
    proportions_group_columns = {
        k: [c for c in ds.coords.keys() if c in group_columns] for k, ds in proportions.items()
    }

    # Calculate the proportions
    # ---- Totals (for normalization)
    group_totals = sum(
        ds["proportion_overall"].sum(
            dim=[c for c in proportions_group_columns[k] if c not in data_group_columns]
        )
        for k, ds in proportions.items()
    )
    # ---- Apply normalization
    proportions_norm = {
        k: (ds["proportion_overall"] / group_totals).fillna(0.0) for k, ds in proportions.items()
    }

    # Distribute the variable over each table
    apportioned_groups = {
        k: proportions_norm[k] * data_array[variable] for k in proportions_norm.keys()
    }
    # ---- Update DataArray names
    for arr in apportioned_groups.values():
        arr.name = variable

    # Return
    if len(apportioned_groups) == 1:
        return next(iter(apportioned_groups.values()))
    else:
        return apportioned_groups
