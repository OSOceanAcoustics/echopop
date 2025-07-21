import numpy as np
import pandas as pd
from typing import Any, Dict, List, Optional, Union
import warnings

from . import utils

warnings.simplefilter("always")

def mesh_biomass_to_nasc(
    mesh_data_df: pd.DataFrame,
    biodata: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    mesh_biodata_link: Dict[str, str],
    stratum_weights_df: pd.DataFrame,
    stratum_sigma_bs_df: pd.DataFrame,
    group_by: List[str] = [],
) -> None:
    """
    Convert kriged biomass density or biomass estimates distributed across a grid or mesh into 
    grouped abundance and NASC
    """

    # Convert to a Dictionary if needed
    if isinstance(biodata, pd.DataFrame):
        biodata = {"": biodata}

    # Get the column indices that will be used for setting shared indices
    column_names = list(set(
        col for df in biodata.values() if isinstance(df.columns, pd.Index) 
        for col in df.columns.names
    ))

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
    biodata_columns = set(col for df in biodata_stk.values() if isinstance(df, pd.DataFrame) 
                        for col in df.columns.tolist())
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
    stacked_proportions = stacked_proportions.add_prefix(f"biomass_")

    # Add to the mesh DataFrame
    mesh_data_df[stacked_proportions.columns.tolist()] = stacked_proportions

    # Get the column names of the resulting biomass values
    biomass_columns = stacked_proportions.columns.tolist() + ["biomass"]

    # Reindex the stratum weights
    stratum_weights_df_idx = stratum_weights_df.reindex(mesh_data_df.index)

    # Calculate abundance
    mesh_data_df[["abundance"] + [f"abundance_{name}" for name in biodata_columns]] = (
        mesh_data_df[biomass_columns].div(stratum_weights_df_idx, axis=0)
    )

    # Index the stratum sigma_bs
    stratum_sigma_bs_df_idx = stratum_sigma_bs_df.reindex(mesh_data_df.index)

    # Calculate NASC
    mesh_data_df["nasc"] = (
        mesh_data_df["abundance"] * stratum_sigma_bs_df_idx["sigma_bs"] * 4. * np.pi
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
) -> Union[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Pass
    """

    # Convert to compatible format
    is_dict = True
    if isinstance(proportions, pd.DataFrame):
        proportions = {"data": proportions}
        is_dict = False

    # Sum variables over indices
    mesh_data_pvt = mesh_data_df.rename(columns=mesh_proportions_link).pivot_table(
        index=stratify_by, 
        values=variable, 
        aggfunc="sum",
        observed=False
    )

    # Parse the additional columns that are required for grouping
    proportions_group_columns = {
            k: [c for c in (list(v.index.names) + list(v.columns)) 
                if c in group_by]
            for k, v in proportions.items()
    }

    # Convert to DataFrame(s) to pivot table(s)
    proportions_grouped_pvt = {
        k: df if utils.is_pivot_table(df)
        else utils.create_pivot_table(
            df, 
            index_cols=proportions_group_columns[k],
            strat_cols=stratify_by,
            value_col="proportion_overall"
        )
        for k, df in proportions.items()
    }

    # Distribute the variable over each table
    apportioned_grouped_pvt = {
        k: df.mul(mesh_data_pvt[variable]).fillna(0.)
        for k, df in proportions_grouped_pvt.items()
    }

    # Return output in the same format as the proportions input
    if is_dict:
        return apportioned_grouped_pvt
    else:
        return apportioned_grouped_pvt["data"]

def impute_kriged_table(
    reference_table_df: pd.DataFrame,
    initial_table_df: pd.DataFrame,
    standardized_table_df: pd.DataFrame,
    group_by: List[str],
    impute_variable: List[str],
    table_reference_indices: List[str],
) -> pd.DataFrame:
    """
    Nearest-neighbor imputation of missing values in population apportionment tables
    """

    # Sum the reference proportions across the impute variable
    reference_stk = (
        reference_table_df.sum(axis=1).unstack(impute_variable).sum(axis=1).unstack(group_by)
    )

    # Unstack the grouping and imputate variable
    reference_group_stk = reference_table_df.sum(axis=1).unstack(group_by + impute_variable)

    # Get the mask for all 0.0's and non-zeros
    ref_zero_mask = reference_stk == 0.
    ref_nonzero_mask = reference_stk != 0.

    # Gather the indices for each column
    ref_zero_indices = {col: reference_stk.index[ref_zero_mask[col]].tolist()
                        for col in reference_stk.columns}
    ref_nonzero_indices = {col: reference_stk.index[ref_nonzero_mask[col]].tolist() 
                        for col in reference_stk.columns}

    # Create translation for row numbers
    interval_to_numeric = {interval: i for i, interval in enumerate(reference_stk.index)}  

    # Convert to row numbers
    ref_zero_rows = {col: np.array([interval_to_numeric[ival] for ival in intervals])
                    for col, intervals in ref_zero_indices.items()}
    ref_nonzero_rows = {col: np.array([interval_to_numeric[ival] for ival in intervals])
                    for col, intervals in ref_nonzero_indices.items()}
    
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
        col: initial_table_df[col].loc[ref_zero_indices[col]] != 0.
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
        standardized_table_copy.columns = standardized_table_copy.columns.reorder_levels(column_indices_rev)   

    # Get the nearest-neighbor rows and recompute the indices
    imputed_rows = {
        col: arr[
            np.argmin(
                np.abs(ref_zero_rows[col][table_nonzeros_mask[col]][:, np.newaxis] - ref_nonzero_rows[col]), 
            axis=1)
        ]
        for col, arr in ref_nonzero_rows.items()
    }

    # Impute to replace these values
    imputed_values = {
        col: (
            initial_table_df[col].loc[nonzero_reference_to_table_indices[col]].to_numpy() * 
            reference_group_stk.iloc[imputed_rows[col]][col].T / 
            reference_stk.iloc[imputed_rows[col]][col]
        ).T
        for col in reference_stk.columns
    }

    # Update the standardized values
    for col in initial_table_df.columns:
        if col in reference_stk.columns:
            standardized_table_copy.iloc[
                nonzero_reference_to_table_rows[col],
                standardized_table_copy.columns.get_loc(col)
            ] = (
                imputed_values[col]
            )
            
    # Return to the original column index order
    standardized_table_copy = standardized_table_copy.reorder_levels(column_indices, axis=1)
        
    # Return the imputed standardized table
    return standardized_table_copy

def standardize_kriged_estimates(
    population_table: Dict[str, pd.DataFrame],
    reference_table: str,
    group_by: List[str],
    impute: bool = True,
    impute_variable: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Pass
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
        population_table_reshape *
        reference_table.sum(axis=1).unstack(reference_stack_indices + group_by) /
        reference_table.unstack(reference_stack_indices).sum(axis=1).unstack(group_by)
    ).fillna(0.)

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
    Combine and sum population estimates across defined tables to yield a single population table
    """
    
    # Subset the table to only include the tables-of-interest
    tables_to_combine = {
        k: v.stack(v.columns.names, future_stack=True).to_frame("value") 
        for k, v in population_table.items() if k in table_names
    }
    
    # Create a pivot table for each table with identical indices
    compatible_tables = {
        k: utils.create_pivot_table(v, index_cols=table_index, 
                                    strat_cols=table_columns, 
                                    value_col="value")
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
    Redistribute population estimates across groups after excluding a specific subset 
    """
    
    # Find any columns that are not in the group_by list
    # ---- Get column names
    column_names = population_table.columns.names
    # ---- Identify extra columns that are not in the group_by list
    extra_columns = [col for col in column_names if col not in group_by]
    # ---- Stack the population table
    stacked_table = population_table.stack(extra_columns, future_stack=True)

    # Apply inverse of exclusion filter to get the values being excluded
    excluded_grouped_table = utils.apply_filters(stacked_table, 
                                                include_filter=exclusion_filter)

    # Replace the excluded values in the full table with 0.
    filtered_grouped_table = utils.apply_filters(stacked_table, 
                                                exclude_filter=exclusion_filter,
                                                replace_value=0.)
    
    # Get the sums for each group across the excluded and filtered tables
    # ---- Excluded
    excluded_grouped_sum = excluded_grouped_table.sum()
    # ---- Filtered/included
    filtered_grouped_sum = filtered_grouped_table.sum()

    # Get the redistributed values that will be added to the filtered table values
    adjustment_table = filtered_grouped_table * excluded_grouped_sum / filtered_grouped_sum

    # Add the adjustments to the filtered table
    filtered_grouped_table += adjustment_table
    
    # Check 
    if np.any(filtered_grouped_table.sum() - stacked_table.sum() > 1e-6):
        # ---- If the sums do not match, raise a warning
        check_sums = filtered_grouped_table.sum() - stacked_table.sum() > 1e-6
        # ---- Raise a warning with the indices where the sums do not match
        warnings.warn(
            f"The sums of the table with the redistributed estimates do not match the original table "
            f"filtered table do not match the original table for indices: "
            f"{', '.join(check_sums[check_sums].index.tolist())}"
        )
        
    # Restore the original column structure
    redistributed_table = (
        filtered_grouped_table.unstack(extra_columns)
        .reorder_levels(column_names, axis=1)
    )
    
    # Return the the redistributed table
    return redistributed_table

    
    
    

# # Overlap with the current `partition_transect_age` function
# def apportion_transect_biomass_abundance(
#     df_nasc: pd.DataFrame,
#     ds_proportions: xr.Dataset,
# ) -> xr.Dataset:
#     """
#     Apportion transect biomass across sex, age_bin, and length bin.

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame containing the apportioned biomass across sex, age_bin, and length bin.

#         TODO: sketch of ds_transect_apportioned
#               this is similar to what you (BL) have in `adult_data`
#               the dimension `location` below are the transect interval locations
#         - dimensions: location, sex, length_bin, age_bin
#         - coorindates: location, sex, length_bin, age_bin
#         - variables:
#            - stratum (location)
#            - transect (location)
#            - latitude (location)
#            - longitude (location)
#            - fraction_hake (location)
#            - biomass_aged (location, sex, length_bin, age_bin)
#            - biomass_unaged (location, sex, length_bin)
#            - abundance_aged (location, sex, length_bin, age_bin)
#            - abundance_unaged (location, sex, length_bin)

#         NOTE: from ds_transect_apportioned you can easily derive
#               the current `biomass_summary_df`
#         NOTE: the current `abundance_unaged_age1_tbl` should be part of `ds_proportions`
#     """
#     ds_transect_apportioned: xr.Dataset
#     return ds_transect_apportioned


# # The biomass parts of the current `apportion_kriged_values` function
# def apportion_kriged_biomass(
#     df_nasc: pd.DataFrame,
#     ds_proportions: xr.Dataset,
# ) -> xr.Dataset:
#     """
#     Apportion kriged biomass across sex, age_bin, and length bin.

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame containing the apportioned biomass across sex, age_bin, and length bin.

#         TODO: sketch of ds_kriged_apportioned
#         - dimensions: sex, length_bin, age_bin
#         - coorindates: sex, length_bin, age_bin
#         - variables:
#            - biomass_aged (stratum, sex, length_bin, age_bin)
#            - biomass_unaged (sex, length_bin)

#     NOTE: Wouldn't it be possible to apportion kriged biomass on a grid-by-grid basis?
#           This way we can have very meaningful maps.
#           The xr.Dataset structure would look like:
#           - dimensions: x, y, sex, length_bin, age_bin
#           - coorindates: lon, lat, sex, length_bin, age_bin
#           - variables:
#              - stratum (lat, lon)
#              - biomass_aged (lat, lon, stratum, sex, length_bin, age_bin)
#              - biomass_unaged (lat, lon, sex, length_bin)

#     """
#     ds_kriged_apportioned: xr.Dataset
#     return ds_kriged_apportioned


# # The current `impute_kriged_values` function
# def fill_missing_aged_from_unaged(
#     ds_kriged_apportioned: xr.Dataset,
#     ds_proportions: xr.Dataset,
# ) -> xr.Dataset:
#     """
#     Fill missing length bins in the aged dataset using unaged data.
#     """
#     pass


# # The current section in biology.py that starts with comment:
# # "# Additional reapportionment if age-1 fish are excluded"
# def reallocate_age1(
#     ds_kriged_apportioned: xr.Dataset,
#     ds_proportions: xr.Dataset,
# ) -> xr.Dataset:
#     """
#     Reallocate age-1 biomass to age-2+ fish.
#     """
#     pass


# # The abundance parts of the current `apportion_kriged_values` function
# def back_calculate_kriged_abundance(
#     ds_kriged_apportioned: xr.Dataset,
#     ds_proportions: xr.Dataset,
# ) -> xr.Dataset:
#     """
#     Back-calculate kriged abundance from apportioned biomass across sex, age_bin, and length bin.

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame containing the apportioned biomass across sex, age_bin, and length bin.

#         TODO: sketch of ds_kriged_apportioned
#         - dimensions: sex, length_bin, age_bin
#         - coorindates: sex, length_bin, age_bin
#         - variables:
#            - biomass_aged (stratum, sex, length_bin, age_bin)
#            - biomass_unaged (sex, length_bin)
#            - abundance_aged (stratum, sex, length_bin, age_bin) -- added in this function
#            - abundance_unaged (sex, length_bin) -- added in this function
#     """
#     ds_kriged_apportioned: xr.Dataset
#     return ds_kriged_apportioned
