import numpy as np
import pandas as pd
from typing import Any, Dict, List, Union
from functools import reduce
from . import utils

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

def standardize_kriged_estimates(
    population_table: Dict[str, pd.DataFrame],
    reference_table: str,
    group_by: List[str],
) -> Union[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Pass
    """

    # Get the table indices
    table_indices = {
        k: list(df.index.names)
        for k, df in population_table.items()
    }    
    
    # Get the shared index names across the reference tables
    reference_index_names = list(
        set.intersection(*(set(v) 
                        for k, v in table_indices.items() 
                        if k in reference_table)) 
    )

    # Get the indices for each reference table required for summation that must be unstacked
    reference_stack_indices = {
        k: list(set(reference_index_names).difference(set(df.index.names)))
        for k, df in population_table.items()
        if k in reference_table
    }
    
    # Consolidate the reference tables into a list
    reference_list = [
        population_table[k].unstack(reference_stack_indices[k]) for k in reference_table
    ]

    # Sum across all of the references
    reference_df = reduce(lambda a, b: a.add(b, fill_value=0), reference_list)

    # Get the indices for each target table
    target_stack_indices = {
        k: list(set(reference_index_names).difference(set(df.index.names)))
        for k, df in population_table.items()
        if k not in reference_table    
    }

    # Standardize the proportions
    proportions_standardized = {
        k: (
            df.sum(axis=1).unstack(group_by) *
            reference_df.sum(axis=1).unstack(target_stack_indices[k] + group_by) /
            reference_df.unstack(target_stack_indices[k]).sum(axis=1).unstack(group_by)
        ).fillna(0.)
        for k, df in population_table.items() if k not in reference_table
    }

    # Return the updated proportions
    return proportions_standardized

    

    

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
