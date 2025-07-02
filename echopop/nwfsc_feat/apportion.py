# import pandas as pd
# import xarray as xr


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
