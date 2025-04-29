from typing import Dict, Union, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

from echopop.nwfsc_feat import get_proportions, ingest_nasc, load_data
from echopop.kriging import Kriging


# ===========================================
# Organize NASC file
nasc_path = "SOME_PATH"
nasc_filename_pattern = "SOME_PATTERN"
region_class_filepath = "SOME_PATH"  # pattern-label mapping under transect_region_mapping/parts
survey_identifier = "YEAR_MONTH"

df_merged = ingest_nasc.merge_echoview_nasc(
    nasc_path, nasc_filename_pattern)

# Optional: only use for years needing this as external resources
df_transect_region_key = ingest_nasc.load_transect_region_key(region_class_filepath)

# Use df.to_csv to save df_transect_region_key, in place of the specialized transect_region_key file
# Keep read_transect_region_file and make sure its output is the same as construct_transect_region_key


# Age-1+
df_nasc_all_ages = ingest_nasc.consolidate_echoview_nasc(
    df_merged,
    region_names=["Age-1 Hake", "Age-1 Hake Mix", "Hake", "Hake Mix"],
    survey_identifier=survey_identifier,
)

# Age-2+ (no age 1)
df_nasc_no_age1 = ingest_nasc.consolidate_echoview_nasc(
    df_merged,
    region_names=["Hake", "Hake Mix"],
    survey_identifier=survey_identifier,
)

# Use df.to_csv to save df_nasc_all_ages and df_nasc_no_age1 if needed

# Use regular pd.read_csv to read df_nasc_*, effectively break up the current load_data()
# -- there is no need to have a one-size-fits-all load_data function
# -- just read them in without validation is fine: these are all files under our control




# ===========================================
# Execute what's in Survey.load_survey_data()
# All *_dict below are a subdict from the original config yaml

root_path = "WHERE_ALL_DATA_ARE"
species_code = "SPECIES_CODE"
df_nasc_no_age1: pd.DataFrame  # extracted nasc data from above, can also be df_nasc_all_ages

bio_path_dict: dict  # the "biological" section of year_config.yml
                     # this will be simplified now that we read from the master spreadsheet
strata_path_dict: dict  # the "stratification" section of year_config.yml

df_bio_dict = load_data.load_biological_data(root_path, bio_path_dict, species_code)
df_strata_dict = load_data.load_stratification(root_path, strata_path_dict)

# Consolidate all input data into df_acoustic_dict
df_nasc_no_age1 = load_data.consolidate_all_data(
    df_nasc=df_nasc_no_age1,
    df_bio_dict=df_bio_dict,
    df_strata_dict=df_strata_dict
)




# ===========================================
# Compute biological composition based on stratum
length_bins: np.array  # length bin specification
df_length_weight, df_regression = get_proportions.length_weight_regression(df_bio_dict["specimen"], length_bins)
# df_regression seems unused afterwards -- good as a record?

# Get counts ----------------
df_aged_counts = get_proportions.fish_count(  # previously "aged_number_distribution"
    df_specimen=df_bio_dict["specimen"],
    df_length=df_bio_dict["length"],
    aged=True,
    sexed=True,
)
df_unaged_counts = get_proportions.fish_count(  # previously "unaged_number_distribution"
    df_specimen=df_bio_dict["specimen"],
    df_length=df_bio_dict["length"],
    aged=False,
    sexed=True,
)
# Previously there was also "aged_number_distribution_filtered" 
# but it is simply df_aged_counts with unsexed fish removed, 
# I think it is better to have that explicitly in the code, 
# so removed OUTSIDE of the get_fish_count function


# Get number proportions ----------------
# TODO: DISCUSS THIS!
# TODO: what does _overall stand for?
da_number_proportion = get_proportions.number_proportion()


# Get weight proportions ----------------
# aged fish - weight distribution
da_sex_length_age: xr.DataArray = get_proportions.weight_distributions(
    df_specimen=df_bio_dict["specimen"],
    df_length=df_bio_dict["length"],
    df_length_weight=df_length_weight,
    aged=True,
)

# unaged fish - weight distribution
da_sex_length: xr.DataArray = get_proportions.weight_distributions(
    df_specimen=df_bio_dict["specimen"],
    df_length=df_bio_dict["length"],
    df_length_weight=df_length_weight,
    aged=False,
)

# Get stratum averaged weight for all sex, male, female
df_averaged_weight = get_proportions.stratum_averaged_weight()


# Get weight proportions ----------------
# TODO: DISCUSS THIS!
# TODO: what does _overall stand for?
da_weight_proportion = get_proportions.weight_proportion()





# ===========================================
# NASC to number density






# ===========================================
# Perform kriging using class Kriging
# TODO:
# put back FEAT-specific kriging files

# Load kriging-related params
kriging_const: dict  # from initalization_config.yaml:
                     # A0, longitude_reference, longitude/latitude_offset
kriging_path_dict: dict  # the "kriging" section of year_config.yml
                         # combined with the "kriging" section of init_config.yml
kriging_param_dict, variogram_param_dict = load_data.load_kriging_variogram_params(
    root_path=root_path,
    file_path_dict=kriging_path_dict,
    kriging_const=kriging_const,
)

kriging = Kriging(
    kriging_param_dict=kriging_param_dict,
    variogram_param_dict=variogram_param_dict,
    mesh_template="PATH_TO_MESH_TEMPLATE",
    isobath_template="PATH_TO_ISOBATH_REFERENCE",
)

# Create kriging mesh including cropping based on transects
# Created mesh is stored in kriging.df_mesh
kriging.create_mesh()

# Perform coordinate transformation based on isobath if needed
# This adds columns x/y to kriging.df_mesh
kriging.latlon_to_xy()

# Perform kriging
# This adds kriging result columns to df_in
df_out = kriging.krige(
    df_in=df_nasc_no_age1,
    variables="biomass"
)
