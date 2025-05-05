from typing import Dict, Union, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

from echopop.nwfsc_feat import get_proportions, ingest_nasc, load_data
from echopop.kriging import Kriging
from echopop.inversion import InversionLengthTS


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

dict_df_bio = load_data.load_biological_data(root_path, bio_path_dict, species_code)
dict_df_strata = load_data.load_stratification(root_path, strata_path_dict)

# Consolidate all input data into df_acoustic_dict
df_nasc_no_age1 = load_data.consolidate_all_data(
    df_nasc=df_nasc_no_age1,
    df_bio_dict=dict_df_bio,
    df_strata_dict=dict_df_strata
)




# ===========================================
# Compute biological composition based on stratum
length_bins: np.array  # length bin specification
df_length_weight, df_regression = get_proportions.length_weight_regression(dict_df_bio["specimen"], length_bins)
# df_regression seems unused afterwards -- good as a record?

# Get counts ----------------
df_aged_counts = get_proportions.fish_count(  # previously "aged_number_distribution"
    df_specimen=dict_df_bio["specimen"],
    df_length=dict_df_bio["length"],
    aged=True,
    sexed=True,
)
df_unaged_counts = get_proportions.fish_count(  # previously "unaged_number_distribution"
    df_specimen=dict_df_bio["specimen"],
    df_length=dict_df_bio["length"],
    aged=False,
    sexed=True,
)
# Previously there was also "aged_number_distribution_filtered" 
# but it is simply df_aged_counts with unsexed fish removed, 
# I think it is better to have that explicitly in the code, 
# so removed OUTSIDE of the get_fish_count function


# Get number proportions ----------------

# in the output dataframes: *_overall = *_aged + *_unaged
# only handle 1 species at a time
dict_df_number_proportion: Dict[pd.DataFrame] = get_proportions.number_proportions(
    df_aged_counts,
    df_unaged_counts,
)



# Get weight proportions ----------------
dict_df_weight_distr: Dict[pd.DataFrame]

# aged fish - weight distribution over sex/length/age
dict_df_weight_distr["aged"] = get_proportions.weight_distributions_over_lenghth_age(
    df_specimen=dict_df_bio["specimen"],
    df_length=dict_df_bio["length"],
    df_length_weight=df_length_weight,
    aged=True,
)

# unaged fish - weight distribution over sex/length
dict_df_weight_distr["unaged"] = get_proportions.weight_distributions_over_lenghth_age(
    df_specimen=dict_df_bio["specimen"],
    df_length=dict_df_bio["length"],
    df_length_weight=df_length_weight,
    aged=False,
)

# Get averaged weight for all sex, male, female for all strata
df_averaged_weight = get_proportions.stratum_averaged_weight(
    df_length_weight=df_length_weight,
    dict_df_bio=dict_df_bio,  # use "specimen" and "length"
)

# Get weight proportion for all sex, male, female for all strata
dict_df_weight_proportion: Dict[pd.DataFrame] = get_proportions.weight_proportions(
    df_catch=dict_df_bio["catch"],
    dict_df_weight_proportion=dict_df_weight_distr,  # weight proportions
    df_length_weight=df_length_weight,  # length-weight regression
)





# ===========================================
# NASC to number density

# Inititate object to perform inversion
# inversion parameters are stored as object attributes
invert_hake = InversionLengthTS(df_model_params=dict_df_bio["model_params"])

# Perform inversion using the supplied df_length
# df_length will be used to compute the mean sigma_bs for each stratum,
# which is then used in .invert() to compute number density on a stratum-by-stratum basis
df_nasc_no_age1 = invert_hake.invert(df_nasc=df_nasc_no_age1, df_length=dict_df_bio["length"])





# ===========================================
# Perform kriging using class Kriging
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
df_nasc_no_age1_kriged = kriging.krige(
    df_in=df_nasc_no_age1,
    variables="biomass"
)
