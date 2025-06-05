from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import xarray as xr

from echopop.inversion import InversionLengthTS
from echopop.kriging import Kriging
from echopop.nwfsc_feat import apportion, get_proportions, ingest_nasc, load_data

# ==================================================================================================
# ==================================================================================================
# Organize NASC file
# ------------------
nasc_path: Path = Path("C:/Users/Brandyn/Documents/GitHub/Data/raw_nasc/")
filename_transect_pattern: str = r"T(\d+)"
default_transect_spacing: float = 10.0  # nmi
default_transect_spacing_latitude: float = 60.0  # deg N

# Outputs:
# ---- 1. Complete intervals regardless of hake region presence/absence [interval_df]
# ---- 2. DataFrame containing the full set of cells-intervals-layers [exports_df]
interval_df, exports_df = ingest_nasc.merge_echoview_nasc(
    nasc_path,
    filename_transect_pattern,
    default_transect_spacing,
    default_transect_spacing_latitude,
)

# ==================================================================================================
# Read in transect-region-haul keys
# ---------------------------------
transect_region_filepath_all_ages: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/Data/Stratification/"
    "US_CAN_2019_transect_region_haul_age1+ auto_final.xlsx"
)
transect_region_sheetname_all_ages: str = "Sheet1"
transect_region_filepath_no_age1: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/Data/Stratification/"
    "US_CAN_2019_transect_region_haul_age2+ auto_20191205.xlsx"
)
transect_region_sheetname_no_age1: str = "Sheet1"
transect_region_file_rename: dict = {
    "tranect": "transect_num",
    "region id": "region_id",
    "trawl #": "haul_num",
}

# Read in the transect-region-haul key files for each group
# Outputs:
# ---- DataFrame containing columns for `transect_num` [float], `region_id` [float], 
# ---- `haul_num` [float]
transect_region_haul_key_all_ages: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
    transect_region_filepath_all_ages,
    transect_region_sheetname_all_ages,
    transect_region_file_rename,
)

transect_region_haul_key_no_age1: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
    transect_region_filepath_no_age1, transect_region_sheetname_no_age1, transect_region_file_rename
)

# ==================================================================================================
# Read in transect-region-haul keys
# ---------------------------------
CAN_haul_offset: float = 200
region_name_expr_dict: Dict[str, dict] = {
    "REGION_CLASS": {
        "Age-1 Hake": "^(?:h1a(?![a-z]|m))",
        "Age-1 Hake Mix": "^(?:h1am(?![a-z]|1a))",
        "Hake": "^(?:h(?![a-z]|1a)|hake(?![_]))",
        "Hake Mix": "^(?:hm(?![a-z]|1a)|hake_mix(?![_]))",
    },
    "HAUL_NUM": {
        "[0-9]+",
    },
    "COUNTRY": {
        "CAN": "^[cC]",
        "US": "^[uU]",
    },
}

# Process the region name codes to define the region classes
# e.g. H5C - Region 2 corresponds to "Hake, Haul #5, Canada"
# Outputs:
# ---- `exports_df` with appended columns representing the updated haul number, region class, and 
# ---- region name
exports_with_regions_df: pd.DataFrame = ingest_nasc.process_region_names(
    exports_df,
    region_name_expr_dict,
    CAN_haul_offset,
)

# ==================================================================================================
# [OPTIONAL] Generate transect-region-haul key from compiled values
# ---------------------------------
region_list_no_age1: List[str] = ["Hake", "Hake Mix"]
region_list_all_ages: List[str] = ["Age-1 Hake", "Age-1", "Hake", "Hake Mix"]

# Outputs:
# ---- DataFrame containing columns for `transect_num` [float], `region_id` [float], 
# ---- `haul_num` [float]
transect_region_haul_key_no_age1: pd.DataFrame = ingest_nasc.generate_transect_region_haul_key(
    exports_with_regions_df, filter_list=region_list_no_age1
)

transect_region_haul_key_all_ages = ingest_nasc.generate_transect_region_haul_key(
    exports_with_regions_df, filter_list=region_list_all_ages
)

# ==================================================================================================
# Consolidate the Echvoiew NASC export files
# ------------------------------------------

# Outputs:
# ---- DataFrame containing columns for `transect_num` [float], etc. required for transect data 
# ---- analysis
df_nasc_no_age1: pd.DataFrame = ingest_nasc.consolidate_echvoiew_nasc(
    df_merged=exports_with_regions_df,
    interval_df=interval_df,
    region_class_names=region_list_no_age1,
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_no_age1,
)

df_nasc_all_ages: pd.DataFrame = ingest_nasc.consolidate_echvoiew_nasc(
    df_merged=exports_with_regions_df,
    interval_df=interval_df,
    region_class_names=region_list_all_ages,
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_all_ages,
)

# ==================================================================================================
# [OPTIONAL] Read in a pre-consolidated NASC data file
# ----------------------------------------------------
nasc_filename: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/Data/Exports/US_CAN_NASC_2019_table_all_ages.xlsx"
)
nasc_sheet: str = "Sheet1"
FEAT_TO_ECHOPOP_COLUMNS: Dict[str, str] = {
    "transect": "transect_num",
    "region id": "region_id",
    "vessel_log_start": "distance_s",
    "vessel_log_end": "distance_e",
    "spacing": "transect_spacing",
    "layer mean depth": "layer_mean_depth",
    "layer height": "layer_height",
    "bottom depth": "bottom_depth",
    "assigned haul": "haul_num",
}

# Outputs:
# ---- DataFrame containing columns for `transect_num` [float], etc. required for transect data 
# ---- analysis
nasc_all_ages_df: pd.DataFrame = ingest_nasc.read_nasc_file(
    filename=nasc_filename, sheetname=nasc_sheet, column_name_map=FEAT_TO_ECHOPOP_COLUMNS
)

# ==================================================================================================
# [OPTIONAL] Filter the transect intervals to account for on- and off-effort
# --------------------------------------------------------------------------
transect_filter_filename: Path = Path("Path/to/file")
# ---- Note: this is only applicable to survey years 2012 and earlier, but this sort of file could
# ---- be generated for any year
transect_filter_sheet: str = "Sheet1"
subset_filter: str = "survey == 201003"

# Outputs:
# ---- DataFrame with filtered intervals representing on-effort
nasc_all_ages_cleaned_df: pd.DataFrame = ingest_nasc.filter_transect_intervals(
    nasc_df=nasc_all_ages_df,
    transect_filter_df=transect_filter_filename,
    subset_filter=subset_filter,
    transect_filter_sheet=transect_filter_sheet,
)

# ==================================================================================================
# Load in the biolodical data
# ---------------------------
ROOT_PATH: Path = Path("C:/Users/Brandyn/Documents/GitHub/Data")
biodata_filepath: Path = ROOT_PATH / "Biological/1995-2023_biodata_redo.xlsx"
biodata_sheet_map: Dict[str, str] = {
    "catch": "biodata_catch", 
    "length": "biodata_length",
    "specimen": "biodata_specimen",
}
subset_dict: Dict[Any, Any] = {
    "ships": {
        160: {
            "survey": 201906
        },
        584: {
            "survey": 2019097,
            "haul_offset": 200
        }
    },
    "species_code": [22500]
}
FEAT_TO_ECHOPOP_BIODATA_COLUMNS = {
    "frequency": "length_count",
    "haul": "haul_num",
    "weight_in_haul": "haul_weight",
}
biodata_label_map: Dict[Any, Dict] = {
    "sex": {
        1: "male",
        2: "female",
        3: "unsexed"
    }
}

# Outputs:
# ---- Dictionary comprising different biological datasets defined via `biodata_sheet_map`
dict_df_bio = load_data.load_biological_data(biodata_filepath, 
                                             biodata_sheet_map, 
                                             FEAT_TO_ECHOPOP_BIODATA_COLUMNS, 
                                             subset_dict, 
                                             biodata_label_map)

# ==================================================================================================
# Load in strata files
# --------------------
strata_filepath = ROOT_PATH / "Stratification/US_CAN strata 2019_final.xlsx"
strata_sheet_map = {
    "inpfc": "INPFC",
    "ks": "Base KS",
}
FEAT_TO_ECHOPOP_STRATA_COLUMNS = {
    "fraction_hake": "nasc_proportion",
    "haul": "haul_num",
    "stratum": "stratum_num",
}

# Outputs:
# ---- Dictionary comprising stratification information
dict_df_strata = load_data.load_strata(strata_filepath, 
                                       strata_sheet_map, 
                                       FEAT_TO_ECHOPOP_STRATA_COLUMNS)

# ==================================================================================================
# Load in geographical strata files
# ---------------------------------
geostrata_filepath = ROOT_PATH / "Stratification/Stratification_geographic_Lat_2019_final.xlsx"
geostrata_sheet_map = {
    "inpfc": "INPFC",
    "ks": "stratification1",
}
FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS = {
    "latitude (upper limit)": "northlimit_latitude",
    "stratum": "stratum_num",
}

# Outputs:
# ---- Dictionary comprising geoigraphical stratification information
dict_df_geostrata = load_data.load_geostrata(geostrata_filepath, 
                                             geostrata_sheet_map, 
                                             FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS)

# ==================================================================================================
# Stratify data based on haul numbers
# -----------------------------------

# Outputs:
# ---- DataFrame merged with stratum number (KS or INPFC)
# -------- INPFC (from strata)
nasc_all_ages_df_inpfc = load_data.join_strata_by_haul(nasc_all_ages_df, 
                                                       dict_df_strata["inpfc"])
dict_df_bio_inpfc = load_data.join_strata_by_haul(dict_df_bio,
                                                  dict_df_strata["inpfc"])
# -------- KS (from strata)
nasc_all_ages_df_ks = load_data.join_strata_by_haul(nasc_all_ages_df, 
                                                    dict_df_strata["ks"])
dict_df_bio_ks = load_data.join_strata_by_haul(dict_df_bio,
                                               dict_df_strata["ks"])

# ==================================================================================================
# Load kriging mesh file
# ----------------------
mesh_filepath = ROOT_PATH / "Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx"
mesh_sheet_name = "krigedgrid2_5nm_forChu"
FEAT_TO_ECHOPOP_MESH_COLUMNS = {
    "centroid_latitude": "latitude",
    "centroid_longitude": "longitude",
    "fraction_cell_in_polygon": "fraction",
}

# Outputs:
# ---- DataFrame containing the (centroid) latitude/longitude coordinates of mesh polygons and the 
# ---- fraction of the polygon included
mesh_df = load_data.load_mesh_data(mesh_filepath, mesh_sheet_name, FEAT_TO_ECHOPOP_MESH_COLUMNS)

# ==================================================================================================
# [OPTIONAL] Stratify data based on latitude intervals
# ----------------------------------------------------

# Outputs:
# NASC
# ---- DataFrame merged with geographically distributed stratum number (KS or INPFC)
# -------- INPFC (from geostrata)
nasc_all_ages_df_inpfc = load_data.join_geostrata_by_latitude(nasc_all_ages_df, 
                                                              dict_df_geostrata["inpfc"])
# -------- KS (from geostrata)
nasc_all_ages_df_ks = load_data.join_geostrata_by_latitude(nasc_all_ages_df, 
                                                           dict_df_geostrata["inpfc"])

# MESH
# ---- DataFrame merged with geographically distributed stratum number (KS or INPFC)
# -------- INPFC (from geostrata)
mesh_df_inpfc = load_data.join_geostrata_by_latitude(mesh_df, dict_df_geostrata["inpfc"])
# -------- KS (from geostrata)
mesh_df_ks = load_data.join_geostrata_by_latitude(mesh_df, dict_df_geostrata["inpfc"])

# ==================================================================================================
# Load kriging and variogram parameters
# -------------------------------------
geostatistic_params_filepath = ROOT_PATH / "Kriging_files/default_vario_krig_settings_2019_US_CAN.xlsx"
geostatistic_params_sheet_name = "Sheet1"
FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS = {
    "hole": "hole_effect_range",
    "lscl": "correlation_range",
    "nugt": "nugget",
    "powr": "decay_power",
    "ratio": "anisotropy",
    "res": "lag_resolution",
    "srad": "search_radius",
}
column_name_map = FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS

# Outputs:
# ---- Dictionaries comprising kriging and variogram model parameterization
kriging_params_dict, variogram_params_dict = load_data.load_kriging_variogram_params(
    geostatistic_params_filepath,
    geostatistic_params_sheet_name,
    FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS
)

# ===========================================
# Compute biological composition based on stratum
length_bins: np.array  # length bin specification
df_length_weight, df_regression = get_proportions.length_weight_regression(
    dict_df_bio["specimen"], length_bins
)
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

# Assemble an xr.Dataset of number and weight proportions
# I (WJ) think what you have in `distribute_length_age` is essentially this step
# This is a temporary step to convert the dataframes to xarray
# NOTE: This is a temporary step to convert the dataframes to xarray.
#       This can be removed once the following functions are implemented
#       to use xarray Dataset/DataArray as input/output directly:
#       -- get_proportions.number_proportions()
#       -- get_proportions.weight_distributions_over_lenghth_age()
#       -- get_proportions.stratum_averaged_weight()
#       -- get_proportions.weight_proportions()
ds_proportions: xr.Dataset = get_proportions.assemble_proportions(
    dict_df_number_proportion=dict_df_number_proportion,
    dict_df_weight_proportion=dict_df_weight_proportion,
)


# ===========================================
# NASC to number density and biomass

# Initiate object to perform inversion
# inversion parameters are stored as object attributes
invert_hake = InversionLengthTS(df_model_params=dict_df_bio["model_params"])

# Perform inversion using the supplied df_length
# df_length will be used to compute the mean sigma_bs for each stratum,
# which is then used in .invert() to compute number density on a stratum-by-stratum basis
df_nasc_no_age1 = invert_hake.invert(df_nasc=df_nasc_no_age1, df_length=dict_df_bio["length"])
df_nasc_all_ages = invert_hake.invert(df_nasc=df_nasc_all_ages, df_length=dict_df_bio["length"])


# Apportion abundance and biomass for transect intervals
# TODO: these apportioned transect results are not used in kriging, is this correct?
ds_nasc_no_age1_apportioned: xr.Dataset = apportion.apportion_transect_biomass_abundance(
    df_nasc=df_nasc_no_age1,
    ds_proportions=ds_proportions,
)
ds_nasc_all_age_apportioned: xr.Dataset = apportion.apportion_transect_biomass_abundance(
    df_nasc=df_nasc_all_ages,
    ds_proportions=ds_proportions,
)


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
df_nasc_no_age1_kriged = kriging.krige(df_in=df_nasc_no_age1, variables="biomass")
df_nasc_all_age_kriged = kriging.krige(df_in=df_nasc_all_ages, variables="biomass")


# ===========================================
# Apportion kriged biomass across sex, length bins, and age bins,
# and from there derive kriged abundance and kriged number density.
# Reference flow diagram: https://docs.google.com/presentation/d/1FOr2-iMQYj21VzVRDC-YUuqpOP0_urtI/edit?slide=id.p1#slide=id.p1  # noqa

# Age 1 kriged biomass -------------
# Apportion biomass
ds_kriged_biomass_age1: xr.Dataset = apportion.apportion_kriged_biomass(
    df_nasc=df_nasc_no_age1_kriged,
    ds_proportions=ds_proportions,
)

# Fill missing length bins of aged fish using length distributions of unaged fish
ds_kriged_biomass_age1: xr.Dataset = apportion.fill_missing_aged_from_unaged(
    ds_kriged_apportioned=ds_kriged_biomass_age1,
    ds_proportions=ds_proportions,
)

# Back-calculate abundance
ds_kriged_biomass_age1: xr.Dataset = apportion.back_calculate_kriged_abundance(
    ds_kriged_apportioned=ds_kriged_biomass_age1,
    ds_proportions=ds_proportions,
)


# All age (age 2+) kriged biomass -------------
# Apportion biomass
ds_kriged_biomass_all_ages: xr.Dataset = apportion.apportion_kriged_biomass(
    df_nasc=df_nasc_all_age_kriged,
    ds_proportions=ds_proportions,
)

# Fill missing length bins of aged fish using length distributions of unaged fish
ds_kriged_biomass_all_ages: xr.Dataset = apportion.fill_missing_aged_from_unaged(
    ds_kriged_apportioned=ds_kriged_biomass_all_ages,
    ds_proportions=ds_proportions,
)

# Reallocate age-1 fish to age-2+ fish
ds_kriged_biomass_all_ages: xr.Dataset = apportion.reallocate_age1(
    ds_kriged_apportioned=ds_kriged_biomass_all_ages,
    ds_proportions=ds_proportions,
)

# Back-calculate abundance
ds_kriged_biomass_all_ages: xr.Dataset = apportion.back_calculate_kriged_abundance(
    ds_kriged_apportioned=ds_kriged_biomass_all_ages,
    ds_proportions=ds_proportions,
)
