from pathlib import Path
from typing import Any, Dict
from functools import partial

import numpy as np
import numpy.typing as npt
import pandas as pd

import copy
from lmfit import Parameters
from echopop import inversion
from echopop.nwfsc_feat.geostatistics import Geostats
from echopop.nwfsc_feat import (
    apportion,
    biology, 
    FEAT,
    ingest_nasc, 
    get_proportions, 
    load_data, 
    mesh,
    spatial,
    stratified,
    transect, 
    utils
)

# ==================================================================================================
# ==================================================================================================
# DEFINE DATA ROOT DIRECTORY
# --------------------------
DATA_ROOT = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019")
# DATA_ROOT = Path("C:/Users/Brandyn Lucca/Documents/Data/echopop_2019")

# ==================================================================================================
# ==================================================================================================
# DATA INGESTION
# ==================================================================================================
# Organize NASC file
# ------------------

# Merge exports
df_intervals, df_exports = ingest_nasc.merge_echoview_nasc(
    nasc_path = DATA_ROOT / "raw_nasc/",
    filename_transect_pattern = r"T(\d+)",
    default_transect_spacing = 10.0,
    default_latitude_threshold = 60.0,
)

# ==================================================================================================
# Read in transect-region-haul keys
# ---------------------------------
TRANSECT_REGION_FILEPATH_ALL_AGES: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/Stratification/"
    "US_CAN_2019_transect_region_haul_age1+ auto_final.xlsx"
)
TRANSECT_REGION_SHEETNAME_ALL_AGES: str = "Sheet1"
TRANSECT_REGION_FILEPATH_NO_AGE1: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/Stratification/"
    "US_CAN_2019_transect_region_haul_age2+ auto_20191205.xlsx"
)
TRANSECT_REGION_SHEETNAME_NO_AGE1: str = "Sheet1"
TRANSECT_REGION_FILE_RENAME: dict = {
    "tranect": "transect_num",
    "region id": "region_id",
    "trawl #": "haul_num",
}

# Read in the transect-region-haul key files for each group
transect_region_haul_key_all_ages: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
    filename=TRANSECT_REGION_FILEPATH_ALL_AGES,
    sheetname=TRANSECT_REGION_SHEETNAME_ALL_AGES,
    rename_dict=TRANSECT_REGION_FILE_RENAME,
)

transect_region_haul_key_no_age1: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
    TRANSECT_REGION_FILEPATH_NO_AGE1, TRANSECT_REGION_SHEETNAME_NO_AGE1, TRANSECT_REGION_FILE_RENAME
)

# ==================================================================================================
# Read in transect-region-haul keys
# ---------------------------------
REGION_NAME_EXPR_DICT: Dict[str, dict] = {
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
df_exports_with_regions: pd.DataFrame = ingest_nasc.process_region_names(
    df=df_exports,
    region_name_expr_dict=REGION_NAME_EXPR_DICT,
    can_haul_offset=200,
)

# ==================================================================================================
# [OPTIONAL] Generate transect-region-haul key from compiled values
# ---------------------------------

# Generate transect-region-haul key from compiled values
df_transect_region_haul_key_no_age1: pd.DataFrame = ingest_nasc.generate_transect_region_haul_key(
    df=df_exports_with_regions, 
    filter_list=["Hake", "Hake Mix"]
)

df_transect_region_haul_key_all_ages = ingest_nasc.generate_transect_region_haul_key(
    df=df_exports_with_regions, 
    filter_list=["Age-1 Hake", "Age-1", "Hake", "Hake Mix"]
)

# ==================================================================================================
# Consolidate the Echvoiew NASC export files
# ------------------------------------------
df_nasc_no_age1: pd.DataFrame = ingest_nasc.consolidate_echvoiew_nasc(
    df_merged=df_exports_with_regions,
    interval_df=df_intervals,
    region_class_names=["Hake", "Hake Mix"],
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_no_age1,
)

df_nasc_all_ages: pd.DataFrame = ingest_nasc.consolidate_echvoiew_nasc(
    df_merged=df_exports_with_regions,
    interval_df=df_intervals,
    region_class_names=["Age-1 Hake", "Age-1", "Hake", "Hake Mix"],
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_all_ages,
)

# ==================================================================================================
# [OPTIONAL] Read in a pre-consolidated NASC data file
# ----------------------------------------------------
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

#
df_nasc_all_ages: pd.DataFrame = ingest_nasc.read_nasc_file(
    filename=DATA_ROOT / "Exports/US_CAN_NASC_2019_table_all_ages.xlsx",
    sheetname="Sheet1",
    column_name_map=FEAT_TO_ECHOPOP_COLUMNS,
)

# ==================================================================================================
# [OPTIONAL] Filter the transect intervals to account for on- and off-effort
# --------------------------------------------------------------------------

# DataFrame with filtered intervals representing on-effort
df_nasc_all_ages_cleaned: pd.DataFrame = ingest_nasc.filter_transect_intervals(
    nasc_df=df_nasc_all_ages,
    transect_filter_df=Path("Path/to/file"),
    subset_filter="survey == 201003",
    transect_filter_sheet="Sheet1",
)

# ==================================================================================================
# Load in the biolodical data
# ---------------------------
BIODATA_SHEET_MAP: Dict[str, str] = {
    "catch": "biodata_catch", 
    "length": "biodata_length",
    "specimen": "biodata_specimen",
}
SUBSET_DICT: Dict[Any, Any] = {
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
    "weight_in_haul": "weight",
}
BIODATA_LABEL_MAP: Dict[Any, Dict] = {
    "sex": {
        1: "male",
        2: "female",
        3: "unsexed"
    }
}

# 
dict_df_bio = load_data.load_biological_data(
    biodata_filepath=DATA_ROOT / "Biological/1995-2023_biodata_redo.xlsx", 
    biodata_sheet_map=BIODATA_SHEET_MAP, 
    column_name_map=FEAT_TO_ECHOPOP_BIODATA_COLUMNS, 
    subset_dict=SUBSET_DICT, 
    biodata_label_map=BIODATA_LABEL_MAP
)

# ==================================================================================================
# Filter the haul data to avoid duplicated weight sums
# ----------------------------------------------------

biology.remove_specimen_hauls(dict_df_bio)

# ==================================================================================================
# Load in strata files
# --------------------
STRATA_SHEET_MAP = {
    "inpfc": "INPFC",
    "ks": "Base KS",
}
FEAT_TO_ECHOPOP_STRATA_COLUMNS = {
    "fraction_hake": "nasc_proportion",
    "haul": "haul_num",
    "stratum": "stratum_num",
}

#
df_dict_strata = load_data.load_strata(
    strata_filepath=DATA_ROOT / "Stratification/US_CAN strata 2019_final.xlsx", 
    strata_sheet_map=STRATA_SHEET_MAP, 
    column_name_map=FEAT_TO_ECHOPOP_STRATA_COLUMNS
)

# ==================================================================================================
# Load in geographical strata files
# ---------------------------------
GEOSTRATA_SHEET_MAP = {
    "inpfc": "INPFC",
    "ks": "stratification1",
}
FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS = {
    "latitude (upper limit)": "northlimit_latitude",
    "stratum": "stratum_num",
}

# 
df_dict_geostrata = load_data.load_geostrata(
    geostrata_filepath=DATA_ROOT / "Stratification/Stratification_geographic_Lat_2019_final.xlsx", 
    geostrata_sheet_map=GEOSTRATA_SHEET_MAP, 
    column_name_map=FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS
)

# ==================================================================================================
# Stratify data based on haul numbers
# -----------------------------------

# Add INPFC
# ---- NASC
df_nasc_all_ages = load_data.join_strata_by_haul(data=df_nasc_all_ages, 
                                                 strata_df=df_dict_strata["inpfc"],
                                                 stratum_name="stratum_inpfc") 
# ---- Biodata
dict_df_bio = load_data.join_strata_by_haul(dict_df_bio,
                                            df_dict_strata["inpfc"],
                                            stratum_name="stratum_inpfc")

# Add KS
# ---- NASC
df_nasc_all_ages = load_data.join_strata_by_haul(df_nasc_all_ages, 
                                                 df_dict_strata["ks"],
                                                 stratum_name="stratum_ks") 
df_nasc_no_age1 = load_data.join_strata_by_haul(df_nasc_no_age1, 
                                                df_dict_strata["ks"],
                                                stratum_name="stratum_ks") 
# ---- Biodata
dict_df_bio = load_data.join_strata_by_haul(dict_df_bio,
                                            df_dict_strata["ks"],
                                            stratum_name="stratum_ks") 

# ==================================================================================================
# Load kriging mesh file
# ----------------------

FEAT_TO_ECHOPOP_MESH_COLUMNS = {
    "centroid_latitude": "latitude",
    "centroid_longitude": "longitude",
    "fraction_cell_in_polygon": "fraction",
}

# 
df_mesh = load_data.load_mesh_data(
    mesh_filepath=DATA_ROOT / "Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx", 
    sheet_name="krigedgrid2_5nm_forChu", 
    column_name_map=FEAT_TO_ECHOPOP_MESH_COLUMNS
)

# ==================================================================================================
# [OPTIONAL] Stratify data based on latitude intervals
# ----------------------------------------------------
# INPFC (from geostrata)
df_nasc_all_ages = load_data.join_geostrata_by_latitude(df_nasc_all_ages, 
                                                        df_dict_geostrata["inpfc"],
                                                        stratum_name="geostratum_inpfc")
df_nasc_no_age1 = load_data.join_geostrata_by_latitude(df_nasc_no_age1, 
                                                       df_dict_geostrata["inpfc"],
                                                       stratum_name="geostratum_inpfc")
# KS (from geostrata)
df_nasc_all_ages = load_data.join_geostrata_by_latitude(df_nasc_all_ages, 
                                                        df_dict_geostrata["ks"],
                                                        stratum_name="geostratum_ks")
df_nasc_no_age1 = load_data.join_geostrata_by_latitude(df_nasc_no_age1, 
                                                       df_dict_geostrata["ks"],
                                                       stratum_name="geostratum_ks")

# MESH
# ---- DataFrame merged with geographically distributed stratum number (KS or INPFC)
# -------- INPFC (from geostrata)
df_mesh = load_data.join_geostrata_by_latitude(df_mesh, 
                                               df_dict_geostrata["inpfc"], 
                                               stratum_name="geostratum_inpfc")
# -------- KS (from geostrata)
df_mesh = load_data.join_geostrata_by_latitude(df_mesh, 
                                               df_dict_geostrata["ks"], 
                                               stratum_name="geostratum_ks")

# ==================================================================================================
# Load kriging and variogram parameters
# -------------------------------------

FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS = {
    "hole": "hole_effect_range",
    "lscl": "correlation_range",
    "nugt": "nugget",
    "powr": "decay_power",
    "ratio": "anisotropy",
    "res": "lag_resolution",
    "srad": "search_radius",
}

# 
dict_kriging_params, dict_variogram_params = load_data.load_kriging_variogram_params(
    geostatistic_params_filepath=(
        DATA_ROOT / "Kriging_files/default_vario_krig_settings_2019_US_CAN.xlsx"
    ),
    sheet_name="Sheet1",
    column_name_map=FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS
)

# ==================================================================================================
# ==================================================================================================
# DATA PROCESSING
# ==================================================================================================
# Generate binned distributions [age, length]
# -------------------------------------------
AGE_BINS: npt.NDArray[np.number] = np.linspace(start=1., stop=22., num=22)
LENGTH_BINS: npt.NDArray[np.number] = np.linspace(start=2., stop=80., num=40)

# 
# ---- Length
utils.binify(
    data=dict_df_bio, bins=LENGTH_BINS, bin_column="length", 
)

# Age
utils.binify(
    data=dict_df_bio, bins=AGE_BINS, bin_column="age",
)

# ==================================================================================================
# Fit length-weight regression to the binned data
# -----------------------------------------------

# Dictionary for length-weight regression coefficients
dict_length_weight_coefs = {}

# For all fish
dict_length_weight_coefs["all"] = dict_df_bio["specimen"].assign(sex="all").groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)

# Sex-specific
dict_length_weight_coefs["sex"] = dict_df_bio["specimen"].groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)

# ==================================================================================================
# Compute the mean weights per length bin
# ---------------------------------------

# Sex-specific (grouped coefficients)
df_binned_weights_sex = biology.length_binned_weights(
    data=dict_df_bio["specimen"],
    length_bins=LENGTH_BINS,
    regression_coefficients=dict_length_weight_coefs["sex"],
    impute_bins=True,
    minimum_count_threshold=5
)

# All fish (single coefficient set)
df_binned_weights_all = biology.length_binned_weights(
    data=dict_df_bio["specimen"].assign(sex="all"),
    length_bins=LENGTH_BINS,
    regression_coefficients=dict_length_weight_coefs["all"],
    impute_bins=True,
    minimum_count_threshold=5,
)

# Combine the pivot tables by adding the "all" column to the sex-specific table
binned_weight_table = pd.concat([df_binned_weights_sex, df_binned_weights_all], axis=1)

# ==================================================================================================
# Compute the count distributions per age- and length-bins
# --------------------------------------------------------

# Dictionary for number counts
dict_df_counts = {}

# Aged
dict_df_counts["aged"] = get_proportions.compute_binned_counts(
    data=dict_df_bio["specimen"].dropna(subset=["age", "length", "weight"]), 
    groupby_cols=["stratum_ks", "length_bin", "age_bin", "sex"], 
    count_col="length",
    agg_func="size"
)

# Unaged
dict_df_counts["unaged"] = get_proportions.compute_binned_counts(
    data=dict_df_bio["length"].copy().dropna(subset=["length"]), 
    groupby_cols=["stratum_ks", "length_bin", "sex"], 
    count_col="length_count",
    agg_func="sum"
)

# ==================================================================================================
# Compute the number proportions
# ------------------------------
dict_df_number_proportion: Dict[str, pd.DataFrame] = get_proportions.number_proportions(
    data=dict_df_counts, 
    group_columns=["stratum_ks"],
    column_aliases=["aged", "unaged"],
    exclude_filters=[{"sex": "unsexed"}, None] 
)

# ==================================================================================================
# Distribute (bin) weight over age, length, and sex
# -------------------------------------------------
# Pre-allocate a dictionary
dict_df_weight_distr: Dict[str, Any] = {}

# Aged
dict_df_weight_distr["aged"] = get_proportions.binned_weights(
    length_dataset=dict_df_bio["specimen"],
    include_filter = {"sex": ["female", "male"]},
    interpolate=False,
    contrast_vars="sex",
    table_cols=["stratum_ks", "sex", "age_bin"]
)

# Unaged
dict_df_weight_distr["unaged"] = get_proportions.binned_weights(
    length_dataset=dict_df_bio["length"],
    length_weight_dataset=binned_weight_table,
    include_filter = {"sex": ["female", "male"]},
    interpolate=True,
    contrast_vars="sex",
    table_cols=["stratum_ks", "sex"]
)

# ==================================================================================================
# Calculate the average weights pre stratum when combining different datasets
# ---------------------------------------------------------------------------
df_averaged_weight = get_proportions.stratum_averaged_weight(
    proportions_dict=dict_df_number_proportion, 
    binned_weight_table=binned_weight_table,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
)

# ==================================================================================================
# Compute the length-binned weight proportions for aged fish
# ----------------------------------------------------------

# Initialize Dictionary container
dict_df_weight_proportion: Dict[str, Any] = {}

# Aged
dict_df_weight_proportion["aged"] = get_proportions.weight_proportions(
    weight_data=dict_df_weight_distr, 
    catch_data=dict_df_bio["catch"], 
    group="aged",
    stratum_col="stratum_ks"
)

# ==================================================================================================
# Compute the standardized haul weights for unaged fish
# -----------------------------------------------------

standardized_sexed_unaged_weights_df = get_proportions.scale_weights_by_stratum(
    weights_df=dict_df_weight_distr["unaged"], 
    reference_weights_df=dict_df_bio["catch"].groupby(["stratum_ks"])["weight"].sum(),
    stratum_col="stratum_ks",
)

# ==================================================================================================
# Compute the standardized weight proportionsfor unaged fish
# ----------------------------------------------------------

dict_df_weight_proportion["unaged"] = get_proportions.scale_weight_proportions(
    weight_data=standardized_sexed_unaged_weights_df, 
    reference_weight_proportions=dict_df_weight_proportion["aged"], 
    catch_data=dict_df_bio["catch"], 
    number_proportions=dict_df_number_proportion,
    binned_weights=binned_weight_table["all"],
    group="unaged",
    group_columns = ["sex"],
    stratum_col = "stratum_ks"
)

# ==================================================================================================
# ==================================================================================================
# NASC TO POPULATION ESTIMATE CONVERSION
# ==================================================================================================
# Initialize the Inversion class
# ------------------------------
model_parameters = {
    "ts_length_regression": {
        "slope": 20.,
        "intercept": -68.
    },
    "stratify_by": ["stratum_ks"],
    "expected_strata": df_dict_strata["ks"].stratum_num.unique(),
    "impute_missing_strata": True,
    "haul_replicates": True,
}

# Initiate object to perform inversion
invert_hake = inversion.InversionLengthTS(MODEL_PARAMETERS)

# ==================================================================================================
# Invert number density
# ---------------------

# If the above haul-averaged `sigma_bs` values were calculated, then the inversion can can 
# completed without calling in additional biodata
df_nasc_all_ages = invert_hake.invert(df_nasc=df_nasc_all_ages,
                                      df_length=[dict_df_bio["length"], dict_df_bio["specimen"]])
df_nasc_no_age1 = invert_hake.invert(df_nasc=df_nasc_no_age1,
                                     df_length=[dict_df_bio["length"], dict_df_bio["specimen"]])
# ---- The average `sigma_bs` for each stratum can be inspected at:
cached_values = invert_hake.sigma_bs_strata

# Alternatively, these can be computed directly by skipping the intermediate averaging applied 
# to hauls first before being averaged for each stratum
invert_hake.invert(df_nasc=df_nasc_no_age1, 
                   df_length=[dict_df_bio["length"], dict_df_bio["specimen"]])
# ---- These yield subtle, but non-zero differences in the average `sigma_bs` per stratum
invert_hake.sigma_bs_strata - cached_values

# ==================================================================================================
# Set transect interval distances
# -------------------------------

# Calculate along-transect interval distances which is required for getting the area-per-interval 
# and therefore going from number density to abundance
transect.set_interval_distance(df_nasc=df_nasc_all_ages, interval_threshold=0.05)
transect.set_interval_distance(df_nasc=df_nasc_no_age1, interval_threshold=0.05)

# ==================================================================================================
# Calculate transect interval areas
# ---------------------------------
df_nasc_all_ages["area_interval"] = (
    df_nasc_all_ages["transect_spacing"] * df_nasc_all_ages["distance_interval"]
)
df_nasc_no_age1["area_interval"] = (
    df_nasc_no_age1["transect_spacing"] * df_nasc_no_age1["distance_interval"]
)

# ==================================================================================================
# Calculate remaining population metrics across all animals 
# ---------------------------------------------------------
biology.set_population_metrics(df_nasc=df_nasc_all_ages, 
                               metrics=["abundance", "biomass", "biomass_density"],
                               stratify_by="stratum_ks",
                               df_average_weight=df_averaged_weight["all"])

biology.set_population_metrics(df_nasc=df_nasc_no_age1, 
                               metrics=["abundance", "biomass", "biomass_density"],
                               stratify_by="stratum_ks",
                               df_average_weight=df_averaged_weight["all"])

# ==================================================================================================
# Apportion age-1 vs age-2+ population estimates 
# ----------------------------------------------
# TODO: This apportionment step is required for kriging

# This step is used in EchoPro
# Otherwise, the mean `sigma_bs` can be computed directly from the data (as shown below), although 
# computing the mean average sigma_bs per haul better accounts for pseudoreplication 
invert_hake.set_haul_sigma_bs(df_length=[dict_df_bio["length"], dict_df_bio["specimen"]])
# ---- This DataFrame can be inspected at:
invert_hake.sigma_bs_haul

# ==================================================================================================
# Invert number density
# ---------------------

# If the above haul-averaged `sigma_bs` values were calculated, then the inversion can can 
# completed without calling in additional biodata
df_nasc_all_ages = invert_hake.invert(df_nasc=df_nasc_all_ages)
df_nasc_no_age1 = invert_hake.invert(df_nasc=df_nasc_no_age1)
# ---- The average `sigma_bs` for each stratum can be inspected at:
invert_hake.sigma_bs_strata

# ==================================================================================================
# Set transect interval distances
# -------------------------------

# Calculate along-transect interval distances which is required for getting the area-per-interval 
# and therefore going from number density to abundance
transect.set_interval_distance(df_nasc=df_nasc_all_ages, interval_threshold=0.05)
transect.set_interval_distance(df_nasc=df_nasc_no_age1, interval_threshold=0.05)

# ==================================================================================================
# Calculate transect interval areas
# ---------------------------------
df_nasc_all_ages["area_interval"] = (
    df_nasc_all_ages["transect_spacing"] * df_nasc_all_ages["distance_interval"]
)
df_nasc_no_age1["area_interval"] = (
    df_nasc_no_age1["transect_spacing"] * df_nasc_no_age1["distance_interval"]
)

# ==================================================================================================
# Calculate (and apportion) number densities to abundance, and number densities/abundance for each 
# sex 
# --------------------------------------------------------------------------------------------------

biology.set_abundance(
    dataset=df_nasc_no_age1,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
    exclude_filter={"sex": "unsexed"},
    number_proportions=dict_df_number_proportion
)

# ==================================================================================================
# Calculate (and apportion) biomass densities and biomass (from number density and abundance, 
# respectively) for the overall transect dataset as well as for each sex
# --------------------------------------------------------------------------------------------------

biology.set_biomass(
    dataset=df_nasc_no_age1,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
    df_average_weight=df_averaged_weight,
)

# ==================================================================================================
# Get proportions for each stratum specific to age-1
# --------------------------------------------------

# Age-1 NASC proportions
age1_nasc_proportions = get_proportions.get_nasc_proportions_slice(
    number_proportions=dict_df_number_proportion["aged"],
    stratify_by=["stratum_ks"],
    ts_length_regression_parameters={"slope": 20., 
                                     "intercept": -68.},
    include_filter = {"age_bin": [1]}
)

# Age-1 number proportions
age1_number_proportions = get_proportions.get_number_proportions_slice(
    number_proportions=dict_df_number_proportion["aged"],
    stratify_by=["stratum_ks"],
    include_filter = {"age_bin": [1]}
)

# Age-1 weight proportions
age1_weight_proportions = get_proportions.get_weight_proportions_slice(
    weight_proportions=dict_df_weight_proportion["aged"],
    stratify_by=["stratum_ks"],
    include_filter={"age_bin": [1]},
    number_proportions=dict_df_number_proportion,
    length_threshold_min=10.0,
    weight_proportion_threshold=1e-10
)

# ==================================================================================================
# Apply the calculated proportions to the abundance, biomass, and NASC estimates
# ------------------------------------------------------------------------------

df_nasc_no_age1_prt = apportion.partition_transect_data(
    dataset=df_nasc_no_age1,
    partition_dict={
        "nasc": age1_nasc_proportions, 
        "abundance": age1_number_proportions,
        "biomass": age1_weight_proportions
    },
)

# ==================================================================================================
# ==================================================================================================
# GEOSTATISTICS
# ==================================================================================================
# Load reference line (isobath)
# -----------------------------

df_isobath = load_data.load_isobath_data(
    isobath_filepath=DATA_ROOT / "Kriging_files/Kriging_grid_files/transformation_isobath_coordinates.xlsx", 
    sheet_name="Smoothing_EasyKrig", 
)

# ==================================================================================================
# Initialize geostatistics class
# ------------------------------

# Define the requisite kriging parameters
kriging_parameters = {
    "search_radius": 0.021,
    "anisotropy": 0.001,
    "k_min": 3,
    "k_max": 10,
}  

# Define the requisite variogram parameters
variogram_parameters = {
    "model": ["bessel", "exponential"], 
    "n_lags": 30, 
    "lag_resolution": 0.002
}

# Initialize
geo = Geostats(
    data_df=df_nasc_no_age1_prt,
    mesh_df=df_mesh,
    kriging_params=kriging_parameters,
    variogram_params=variogram_parameters,    
)

# ==================================================================================================
# Standardize coordinates [transect & mesh]
# -----------------------------------------
geo.project_coordinates(
    reference_df=df_isobath,
    x_offset=-124.78338,
    y_offset=45.,
    normalize=True,
)

# ==================================================================================================
# Compute the empirical variogram
# -------------------------------
geo.calculate_empirical_variogram(
    variable="biomass_density",
    azimuth_filter=True,
    azimuth_angle_threshold=180.,
    force_lag_zero=True,
)

# ==================================================================================================
# Fit theoretical/modeled variogram to the transect data
# ------------------------------------------------------

# Set up `lmfit` parameters
variogram_parameters_lmfit = Parameters()
variogram_parameters_lmfit.add_many(
    ("nugget", dict_variogram_params["nugget"], True, 0., None),
    ("sill", dict_variogram_params["sill"], True, 0., None),
    ("correlation_range", dict_variogram_params["correlation_range"], True, 0., None),
    ("hole_effect_range", dict_variogram_params["hole_effect_range"], True, 0., None),
    ("decay_power", dict_variogram_params["decay_power"], True, 0., None),
)

# Set up optimization parameters used for fitting the variogram
dict_optimization = {"max_nfev": 500, "ftol": 1e-06, "gtol": 0.0001, "xtol": 1e-06, 
                     "diff_step": 1e-08, "tr_solver": "exact", "x_scale": "jac", 
                     "jac": "3-point"}

# Get the best-fit variogram parameters
geo.fit_variogram_model(
    variogram_parameters_lmfit, dict_optimization,
)

# ==================================================================================================
# Mesh cropping using the FEAT methods
# ------------------------------------
geo.crop_mesh(
    crop_function=mesh.transect_ends_crop,
    # cropping_method = ""
    latitude_resolution=1.25/60.,
    transect_mesh_region_function=FEAT.transect_mesh_region_2019,
)

# ==================================================================================================
# [OPTIONAL] Mesh cropping using the hull convex
# ----------------------------------------------
GEO_COPY = copy.deepcopy(geo)

GEO_COPY.crop_mesh(
    crop_function=mesh.hull_crop,
    num_nearest_transects=3,
    mesh_buffer_distance=2.5,
)

# ==================================================================================================
# Get the western extent of the transect bounds
# ---------------------------------------------
transect_western_extents = spatial.get_survey_western_extents(
    transect_df=geo.data_df,
    coordinate_names=("x", "y"),
    latitude_threshold=51.
)

# ==================================================================================================
# Krige the biomass density to get kriged biomass
# -----------------------------------------------

# Pre-define arguments within a partial function defining the western boundary search strategy
boundary_search_strategy = partial(spatial.western_boundary_search_strategy, 
                                   western_extent=transect_western_extents,
                                   kriging_mesh=geo.mesh_df,
                                   coordinate_names=("x", "y"))

# Krige
df_kriged_results = geo.krige(
    default_mesh_cell_area=6.25,
    adaptive_search_strategy=boundary_search_strategy,
)

# ##################################################################################################
# Back-calculate sex-specific biomass and abundance, and total NASC from the kriged biomass 
# density estimates
# -----------------

apportion.mesh_biomass_to_nasc(
    mesh_data_df=df_kriged_results,
    biodata=dict_df_weight_proportion,
    group_by=["sex"],
    mesh_biodata_link={"geostratum_ks": "stratum_ks"},
    stratum_weights_df=df_averaged_weight["all"],
    stratum_sigma_bs_df=invert_hake.sigma_bs_strata,    
)

# ##################################################################################################
# Distribute kriged abundance estimates over length and age/length
# ----------------------------------------------------------------

dict_kriged_abundance_table = apportion.distribute_kriged_estimates(
    mesh_data_df=df_kriged_results,
    proportions=dict_df_number_proportion,
    variable="abundance",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
    mesh_proportions_link={"geostratum_ks": "stratum_ks"},
)

# ##################################################################################################
# Distribute kriged biomass estimates over length and age/length
# --------------------------------------------------------------

dict_kriged_biomass_table = apportion.distribute_kriged_estimates(
    mesh_data_df=df_kriged_results,
    proportions=dict_df_weight_proportion,
    variable="biomass",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
    mesh_proportions_link={"geostratum_ks": "stratum_ks"},
)

# ##################################################################################################
# Standardize the unaged abundance estimates to be distributed over age
# ---------------------------------------------------------------------

dict_kriged_abundance_table["standardized_unaged"] = apportion.standardize_kriged_estimates(
    population_table=dict_kriged_abundance_table["unaged"],
    reference_table=dict_kriged_abundance_table["aged"],
    group_by=["sex"],
    impute=False,    
)

# ##################################################################################################
# Standardize the unaged abundance estimates to be distributed over age
# ---------------------------------------------------------------------
# THIS NEEDS TO BE FULLY IMPLEMENTED WITH TESTS, ETC.

dict_kriged_biomass_table["standardized_unaged"] = apportion.standardize_kriged_estimates(
    population_table=dict_kriged_biomass_table["unaged"],
    reference_table=dict_kriged_biomass_table["aged"],
    group_by=["sex"],
    impute=True,
    impute_variable=["age_bin"],
)

# ##################################################################################################
# Consolidate the kriged abundance estimates into a single DataFrame table
# ------------------------------------------------------------------------

df_kriged_abundance_table = apportion.combine_population_tables(
    population_table=dict_kriged_abundance_table,
    table_names=["aged", "standardized_unaged"],
    table_index=["length_bin"],
    table_columns=["age_bin", "sex"],
)

# ##################################################################################################
# Consolidate the kriged biomass estimates into a single DataFrame table
# -----------------------------------------------------------------------

df_kriged_biomass_table = apportion.combine_population_tables(
    population_table=dict_kriged_biomass_table,
    table_names=["aged", "standardized_unaged"],
    table_index=["length_bin"],
    table_columns=["age_bin", "sex"],
)

# ##################################################################################################
# Redistribute the kriged abundance estimates
# -------------------------------------------
# THIS NEEDS TO BE FULLY IMPLEMENTED WITH TESTS, ETC.

# Re-allocate the age-1 abundance estimates 
df_kriged_abundance_table_noage1 = apportion.redistribute_population_table(
    population_table=df_kriged_abundance_table,
    exclusion_filter={"age_bin": [1]},
    group_by=["sex"],
)

####################################################################################################
# Redistribute the kriged biomass estimates
# -----------------------------------------

# Re-allocate the age-1 abundance estimates 
df_kriged_biomass_table_noage1 = apportion.redistribute_population_table(
    population_table=df_kriged_biomass_table,
    exclusion_filter={"age_bin": [1]},
    group_by=["sex"],
)

# ##################################################################################################
# Instantiate stratified analysis
# -------------------------------
# Model parameters for stratified analysis initialization
JOLLYHAMPTON_PARAMETERS = {
    "transects_per_latitude": 5,
    "strata_transect_proportion": 0.75,
    "num_replicates": 1000,
}

jh = stratified.JollyHampton(JOLLYHAMPTON_PARAMETERS)

# ##################################################################################################
# Compute coefficient of variation (CV) and other estimators for transect biomass
# -------------------------------------------------------------------------------

# Run bootstrapping procedure
jh.stratified_bootstrap(data_df=df_nasc_no_age1, 
                        stratify_by=["geostratum_inpfc"], 
                        variable="biomass")

# Compute summary statistics for each stratum and overall survey
transect_results = jh.summarize(ci_percentile=0.95, ci_method="t-jackknife")
print(transect_results)

# ##################################################################################################
# Create virtual transects for the gridded kriged biomass results
# ---------------------------------------------------------------

kriged_transects = jh.create_virtual_transects(
    data_df=df_kriged_results,
    geostrata_df=df_dict_geostrata["inpfc"], 
    stratify_by=["geostratum_inpfc"],
    variable="biomass",
)

# ##################################################################################################
# Compute CV and other estimators for kriged biomass
# --------------------------------------------------

# Run bootstrapping procedure
jh.stratified_bootstrap(data_df=kriged_transects, 
                        stratify_by=["geostratum_inpfc"], 
                        variable="biomass")

# Compute summary statistics for each stratum and overall survey
kriged_results = jh.summarize(ci_percentile=0.95, ci_method="t-jackknife")
print(kriged_results)

# ##################################################################################################
# Compare transect and kriged outputs
# -----------------------------------

print(
    kriged_results.xs("mean", axis=1, level="metric") - 
    transect_results.xs("mean", axis=1, level="metric")
)
