from pathlib import Path
from typing import Any, Dict, List
from typing import Any, Dict, List

import numpy as np
import numpy.typing as npt
import pandas as pd
import xarray as xr

from echopop.inversion import InversionLengthTS
from echopop.kriging import Kriging
from echopop.nwfsc_feat import biology, ingest_nasc, get_proportions, load_data, utils
# from echopop.nwfsc_feat import apportion, get_proportions, ingest_nasc, load_data

# ==================================================================================================
# ==================================================================================================
# DATA INGESTION
# ==================================================================================================
# Organize NASC file
# ------------------

# Merge exports
df_intervals, df_exports = ingest_nasc.merge_echoview_nasc(
    nasc_path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/raw_nasc/"),
    filename_transect_pattern = r"T(\d+)",
    default_transect_spacing = 10.0,
    default_latitude_threshold = 60.0,
)

# ==================================================================================================
# Read in transect-region-haul keys
# ---------------------------------
transect_region_filepath_all_ages: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/Stratification/"
    "US_CAN_2019_transect_region_haul_age1+ auto_final.xlsx"
)
transect_region_sheetname_all_ages: str = "Sheet1"
transect_region_filepath_no_age1: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/Stratification/"
    "US_CAN_2019_transect_region_haul_age2+ auto_20191205.xlsx"
)
transect_region_sheetname_no_age1: str = "Sheet1"
transect_region_file_rename: dict = {
    "tranect": "transect_num",
    "region id": "region_id",
    "trawl #": "haul_num",
}

# Read in the transect-region-haul key files for each group
transect_region_haul_key_all_ages: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
transect_region_haul_key_all_ages: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
    transect_region_filepath_all_ages,
    transect_region_sheetname_all_ages,
    transect_region_file_rename,
)

transect_region_haul_key_no_age1: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
transect_region_haul_key_no_age1: pd.DataFrame = ingest_nasc.read_transect_region_haul_key(
    transect_region_filepath_no_age1, transect_region_sheetname_no_age1, transect_region_file_rename
)

# ==================================================================================================
# Read in transect-region-haul keys
# ---------------------------------
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
df_exports_with_regions: pd.DataFrame = ingest_nasc.process_region_names(
    df_exports,
df_exports_with_regions: pd.DataFrame = ingest_nasc.process_region_names(
    df_exports,
    region_name_expr_dict,
    can_haul_offset = 200,
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
    filename=Path(
    "C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019//Exports/"
    "US_CAN_NASC_2019_table_all_ages.xlsx"
    ), 
    sheetname="Sheet1", 
    column_name_map=FEAT_TO_ECHOPOP_COLUMNS
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
ROOT_PATH: Path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/")

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
    "weight_in_haul": "weight",
}
biodata_label_map: Dict[Any, Dict] = {
    "sex": {
        1: "male",
        2: "female",
        3: "unsexed"
    }
}

# 
dict_df_bio = load_data.load_biological_data(
    biodata_filepath=ROOT_PATH / "Biological/1995-2023_biodata_redo.xlsx", 
    biodata_sheet_map=biodata_sheet_map, 
    column_name_map=FEAT_TO_ECHOPOP_BIODATA_COLUMNS, 
    subset_dict=subset_dict, 
    biodata_label_map=biodata_label_map
)

# ==================================================================================================
# Load in strata files
# --------------------
strata_sheet_map = {
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
    strata_filepath=ROOT_PATH / "Stratification/US_CAN strata 2019_final.xlsx", 
    strata_sheet_map=strata_sheet_map, 
    column_name_map=FEAT_TO_ECHOPOP_STRATA_COLUMNS
)

# ==================================================================================================
# Load in geographical strata files
# ---------------------------------
geostrata_sheet_map = {
    "inpfc": "INPFC",
    "ks": "stratification1",
}
FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS = {
    "latitude (upper limit)": "northlimit_latitude",
    "stratum": "stratum_num",
}

# 
df_dict_geostrata = load_data.load_geostrata(
    geostrata_filepath=ROOT_PATH / "Stratification/Stratification_geographic_Lat_2019_final.xlsx", 
    geostrata_sheet_map=geostrata_sheet_map, 
    column_name_map=FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS
)

# ==================================================================================================
# Stratify data based on haul numbers
# -----------------------------------

# Add INPFC
# ---- NASC
df_nasc_all_ages = load_data.join_strata_by_haul(df_nasc_all_ages, 
                                                 df_dict_strata["inpfc"],
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
    mesh_filepath=ROOT_PATH / "Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx", 
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
# KS (from geostrata)
df_nasc_all_ages = load_data.join_geostrata_by_latitude(df_nasc_all_ages, 
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
        ROOT_PATH / "Kriging_files/default_vario_krig_settings_2019_US_CAN.xlsx"
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
age_bins: npt.NDArray[np.number] = np.linspace(start=1., stop=22., num=22)
length_bins: npt.NDArray[np.number] = np.linspace(start=2., stop=80., num=40)

# 
# ---- Length
utils.binify(
    data=dict_df_bio, bins=length_bins, bin_column="length", 
)

# Age
utils.binify(
    data=dict_df_bio, bins=age_bins, bin_column="age",
)

# ==================================================================================================
# Fit length-weight regression to the binned data
# -----------------------------------------------

# Dictionary for length-weight regression coefficients
dict_length_weight_coefs = {}

# For all fish
dict_length_weight_coefs["all"] = biology.fit_length_weight_regression(
    data=dict_df_bio["specimen"]
)

# Sex-specific
dict_length_weight_coefs["sex"] = dict_df_bio["specimen"].groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)

# ==================================================================================================
# Compute the mean weights per length bin
# ---------------------------------------

# All fish (single coefficient set)
df_binned_weights_df_all = biology.length_binned_weights(
    data=dict_df_bio["specimen"],
    length_bins=length_bins,
    regression_coefficients=dict_length_weight_coefs["all"],
    impute_bins=True,
    minimum_count_threshold=5
)

# Sex-specific (grouped coefficients)
df_binned_weights_df_sexed = biology.length_binned_weights(
    data=dict_df_bio["specimen"],
    length_bins=length_bins,
    regression_coefficients=dict_length_weight_coefs["sex"],
    impute_bins=True,
    minimum_count_threshold=5
)

# Combine the pivot tables by adding the "all" column to the sex-specific table
binned_weight_table = pd.concat([df_binned_weights_df_all.assign(sex="all"), 
                                 df_binned_weights_df_sexed], 
                                axis=0,
                                ignore_index=True)

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
proportions_dict: Dict[str, pd.DataFrame] = dict_df_number_proportion
binned_weight_table: pd.DataFrame = binned_weight_table

#
df_averaged_weight = get_proportions.stratum_averaged_weight(
    proportions_dict=dict_df_number_proportion, 
    binned_weight_table=binned_weight_table,
    stratum_col="stratum_ks",
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

standardized_sexed_unaged_weights_df = get_proportions.standardize_weights_by_stratum(
    weights_df=dict_df_weight_distr["unaged"], 
    reference_weights_df=dict_df_bio["catch"].groupby(["stratum_ks"])["weight"].sum(),
    stratum_col="stratum_ks",
)

# ==================================================================================================
# Compute the standardized weight proportionsfor unaged fish
# ----------------------------------------------------------

dict_df_weight_proportion["unaged"] = get_proportions.standardize_weight_proportions(
    weight_data=standardized_sexed_unaged_weights_df, 
    reference_weight_proportions=dict_df_weight_proportion["aged"], 
    catch_data=dict_df_bio["catch"], 
    number_proportions=dict_df_number_proportion,
    binned_weights=df_binned_weights_df_all,
    group="unaged",
    group_columns = ["sex"],
    stratum_col = "stratum_ks"
)

# ==================================================================================================
# ==================================================================================================
# NASC TO POPULATION ESTIMATE CONVERSION
# ==================================================================================================
# ===========================================


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
