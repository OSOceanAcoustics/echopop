from pathlib import Path
from typing import Any, Dict, List, Union

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
nasc_path: Path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/raw_nasc/")
filename_transect_pattern: str = r"T(\d+)"
default_transect_spacing: float = 10.0  # nmi
default_transect_spacing_latitude: float = 60.0  # deg N

# Merge exports
df_intervals, df_exports = ingest_nasc.merge_echoview_nasc(
    nasc_path,
    filename_transect_pattern,
    default_transect_spacing,
    default_transect_spacing_latitude,
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
df_exports_with_regions: pd.DataFrame = ingest_nasc.process_region_names(
    df_exports,
    region_name_expr_dict,
    CAN_haul_offset,
)

# ==================================================================================================
# [OPTIONAL] Generate transect-region-haul key from compiled values
# ---------------------------------
region_list_no_age1: List[str] = ["Hake", "Hake Mix"]
region_list_all_ages: List[str] = ["Age-1 Hake", "Age-1", "Hake", "Hake Mix"]

# Generate transect-region-haul key from compiled values
df_transect_region_haul_key_no_age1: pd.DataFrame = ingest_nasc.generate_transect_region_haul_key(
    df_exports_with_regions, filter_list=region_list_no_age1
)

df_transect_region_haul_key_all_ages = ingest_nasc.generate_transect_region_haul_key(
    df_exports_with_regions, filter_list=region_list_all_ages
)

# ==================================================================================================
# Consolidate the Echvoiew NASC export files
# ------------------------------------------
df_nasc_no_age1: pd.DataFrame = ingest_nasc.consolidate_echvoiew_nasc(
    df_merged=df_exports_with_regions,
    interval_df=df_intervals,
    region_class_names=region_list_no_age1,
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_no_age1,
)

df_nasc_all_ages: pd.DataFrame = ingest_nasc.consolidate_echvoiew_nasc(
    df_merged=df_exports_with_regions,
    interval_df=df_intervals,
    region_class_names=region_list_all_ages,
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_all_ages,
)

# ==================================================================================================
# [OPTIONAL] Read in a pre-consolidated NASC data file
# ----------------------------------------------------
nasc_filename: Path = Path(
    "C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019//Exports/US_CAN_NASC_2019_table_all_ages.xlsx"
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

#
df_nasc_all_ages: pd.DataFrame = ingest_nasc.read_nasc_file(
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
df_nasc_all_ages_cleaned: pd.DataFrame = ingest_nasc.filter_transect_intervals(
    nasc_df=df_nasc_all_ages,
    transect_filter_df=transect_filter_filename,
    subset_filter=subset_filter,
    transect_filter_sheet=transect_filter_sheet,
)

# ==================================================================================================
# Load in the biolodical data
# ---------------------------
ROOT_PATH: Path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/")
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

#
df_dict_strata = load_data.load_strata(strata_filepath, 
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

# 
df_dict_geostrata = load_data.load_geostrata(geostrata_filepath, 
                                             geostrata_sheet_map, 
                                             FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS)

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
mesh_filepath = ROOT_PATH / "Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx"
mesh_sheet_name = "krigedgrid2_5nm_forChu"
FEAT_TO_ECHOPOP_MESH_COLUMNS = {
    "centroid_latitude": "latitude",
    "centroid_longitude": "longitude",
    "fraction_cell_in_polygon": "fraction",
}

# 
df_mesh = load_data.load_mesh_data(mesh_filepath, mesh_sheet_name, FEAT_TO_ECHOPOP_MESH_COLUMNS)

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
length_weight_dataset: pd.DataFrame = dict_df_bio_binned_ks["specimen"].copy()

# Outputs: 
# ---- `pandas.Series` that includes the fitted slope and intercepts for the length-weight 
# ---- log-linear regressions. This can either be ran directly on a `pandas.DataFrame` or via 
# ---- `pandas.DataFrame.groupby`
length_weight_coefs_all: pd.Series = biology.fit_length_weight_regression(length_weight_dataset)
length_weight_coefs_sex = length_weight_dataset.groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)

# ==================================================================================================
# Compute the mean weights per length bin
# ---------------------------------------
data: pd.DataFrame = dict_df_bio_binned_ks["specimen"].copy()
regression_coefficients: Union[pd.Series, pd.DataFrame] = length_weight_coefs_all
impute_bins: bool = True
minimum_count_threshold: int = 5

# Outputs: 
# ---- `pandas.DataFrame` with the fitted mean weights per length bin with the corresponding 
# ---- index column, if it exists (e.g. `"sex"`)

# All fish
binned_weights_df_all = biology.length_binned_weights(
    data, length_distribution, regression_coefficients, impute_bins, minimum_count_threshold
)

# Sex-specific
binned_weights_df_sexed = biology.length_binned_weights(
    data=dict_df_bio_binned_ks["specimen"], 
    length_distribution=length_distribution,
    regression_coefficients=length_weight_coefs_sex, 
    impute_bins=True, 
    minimum_count_threshold=5
)

# Concatenate the sex-specific and all-specimen dataframes
binned_weight_table = pd.concat(
    [binned_weights_df_sexed.copy(),
    binned_weights_df_all.assign(sex="all")],
    ignore_index=True
)

# ==================================================================================================
# Compute the count distributions per age- and length-bins
# --------------------------------------------------------
# !!! When it comes down to this operation, it seems like overkill to have this set operation  
# !!! defined as a totally separate function when the operation is pretty straightforward. It is 
# !!! perhaps just a tad repetitive.

# Get the binned aged (age-length-bins) and unaged (length-bins) counts
# ---- Aged (all fish) -- inclusive of unsexed
# aged_counts_all = (
#     dict_df_bio_binned_ks["specimen"]
#     .groupby(["stratum_num", "length_bin", "age_bin"], observed=False)["length"].size()
#     .to_frame("count")
# )
# # ---- Aged (sexed fish) -- in this case exclude "unsexed" to leave only "female"/"male"
# aged_counts_sexed = (
#     dict_df_bio_binned_ks["specimen"]
#     .loc[dict_df_bio_binned_ks["specimen"].sex != "unsexed"]
#     .groupby(["stratum_num", "length_bin", "age_bin", "sex"], observed=False)["length"].size()
#     .to_frame("count")
# )
# # ---- Aged (all fish) -- no unsexed [equivalent to the `binned_aged_counts_filtered_df`]
# aged_counts_no_unsexed_all = (
#     dict_df_bio_binned_ks["specimen"]
#     .loc[dict_df_bio_binned_ks["specimen"].sex != "unsexed"]
#     .groupby(["stratum_num", "length_bin", "age_bin"], observed=False)["length"].size()
#     .to_frame("count")
# )
# # ---- Unaged (all fish)
# unaged_counts_all = (
#     dict_df_bio_binned_ks["length"]
#     .loc[dict_df_bio_binned_ks["length"].sex != "unsexed"]
#     .groupby(["stratum_num", "length_bin"], observed=False)["length"].sum()
#     .to_frame("count")
# )
# # ---- Unaged (sexed fish)
# unaged_counts_sexed = (
#     dict_df_bio_binned_ks["length"]
#     .loc[dict_df_bio_binned_ks["length"].sex != "unsexed"]
#     .groupby(["stratum_num", "length_bin", "sex"], observed=False)["length_count"].sum()
#     .to_frame("count")
# )

base_groupby_cols: List[str] = ["stratum_num", "length_bin"]
aged_df: pd.DataFrame = dict_df_bio_binned_ks["specimen"].copy().dropna(subset=["age", "length", "weight"])
unaged_df: pd.DataFrame = dict_df_bio_binned_ks["length"].copy().dropna(subset=["length"])

# Get the binned aged (age-length-bins) and unaged (length-bins) counts
# ---- Aged (all fish) -- inclusive of unsexed
# aged_counts_all = get_proportions.compute_binned_counts(
#     data=aged_df, 
#     groupby_cols=base_groupby_cols + ["age_bin"], 
#     count_col="length",
#     agg_func="size"
# )
# ---- Aged (sexed fish) -- in this case exclude "unsexed" to leave only "female"/"male"
aged_counts_sexed = get_proportions.compute_binned_counts(
    data=aged_df, 
    groupby_cols=base_groupby_cols + ["age_bin", "sex"], 
    count_col="length",
    agg_func="size"
)

# # ---- Aged (all fish) -- no unsexed [equivalent to the `binned_aged_counts_filtered_df`]
# aged_counts_no_unsexed_all = get_proportions.compute_binned_counts(
#     data=aged_df, 
#     groupby_cols=base_groupby_cols + ["age_bin"], 
#     count_col="length",
#     agg_func="size",
#     exclude_filters=exclude_filter
# )
# # ---- Aged (sexed fish) -- no unsexed (just "female"/"male")
# aged_counts_no_unsexed_sexed = get_proportions.compute_binned_counts(
#     data=aged_df, 
#     groupby_cols=base_groupby_cols + ["age_bin"], 
#     count_col="length",
#     agg_func="size",
#     exclude_filters=exclude_filter
# )
# # ---- Unaged (all fish)
# unaged_counts_all = get_proportions.compute_binned_counts(
#     data=unaged_df, 
#     groupby_cols=base_groupby_cols, 
#     count_col="length_count",
#     agg_func="sum",
#     exclude_filters=exclude_filter
# )
# ---- Unaged (sexed fish)
unaged_counts_sexed = get_proportions.compute_binned_counts(
    data=unaged_df, 
    groupby_cols=base_groupby_cols + ["sex"], 
    count_col="length_count",
    agg_func="sum"
)

# ==================================================================================================
# Compute the number proportions
# ------------------------------
group_columns = ["stratum_num"]
column_aliases = ["aged", "unaged"]
exclude_filters: Dict[str, Any] = [{"sex": "unsexed"}, None] # Positional filter match dataframe order

# Outputs: 
# ---- `Dict[str, pandas.DataFrame]` with the number proportions computed for aged and unaged fish. 
# ----- This includes columns for `proportion` (within group) and `proportion_overall` 
# ---- (across groups)
dict_df_number_proportion: Dict[str, pd.DataFrame] = get_proportions.number_proportions(
    aged_counts_sexed, unaged_counts_sexed, 
    group_columns=group_columns,
    column_aliases=column_aliases,
    exclude_filters=exclude_filters
)

# ==================================================================================================
# Distribute (bin) weight over age, length, and sex
# -------------------------------------------------

# Outputs: 
# ---- `pandas.DataFrame` with the binned weights distributed across defined bins/strata/sex/etc.

# Pre-allocate a dictionary
dict_df_weight_distr: Dict[str, Any] = {}

# Aged
dict_df_weight_distr["aged"] = get_proportions.binned_weights(
    dict_df_bio_binned_ks["specimen"].copy(),
    include_filter = {"sex": ["female", "male"]},
    interpolate=False,
    contrast_vars="sex",
    table_cols=["stratum_num", "sex", "age_bin"]
)

# Unaged
dict_df_weight_distr["unaged"] = get_proportions.binned_weights(
    dict_df_bio_binned_ks["length"].copy(),
    length_weight_dataset=binned_weights_df_sexed.copy(),
    include_filter = {"sex": ["female", "male"]},
    interpolate=True,
    contrast_vars="sex",
    table_cols=["stratum_num", "sex"]
)

# ==================================================================================================
# Calculate the average weights pre stratum when combining different datasets
# ---------------------------------------------------------------------------
proportion_dict: Dict[str, pd.DataFrame] = dict_df_number_proportion
binned_weight_table: pd.DataFrame = binned_weight_table

# Output:
# ---- `pandas.DataFrame` indexed by stratum with the average weights per stratum for male, female, 
# ---- and all fish
df_averaged_weight = get_proportions.stratum_averaged_weight(
    proportion_dict, binned_weight_table
)

# ==================================================================================================
# Compute the length-binned weight proportions for aged fish
# ----------------------------------------------------------
weight_data: Dict[str, pd.DataFrame] = dict_df_weight_distr
catch_data: pd.DataFrame = dict_df_bio_binned_ks["catch"]
group: str = "aged"

# Initialize Dictionary container
dict_df_weight_proportion: Dict[str, Any] = {}

# Output:
# ---- `pandas.DataFrame` indexed by, in this case, length, age, and sex across all strata 
# ---- containing weight proportions
dict_df_weight_proportion["aged"] = get_proportions.weight_proportions(
    weight_data, catch_data, group
)

# ==================================================================================================
# Compute the standardized haul weights for unaged fish
# -----------------------------------------------------
weight_data: pd.DataFrame = dict_df_weight_distr["unaged"]
stratum_weights: pd.DataFrame = (
    dict_df_bio_binned_ks["catch"].groupby(["stratum_num"])["haul_weight"]
    .sum().reset_index(name="weight")
)

# Output:
# ---- `pandas.DataFrame` indexed by, in this case, sex across all strata containing the 
# ---- standardized haul weights
standardized_sexed_unaged_weights_df = get_proportions.standardize_weights_by_stratum(
    weight_data, stratum_weights
)

# ==================================================================================================
# Compute the standardized weight proportionsfor unaged fish
# ----------------------------------------------------------
weight_data: pd.DataFrame = standardized_sexed_unaged_weights_df
catch_data: Dict[str, pd.DataFrame] = dict_df_bio_binned_ks["catch"]
reference_data: pd.DataFrame = dict_df_weight_proportion["aged"]
proportion_dict: Dict[str, pd.DataFrame] = proportion_dict
binned_weight_table: pd.DataFrame = binned_weights_df_all
group: str = "unaged"
group_columns: List[str] = ["sex"]

# Output:
# ---- `pandas.DataFrame` indexed by, in this case, length, and sex across all strata containing 
# ---- the standardized weight proportions based on the sexed ratios and fitted length-weight 
# ---- distributions computed using the aged measurements
dict_df_weight_proportion["unaged"] = get_proportions.standardize_weight_proportions(
    weight_data, reference_data, catch_data, 
    proportion_dict, binned_weight_table, group, group_columns
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
