####################################################################################################
# 2012
# ----
from pathlib import Path
from echopop.workflows.nwfsc_feat import cli_utils
####################################################################################################
# PARAMETER ENTRY
# ---------------
# ** ENTER FILE INFORMATION FOR ALL INGESTED DATASETS.
# ** ADDITIONAL PARAMETERIZATIONS THROUGHOUT THE SCRIPT SHOULD BE EDITED BASED ON SPECIFIC NEEDS. 
# ** MAKE SURE TO EDIT WITH CARE.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRINT CONSOLE LOGGING MESSAGES
# ---- When set to `True`, logging information will be printed in the terminal/console as the 
# ---- script progresses
try: 
    # ---- FOR CLI USE
    VERBOSE = cli_utils.get_verbose()
except Exception:
    # ---- FOR INTERACTIVE REPL USE
    VERBOSE = True
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA ROOT DIRECTORY
DATA_ROOT = Path("C:/Data/EchopopData/echopop_2012")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REPORTS SAVE DIRECTORY
REPORTS_DIR = DATA_ROOT / "reports"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ALREADY PROCESSED NASC FILE ? 
# ---- When False, the raw NASC exports will be processed. When True, the pre-formatted NASC 
# ---- spreadsheet will be read in. This also requires defining `NASC_EXPORTS_SHEET`
NASC_PREPROCESSED = False
# NASC EXPORTS FILE(S)
NASC_EXPORTS_FILES = DATA_ROOT / "raw_nasc/"
# NASC EXPORTS SHEET
NASC_EXPORTS_SHEET = "Sheet1"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TRANSECT REGION HAUL MAPPING FILE
TRANSECT_REGION_HAUL_FILE = (
    DATA_ROOT / "Stratification/US&CAN_2012_transect_region_haul_age1+ auto final_new.xlsx"
)
# TRANSECT REGION HAUL MAPPING SHEET
TRANSECT_REGION_HAUL_SHEET = "Sheet1"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BIODATA FILE
BIODATA_FILE = DATA_ROOT / "Biological/1995-2023_biodata_redo.xlsx"
# BIODATA SHEETS
# ---- Assign the sheetnames to 'catch', 'length', 'specimen'
BIODATA_SHEETS = {
    "catch": "biodata_catch",
    "length": "biodata_length",
    "specimen": "biodata_specimen",
}
# BIODATA PROCESSING
# ---- This is used to parse the biodata master spreadsheet, which is required for aligning the 
# ---- biodata with ancillary files such as stratification and transect-haul mappings. This should 
# ---- define the "ships" based on their IDs with the associated survey IDs. If an offset should be 
# ---- added to the haul numbers, that must also be defined here. The target species should also be 
# ---- defined here. 
BIODATA_SHIP_SPECIES = {
    "ships": {
        19: {
            "survey": 201201
        },
        499: {
            "survey": 201205,
            "haul_offset": 100
        }
    },
    "species_code": [22500]
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HAUL STRATIFICATION FILE
HAUL_STRATA_FILE = (
    DATA_ROOT / 
    "Stratification/US&CAN strata 2012 11-30-2012.xlsx"
)
# HAUL STRATIFICATION SHEET MAP
# ---- Valid keys are limited to "ks" and "inpfc"
HAUL_STRATA_SHEETS = {
    "inpfc": "INPFC",
    "ks": "length strata byhaul_stratum1",
}
# GEOGRAPHIC STRATIFICATION FILE
GEOSTRATA_FILE = (
    DATA_ROOT / 
    "Stratification/Stratification_geographic_Lat.xlsx"
)
# GEOGRAPHIC STRATIFICATION SHEET MAP
# ---- Valid keys are limited to "ks" and "inpfc"
GEOSTRATA_SHEETS = {
    "inpfc": "INPFC",
    "ks": "stratification1",
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KRIGING MESH FILE 
KRIGING_MESH_FILE = (
    DATA_ROOT / "Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx"
)
# KRIGING MESH SHEET
KRIGING_MESH_SHEET = "krigedgrid2_5nm_forChu"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KRIGING AND VARIOGRAM PARAMETERS FILE
KRIGING_VARIOGRAM_PARAMETERS_FILE = (
    DATA_ROOT / "Kriging_files/default_vario_krig_settings_orig.xlsx"
)
# KRIGING AND VARIOGRAM PARAMETERS SHEET
KRIGING_VARIGORAM_PARAMETERS_SHEET = "Sheet1"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 200m ISOBATH FILE
ISOBATH_FILE = (
    DATA_ROOT / "Kriging_files/Kriging_grid_files/transformation_isobath_coordinates.xlsx"
)
# 200m ISOBATH SHEET
ISOBATH_SHEET = "Smoothing_EasyKrig"
####################################################################################################
####################################################################################################
# !!! START OF PROCESSING SCRIPT !!
# !!! EDIT CODE BELOW WITH CARE !!
####################################################################################################
####################################################################################################
import logging
import numpy as np
import pandas as pd
from lmfit import Parameters
from echopop.workflows.nwfsc_feat import functions as feat, parameters as feat_parameters, Reporter
import echopop.ingest as ingestion
from echopop import geostatistics, inversion, utils
from echopop.survey import biology, proportions, stratified, transect
from echopop.workflows.nwfsc_feat import apportionment as feat_apportion, biology as feat_biology
####################################################################################################
# FORMAT LOGGER
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
    
logging.basicConfig(
    level=logging.INFO if VERBOSE else logging.WARNING,
    format="%(message)s")
# ==================================================================================================
# DATA INGESTION 
# ==================================================================================================
# INGEST NASC DATA 
if NASC_PREPROCESSED:
    logging.info(f"Reading pre-generated NASC export file: '{NASC_EXPORTS_FILES.as_posix()}'.")

    # DEFINE COLUMN MAPPING
    FEAT_TO_ECHOPOP_COLUMNS = {
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

    # Read file
    df_nasc = ingestion.nasc.read_nasc_file(
        filename=NASC_EXPORTS_FILES,
        sheetname=NASC_EXPORTS_SHEET,
        column_name_map=FEAT_TO_ECHOPOP_COLUMNS
    )
else:
    logging.info(
        f"Beginning NASC export ingestion for files in: '{NASC_EXPORTS_FILES.as_posix()}'."
    )

    # MERGE EXPORTS
    logging.info(
        "---- Merging NASC exports...\n"
        "     Filename transect pattern: 'T(\\d+)'\n"
        "     Default transect spacing: 10.0 nmi\n"
        "     Default latitude threshold: 60.0 deg."
    )    
    df_intervals, df_exports = ingestion.nasc.merge_echoview_nasc(
        nasc_path = NASC_EXPORTS_FILES,
        filename_transect_pattern = r"T(\d+)",
        default_transect_spacing = 10.0,
        default_latitude_threshold = 60.0,
    )

    # GENERATE TRANSECT-REGION-HAUL KEY
    logging.info(
        "---- Loading transect-region-haul key mapping\n"
    )

    # TRANSECT REGION HAUL KEY NAME MAPPING
    TRANSECT_REGION_FILE_RENAME = {
        "tranect": "transect_num",
        "region id": "region_id",
        "trawl #": "haul_num",
    }

    # LOAD
    df_transect_region_haul_key = ingestion.nasc.read_transect_region_haul_key(
        filename=TRANSECT_REGION_HAUL_FILE,
        sheetname=TRANSECT_REGION_HAUL_SHEET,
        rename_dict=TRANSECT_REGION_FILE_RENAME
    )

    # CONSOLIDATE THE EXPORTS WITH TRANSECT-REGION-HAUL MAPPINGS
    logging.info(
        "---- Finalizing NASC export ingestion\n"
        "     Searching for the export regions: 'Age-1 Hake', 'Age-1 Hake Mix', 'Hake', 'Hake Mix'\n"
        "     Imputing overlapping region IDs within each interval: True"
    )
    df_nasc = ingestion.nasc.consolidate_echvoiew_nasc(
        df_merged=df_exports,
        interval_df=df_intervals,
        region_class_names=["Age-1 Hake", "Age-1 Hake Mix", "Hake", "Hake Mix"],
        impute_region_ids=True,
        transect_region_haul_key_df=df_transect_region_haul_key
    )
logging.info(
    "NASC ingestion complete\n"
    "'df_nasc' created."
)
# ==================================================================================================
# INGEST BIODATA
logging.info(
    f"Beginning biodata ingestion for: '{BIODATA_FILE.as_posix()}'."
)

# BIODATA DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_BIODATA_COLUMNS = {
    "frequency": "length_count",
    "haul": "haul_num",
    "weight_in_haul": "weight",
}

# BIODATA LABEL MAPPING
BIODATA_SEX = {
    "sex": {
        1: "male",
        2: "female",
        3: "unsexed"
    }
}

# READ IN DATA
dict_df_bio = ingestion.load_biological_data(
    biodata_filepath=BIODATA_FILE, 
    BIODATA_SHEETS=BIODATA_SHEETS, 
    column_name_map=FEAT_TO_ECHOPOP_BIODATA_COLUMNS, 
    subset_dict=BIODATA_SHIP_SPECIES, 
   biodata_label_map=BIODATA_SEX
)
# ---- Remove specimen hauls
feat_biology.remove_specimen_hauls(dict_df_bio)
logging.info(
    "Biodata ingestion complete\n"
    "'dict_df_bio' created."
)
# ==================================================================================================
# INGEST STRATIFICATION DATA
logging.info(
    "Loading stratification files..."
)

# HAUL-BASED STRATIFICATION DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_STRATA_COLUMNS = {
    "wt": "nasc_proportion",
    "haul": "haul_num",
    "strata": "stratum_num"
}

# READ IN STRATA FILE 
logging.info(
    f"Load in haul-based stratification: '{HAUL_STRATA_FILE.as_posix()}'."
)
df_dict_strata = ingestion.load_strata(
    strata_filepath=HAUL_STRATA_FILE, 
    strata_sheet_map=HAUL_STRATA_SHEETS, 
    column_name_map=FEAT_TO_ECHOPOP_STRATA_COLUMNS
)
logging.info(
    "Haul-based stratification loading complete\n"
    "'df_dict_strata' created."
)

# GEOGRAPHIC-BASED STRATIFICATION DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS = {
    "latitude (upper limit)": "northlimit_latitude",
    "strata": "stratum_num",
    "strata index": "stratum_num",
}

# READ IN GEOSTRATA FILE
logging.info(
    f"Load in geographic-based stratification: '{GEOSTRATA_FILE.as_posix()}'."
)
df_dict_geostrata = ingestion.load_geostrata(
    geostrata_filepath=GEOSTRATA_FILE, 
    geostrata_sheet_map=GEOSTRATA_SHEETS, 
    column_name_map=FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS
)
logging.info(
    "Geographic-based stratification loading complete\n"
    "'df_dict_geostrata' created."
)
# ==================================================================================================
# LOAD KRIGING MESH FILE
logging.info(
    f"Loading kriging mesh file: '{KRIGING_MESH_FILE.as_posix()}'."
)

# KRIGING MESH DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_MESH_COLUMNS = {
    "latitude of centroid": "latitude",
    "longitude of centroid": "longitude",
    "cell portion": "fraction",
}

# LOAD MESH
df_mesh = ingestion.load_mesh_data(
    mesh_filepath=KRIGING_MESH_FILE, 
    sheet_name=KRIGING_MESH_SHEET, 
    column_name_map=FEAT_TO_ECHOPOP_MESH_COLUMNS
)
logging.info(
    "Kriging mesh loading complete\n"
    "'df_mesh' created."
)
# ==================================================================================================
# LOAD ISOBATH FILE
logging.info(
    f"Loading isobath file: '{ISOBATH_FILE}'."
)
df_isobath = ingestion.load_isobath_data(
    isobath_filepath=ISOBATH_FILE,
    sheet_name=ISOBATH_SHEET
)
logging.info(
    f"Iosbath loading complete\n"
    "'df_isobath' created."
)
# ==================================================================================================
# LOAD KRIGING AND VARIOGRAM PARAMETERS
logging.info(
    f"Loading variogram and kriging parameters: '{KRIGING_VARIOGRAM_PARAMETERS_FILE.as_posix()}'."
)

# PARAMETERS DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS = {
    "hole": "hole_effect_range",
    "lscl": "correlation_range",
    "nugt": "nugget",
    "powr": "decay_power",
    "ratio": "aspect_ratio",
    "res": "lag_resolution",
    "srad": "search_radius",
}

# LOAD IN PARAMETERS
dict_kriging_params, dict_variogram_params = ingestion.load_kriging_variogram_params(
    geostatistic_params_filepath=KRIGING_VARIOGRAM_PARAMETERS_FILE,
    sheet_name=KRIGING_VARIGORAM_PARAMETERS_SHEET,
    column_name_map=FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS
)
logging.info(
    "Variogram and kriging parameter loading complete\n"
    "---- 'dict_variogram_params' created [variogram]\n"
    "---- 'dict_kriging_params' created [kriging]"
)
# ==================================================================================================
# INITIAL DATA PROCESSING
# ==================================================================================================
# APPLY STRATIFICATION DEFINITIONS TO BIODATA AND NASC
logging.info("Applying strata to datasets...")

# HAUL-BASED STRATA
logging.info(
    "Applying haul-based strata to 'dict_df_bio' and 'df_nasc'.\n"
    "     Default stratum: 0\n"
    "     New columns:\n"
    "         INPFC: 'stratum_inpfc'\n"
    "         KS: 'stratum_ks'"
)
# ---- BIODATA [INPFC]
dict_df_bio = ingestion.join_strata_by_haul(data=dict_df_bio,
                                            strata_df=df_dict_strata["inpfc"],
                                            default_stratum=0,
                                            stratum_name="stratum_inpfc")
# ---- BIODATA [KS]
dict_df_bio = ingestion.join_strata_by_haul(data=dict_df_bio,
                                            strata_df=df_dict_strata["ks"],
                                            default_stratum=0,
                                            stratum_name="stratum_ks")
# ---- NASC [INPFC]
df_nasc = ingestion.join_strata_by_haul(data=df_nasc,
                                        strata_df=df_dict_strata["inpfc"],
                                        default_stratum=0,
                                        stratum_name="stratum_inpfc")
# ---- NASC [KS]
df_nasc = ingestion.join_strata_by_haul(data=df_nasc,
                                        strata_df=df_dict_strata["ks"],
                                        default_stratum=0,
                                        stratum_name="stratum_ks")

# GEOGRAPHIC-BASED STRATA
logging.info(
    "Applying geographic-based strata to 'df_nasc' and 'df_mesh.\n"
    "     New columns:\n"
    "         INPFC: 'geostratum_inpfc'\n"
    "         KS: 'geostratum_ks'"
)
# ---- NASC [INPFC]
df_nasc = ingestion.join_geostrata_by_latitude(data=df_nasc,
                                               geostrata_df=df_dict_geostrata["inpfc"],
                                               stratum_name="geostratum_inpfc")
# ---- NASC [KS]
df_nasc = ingestion.join_geostrata_by_latitude(data=df_nasc,
                                               geostrata_df=df_dict_geostrata["ks"],
                                               stratum_name="geostratum_ks")
# ---- MESH [INPFC]
df_mesh = ingestion.join_geostrata_by_latitude(data=df_mesh, 
                                               geostrata_df=df_dict_geostrata["inpfc"], 
                                               stratum_name="geostratum_inpfc")
# ---- MESH [KS]
df_mesh = ingestion.join_geostrata_by_latitude(data=df_mesh, 
                                               geostrata_df=df_dict_geostrata["ks"], 
                                               stratum_name="geostratum_ks")
logging.info("Strata application complete!")
# ==================================================================================================
# BINIFY DATA 
logging.info(
    "Binning biodata ['dict_df_bio'] into discrete age and length bins\n"
    "     Age bins: [1, 2, 3, ..., 20, 21, 22]\n"
    "     Length bins: [2.0, 4.0, 6.0, ... 76.0, 78.0, 80.0]\n"
    "     New columns:\n"
    "         Age: 'age_bin'\n"
    "         Length: 'length_bin'"
)

# AGE-BINS
AGE_BINS = np.linspace(start=1., stop=22, num=22)
utils.binify(
    data=dict_df_bio, bins=AGE_BINS, bin_column="age",
)

# LENGTH-BINS
LENGTH_BINS = np.linspace(start=2., stop=80., num=40)
utils.binify(
    data=dict_df_bio, bins=LENGTH_BINS, bin_column="length", 
)
logging.info("Age and length binning complete!")
# ==================================================================================================
# FIT LENGTH-WEIGHT REGRESSION
logging.info("Fitting length-weight regression for each sex and for all fish.")

# CREATE DICTIONARY CONTAINER
dict_length_weight_coefs = {}

# ALL FISH
dict_length_weight_coefs["all"] = dict_df_bio["specimen"].assign(sex="all").groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)

# SEX-SPECIFIC
dict_length_weight_coefs["sex"] = dict_df_bio["specimen"].groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)
logging.info("Fitting length-weight regression for each sex and for all fish complete!")
# ==================================================================================================
# COMPUTE MEAN WEIGHTS PER LENGTH BIN
logging.info(
    "Computing the mean weight per length bin for each sex and for all fish.\n"
    "     Impute missing length bins using modeled weights: True"
    "     Minimum specimen count per bin: 5"
    )

# SEX-SPECIFIC
df_binned_weights_sex = feat_biology.length_binned_weights(
    data=dict_df_bio["specimen"],
    length_bins=LENGTH_BINS,
    regression_coefficients=dict_length_weight_coefs["sex"],
    impute_bins=True,
    minimum_count_threshold=5
)

# ALL FISH
df_binned_weights_all = feat_biology.length_binned_weights(
    data=dict_df_bio["specimen"].assign(sex="all"),
    length_bins=LENGTH_BINS,
    regression_coefficients=dict_length_weight_coefs["all"],
    impute_bins=True,
    minimum_count_threshold=5,
)

# COMBINE
binned_weight_table = pd.concat([df_binned_weights_sex, df_binned_weights_all], axis=1)
logging.info(
    "Length-binned mean weight calculations complete\n"
    "'binned_weight_table' created."
    )
# ==================================================================================================
# COMPUTE COUNT DISTRIBUTIONS PER AGE- AND LENGTH-BINS
logging.info(
    "Computing the counts per age- and length-bins across sex.\n"
    "     Stratifying by: 'stratum_ks'"
    "     Grouping by: 'sex'"
    )

# DICTIONARY CONTAINER
dict_df_counts = {}

# AGED
dict_df_counts["aged"] = proportions.compute_binned_counts(
    data=dict_df_bio["specimen"].dropna(subset=["age", "length", "weight"]), 
    groupby_cols=["stratum_ks", "length_bin", "age_bin", "sex"], 
    count_col="length",
    agg_func="size"
)

# UNAGED
dict_df_counts["unaged"] = proportions.compute_binned_counts(
    data=dict_df_bio["length"].copy().dropna(subset=["length"]), 
    groupby_cols=["stratum_ks", "length_bin", "sex"], 
    count_col="length_count",
    agg_func="sum"
)
logging.info(
    "Count distributions across age, length, and sex complete\n"
    "'dict_df_counts' created."
    )
# ==================================================================================================
# COMPUTE NUMBER PROPORTIONS
logging.info(
    "Computing number proportions across age and length bins\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Excluding: 'sex'='unsexed' from 'dict_df_counts['aged']'"
    )
dict_df_number_proportions = proportions.number_proportions(
    data=dict_df_counts, 
    group_columns=["stratum_ks"],
    exclude_filters={"aged": {"sex": "unsexed"}},
)
logging.info(
    "Number proportions calculation complete\n"
    "'dict_df_number_proportions' created\n"
    )
# ==================================================================================================
# COMPUTE BINNED WEIGHTS
logging.info(
    "Computing the summed weights per age- and length-bins across sex.\n"
    "     Stratifying by: 'stratum_ks'"
    "     Grouping by: 'sex'"
    "     Excluding: 'sex'='unsexed'"
    )

# DICTIONARY CONTAINER
dict_df_weight_distr = {}

# AGED
dict_df_weight_distr["aged"] = proportions.binned_weights(
    length_dataset=dict_df_bio["specimen"],
    include_filter = {"sex": ["female", "male"]},
    interpolate_regression=False,
    contrast_vars="sex",
    table_cols=["stratum_ks", "sex", "age_bin"]
)

# UNAGED
logging.info(
    "Unaged binned weights require additional processing steps.\n"
    "     Interpolating binned length-weight regression estimates: True"
    )
dict_df_weight_distr["unaged"] = proportions.binned_weights(
    length_dataset=dict_df_bio["length"],
    length_weight_dataset=binned_weight_table,
    include_filter = {"sex": ["female", "male"]},
    interpolate_regression=True,
    contrast_vars="sex",
    table_cols=["stratum_ks", "sex"]
)
logging.info(
    "Summed weights per age- and length-bins across sex computation complete\n"
    "'dict_df_weight_distr' created."
    )
# ==================================================================================================
# COMPUTE WEIGHT PROPORTIONS
logging.info(
    "Computing weight proportions across age and length bins\n"
    "     Stratifying by: 'stratum_ks'"
    "     Grouping by: 'sex'"
    )

# DICTIONARY CONTAINER
dict_df_weight_proportions = {}

# AGED WEIGHT PROPORTIONS
logging.info("Computing aged weight proportions...")
dict_df_weight_proportions["aged"] = proportions.weight_proportions(
    weight_data=dict_df_weight_distr, 
    catch_data=dict_df_bio["catch"], 
    group="aged",
    stratum_col="stratum_ks"
)

# UNAGED SCALING
logging.info("Scaling unaged binned weights...")
standardized_sexed_unaged_weights_df = proportions.scale_weights_by_stratum(
    weights_df=dict_df_weight_distr["unaged"], 
    reference_weights_df=dict_df_bio["catch"].groupby(["stratum_ks"])["weight"].sum(),
    stratum_col="stratum_ks",
)

# UNAGED WEIGHT PROPORTIONS
logging.info(
    "Computing unaged weight proportions\n"
    "     Scaling weight proportions in reference to the aged estimates"
    )
dict_df_weight_proportions["unaged"] = proportions.scale_weight_proportions(
    weight_data=standardized_sexed_unaged_weights_df, 
    reference_weight_proportions=dict_df_weight_proportions["aged"], 
    catch_data=dict_df_bio["catch"], 
    number_proportions=dict_df_number_proportions,
    binned_weights=binned_weight_table["all"],
    group="unaged",
    group_columns = ["sex"],
    stratum_col = "stratum_ks"
)
logging.info(
    "Weight proportions calculation complete\n"
    "'dict_df_weight_proportions' created."
    )
# ==================================================================================================
# NASC TO BIOMASS CONVERSION
# ==================================================================================================
# INVERSION
logging.info(
    "Beginning inversion based on hake-specific TS-length regression coefficients\n"
    "     Model: 20.0 x log[10](L) + -68.0\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Imputing missing strata: True\n"
    "     Treating hauls as replicates: True"
    )

# DEFINE INVERSION MODEL PARAMETERS
MODEL_PARAMETERS = {
    "ts_length_regression": {
        "slope": 20.,
        "intercept": -68.
    },
    "stratify_by": ["stratum_ks"],
    "expected_strata": df_dict_strata["ks"].stratum_num.unique(),
    "impute_missing_strata": True,
    "haul_replicates": True,
}

# INITIALIZE INVERSION OBJECT
invert_hake = inversion.InversionLengthTS(MODEL_PARAMETERS)
logging.info("Inversion-class object 'invert_hake' created...")

# INVERT NUMBER DENSITY
df_nasc = invert_hake.invert(
    df_nasc=df_nasc, df_length=[dict_df_bio["length"], dict_df_bio["specimen"]]
)
logging.info(
    "Number density inversion complete\n"
    "     New column in 'df_nasc':\n"
    "         Number density (animals nm^-2): 'number_density'"
    )
# ==================================================================================================
# CONVERT TO BIOMASS
logging.info("Converting number density estimates into biomass estimates")

# SET TRANSECT INTERVAL DISTANCES
logging.info(
    "Defining transect interval distances...\n"
    "     Along-transect interval distance (nmi) threshold: 0.5 nmi"
)
transect.compute_interval_distance(df_nasc=df_nasc, interval_threshold=0.05)

# SET TRANSECT INTERVAL AREAS
logging.info("Defining transect interval areas...")
df_nasc["area_interval"] = (
    df_nasc["transect_spacing"] * df_nasc["distance_interval"]
)

# COMPUTE ABUNDANCE
logging.info(
    "Compute interval abundances...\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Grouping by: 'sex'\n"
    "     Excluding: 'sex'='unsexed' from 'dict_df_number_proportions'"    
)
feat_biology.compute_abundance(
    dataset=df_nasc,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
    exclude_filter={"sex": "unsexed"},
    number_proportions=dict_df_number_proportions
)

# COMPUTE STRATUM-AVERAGED WEIGHTS
df_averaged_weight = proportions.stratum_averaged_weight(
    proportions_dict=dict_df_number_proportions, 
    binned_weight_table=binned_weight_table,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
)

# COMPUTE BIOMASS
logging.info(
    "Compute interval biomass...\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Grouping by: 'sex'\n"  
)
feat_biology.compute_biomass(
    dataset=df_nasc,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
    df_average_weight=df_averaged_weight,
)
logging.info(
    "NASC to biomass conversion complete\n"
    "     New columns in 'df_nasc':\n"
    "         Sex-specific number densities (animals nmi^-2): "
    "'number_density_female'/'number_density_male'\n"
    "         Abundance (animals): 'abundance'/'abundance_female'/'abundance_male'\n"
    "         Biomass density (kg nmi^-2): 'biomass_density'/'biomass_density_female'/"
    "'biomass_density_male'\n"
    "         Biomass (kg): 'biomass'/'biomass_female'/'biomass_male'"
    )
# ==================================================================================================
# DISTRIBUTE POPULATION ESTIMATES ACROSS AGE AND LENGTH BINS
logging.info(
    "Distribute population estimates across age- and length-bins\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Grouping by: 'sex'"
)

# ABUNDANCE
logging.info("Distributing abundances...")
dict_transect_abundance_table = feat_apportion.distribute_population_estimates(
    data=df_nasc,
    proportions=dict_df_number_proportions,
    variable="abundance",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
)
logging.info("Abundance distributions complete\n'dict_transect_abundance_table' created.")
# BIOMASS [ALL]
logging.info("Distributing biomass...")
dict_transect_biomass_table = feat_apportion.distribute_population_estimates(
    data=df_nasc,
    proportions=dict_df_weight_proportions,
    variable="biomass",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
)
logging.info("Biomass distribution complete\n'dict_transect_biomass_table' created.")
# BIOMASS [AGED-ONLY]
logging.info("Distributing biomass...\n     Aged-only weight proportions: True")
df_transect_aged_biomass_table = feat_apportion.distribute_population_estimates(
    data=df_nasc,
    proportions=dict_df_weight_proportions["aged"],
    variable="biomass",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
)
logging.info("Aged-biomass distribution complete\n'df_transect_aged_biomass_table' created.")
# ==================================================================================================
# GEOSTATISTICS
# ==================================================================================================
# INITIAL
logging.info("Beginning geostatistical analysis...")
# ==================================================================================================
# COORDINATE TRANSFORMATION
logging.info(
    "Transform spatial coordinates for 'df_nasc' and 'df_mesh'\n"
    "     Reference coordinates: 'df_isobath'\n"
    "     Longitudinal offset: -124.78338 deg.E\n"
    "     Latitudinal offset: 45.0 deg.N"
)

# NASC
df_nasc, delta_longitude, delta_latitude = geostatistics.transform_coordinates(
    data = df_nasc,
    reference = df_isobath,
    x_offset = -124.78338,
    y_offset = 45.,   
)

# MESH
df_mesh, _, _ = geostatistics.transform_coordinates(
    data = df_mesh,
    reference = df_isobath,
    x_offset = -124.78338,
    y_offset = 45.,   
    delta_x=delta_longitude,
    delta_y=delta_latitude
)
logging.info(
    "Coordinate transformation complete\n"
    "     New columns:\n"
    "          Transformed longitude: 'x'\n"
    "          Transformed latitude: 'y'\n"
)
# ==================================================================================================
# VARIOGRAM ANALYSIS
logging.info(
    "Beginning variogram analysis\n"
    "     Normalized lag resolution: 0.002\n"
    "     Number of lags: 30\n"
)

# INITIALIZE VARIOGRAM-CLASS OBJECT
vgm = geostatistics.Variogram(
    lag_resolution=0.002,
    n_lags=30,
    coordinate_names=("x", "y"),
)
logging.info("Variogram-class object 'vgm' created...")

# EMPIRICAL VARIOGRAM
logging.info(
    "Computing the empirical variogram\n"
    "     Variable: 'biomass_density'\n"
    "     Applying azimuth angle filter: True\n"
    "     Azimuth angle filter: 180.0 deg.\n"
)
vgm.calculate_empirical_variogram(
    data=df_nasc,
    variable="biomass_density",
    azimuth_filter=True,
    azimuth_angle_threshold=180.,
)

# SET UP FITTING PARAMETERS
# ----- lmfit.Parameters tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
logging.info(
    f"Optimizing variogram parameters using non-linear least-squares\n"
    f"     Model: Exponential-Bessel (['exponential', 'bessel'])\n"
    f"     Initial values:\n"
    f"          Nugget: {dict_variogram_params["nugget"]}\n"
    f"          Sill: {dict_variogram_params["sill"]}\n"
    f"          Correlation range: {dict_variogram_params["correlation_range"]}\n"
    f"          Hole effect range: {dict_variogram_params["hole_effect_range"]}\n"
    f"          Decay power exponent: {dict_variogram_params["decay_power"]}"
)
variogram_parameters_lmfit = Parameters()
variogram_parameters_lmfit.add_many(
    ("nugget", dict_variogram_params["nugget"], True, 0.),
    ("sill", dict_variogram_params["sill"], True, 0.),
    ("correlation_range", dict_variogram_params["correlation_range"], True, 0.),
    ("hole_effect_range", dict_variogram_params["hole_effect_range"], True, 0.),
    ("decay_power", dict_variogram_params["decay_power"], True, 1.25, 1.75),
)

# OPTIMIZATION PARAMETERS
OPTIM_ARGS = {
    "max_nfev": None, "ftol": 1e-08, "gtol": 1e-8, "xtol": 1e-8, "diff_step": None, 
    "tr_solver": "exact", "x_scale": 1., "jac": "2-point"
}
logging.info(
    f"Optimization arguments:\n"
    f"{OPTIM_ARGS}"
)

# RUN MINIMIZER
best_fit_parameters = vgm.fit_variogram_model(
    model=["exponential", "bessel"],
    model_parameters=variogram_parameters_lmfit,
    optimizer_kwargs=OPTIM_ARGS,
)
logging.info(
    f"Variogram parameter fitting complete\n"
    f"     Best-fit parameters:\n"
    f"     {best_fit_parameters}"
)
# ==================================================================================================
# KRIGING ANALYSIS
logging.info(
    f"Beginning kriging analysis\n"
    f"     Using best-fit variogram model and parameters: True\n"
    f"     Kriging parameters:\n"
    f"          Normalized search radius: {best_fit_parameters["correlation_range"] * 3}\n"
    f"          Minimum nearest-neighbor count: 3\n"
    f"          Maximum nearest-neighbor count: 10\n"
    f"          Anisotropic aspect ratio: 0.001"
)

# KRIGING PARAMETERS CONTAINER
KRIGING_PARAMETERS = {
    "search_radius": best_fit_parameters["correlation_range"] * 3,
    "aspect_ratio": 0.001,
    "k_min": 3,
    "k_max": 10,
}  

# VARIOGRAM PARAMETERS CONTAINER
VARIOGRAM_PARAMETERS = {
    "model": ["exponential", "bessel"],
    **best_fit_parameters
}

# INITIALIZE CLASS OBJECT
krg = geostatistics.Kriging(
    mesh=df_mesh,
    kriging_params=KRIGING_PARAMETERS,
    variogram_params=VARIOGRAM_PARAMETERS,
    coordinate_names=("x", "y"),
)
logging.info("Kriging-class object 'krg' created...")

# RUN KRIGING
logging.info(
    "Interpolating population estimates using ordinary kriging\n"
    "     Variable: 'biomass_density'\n"
    "     Extrapolation (full mesh): True\n"
    "     Default mesh cell area: 6.25 nmi^2\n"
)
df_kriged_results = krg.krige(
    transects=df_nasc,
    variable="biomass_density",
    extrapolate=True,
    default_mesh_cell_area=6.25,
)
logging.info(
    f"Kriging complete\n"
    f"'df_kriged_results' created."
)
# ==================================================================================================
# CONVERT BIOMASS DENSITY TO NASC
logging.info(
    "Converting biomass density estimates into NASC\n"
    "     Stratifying by: 'geostratum_ks'/'stratum_ks'\n"
    "     Grouping by: 'sex'\n"
    "     Using stratum weights for all fish: True"
)

# CONVERT TO BIOMASS
df_kriged_results["biomass"] = df_kriged_results["biomass_density"] * df_kriged_results["area"]
logging.info("New column in 'df_kriged_results': 'biomass'")

# BIOMASS TO NASC
feat_apportion.mesh_biomass_to_nasc(
    mesh_data_df=df_kriged_results,
    biodata=dict_df_weight_proportions,
    group_by=["sex"],
    mesh_biodata_link={"geostratum_ks": "stratum_ks"},
    stratum_weights_df=df_averaged_weight["all"],
    stratum_sigma_bs_df=invert_hake.sigma_bs_strata,    
)
logging.info(
    "Biomass density to NASC conversion complete\n"
    "    New columns in `df_kriged_results`\n"
    "        Biomass (kg): 'biomass'/'biomass_female'/'biomass_male'\n"
    "        Abundance (animals): 'abundance'/'abundance_female'/'abundance_male'\n"
    "        NASC (m^2 nmi^-2): 'nasc'"
)
# ==================================================================================================
# DISTRIBUTE POPULATION ESTIMATES ACROSS AGE AND LENGTH BINS
logging.info(
    "Distribute kriged population estimates across age- and length-bins\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Grouping by: 'sex'"
)

# ABUNDANCE [ALL]
logging.info("Distributing abundances...")
dict_kriged_abundance_table = feat_apportion.distribute_population_estimates(
    data=df_kriged_results,
    proportions=dict_df_number_proportions,
    variable="abundance",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
    data_proportions_link={"geostratum_ks": "stratum_ks"},
)
logging.info("Abundance distributions complete\n'dict_kriged_abundance_table' created.")

# SCALE UNAGED ABUNDANCE
logging.info(
    "Scaling unaged abundance...\n"     
    "     Reference: Aged abundances\n"
    "     Imputing missing bins: False"
)
dict_kriged_abundance_table["standardized_unaged"] = feat_apportion.distribute_unaged_from_aged(
    population_table=dict_kriged_abundance_table["unaged"],
    reference_table=dict_kriged_abundance_table["aged"],
    group_by=["sex"],
    impute=False,    
)

# BIOMASS [ALL]
logging.info("Distributing biomass...")
dict_kriged_biomass_table = feat_apportion.distribute_population_estimates(
    data=df_kriged_results,
    proportions=dict_df_weight_proportions,
    variable="biomass",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
    data_proportions_link={"geostratum_ks": "stratum_ks"},
)
logging.info("Biomass distribution complete\n'dict_kriged_biomass_table' created.")

# SCALE UNAGED BIOMASS
logging.info(
    "Scaling unaged biomass...\n"     
    "     Reference: Aged biomass\n"
    "     Imputing missing bins: True"
)
dict_kriged_biomass_table["standardized_unaged"] = feat_apportion.distribute_unaged_from_aged(
    population_table=dict_kriged_biomass_table["unaged"],
    reference_table=dict_kriged_biomass_table["aged"],
    group_by=["sex"],
    impute=True,
    impute_variable=["age_bin"],
)

# CONSOLIDATE
# ---- ABUNDANCE
logging.info("Consolidating abundance tables...")
df_kriged_abundance_table = feat_apportion.sum_population_tables(
    population_table=dict_kriged_abundance_table,
    table_names=["aged", "standardized_unaged"],
    table_index=["length_bin"],
    table_columns=["age_bin", "sex"],
)
logging.info("Abundance table complete\n'df_kriged_abundance_table' created.")
# ---- Biomass
logging.info("Consolidating biomass tables...")
df_kriged_biomass_table = feat_apportion.sum_population_tables(
    population_table=dict_kriged_biomass_table,
    table_names=["aged", "standardized_unaged"],
    table_index=["length_bin"],
    table_columns=["age_bin", "sex"],
)
logging.info("Biomass table complete\n'df_kriged_biomass_table' created.")
# ==================================================================================================
# JOLLY AND HAMPTON (1990) ANALYSIS
# ==================================================================================================
logging.info("Beginning stratified analysis to estimate uncertainties (Jolly and Hampton, 1990)...")

# ANALYSIS PARAMETERS CONTAINER
JOLLYHAMPTON_PARAMETERS = {
    "transects_per_latitude": 5,
    "strata_transect_proportion": 0.75,
    "num_replicates": 1000,
}

# INITIALIZE JOLLYHAMPTON CLASS OBJECT
jh = stratified.JollyHampton(JOLLYHAMPTON_PARAMETERS)
logging.info("Stratified-analysis-class object 'jh' created...")

# RUN ON TRANSECT DATA
logging.info(
    "Running Jolly and Hampton (1990) algorithm for transect data\n"
    "     Variable: 'biomass'\n"
    "     Number of bootstrap replicates: 1000\n"
    "     Stratum transect sampling proportion: 0.75\n"
    "     Stratifying by: 'geostratum_ks'"
)
jh.stratified_bootstrap(data_df=df_nasc, 
                        stratify_by=["geostratum_inpfc"], 
                        variable="biomass")
logging.info(
    "Summarizing results....\n"
    "     Confindence interval percentile: 0.95\n"
    "     Confidence interval method: Jackknife studentized interval ('t-jackknife')"
)
df_jh_transect_results = jh.summarize(ci_percentile=0.95, ci_method="t-jackknife")
logging.info("Stratified transect analysis results complete\n'df_jh_transect_results' created.")

# RUN ON KRIGED DATA
# ---- Create virtual transects
logging.info(
    "Creating virtual transects for kriged mesh data\n"
    "     Stratifying by: 'geostratum_inpfc'\n"
    "     Virtual transects per latitude: 5"
)
kriged_transects = jh.create_virtual_transects(
    data_df=df_kriged_results,
    geostrata_df=df_dict_geostrata["inpfc"], 
    stratify_by=["geostratum_inpfc"],
    variable="biomass",
)
# ---- Run rest of flow
logging.info(
    "Running Jolly and Hampton (1990) algorithm for kriged data\n"
    "     Variable: 'biomass'\n"
    "     Number of bootstrap replicates: 1000\n"
    "     Virtual transect sampling proportion: 0.75\n"
    "     Stratifying by: 'geostratum_ks'"
)
jh.stratified_bootstrap(data_df=kriged_transects, 
                        stratify_by=["geostratum_inpfc"], 
                        variable="biomass")
logging.info(
    "Summarizing results....\n"
    "     Confindence interval percentile: 0.95\n"
    "     Confidence interval method: Jackknife studentized interval ('t-jackknife')"
)
df_jh_kriged_results = jh.summarize(ci_percentile=0.95, ci_method="t-jackknife")
logging.info("Stratified kriged analysis results complete\n'df_jh_kriged_results' created.")
# ==================================================================================================
# REPORT GENERATION
# ==================================================================================================
# CREATE REPORTS
logging.info(
    f"Writing reports to: '{REPORTS_DIR.as_posix()}'."
)
reporter = Reporter(REPORTS_DIR, verbose=VERBOSE)

# AGED-LENGTH HAUL
reporter.aged_length_haul_counts_report(
    filename="aged_length_haul_counts.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    bio_data=dict_df_bio["specimen"].dropna(subset=["age", "length", "weight"])
)

# TOTAL LENGTH HAUL COUNTS
reporter.total_length_haul_counts_report(
    filename="total_length_haul_counts.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    bio_data=dict_df_bio
)

# KRIGED AGED BIOMASS

# All values
reporter.kriged_aged_biomass_mesh_report(
    filename="kriged_aged_biomass_mesh_full.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    kriged_data=df_kriged_results,
    weight_data=dict_df_weight_distr["aged"],
    kriged_stratum_link={"geostratum_ks": "stratum_ks"},
)

# Nonzero values
reporter.kriged_aged_biomass_mesh_report(
    filename="kriged_aged_biomass_mesh_nonzero.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    kriged_data=df_kriged_results[df_kriged_results["biomass"] > 0.],
    weight_data=dict_df_weight_distr["aged"],
    kriged_stratum_link={"geostratum_ks": "stratum_ks"},
)

# KRIGERD MESH RESULTS

# All values
reporter.kriged_mesh_results_report(
    filename="kriged_biomass_mesh_full.xlsx",
    sheetname="Sheet1",
    kriged_data=df_kriged_results,
    kriged_stratum="geostratum_ks",
    kriged_variable="biomass",
    sigma_bs_data=invert_hake.sigma_bs_strata,
    sigma_bs_stratum="stratum_ks",
)

# Nonzero values
reporter.kriged_mesh_results_report(
    filename="kriged_biomass_mesh_nonzero.xlsx",
    sheetname="Sheet1",
    kriged_data=df_kriged_results[df_kriged_results["biomass"] > 0.],
    kriged_stratum="geostratum_ks",
    kriged_variable="biomass",
    sigma_bs_data=invert_hake.sigma_bs_strata,
    sigma_bs_stratum="stratum_ks",
)

# KRIGED LENGTH-AGE ABUNDANCES
reporter.kriged_length_age_abundance_report(
    filename="kriged_length_age_abundance_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatables=dict_kriged_abundance_table,
)

# KRIGED LENGTH-AGE BIOMASS
reporter.kriged_length_age_biomass_report(
    filename="kriged_length_age_biomass_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatables=dict_kriged_biomass_table,
)

# KRIGING INPUT
reporter.kriging_input_report(
    filename="kriging_input_report.xlsx",
    sheetname="Sheet1",
    transect_data=df_nasc,
)

# TRANSECT LENGTH-AGE ABUNDANCES
reporter.transect_length_age_abundance_report(
    filename="transect_length_age_abundance_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatables=dict_transect_abundance_table,
)

# TRANSECT LENGTH-AGE BIOMASS
reporter.transect_length_age_biomass_report(
    filename="transect_length_age_biomass_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatable=dict_transect_biomass_table["aged"],
)

# TRANSECT AGED BIOMASS

# Full values
reporter.transect_aged_biomass_report(
    filename="transect_aged_biomass_report_full.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    transect_data=df_nasc,
    weight_data=dict_df_weight_distr["aged"],
)

# Nonzero values
reporter.transect_aged_biomass_report(
    filename="transect_aged_biomass_report_nonzero.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    transect_data=df_nasc[df_nasc["biomass"] > 0.],
    weight_data=dict_df_weight_distr["aged"],
)

# TRANSECT RESULTS

# Full values
reporter.transect_population_results_report(
    filename="transect_population_results_full.xlsx",
    sheetname="Sheet1",
    transect_data=df_nasc,
    weight_strata_data=df_averaged_weight,
    sigma_bs_stratum=invert_hake.sigma_bs_strata,
    stratum_name="stratum_ks",
)

# Nonzero values
reporter.transect_population_results_report(
    filename="transect_population_results_nonzero.xlsx",
    sheetname="Sheet1",
    transect_data=df_nasc[df_nasc["nasc"] > 0.],
    weight_strata_data=df_averaged_weight,
    sigma_bs_stratum=invert_hake.sigma_bs_strata,
    stratum_name="stratum_ks",
)

