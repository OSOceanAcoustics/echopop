from pathlib import Path
from echopop.nwfsc_feat.ingest_nasc import consolidate_echvoiew_nasc, generate_transect_region_haul_key, process_region_names, merge_echoview_nasc, read_transect_region_haul_key, read_echoview_nasc, map_transect_num, validate_transect_exports, clean_echoview_cells_df
# ===========================================
# Organize NASC file
nasc_path: Path = Path("C:/Users/Brandyn Lucca/Documents/Data/echopop_2019/Raw_NASC")
filename_transect_pattern: str = r"T(\d+)"
default_transect_spacing = 10. # nmi
default_transect_spacing_latitude = 60. # deg N

transect_region_filepath_all_ages: Path = Path("C:/Users/Brandyn Lucca/Documents/Data/echopop_2019/Stratification/US&CAN_2019_transect_region_haul_age1+ auto_final.xlsx")
transect_region_sheetname_all_ages: str = "Sheet1"
transect_region_filepath_no_age1: Path = Path("C:/Users/Brandyn Lucca/Documents/Data/echopop_2019/Stratification/US&CAN_2019_transect_region_haul_age2+ auto_20191205.xlsx")
transect_region_sheetname_no_age1: str = "Sheet1"
transect_region_file_rename: dict = {
    "tranect": "transect_num",
    "region id": "region_id",
    "trawl #": "haul_num",
}

CAN_haul_offset = 200
# region_name_pattern = "{REGION_CLASS}{HAUL_NUM}{COUNTRY}"
region_name_expr_dict = {
    "REGION_CLASS": {
        # "Age-1 Hake": "^(?:h1a(?![a-z]|1a))",
        # "Age-1 Hake": "^(?:h1a$)",
        "Age-1 Hake": "^(?:h1a(?![a-z]|m))",
        "Age-1 Hake Mix": "^(?:h1am(?![a-z]|1a))",
        "Hake": "^(?:h(?![a-z]|1a)|hake(?![_]))",
        "Hake Mix": "^(?:hm(?![a-z]|1a)|hake_mix(?![_]))"
    },
    "HAUL_NUM": {
        "[0-9]+",
    },
    "COUNTRY": {
        "CAN": "^[cC]",
        "US": "^[uU]",
    }
}
filter_list_no_age1 = ["Hake", "Hake Mix"]
filter_list_all_ages = ["Age-1 Hake", "Age-1", "Hake", "Hake Mix"]

interval_template, exports_df = merge_echoview_nasc(nasc_path, 
                                                    filename_transect_pattern, 
                                                    default_transect_spacing, 
                                                    default_transect_spacing_latitude)

transect_region_haul_key_all_ages = read_transect_region_haul_key(
    transect_region_filepath_all_ages,
    transect_region_sheetname_all_ages,
    transect_region_file_rename
)

transect_region_haul_key_no_age1 = read_transect_region_haul_key(
    transect_region_filepath_no_age1,
    transect_region_sheetname_no_age1,
    transect_region_file_rename
)

# TODO: Why include ALL of the empty rows
exports_with_regions_df = process_region_names(
    exports_df, 
    region_name_expr_dict, 
    CAN_haul_offset, 
)

# When there is no transect-region-haul key file 
transect_region_haul_key_no_age1 = generate_transect_region_haul_key(
    exports_with_regions_df, 
    filter_list=filter_list_no_age1
)

transect_region_haul_key_all_ages = generate_transect_region_haul_key(
    exports_with_regions_df, 
    filter_list=filter_list_all_ages
)

# Consolidate
df_nasc_no_age1 = consolidate_echvoiew_nasc(
    df_merged=exports_with_regions_df,
    interval_df=interval_template,
    region_class_names=filter_list_no_age1,
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_no_age1,
)

df_nasc_all_ages = consolidate_echvoiew_nasc(
    df_merged=exports_with_regions_df,
    interval_df=interval_template,
    region_class_names=filter_list_all_ages,
    impute_region_ids=True,
    transect_region_haul_key_df=transect_region_haul_key_all_ages,
)

# TODO: 
# --- Transect interval filtering 
# --- Stratify (KS)
# --- Reading in pre-consolidated files
pre_consolidated_export_all_ages = Path("C:/Users/Brandyn Lucca/Documents/Data/echopop_2019/Exports/US_CAN_detailsa_2019_table1y+_ALL_final - updated.xlsx")
file_sheetname = "Sheet1"
pre_consolidated_export_all_ages.exists()

# --- Column mapping [optional input; otherwise need to manually change column names -- NOT built-in]
FEAT_TO_ECHOPOP_COLUMNS = {
    "transect": "transect_num",
    "region id": "region_id",
    # "vl start": "distance_s",
    # "vl end": "distance_e",
    "vessel_log_start": "distance_s",
    "vessel_log_end": "distance_e",
    "spacing": "transect_spacing",
    "layer mean depth": "layer_mean_depth",
    "layer height": "layer_height",
    "bottom depth": "bottom_depth",
    "assigned haul": "haul_num"
}

# 
import pandas as pd
from typing import Optional
from echopop.nwfsc_feat.ingest_nasc import read_echoview_export, impute_bad_coordinates

# INPUTS
filename = pre_consolidated_export_all_ages
sheetname = file_sheetname
column_name_map: Optional[dict] = FEAT_TO_ECHOPOP_COLUMNS
validator = None
# -----

def read_nasc_file(filename, sheetname, column_name_map=None, validator=None):
    """
    Read NASC data from a consolidated XLSX file
    
    Parameters
    ----------
    filename : str or Path
        Path to the Excel file
    sheetname : str
        Name of the sheet to read
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names
    validator : callable, optional
        Function to validate the dataframe

    Examples
    --------
    >>> column_map = {"transect": "transect_num", "region id": "region_id"}
    >>> df = read_nasc_file("data.xlsx", "Sheet1", column_map)
    
    >>> # Without column mapping
    >>> df = read_nasc_file("data.xlsx", "Sheet1")
    
    Returns
    -------
    pandas.DataFrame
        Cleaned DataFrame with renamed columns and imputed coordinates
    """
    # Read in the defined file
    consolidated_file = pd.read_excel(filename, sheet_name=sheetname, index_col=None, header=0) 

    # Set column names to lowercase
    consolidated_file.columns = consolidated_file.columns.str.lower()

    # Rename columns
    if column_name_map:
        consolidated_file.rename(columns=column_name_map, inplace=True)

    # Fix latitude and longitude
    # ---- Latitude
    if "latitude" in consolidated_file.columns:
        impute_bad_coordinates(consolidated_file, "latitude")
    # ---- Longitude
    if "longitude" in consolidated_file.columns:
        impute_bad_coordinates(consolidated_file, "longitude")

    # Return the cleaned DataFrame
    return consolidated_file

pre_consolidated_export_all_ages = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/Exports/US_CAN_detailsa_2019_table1y+_ALL_final - updated.xlsx")
file_sheetname = "Sheet1"
read_nasc_file(pre_consolidated_export_all_ages, file_sheetname, column_name_map)
