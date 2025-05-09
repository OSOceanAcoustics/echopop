from pathlib import Path
from typing import Dict, Generator, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
import glob
import re

from echopop.kriging import Kriging
from echopop.nwfsc_feat import get_proportions, ingest_nasc, load_data
from echopop.nwfsc_feat.ingest_nasc import read_echoview_nasc, map_transect_num, validate_transect_exports, clean_echoview_cells_df
from echopop.ingest.common import read_csv_file
from echopop.core.echopop_columns import ECHOVIEW_TO_ECHOPOP
import functools
# ===========================================
# Organize NASC file
nasc_path: Path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/raw_nasc")
filename_transect_pattern: str = r"T(\d+)"
default_transect_spacing = 10. # nmi
default_transect_spacing_latitude = 60. # deg N

transect_region_filepath_all_ages: Path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/Stratification/US_CAN_2019_transect_region_haul_age1+ auto_final.xlsx")
transect_region_sheetname_all_ages: str = "Sheet1"
transect_region_filepath_no_age1: Path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/Stratification/US_CAN_2019_transect_region_haul_age2+ auto_20191205.xlsx")
transect_region_sheetname_no_age1: str = "Sheet1"
transect_region_file_rename: dict = {
    "tranect": "transect_num",
    "region id": "region_id",
    "trawl #": "haul_num",
}

CAN_haul_offset = 200
region_name_pattern = "{REGION_CLASS}{HAUL_NUM}{COUNTRY}"
region_name_expr_dict = {
    "REGION_CLASS": {
        "Age-1 Hake": "^(?:h1a(?![a-z]|1a))",
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

process_region_names(
    exports_df, 
    region_name_expr_dict, 
    CAN_haul_offset, 
)


pattern_dict = region_name_expr_dict
df = exports_df.copy()
# Get unique region names
unique_regions = pd.DataFrame({"region_name": df["region_name"].unique()})

# Compile patterns directly from the provided dictionary
regex_patterns = compile_patterns(pattern_dict)


# Apply valid types
valid_dtypes = {
    "region_class": str,
    "haul_num": float,
    "country": str,
}
# ---- Replace missing `haul_num` with 0
# if "haul_num" in unique_regions_coded:
#     unique_regions_coded["haul_num"] = unique_regions_coded["haul_num"]
# ---- Apply conversion
unique_regions_filtered = unique_regions_coded.apply(
    lambda col: col.astype(valid_dtypes.get(col.name, type(col.iloc[0])))
)
# ---- Apply haul number offset if so defined
if "country" in unique_regions_coded:
    unique_regions_filtered.loc[unique_regions_filtered["country"] == "CAN", "haul_num"] = (
        unique_regions_filtered.loc[unique_regions_filtered["country"] == "CAN", "haul_num"]
        + CAN_haul_offset
    )

# Map the regions-hauls
# ---- Index the data based on region name
transect_data_filtered = exports_df.set_index("region_name")
# ---- Map the values
transect_data_filtered.loc[:, unique_regions_filtered.columns] = unique_regions_filtered
# ---- Reset the index
transect_data_filtered.reset_index(inplace=True)
# ---- Drop NA's
# transect_data_filtered_test.dropna(subset=["region_class"], inplace=True)

filter_list = ["Age-1 Hake", "Age-1 Hake Mix", "Hake", "Hake Mix"]
# filter_list = ["Hake", "Hake Mix"]
# ---- Format region pattern
region_pattern = (
    rf"^(?:{'|'.join([re.escape(name.lower()) for name in filter_list])})"
)
# ---- Drop the NA's
# ---- Filter the dataframe
transect_region = transect_data_filtered[
    transect_data_filtered["region_class"].str.contains(
        region_pattern, case=False, regex=True
    )
]
# ---- Ensure that `region_class` is lowercase
transect_region.loc[:, "region_class"] = transect_region.loc[:, "region_class"].str.lower()
# ---- Find the unique keys
unique_regions_map = (
    transect_region.groupby(["transect_num", "haul_num", "region_id"])[
        ["region_class", "region_name"]
    ]
    .first()
    .reset_index()
    .sort_values(["haul_num"])
)
