from pathlib import Path
from typing import Dict, Generator, List, Optional, Tuple, Union

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
