# import copy
# import glob
# import os
# import re
# from pathlib import Path
# from typing import List, Literal, Optional, Tuple, Union

# import numpy as np
# import pandas as pd

# from echopop.analysis import (
#     acoustics_to_biology,
#     apportion_kriged_values,
#     krige,
#     process_transect_data,
#     stratified_summary,
# )
# from echopop.core import CONFIG_MAP, DATA_STRUCTURE, ECHOVIEW_EXPORT_MAP, REGION_EXPORT_MAP
# from echopop.core import BIODATA_HAUL_MAP
# from echopop.spatial.transect import export_transect_layers, export_transect_spacing
# from echopop.survey import Survey
# from echopop.utils import batch_load as ebl, load as el, message as em
# from echopop.utils.batch_load import (
#     consolidate_exports,
#     get_haul_transect_key,
#     get_transect_numbers,
#     load_export_regions,
#     write_haul_to_transect_key,
#     write_transect_region_key,
#     validate_echoview_exports,
#     validate_export_directories
# )
# from echopop.utils.operations import group_merge

# #
# init_config_path = (
#     "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
# )
# # NOTE: File configuration
# survey_year_config_path = (
#     "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
# )

# self = Survey

# construct_nasc: bool = True
# transect_pattern: Optional[str] = r"T(\d+)"
# index_variable: Union[str, List[str]] = ["transect_num", "interval"]
# unique_region_id: str = "region_id"
# region_class_column: str = "region_class"
# # # !: DELETE THESE
# # export_file_directory: Optional[Union[str, Path]] = None
# # export_save_directory: Optional[Union[str, Path]] = None
# # # !: ----

# # ###
# # # Initialize `meta` attribute
# self.meta = copy.deepcopy(DATA_STRUCTURE["meta"])

# # # Loading the configuration settings and definitions that are used to initialize the Survey
# # # class object
# self.config = el.load_configuration(Path(init_config_path), Path(survey_year_config_path))

# # NASC export file batch processing
# # if construct_nasc:
# #     print("Constructing consolidated NASC export files.")
# #     # ---- Batch processing
# #     ebl.batch_read_echoview_exports(
# #         self.config,
# #         transect_pattern,
# #         index_variable,
# #         unique_region_id,
# #         region_class_column,
# #         export_file_directory,
# #         export_save_directory,
# #     )
# configuration_dict = self.config
# transect_pattern = transect_pattern
# index_variable = index_variable
# unique_region_id = unique_region_id
# region_class_column = region_class_column
# # !: DELETE THESE
# file_directory = export_file_directory
# save_directory = export_save_directory
# # !: ----
# configuration_dict = self.config


# # Extract the placeholder names from a pattern string
# def extract_placeholders(pattern_template):
#     placeholder_pattern = re.compile(r"\{(\w+)\}")
#     matches = placeholder_pattern.findall(pattern_template)
#     return matches


# pattern_config = configuration_dict["transect_region_mapping"]["parts"]
# template = configuration_dict["transect_region_mapping"]["pattern"]

# template_parts = extract_placeholders(template)
