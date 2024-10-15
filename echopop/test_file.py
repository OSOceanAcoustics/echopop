import copy
import os
import re
from pathlib import Path
from typing import List, Literal, Optional, Union

import numpy as np
import pandas as pd

from echopop.core import (
    BIODATA_HAUL_MAP,
    DATA_STRUCTURE,
    ECHOVIEW_EXPORT_MAP,
    LAYER_NAME_MAP,
    NAME_CONFIG,
    REGION_EXPORT_MAP,
)
from echopop.spatial.transect import export_transect_layers, export_transect_spacing
from echopop.survey import Survey
from echopop.utils.data_structure_utils import map_imported_datasets
from echopop.utils.load import prepare_input_data, read_validated_data
from echopop.utils.load_nasc import (
    compile_patterns,
    consolidate_exports,
    construct_transect_region_key,
    filter_export_regions,
    get_haul_transect_key,
    get_transect_numbers,
    validate_echoview_exports,
    validate_export_directories,
)
from echopop.utils.operations import compile_patterns, extract_parts_and_labels, group_merge
from echopop.utils.validate_df import DATASET_DF_MODEL, KSStrata

survey = Survey(
    init_config_path="C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config.yml",
    survey_year_config_path="C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml",
)
survey.load_survey_data()
survey.transect_analysis()
survey.stratified_analysis()

self = survey
index_variable: Union[str, List[str]] = ["transect_num", "interval"]
ingest_exports: Optional[Literal["echoview", "echopype"]] = "echoview"
region_class_column: str = "region_class"
transect_pattern: str = r"T(\d+)"
unique_region_id: str = "region_id"
verbose: bool = True
configuration_dict = self.config
input_dict = self.input
dataset_type = ["biological", "kriging", "stratification"]


def dataset_integrity(
    input_dict: dict, analysis: Literal["transect", "stratified", "variogram", "kriging"]
):
    """
    Determine whether all of the necessary datasets are contained within the `Survey`-class object
    for each analysis

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    analysis: Literal["transect", "stratified", "variogram", "kriging"]
        The name of the analysis to be performed
    """

    # Map the imported datasets
    imported_data = map_imported_datasets(input_dict)

    # Initialize `missing`
    missing = None

    # Transect analysis
    if analysis == "transect":
        try:
            assert set(imported_data).issubset(["acoustics", "biology", "spatial"])
        except AssertionError as e:
            # ---- Collect missing values
            missing = list(set(["acoustics", "biology", "spatial"]) - set(imported_data))
            # ---- Missing analysis str
            missing_analysis_method = "'transect_analysis'"

    # Stratified analysis (transect)
    if analysis == "stratified:transect":
        try:
            assert set(imported_data).issubset(["acoustics", "spatial"])
        except AssertionError as e:
            # ---- Collect missing values
            missing = list(set(["acoustics", "spatial"]) - set(imported_data))
            # ---- Missing analysis str
            missing_analysis_method = "'stratified_analysis' (for transect data)"

    # Stratified analysis (kriging)
    if analysis == "stratified:kriging":
        try:
            assert set(imported_data).issubset(["acoustics", "spatial", "statistics"])
        except AssertionError as e:
            # ---- Collect missing values
            missing = list(set(["acoustics", "spatial", "statistics"]) - set(imported_data))
            # ---- Missing analysis str
            missing_analysis_method = "'stratified_analysis' (for kriged data)"

    # Kriging analysis
    if analysis == "kriging":
        try:
            assert set(imported_data).issubset(["acoustics", "biology", "spatial", "statistics"])
        except AssertionError as e:
            # ---- Collect missing values
            missing = list(
                set(["acoustics", "biology", "spatial", "statistics"]) - set(imported_data)
            )
            # ---- Missing analysis str
            missing_analysis_method = "'kriging_analysis'"

    # Variogram analysis
    if analysis == "variogram":
        try:
            assert set(imported_data).issubset(["acoustics", "spatial", "statistics"])
        except AssertionError as e:
            # ---- Collect missing values
            missing = list(set(["acoustics", "spatial", "statistics"]) - set(imported_data))
            # ---- Missing analysis str
            missing_analysis_method = "'fit_variogram'/'variogram_gui'"

    if missing:
        # ---- Format string
        missing_str = ", ".join([f"'{i}'" for i in missing])
        # ---- Raise error
        raise ValueError(
            f"The following datasets are missing for {missing_analysis_method}: {missing_str}."
        )
