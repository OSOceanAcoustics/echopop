import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union
from itertools import product 
import numpy as np
import pandas as pd
import yaml
import re
import glob
import os
import warnings
import pytest
from lmfit import Minimizer, Parameters
from scipy import special
from IPython.display import display
from pandera import DataFrameModel, Field, check
from pandera.errors import SchemaError, SchemaErrors
from pandera.typing import Series
from pydantic import BaseModel, Field, RootModel, ValidationError, field_validator, model_validator
import json 

from echopop.core import (
    BIODATA_HAUL_MAP, 
    DATA_STRUCTURE, 
    LAYER_NAME_MAP, 
    NAME_CONFIG, 
    REGION_EXPORT_MAP,
    ECHOVIEW_TO_ECHOPOP_NAMES, 
)
from echopop.utils.validate_df import (
    DATASET_DF_MODEL,
    AcousticData,
    BaseDataFrame,
    CatchBiodata,
    GeoStrata,
    IsobathData,
    KrigedMesh,
    KSStrata,
    LengthBiodata,
    SpecimenBiodata,
    VarioKrigingPara,
)
from echopop.utils.validate_dict import (
    BiologicalFiles,
    CONFIG_DATA_MODEL,
    CONFIG_INIT_MODEL,
    CSVFile,
    XLSXFile,
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
    StratificationFiles,
    SpeciesDefinition,
    KrigingFiles,
    TransectRegionMap, 
    PatternParts, 
    INPFCRegionMap, 
    InputModel
)
from echopop.analysis import (
    acoustics_to_biology,
    apportion_kriged_values,
    krige,
    process_transect_data,
    stratified_summary,
    variogram_analysis,
    back_calculate_abundance_nasc,
)
from echopop.utils.load import (
    dataset_integrity,
    map_imported_datasets, 
    read_validated_data, 
    prepare_input_data, 
    preprocess_acoustic_biology_spatial, 
    preprocess_acoustic_spatial, 
    preprocess_biodata, 
    preprocess_biology_spatial, 
    preprocess_spatial, 
    preprocess_statistics
)
from echopop.acoustics import (
    aggregate_sigma_bs, 
    nasc_to_biomass
)
from echopop.biology import (
    age1_metric_proportions,
    distribute_length_age,
    filter_species,
    fit_length_weight_relationship,
    fit_length_weights,
    impute_kriged_values,
    number_proportions,
    partition_transect_age,
    quantize_number_counts,
    quantize_weights,
    reallocate_kriged_age1,
    weight_proportions,
    age1_metric_proportions
)
from echopop.spatial.transect import (
    correct_transect_intervals,
    edit_transect_columns,
    export_transect_layers, 
    export_transect_spacing,
    save_transect_coordinates,
    summarize_transect_strata,
    transect_spatial_features,
    transect_bearing,
    transect_array,
    transect_coordinate_centroid,
    transect_extent
)
from echopop.spatial.variogram import (
    empirical_variogram,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
    VARIOGRAM_MODELS,
    prepare_variogram_matrices,
    quantize_lags,
    semivariance,
    variogram_matrix_filter,
)
from echopop.utils.operations import (
    compile_patterns, 
    extract_parts_and_labels, 
    group_merge,
)
from echopop.spatial.mesh import (
    crop_mesh, 
    griddify_lag_distances,
    mesh_to_transects, 
    stratify_mesh,
    interpolate_survey_extent
)
from echopop.utils.validate import (
    posfloat, 
    posint, 
    realcircle, 
    realposfloat
)
from echopop.utils.data_structure_utils import map_imported_datasets
from echopop.spatial.projection import transform_geometry
from echopop.spatial.krige import kriging

from echopop.spatial.projection import transform_geometry
from echopop.statistics import stratified_transect_statistic

from echopop.graphics import plotting as egp, variogram_interactive as egv
from echopop.utils import load as el, load_nasc as eln, message as em
import echopop.spatial.variogram as esv
from echopop.tests.conftest import assert_dataframe_equal, load_json_data
from echopop.utils.data_structure_utils import map_imported_datasets
from echopop.utils.load_nasc import (
    read_echoview_export, 
    export_transect_layers,
    export_transect_spacing,
    load_export_regions,
    filter_export_regions,
    ingest_echoview_exports,
    validate_export_directories,
    compile_patterns,
    consolidate_exports,
    construct_transect_region_key,
    get_haul_strata_key,
    get_transect_numbers,
)
from echopop.extensions.feat_report import (
    format_file_sheet,
    initialize_workbook,
    format_table_headers,
    append_datatable_rows,
    append_sheet_label,
    append_table_aggregates,
    pivot_aged_weight_proportions,
    pivot_dataframe_data,
    pivot_haul_tables,
    repivot_table,
    write_age_length_table_report,
    write_aged_dataframe_report,
    write_haul_report,
    FEATReports
)
from echopop.survey import Survey

SURVEY_YEAR = 2021

# Initialization configuration
init_config_path = f"C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_\
config_{SURVEY_YEAR}.yml"

# Filepath/dataset configuration
survey_year_config_path = f"C:/Users/Brandyn/Documents/GitHub/echopop/config_files\
/survey_year_{SURVEY_YEAR}_config.yml"

# Load json settings
# ---- File open
with open(Path(os.getcwd() + "\\echopop\\compatibility_parameters_test.json").as_posix()) as f:
    json_dict = json.load(f)
# ---- Load
parameters = json_dict[f"{SURVEY_YEAR}"]

####################################################################################################
# Run
# ---- Create object
survey = Survey(init_config_path, survey_year_config_path)
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis()
survey.kriging_analysis()
self = survey
cropping_parameters: Dict[str, Any] = {}
kriging_parameters: Dict[str, Any] = {}
coordinate_transform: bool = True
extrapolate: bool = False
best_fit_variogram: bool = False
variable: Literal["biomass"] = "biomass"
variogram_parameters: Optional[Dict[str, Any]] = None
verbose: bool = True

aged_apportioned["abundance_apportioned"].sum().sum() + unaged_apportioned["abundance_apportioned_unaged"].sum().sum()

survey.results["kriging"]

input_dict, analysis_dict, settings_dict = self.input, self.analysis, self.analysis["settings"]["kriging"]

mesh_results_df, transect_dict, settings_dict = mesh_results.reset_index(), analysis_dict["transect"], settings_dict

idx = np.argmax(mesh_results_df.biomass)
mesh_results_df.reset_index().loc[idx]
abundance = survey.results["kriging"]["mesh_results_df"].loc[idx, "abundance"]
stratum = survey.results["kriging"]["mesh_results_df"].loc[idx, "stratum_num"]
survey.results["kriging"]["mesh_results_df"].abundance.sum()

N_len = proportions_dict["weight"]["aged_unaged_sex_weight_proportions_df"].set_index(["stratum_num", "sex"])
Aged_unaged = proportions_dict["weight"]["aged_unaged_weight_proportions_df"].set_index(["stratum_num"])

N_len.loc[stratum, "weight_proportion_unaged"] * Aged_unaged.loc[stratum, "unaged_proportions"]

proportions_dict["weight"]["unaged_weight_proportions_df"].loc[lambda x: x["stratum_num"] == stratum]
aged_apportioned_abundance.pivot_table(columns=["age_bin"], index=["sex", "length_bin"], observed=False)


analysis_dict, kriged_mesh, settings_dict = self.analysis, kriged_results["mesh_results_df"], self.analysis["settings"]["kriging"]
((summed_abundance * aged_abundance_proportions_pvt.transpose()).fillna(0.0))
aged_apportioned_abundance

# Sum the kriged weights for each stratum
# ---- Extract stratum column name
stratum_col = settings_dict["stratum_name"]
# ---- Extract the biological variable (independent of area)
# biology_col = settings_dict["variable"].replace("_density", "")
# ---- Sum biomass for each stratum
summed_biomass = kriged_mesh.groupby([stratum_col], observed=False)["biomass"].sum()
# ---- Sum abundance for each stratum
summed_abundance = kriged_mesh.groupby([stratum_col], observed=False)["abundance"].sum()

# Extract the weight proportions from the analysis object
proportions_dict = analysis_dict["transect"]["biology"]["proportions"]
# ---- Aged (abundance and biomass)
aged_abundance_proportions = proportions_dict["number"]["aged_length_proportions_df"].copy()
aged_biomass_proportions = proportions_dict["weight"]["aged_weight_proportions_df"].copy()
# ---- Unaged (abundance and biomass)
unaged_abundance_proportions = proportions_dict["number"]["unaged_length_proportions_df"].copy()
unaged_biomass_proportions = proportions_dict["weight"]["unaged_weight_proportions_df"].copy()
# ---- Aged-unaged sexed weight proportions
aged_unaged_sex_biomass_proportions = proportions_dict["weight"][
    "aged_unaged_sex_weight_proportions_df"
].copy()[[stratum_col, "sex", "weight_proportion_overall_unaged"]]
# ---- Aged-unaged sex number proportions
aged_unaged_sex_abundance_proportions = proportions_dict["number"]["sex_proportions_df"].copy()

# Compute the apportioned unaged kriged biological values per stratum
# ---- Merge the unaged abundance proportions

df = unaged_abundance_proportions.pivot_table(columns=["stratum_num"], 
                                              index=["sex", "length_bin"], 
                                              values=["proportion_number_unaged", 
                                                      "proportion_number_overall_unaged"],
                                              observed=False)

a = unaged_biomass_proportions.pivot_table(index=[stratum_col], columns=["length_bin"], values="weight_proportion", observed=False)
b = unaged_sex_biomass_proportions.pivot_table(columns=[stratum_col, "sex"], values="weight_proportion_overall_unaged", observed=False)
b["male"] * a.transpose()
b.stack().transpose()
(a.stack() * b.stack()).stack().sum(axis=1)
b * a
b
kriged_output[kriged_output.sex=="male"].biomass_apportioned.sum() + kriged_output[kriged_output.sex=="female"].biomass_apportioned.sum()
survey.results["kriging"]["tables"]["aged_tbl"]["abundance_apportioned"].sum().sum()
aged_apportioned["abundance_apportioned"].sum().sum()
survey.results["kriging"]["tables"]["unaged_tbl"]["abundance_apportioned_unaged"].sum().sum()
unaged_apportioned["abundance_apportioned_unaged"].sum().sum()
unaged_biomass_proportions.set_index([stratum_col], inplace=True)
unaged_sex_biomass_proportions.set_index([stratum_col], inplace=True)
kriged_apportioned_table["biomass_apportioned"].sum().sum()
aged_pivot.loc["female"].transpose().sum(axis=1)
unaged_biomass_proportions

summed_abundance
df3 = unaged_abundance_proportions.pivot_table(columns=["sex", "length_bin"], 
                                               index=[stratum_col], 
                                               values="proportion_number_overall_unaged", 
                                               observed=False)
(summed_abundance * df3.transpose()).fillna(0.0).sum(axis=1).loc["all"]


df3 = aged_abundance_proportions.pivot_table(columns=[stratum_col, "age_bin"],
                                       index=["sex", "length_bin"],
                                       values="proportion_number_overall_aged",
                                       observed=False).loc[("male", ), 7]

(abundance * df3) * 1e-5


df3 = aged_abundance_proportions.pivot_table(
    columns=["sex", "age_bin", "length_bin"],
    index=[stratum_col],
    values = "proportion_number_overall_aged",
    observed=False,
)

aged_apportioned_abundance.loc[aged_apportioned_abundance.sex == "all"].abundance.sum() + unaged_apportioned_abundance.loc[unaged_apportioned_abundance.sex == "all"].abundance.sum() 

(summed_abundance * df3.transpose()).fillna(0.0).sum(axis=1).loc["all"].sum()

summed_abundance.sum() - ant.loc["all"].sum()

ant = (df3.transpose() * summed_abundance).fillna(0.0).sum(axis=1)
ant.loc["all"].sum()
ant = (summed_abundance * df3.transpose()).fillna(0.0).sum(axis=1).loc["all"].reset_index(name="abundance")
ant.abundance.sum()

ant.pivot_table(columns=["age_bin"], index=["length_bin"], aggfunc="sum")
7.057882e+05 * 1e-5
7012856 * 1e-5
(summed_abundance * df3.transpose()).fillna(0.0).sum(axis=1).loc["all"].reset_index(name="abundance").pivot(columns=["age_bin"], index=["length_bin"]).sum(axis=0) * 1e-9
ant = (summed_abundance * df3.transpose()).fillna(0.0)[7].reset_index(name="abundance").pivot(columns=["age_bin"], index=["sex", "length_bin"], values="abundance")
ant.loc["male"]
df3.transpose() * summed_abundance

summed_abundance * df3.transpose()

df1 = aged_unaged_sex_abundance_proportions.set_index([stratum_col, "sex"])
N_len = abundance * df1.loc[stratum, "proportion_number_overall_unaged"]
N_len_age = abundance * df1.loc[stratum, "proportion_number_overall_aged"]
df1 = aged_abundance_proportions.pivot_table(
    index=["sex", "length_bin"],
    columns=[stratum_col],
    values=["proportion_number_aged",
            "proportion_number_overall_aged"],
    observed=False    
)

N_len

df2 = df.stack(future_stack=True).reset_index().pivot_table(
    columns=[stratum_col, "length_bin"], 
    index=["sex"],
    values=["proportion_number_overall_unaged", "proportion_number_unaged"],
    observed=True)

N_len * df2["proportion_number_unaged"].loc[:, stratum].transpose()
(abundance * df2["proportion_number_overall_unaged"].loc[:, stratum]).transpose()


N_len.index
df2["proportion_number_unaged"].loc[:, stratum].index
df1

df.loc["male"].max()


abundance_mat.loc["male"].max()


df.loc["all"].sum()


df.loc[("male", ), ("proportion_number_overall_unaged", 7)]

df1 = unaged_sexed_biomass_apportioned.pivot_table(columns=["stratum_num"], 
                                                   index=["sex", "length_bin"], 
                                                   values=["weight_proportion", 
                                                           "weight_proportion_overall_unaged"],
                                                   observed=False)

df1.loc["male"].sum()

unaged_sexed_abundance_apportioned = unaged_abundance_proportions.merge(
    aged_unaged_sex_abundance_proportions
)
# ---- Merge the unaged biomass proportions
unaged_sexed_biomass_apportioned = unaged_biomass_proportions.merge(
    aged_unaged_sex_biomass_proportions
)
# ---- Set index to stratum column
unaged_sexed_biomass_apportioned.set_index([stratum_col], inplace=True)
# ---- Set the index based on `summed_abundance`
summed_abundance_indexed = summed_abundance.reindex(unaged_sexed_apportioned.index)
# ---- Append the stratum-aggregated abundance values
unaged_sexed_apportioned["abundance_apportioned_unaged"] = (
    unaged_sexed_apportioned["weight_proportion"]
    * unaged_sexed_apportioned["weight_proportion_overall_unaged"]
    * summed_abundance_indexed
)
# ---- Set the index based on `summed_biomass`
summed_biomass_indexed = summed_biomass.reindex(unaged_sexed_apportioned.index)
# ---- Append the stratum-aggregated biomass values
unaged_sexed_apportioned["biomass_apportioned_unaged"] = (
    unaged_sexed_apportioned["weight_proportion"]
    * unaged_sexed_apportioned["weight_proportion_overall_unaged"]
    * summed_biomass_indexed
)

# Distribute biological values over the overall proportions (i.e. relative to aged and unaged
# fish) for aged fish
# ---- Set index to stratum column
aged_proportions.set_index([stratum_col], inplace=True)
# ---- Compute the distributed abundance values
aged_proportions["abundance_apportioned"] = (
    aged_proportions["weight_proportion_overall"] * summed_abundance
).fillna(0.0)
# ---- Compute the distributed biomass values
aged_proportions["biomass_apportioned"] = (
    aged_proportions["weight_proportion_overall"] * summed_biomass
).fillna(0.0)

# Distribute the aged biological distributions over unaged length distributions to estimate
# aged distributions
# ---- Pivot aged data
aged_pivot = aged_proportions.reset_index().pivot_table(
    index=["sex", "length_bin"],
    columns=["age_bin"],
    values=["abundance_apportioned", "biomass_apportioned"],
    aggfunc="sum",
    observed=False,
)
# ---- Calculate the total biomass values for each sex per length bin
aged_length_biomass_totals = aged_pivot["biomass_apportioned"].sum(axis=1).unstack("sex")
# ---- Pivot unaged data
unaged_pivot = unaged_sexed_apportioned.reset_index().pivot_table(
    index=["length_bin"],
    columns=["sex"],
    values=["abundance_apportioned_unaged", "biomass_apportioned_unaged"],
    aggfunc="sum",
    observed=False,
)
# ---- Calculate the new unaged biomass values distributed over age
unaged_apportioned_biomass_values = (
    unaged_pivot["biomass_apportioned_unaged"]
    * aged_pivot.unstack("sex")["biomass_apportioned"]
    / aged_length_biomass_totals
).fillna(0)

# Imputation is required when unaged values are present but aged values are absent at shared
# length bins! This requires an augmented implementation to address this accordingly
# ---- Biomass
kriged_full_table = impute_kriged_values(
    aged_pivot["biomass_apportioned"],
    unaged_pivot["biomass_apportioned_unaged"],
    aged_length_biomass_totals,
    unaged_apportioned_biomass_values,
    settings_dict,
    variable="biomass",
)

# Additional reapportionment if age-1 fish are excluded
if settings_dict["exclude_age1"]:
    # ---- Re-allocate biomass
    kriging_full_table = reallocate_kriged_age1(
        kriged_full_table, settings_dict, variable="biomass_apportioned"
    )
    # ---- Stack the aged-pivot table
    aged_data = (
        aged_pivot["abundance_apportioned"].stack().reset_index(name="abundance_apportioned")
    )
    # ---- Re-allocate abundance
    aged_table = reallocate_kriged_age1(
        aged_data, settings_dict, variable="abundance_apportioned"
    )
    # ---- Re-pivot
    aged_pivot["abundance_apportioned"] = aged_table.pivot_table(
        index=["sex", "length_bin"],
        columns=["age_bin"],
        values="abundance_apportioned",
        observed=False,
    )
    # ---- Validate that apportioning age-1 values over all adult values did not 'leak'
    # -------- Previous apportioned totals by sex
    previous_totals = kriged_full_table.groupby(["sex"])["biomass_apportioned"].sum()
    # -------- New apportioned totals by sex
    new_totals = kriging_full_table.groupby(["sex"])["biomass_apportioned"].sum()
    # -------- Check (1 kg tolerance)
    if np.any((previous_totals - new_totals) > 1e-6):
        warnings.warn(
            "Apportioned kriged apportioned biomass for age-1 not fully distributed over all "
            "age-2+ age bins."
        )
    # ----- Return
    kriged_output = kriging_full_table.copy()
else:
    kriged_output = kriged_full_table.copy()

# Check equality between original kriged estimates and (imputed) apportioned estimates
if (kriged_output["biomass_apportioned"].sum() - summed_biomass.sum()) > 1e-6:
    # ---- If not equal, generate warning
    warnings.warn(
        "Apportioned kriged apportioned biomass does not equal the total kriged mesh "
        "apportioned biomass! Check for cases where kriged values may only be present in aged "
        "(`self.results['kriging']['tables']['aged_tbl']`) or unaged ("
        "(`self.results['kriging']['tables']['unaged_tbl']`) distributions for each sex."
    )

    # Return output
    # return aged_pivot, unaged_pivot, kriged_output

survey.generate_reports(reports=["kriged_length_age_abundance", "kriged_length_age_biomass"])
self = FEATReports(survey, reports=["kriged_length_age_abundance", "kriged_length_age_biomass"])
title = "Kriged Acoustically Weighted Abundance ({SEX})"
filename = "kriged_len_age_abundance_table.xlsx"

# Get the dataset
dataset = self.data.results["kriging"]["tables"]
# ---- Aged data
aged_data = dataset["aged_tbl"]["abundance_apportioned"].copy()
# ---- Unaged data
unaged_data = dataset["unaged_tbl"]["abundance_apportioned_unaged"].copy()

# Process the aged data
# ---- Stack
aged_stk = aged_data.stack(future_stack=True).reset_index(name="abundance")
# ---- Expand to include 'all'
full_aged_stk = pd.concat(
    [
        aged_stk,
        aged_stk.groupby(["length_bin", "age_bin"], observed=False)["abundance"]
        .sum()
        .reset_index()
        .assign(sex="all"),
    ],
    ignore_index=True,
)
# ---- Re-pivot
full_aged_pvt = full_aged_stk.pivot_table(
    index=["sex", "length_bin"], columns=["age_bin"], values="abundance", observed=False
)

# Process the unaged data
# ---- Stack
unaged_stk = unaged_data.stack(future_stack=True).reset_index(name="abundance")
# ---- Expand to include 'all'
full_unaged_stk = pd.concat(
    [
        unaged_stk,
        unaged_stk.groupby(["length_bin"], observed=False)["abundance"]
        .sum()
        .reset_index()
        .assign(sex="all"),
    ],
    ignore_index=True,
)
# ---- Re-pivot
full_unaged_pvt = full_unaged_stk.pivot_table(
    index=["sex", "length_bin"], values="abundance", observed=False
)

# Repivot the datasets for each sex
tables = {
    sex: repivot_table(
        age_length_dataframe=full_aged_pvt.loc[sex, :],
        length_dataframe=full_unaged_pvt.loc[sex, :],
        variable="abundance",
    )
    for sex in ["all", "male", "female"]
}

# Get filepath
filepath = self.save_directory / filename

# Write the *.xlsx sheet
# write_age_length_table_report(tables, title, filepath)
tables_dict = tables
table_title = title

# Get sheetnames
sheetnames = {
    "male": "Sheet1",
    "female": "Sheet2",
    "all": "Sheet3",
}

sex = "all"
# Subset the data dictionary for the particular sheet
sheet_data = tables_dict[sex]

survey.results["kriging"]["mesh_results_df"]["abundance"].sum()

sheet_data.loc["Subtotal"].sum() * 1e-9


ingest_exports="echoview"
read_transect_region_file=parameters["read_transect_region_file"]
write_transect_region_file=parameters["write_transect_region_file"]
index_variable: Union[str, List[str]] = ["transect_num", "interval"]
region_class_column: str = "region_class"
transect_pattern: str = r"T(\d+)"
unique_region_id: str = "region_id"
verbose: bool = True
configuration_dict = survey.config

# Get NASC export settings
export_settings = configuration_dict["nasc_exports"]

# Validate relevant directories and file existence
save_folder, file_directory, export_files = validate_export_directories(configuration_dict)

# Get the transect numbers
transect_reference = get_transect_numbers(export_files, transect_pattern, file_directory)

# Read in and consolidate all export files (intervals, layers, cells)
transect_data, interval_template = consolidate_exports(
    transect_reference, export_settings["max_transect_spacing"]
)
# ---- Rename a column
interval_template = interval_template.rename(columns={"max_depth": "bottom_depth"})

# ---- Transect-region mapping file metadata
region_files = configuration_dict["export_regions"]
# ---- Region names
region_names = export_settings["regions"]
# ---- Root directory
root_dir = configuration_dict["data_root_dir"]
# ---- Read in the file
# transect_region_key = load_export_regions(region_files, region_names, root_dir)
root_directory = root_dir
# ---- Read in the file
transect_region_key = load_export_regions(region_files, region_names, root_dir)
# ---- Get the country-code information, if parameterized
if "inpfc_strata_region" in configuration_dict["transect_region_mapping"]:
    # ---- Assign country
    country_codes = get_haul_strata_key(configuration_dict, root_dir, transect_region_key)
    # ---- Get columns present that can be used for indexing
    cols_idx = list(
        set(transect_region_key.columns).intersection(
            set(["haul_num", "transect_num", "region_id", "region_name", "region_class"])
        )
    )
    # ---- Set index
    transect_region_key.set_index(cols_idx, inplace=True)
    # ---- Assign
    transect_region_key["country"] = country_codes
    # ---- Reset index
    transect_region_key.reset_index(inplace=True)
    
b = output_nasc.groupby(["transect_num", "region_id"])["nasc"].sum()
b
output_nasc
transect_data[transect_data.transect_num == 0].region_id.unique()

survey.results

transect_pattern

re.search("T(\\d+(?:\\.\\d+)?)(?=-)", "V160-S201103-X2-F38-T0.3-Z0- (analysis)").group(1)
re.search(parameters["transect_pattern"].replace(r"\\", "\\"), "V160-S201103-X2-F38-T0.3-Z0- (analysis)").group(1)

survey.kriging_analysis(extrapolate=False, 
                        variogram_parameters={"model": ["bessel", "exponential"],
                                              "n_lags": 30,
                                              "sill": 0.94,
                                              "correlation_range": 0.008,
                                              "search_radius": 0.024,
                                              "anisotropy": 0.001},
                        cropping_parameters={"crop_method": "convex_hull",
                                             "num_nearest_transect": 2})

survey.results["kriging"]["mesh_results_df"].nasc.sum()
survey.results["kriging"]["mesh_results_df"].biomass.sum()
survey.analysis["transect"]["acoustics"]["adult_transect_df"].nasc.sum()

survey.input["acoustics"]["nasc_df"].haul_all_ages.unique()

survey.input["acoustics"]["nasc_df"].loc[lambda x: ~x.haul_no_age1.isin([1, 7, 9, 11, 12, 13, 15, 26])]
survey.analysis["transect"]["acoustics"]["adult_transect_df"].loc[lambda x: x.haul_num not in [1, 7, 9, 11, 12, 13, 15, 26]]

survey.summary("transect")

self = survey
cropping_parameters: Dict[str, Any] = {"crop_method": "convex_hull"}
kriging_parameters: Dict[str, Any] = {}
coordinate_transform: bool = True
extrapolate: bool = False
best_fit_variogram: bool = False
variable: Literal["biomass"] = "biomass"
variogram_parameters: Optional[Dict[str, Any]] = None
verbose: bool = True

input_dict, analysis_dict, settings_dict = self.input, self.analysis, self.analysis["settings"]["kriging"]

transect_data, mesh_data, cropping_parameters = transect_data, mesh_data, validated_cropping_methods

transect_data, mesh_data, cropping_parameters = transect_data.copy(), mesh, cropping_parameters