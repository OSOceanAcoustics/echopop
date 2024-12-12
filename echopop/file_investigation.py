import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union
from itertools import product 

import numpy as np
import pandas as pd
import yaml
import re

from echopop.core import BIODATA_HAUL_MAP, DATA_STRUCTURE, LAYER_NAME_MAP, NAME_CONFIG, REGION_EXPORT_MAP
from echopop.utils.data_structure_utils import map_imported_datasets
from echopop.utils.validate_df import DATASET_DF_MODEL
from echopop.utils.validate_dict import CONFIG_DATA_MODEL, CONFIG_INIT_MODEL
from IPython.display import display

from echopop.analysis import (
    acoustics_to_biology,
    apportion_kriged_values,
    krige,
    process_transect_data,
    stratified_summary,
    variogram_analysis,
)
from echopop.graphics import plotting as egp, variogram_interactive as egv
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import edit_transect_columns
from echopop.utils import load as el, load_nasc as eln, message as em
from echopop.utils.load import map_imported_datasets, read_validated_data, prepare_input_data, preprocess_acoustic_biology_spatial, preprocess_acoustic_spatial, preprocess_biodata, preprocess_biology_spatial, preprocess_spatial, preprocess_statistics
from echopop.utils.load import dataset_integrity
import numpy as np
import pandas as pd
from echopop.acoustics import aggregate_sigma_bs, nasc_to_biomass
from echopop.biology import (
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
from echopop.spatial.krige import kriging
from echopop.spatial.mesh import crop_mesh, mesh_to_transects, stratify_mesh
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import (
    edit_transect_columns,
    save_transect_coordinates,
    summarize_transect_strata,
    transect_spatial_features,
)
from echopop.spatial.variogram import (
    empirical_variogram,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
)
from echopop.statistics import stratified_transect_statistic
from echopop.utils.validate_dict import (
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
)
from echopop.biology import age1_metric_proportions
from echopop.spatial.transect import correct_transect_intervals

from pandera import DataFrameModel, Field, check
from pandera.errors import SchemaError, SchemaErrors
from pandera.typing import Series
from echopop.utils.validate_df import (
    extract_errors,
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
import glob
import os
import re
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd

from echopop.spatial.transect import export_transect_layers, export_transect_spacing
from echopop.utils.validate_df import KSStrata
from echopop.utils.operations import compile_patterns, extract_parts_and_labels, group_merge
from echopop.utils.load_nasc import (
    load_export_regions,
    validate_export_directories,
    consolidate_exports,
    construct_transect_region_key,
    filter_export_regions,
    compile_patterns,
    get_transect_numbers,
    read_echoview_export,
    load_export_regions,
    export_transect_layers,
    extract_parts_and_labels,
    get_haul_strata_key,
)

import inspect
import warnings
from typing import Any, Dict, List, Optional, Union

from lmfit import Minimizer, Parameters
from scipy import special
from echopop.utils.validate_dict import VariogramBase, VariogramInitial, VariogramOptimize
from echopop.spatial.mesh import griddify_lag_distances
from pydantic import BaseModel, Field, RootModel, ValidationError, field_validator, model_validator
from echopop.utils.validate import posfloat, posint, realcircle, realposfloat
from echopop.utils.validate_dict import *
FutureWarning: The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. To retain the old behavior, exclude the relevant entries before the concat operation.
# LOAD =============================================================================================
new_datasets = ["biological", "kriging", "stratification"]
dataset_type = new_datasets

a = pd.cut(
    np.unique(input_dict["spatial"]["inpfc_strata_df"]["northlimit_latitude"]) * 0.99, 
    latitude_bins,
)
len(a)
len(input_dict["spatial"]["inpfc_strata_df"]["stratum_inpfc"].unique())
pd.cut(
    input_dict["spatial"]["inpfc_strata_df"]["northlimit_latitude"].unique() * 0.99, 
    latitude_bins,
    # labels=input_dict["spatial"]["inpfc_strata_df"]["stratum_inpfc"].unique()
)

self = survey
input_dict = self.input
configuration_dict = self.config

init_config_path = Path(init_config_path)
survey_year_config_path = Path(survey_year_config_path)

dummy = {
    "stratified_survey_mean_parameters": {"strata_transect_proportion": 0.85, 
                                          "num_replicates": 10, 
                                          "mesh_transects_per_latitude": 4},
    "kriging_parameters": {"A0": 6.25, 
                           "longitude_reference": 0.0,
                           "longitude_offset": 5.0,
                           "latitude_offset": 5.0},
    "bio_hake_age_bin": [1, 20, 1],
    "bio_hake_len_bin": [1, 20, 1],
    "TS_length_regression_parameters": {
        "pacific_hake": {"number_code": 22500,
                         "TS_L_slope": 20.0,
                         "TS_L_intercept": -70.0,
                         "length_units": "cm"}},
    "geospatial": {"init": "epsg:4326"}
}            

input_dict["biology"][keys].set_index(["haul_num"], inplace=True)
input_dict["biology"][keys]["stratum_num"] = strata_df["stratum_num"]
input_dict["biology"][keys].loc[1]
input_dict["biology"][keys]["stratum_num"] = strata_df["stratum_num"]
TSLRegressionParameters.create()
CONFIG_INIT_MODEL(init_config_path, **dummy)
input_dict["biology"][keys]["haul_bin"].replace(np.nan, pd.Categorical([0, 0]))
np.isnan(input_dict["biology"][keys].haul_bin[0])
input_dict["spatial"]["inpfc_strata_df"] = input_dict["spatial"]["inpfc_strata_df"].filter(["haul_start", "haul_end", "northlimit_latitude", "stratum_inpfc", "latitude_interval"])
survey.input["biology"]

nan_mask = input_dict["biology"][keys]['haul_bin'].isna()

# Create a placeholder for out-of-range values
out_of_range_label = "Out of Range"

# Add 'Out of Range' category if it's not already present
current_categories = input_dict["biology"][keys]['haul_bin'].cat.categories
if out_of_range_label not in current_categories:
    input_dict["biology"][keys]['haul_bin'] = input_dict["biology"][keys]['haul_bin'].cat.add_categories([out_of_range_label])


keys = "length_df"
values = input_dict["biology"][keys]

input_dict["biology"]["length_df"]
input_dict["biology"][keys]['haul_bin'] = input_dict["biology"][keys]['haul_bin'].cat.add_categories([0])

input_dict["biology"][keys]['haul_bin'].combine_first(
    pd.cut(
        input_dict["biology"][keys]['haul_weight'],
        bins=input_dict["biology"][keys]['haul_bin'].cat.categories
    )
)
input_dict["biology"][keys] = input_dict["biology"][keys].reset_index().filter(["haul_num", "haul_weight", "species_id", "region"])
input_dict["biology"][keys]["haul_bin"].fillna(0)
strata_df["stratum_num"].reindex(input_dict["biology"][keys].index)

inpfc_df["stratum_inpfc"].reindex(input_dict["biology"][keys].index)

self = survey
configuration_dict = self.config
index_variable: Union[str, List[str]] = ["transect_num", "interval"]
ingest_exports: Optional[Literal["echoview", "echopype"]] = "echoview"
read_transect_region_file: bool = False
region_class_column: str = "region_class"
transect_pattern: str = r"T(\d+)"
unique_region_id: str = "region_id"
verbose: bool = True
write_transect_region_file: bool = False

transect_data["region_name"].unique()

key = "all_ages"
values = region_names[key]

part_name = "COUNTRY"
patterns = compiled_patterns[part_name]
pattern = patterns[0]
transect_regions.loc[lambda x: x.region_id == 63]
match = pattern.search(remaining_name)
grouped_region["country"][0]
transect_region_key.loc[lambda x: ~x.country.isin(["US", "CAN"])]

row = unique_regions.loc[94, :]

unique_regions_coded["country"].isna()

unique_regions.loc[94]
unique_regions_coded.country[0]
transect_region_key.loc[lambda x: not x.country]
# eln.ingest_echoview_exports(
#     self.config,
#     transect_pattern,
#     index_variable,
#     unique_region_id,
#     region_class_column,
#     verbose,
# )
configuration_dict = self.config
default_transect_spacing = 10.0
Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2023/Biological/CAN/2023061_DFO_biodata_length.xlsx").exists()
os.listdir()
# TRANSECT ANALYSIS ================================================================================
self = survey
species_id: Union[float, list[float]] = 22500
exclude_age1: bool = True
stratum: Literal["inpfc", "ks"] = "inpfc"
verbose: bool = True

# biomass_summary, self.analysis["transect"] = acoustics_to_biology(
#     self.input, self.analysis["transect"], self.config, self.analysis["settings"]
# )
input_dict = self.input
analysis_dict = self.analysis["transect"]
configuration_dict = self.config
settings_dict = self.analysis["settings"]

nasc_biology_df = nasc_to_biology
fitted_weight_dict = analysis_dict["biology"]["weight"]
population_dict = analysis_dict["biology"]["population"]
strata_adult_proportions_df = strata_adult_proportions

distributions_dict = input_dict["biology"]["distributions"]
proportions_dict = analysis_dict["biology"]["proportions"]
TS_L_parameters = configuration_dict["TS_length_regression_parameters"]["pacific_hake"]

nasc_interval_df["fraction_hake"].max()
nasc_interval_df.loc[lambda x: x.fraction_hake > 0]
nasc_interval_df.merge(sigma_bs_strata, on=[stratum_col], how="outer")

from itertools import product 

list(product(unique_strata, sigma_bs_strata["species_id"].unique()))

nasc_biology["stratum_num"].unique()
nasc_interval_df["stratum_num"].unique()
nasc_interval_df.merge(sigma_bs_strata_bfill.reset_index(), on=[stratum_col])

file = "C:/Users/Brandyn/Documents/GitHub/EchoPro_data/Data/reports/EchoPro_un-kriged_output_0.xlsx"
filepath = Path(file)

import os

if filepath.exists():
    os.remove(filepath)
transect_region_key.dropna
proportions_dict = analysis_dict["biology"]["proportions"]["number"]
length_weight_dict = analysis_dict["biology"]["weight"]
stratum_col = settings_dict["transect"]["stratum_name"]
# analysis_dict["biology"]["proportions"].update(
#     {
#         "weight": weight_proportions(
#             catch_data,

#             analysis_dict["biology"]["proportions"]["number"],
#             length_weight_df,
#             analysis_dict["biology"]["distributions"]["weight"],
#             settings_dict["transect"]["stratum_name"],
#         )
#     }
# )
proportions_dict = analysis_dict["biology"]["proportions"]["number"]
distributions_dict = analysis_dict["biology"]["distributions"]["weight"]
stratum_column = settings_dict["transect"]["stratum_name"]

# age1_proportions = age1_metric_proportions(
#     input_dict["biology"]["distributions"],
#     analysis_dict["biology"]["proportions"],
#     configuration_dict["TS_length_regression_parameters"]["pacific_hake"],
#     settings_dict,
# )
distributions_dict = input_dict["biology"]["distributions"]
proportions_dict = analysis_dict["biology"]["proportions"]
TS_L_parameters = configuration_dict["TS_length_regression_parameters"]["pacific_hake"]

age_weight_proportions_table.loc[:, 1]
age_proportions_table.loc[:, 1]
aged_weights_binned_pvt.sum(axis=0)

from echopop.extensions.feat_report import *

self = FEATReports(survey, reports=["kriged_length_age_abundance"])
v = report_methods["kriged_length_age_abundance"]
title = v["title"]
filename = v["filename"]
tables_dict = tables
table_title = title
filepath
# STRATIFIED ANALYSIUS =============================================================================
self = survey
dataset: Literal["transect", "kriging"] = "kriging"
stratum: Literal["inpfc", "ks"] = "inpfc"
variable: Literal["abundance", "biomass", "nasc"] = "biomass"
mesh_transects_per_latitude: Optional[int] = None
transect_sample: Optional[float] = None
transect_replicates: Optional[int] = None
bootstrap_ci: float = 0.95
bootstrap_ci_method: Literal[
    "BC", "BCa", "empirical", "percentile", "standard", "t-jackknife", "t-standard"
] = "t-jackknife"
bootstrap_ci_method_alt: Optional[
    Literal["empirical", "percentile", "standard", "t-jackknife", "t-standard"]
] = "t-standard"
bootstrap_adjust_bias: bool = True
verbose=True

import warnings
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
import scipy.stats as st

from echopop.spatial.transect import transect_array
# stratified_results, self.analysis = stratified_summary(
#     self.analysis,
#     self.results,
#     self.input["spatial"],
#     self.analysis["settings"]["stratified"],
# )
analysis_dict = self.analysis
results_dict = self.results
spatial_dict = self.input["spatial"]
settings_dict = self.analysis["settings"]["stratified"]

unique_sexes = sex_stratum_proportions["sex"].unique()
missing_rows = idx_template.reset_index()[~idx_template.reset_index()["stratum_num"].isin(sex_stratum_proportions["stratum_num"])]

pd.DataFrame(
    [(row["stratum_num"], row["species_id"], sex) for _, row in missing_rows.iterrows() for sex in unique_sexes],
    columns=["stratum_num", "species_id", "sex"]
)
pd.concat([length_weight_strata.set_index([stratum_col, "sex"]).reindex(idx_template_lw_sex.index).fillna(0.0).reset_index(), length_weight_strata], ignore_index=True)
            
# VARIOGRAM ANALYSIS ===============================================================================
from echopop.spatial.variogram import prepare_variogram_matrices, quantize_lags, semivariance

self = survey
variogram_parameters: Dict[str, Any] = {}
optimization_parameters: Dict[str, Any] = {}
model: Union[str, List[str]] = ["bessel", "exponential"]
n_lags: int = 30
azimuth_range: float = 360.0
standardize_coordinates: bool = True
force_lag_zero: bool = True
initialize_variogram: Union[List[str], Dict[str, Any]] = [
    "nugget",
    "sill",
    "correlation_range",
    "hole_effect_range",
    "decay_power",
]
variable: Literal["biomass"] = "biomass"
verbose: bool = True

self.analysis["transect"]["acoustics"]["adult_transect_df"]
self.results["transect"]

lag_resolution = variogram_parameters["lag_resolution"]
angles=True
coordinates_1 = transect_data.copy()
coordinates_2 = transect_data.copy()

transect_dict = self.analysis["transect"]
settings_dict = self.analysis["settings"]["variogram"]
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

variogram_parameters = {**valid_variogram_params, **empirical_variogram_params}

# KRIGING ANALYSIS =================================================================================
cropping_parameters: Dict[str, Any] = {}
kriging_parameters: Dict[str, Any] = {}
coordinate_transform: bool = True
extrapolate: bool = False
best_fit_variogram: bool = False
variable: Literal["biomass"] = "biomass"
variogram_parameters: Optional[Dict[str, Any]] = None
verbose: bool = True

a = transect_info.reset_index()
x = a["longitude"][0]
y = a["latitude"][0]

b = transect_data.reset_index()
b.loc[lambda m: m["longitude"] == x]

b.drop_duplicates()
# aged_apportioned, unaged_apportioned, kriged_apportioned_table = apportion_kriged_values(
#     self.analysis, kriged_results["mesh_results_df"], self.analysis["settings"]["kriging"]
# )

analysis_dict = self.analysis
kriged_mesh = kriged_results["mesh_results_df"]
settings_dict = self.analysis["settings"]["kriging"]

# Get stratum column name
stratum_col = settings_dict["transect"]["stratum_name"]

# Get the aged and unaged proportions
# ---- Aged proportions
age_proportions = proportions_dict["number"]["aged_length_proportions_df"].copy()
# ---- Unaged proportions
unage_proportions = proportions_dict["number"]["unaged_length_proportions_df"].copy()
# ---- Aged weight proportions
age_weight_proportions = proportions_dict["weight"]["aged_weight_proportions_df"].copy()

# Match table shapes
# ---- All aged strata
all_aged_strata = age_proportions["stratum_num"].unique().tolist()
# ---- All unaged strata
all_unaged_strata = unage_proportions["stratum_num"].unique().tolist()
# # ---- Back-fill missing strata
all_strata = list(set(all_aged_strata + all_unaged_strata))

# Consider the proportions only for `sex="all"` 
# ---- Aged
age_proportions_all = age_proportions[age_proportions["sex"] == "all"].reset_index(drop=True)   
# ---- Unaged
unage_proportions_all = unage_proportions[unage_proportions["sex"] == "all"]

# Calculate the new length-averaged sigma_bs for each stratum
# ---- Extract the length-binned values
length_bins = distributions_dict["length_bins_df"].copy()
# ---- Square the length values
length_bins["length_sq"] = length_bins["length_bins"] ** 2.0
# ---- Multiply by the TS-length regression coefficient (in the linear domain)
length_bins["length_sq"] = length_bins["length_sq"] * 10 ** (
    TS_L_parameters["TS_L_intercept"] / 10.0
)
# ---- Create temporary pivot table for the age proportions to compute the dot product
age_proportions_dot_table = age_proportions_all.pivot_table(
    index=["length_bin"],
    columns=[stratum_col],
    values="proportion_number_aged",
    aggfunc="sum",
    observed=False,
)
# ---- Dot product to calculate the new average sigma_bs for all ages
updated_sigma_bs = length_bins["length_sq"].values.dot(age_proportions_dot_table)

# Compute the age-1 proportions
# ---- Index the age-1 proportions
age1_proportions_all = age_proportions_all.set_index(["age_bin"]).loc[1]
# ---- Pivot into a table for the dot product
age1_proportions_table = age1_proportions_all.pivot_table(
    index=["length_bin"],
    columns=[stratum_col],
    values="proportion_number_aged",
    aggfunc="sum",
    observed=False,
)
# ---- Dot product to calculate the average sigma_bs for age-1 fish
age1_sigma_bs = length_bins["length_sq"].values.dot(age1_proportions_table)
# ---- Calculate age-1 NASC proportioon per stratum
age1_nasc_proportions = age1_sigma_bs / updated_sigma_bs
# ---- Sum the number proportions for each age bin within each stratum
age1_strata_proportions = age1_proportions_table.sum(axis=0).to_numpy()

# Convert the primary tables into pivot tables
# ---- Aged
age_proportions_table = age_proportions_all.pivot_table(
    index=["length_bin"],
    columns=[stratum_col, "age_bin"],
    values="proportion_number_aged",
    aggfunc="sum",
    observed=False,
)
# ---- Unaged
unage_proportions_table = unage_proportions_all.pivot_table(
    columns=[stratum_col],
    index=["length_bin"],
    values="proportion_number_unaged",
    aggfunc="sum",
    observed=False,
)

# Compute the length index thresholds
min_index = np.where(length_bins["length_bins"] == 10.0)[0]
if len(min_index) == 0:
    min_index = 0
else:
    min_index = min_index[0]
max_index = length_bins["length_bins"].size

# Calculate thresholds derived from the summed length distributions of age-1 fish
# ---- General length distribution threshold
age1_length_distribution_threshold = (
    unage_proportions_table.iloc[np.arange(min_index, max_index), :]
    * age1_proportions_table.iloc[np.arange(min_index, max_index), :]
).sum()
# ---- Just aged length distribution (age-1)
age1_specific_length_distribution_threshold = age1_proportions_table.sum(axis=0)

# Pivot the aged weight proportions table
age_weight_proportions_table = age_weight_proportions.pivot_table(
    index=["length_bin"],
    columns=["age_bin", stratum_col],
    values="weight_proportions",
    aggfunc="sum",
    observed=False,
)
# ---- Repivot a separate table for calculating the aged weight proportions for age-1 fish
summed_age_weight_proportions = age_weight_proportions.pivot_table(
    index=["age_bin"],
    columns=[stratum_col],
    values="weight_proportions",
    aggfunc="sum",
    observed=False,
)

# Calculate the age-1 weight proportions
age1_weight_proportions = np.where(
    (age1_length_distribution_threshold <= 1e-10)
    & (age1_specific_length_distribution_threshold <= 1e-10),
    0.0,
    age_weight_proportions_table[1].sum() / summed_age_weight_proportions.sum(),
)

# Return output
# ---- Create DataFrame
apportioned_age1 = pd.DataFrame({f"{stratum_col}": np.unique(age_proportions[stratum_col])})
# ---- Number proportions
apportioned_age1["number_proportion"] = age1_strata_proportions
# ---- Weight proportions
apportioned_age1["weight_proportion"] = age1_weight_proportions
# ---- NASC
apportioned_age1["nasc_proportion"] = age1_nasc_proportions