import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
import pandas as pd
import yaml
import re

from echopop.core import BIODATA_HAUL_MAP, DATA_STRUCTURE, LAYER_NAME_MAP, NAME_CONFIG
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
from echopop.core import DATA_STRUCTURE
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

from echopop.core import ECHOVIEW_EXPORT_MAP, REGION_EXPORT_MAP
from echopop.spatial.transect import export_transect_layers, export_transect_spacing
from echopop.utils.validate_df import KSStrata
from echopop.utils.operations import compile_patterns, extract_parts_and_labels, group_merge
from echopop.utils.load_nasc import (
    validate_echoview_exports,
    validate_export_directories,
    consolidate_exports,
    construct_transect_region_key,
    filter_export_regions,
    compile_patterns,
    get_haul_transect_key,
    get_transect_numbers,
    read_echoview_export,
    load_export_regions,
    export_transect_layers,
    extract_parts_and_labels,
)

# LOAD =============================================================================================
new_datasets = ["biological", "kriging", "stratification"]
dataset_type = new_datasets

self = survey
index_variable: Union[str, List[str]] = ["transect_num", "interval"]
ingest_exports: Optional[Literal["echoview", "echopype"]] = "echoview"
region_class_column: str = "region_class"
transect_pattern: str = r"T(\d+)"
unique_region_id: str = "region_id"
verbose: bool = True

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

filename = cells_files[63]
transect_num = 49
pd.read_csv(filename)
read_echoview_export(cells_files[63], 5)
validate_echoview_exports(read_echoview_export(cells_files[64], 5))

cells_files[64]

pd.concat(generate_dataframes(cells_files[61:64], transect_reference), axis=0, ignore_index=True)

pd.concat(
    [
        validate_echoview_exports(read_echoview_export(interval_files[85], transect_reference.get(interval_files[85], None))),
        validate_echoview_exports(read_echoview_export(interval_files[86], transect_reference.get(interval_files[86], None))),
        validate_echoview_exports(read_echoview_export(interval_files[87], transect_reference.get(interval_files[87], None))),   
    ],
    axis=0,
    ignore_index=True,
)

export_dataframe = read_echoview_export(interval_files[87], transect_reference.get(interval_files[87], None))

pd.concat(generate_dataframes(interval_files[85:89], transect_reference), axis=0, ignore_index=True)

len(interval_files)
len(np.unique(interval_files)

dataset = "biological"
datalayer = "specimen"
region_id = "US"
input_dict = self.input
configuration_dict = self.config
# sheet_name = sheets
sheet_name = sheet_name[0]

column = "haul_start"
dtype = valid_cols[column]

sheet_name = sheet_name[0]
# cls = SpecimenBiodata
cls = copy.deepcopy(validation_settings)
cls = SpecimenBiodata
# cls = VarioKrigingPara
dff = df.copy()
# dff = df_initial.copy()
data = dff.copy()

column_name = 'weight'
col = 'weight'



[isinstance(c, (np.integer, np.floating)) if not isinstance(c, list) else c in valid_cols[col] for c in valid_cols[col]]


c in valid_cols[col] or 
[c for c in valid_cols[col] if c in [int, float] else ]

valid_cols[col].type

df.loc[df[col].astype(str).str.contains(r"\s"), col]

da = pd.DataFrame({
    "x": [1, 2, 3, "", " ", "  ", "   ", None, 4, 5]
})
da.loc[da["x"].str.isnumeric()]

da.loc[da.x.astype(str).str.contains(r"\s"), :]

len(list(cls.to_schema().columns))
cls.to_schema().columns["krig.dx"]
len(list(default_annotations))

column = "fraction_hake"
annotation = valid_cols[column]

column_name = "length"
col_types = column_types

annotation_test = [re.match(annotation, col).group(0) for annotation in col_types 
                    if re.match(annotation, col)][0]

col_types

for col in valid_cols:
    # ---- Match the underlying annotation
    annotation_test = [re.match(annotation, col).group(0) for annotation in col_types 
                       if re.match(annotation, col)][0]
    # ---- Get the `dtype`
    dtype = col_types[annotation_test]
    # ---- Test whether the Series is coercible
    if isinstance(dtype, list):
        # ---- Iterate through the list
        for typing in dtype:
            # ---- Initialize the `dtype` of the validator annotation
            cls.__annotations__[annotation_test] = Series[typing]
            # ---- Get the `test`
            test = cls._DTYPE_TESTS.get(typing, None)
            # ---- Apply test
            if test and test(df[col]):
                # ---- Adjust typing annotation
                cls.__annotations__[col] = Series[typing]
            # ---- Coerce the datatype (or attempt to)
            try:
                df[col] = cls._DTYPE_COERCION.get(typing)(df[col])
                # ---- If successful, break
                break
            except Exception as e:
                # ---- Drop traceback
                e.__traceback__ = None
                # ---- Format message
                message = (
                    f"{col.capitalize()} column must be a Series of '{str(dtype)}' "
                    f"values. Series values could not be automatically coerced."
                )
                # ---- Add error to collector
                errors_coerce = pd.concat(
                    [errors_coerce, pd.DataFrame(dict(Column=col, error=message))]
                )      
    # ---- If not a List supplied by the metadata attribute
    else:
        # ---- Attempt coercion
        try:
            df[col] = df[col].astype(str(dtype))
        except Exception as e:
            # ---- Drop traceback
            e.__traceback__ = None            
            # ---- Format message
            message = (
                f"{col.capitalize()} column must be a Series of '{str(dtype)}' "
                f"values. Series values could not be automatically coerced."
            )
            # ---- Add error to collector
            errors_coerce = pd.concat(
                [errors_coerce, pd.DataFrame(dict(Column=col, error=message))]
            )      
    

# TRANSECT ANALYSIS ================================================================================
self = survey
species_id: Union[float, list[float]] = 22500
exclude_age1: bool = True
stratum: Literal["inpfc", "ks"] = "ks"
verbose: bool = True

# biomass_summary, self.analysis["transect"] = acoustics_to_biology(
#     self.input, self.analysis["transect"], self.config, self.analysis["settings"]
# )
input_dict = self.input
analysis_dict = self.analysis["transect"]
configuration_dict = self.config
settings_dict = self.analysis["settings"]

self.analysis["transect"]["acoustics"]["adult_transect_df"]

nasc_interval_df["stratum_num"].unique()

proportions_dict = analysis_dict["biology"]["proportions"]["number"]
length_weight_dict = analysis_dict["biology"]["weight"]
np.concatenate([weight_all, weight_male, weight_female]).size
pd.DataFrame(
    {
        f"{stratum_col}": np.tile(
                np.unique(station_proportions[stratum_col]), len(np.unique(station_proportions.sex))
            ),
        "sex": np.repeat(
            ["all", "male", "female"], len(np.unique(station_proportions[stratum_col]))
        ),
    }
)

# strata_adult_proportions, nasc_to_biology = nasc_to_biomass(
#     input_dict, analysis_dict, configuration_dict, settings_dict
# )
self.input["biology"]["specimen_df"]["stratum_num"].unique()
self.input["biology"]["length_df"]["stratum_num"].unique()
self.input["biology"]["catch_df"]["stratum_num"].unique()
self.input["acoustics"]["nasc_df"]["NASC_all_ages"].unique()
self.input["spatial"]
# age1_proportions = age1_metric_proportions(
#     input_dict["biology"]["distributions"],
#     analysis_dict["biology"]["proportions"],
#     configuration_dict["TS_length_regression_parameters"]["pacific_hake"],
#     settings_dict,
# )
distributions_dict = input_dict["biology"]["distributions"]
proportions_dict = analysis_dict["biology"]["proportions"]
TS_L_parameters = configuration_dict["TS_length_regression_parameters"]["pacific_hake"]

# STRATIFIED ANALYSIUS =============================================================================
self = survey
dataset: Literal["transect", "kriging"] = "transect"
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

data_series = transect_distances
index_matrix = transect_numbers_arr
[transect_distances[idx] for idx in transect_numbers_arr]
transect_distances.index.unique()
transect_summary[transect_summary.duplicated()]
transect_distances.loc[1029]

transect_numbers_arr.shape
transect_summary[transect_summary["transect_num"] == 1029]
transect_distances.loc[transect_distances.index.duplicated()]

transect_numbers_arr = [
    np.random.choice(
        transect_numbers[j].values, num_transects_to_sample[j], replace=False
    )
    for i in range(transect_replicates)
]
print([x.shape for x in transect_numbers_arr])  

[transect_distances[idx] for idx in transect_numbers_arr]

self.input["acoustics"]["nasc_df"]

self.analysis["transect"]["acoustics"]["adult_transect_df"]

print([x.shape for x in transect_numbers_arr])


# KRIGING ANALYSIS =================================================================================
cropping_parameters: Dict[str, Any] = {}
kriging_parameters: Dict[str, Any] = {}
coordinate_transform: bool = True
extrapolate: bool = False
best_fit_variogram: bool = False
variable: Literal["biomass"] = "biomass"
variogram_parameters: Optional[Dict[str, Any]] = None
verbose: bool = True

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