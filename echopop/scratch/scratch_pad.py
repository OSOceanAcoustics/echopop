import abc
import numpy as np
import pandas as pd
from typing import Union, Dict, List, Optional, Any
from functools import reduce
import pytest

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs
from echopop.nwfsc_feat import utils
from typing import Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate

#######

########
from functools import partial
from echopop.spatial.variogram import variogram
variable: str = "biomass_density"
anisotropy: float = dict_kriging_params["anisotropy"]
search_radius: float = dict_kriging_params["search_radius"]
# OR: 
correlation_range: float = dict_best_fit_variogram_params["correlation_range"]
k_min: int = int(dict_kriging_params["kmin"])
k_max: int = int(dict_kriging_params["kmax"])
kriging_mesh: pd.DataFrame = df_mesh.copy()
transect_df: pd.DataFrame = df_nasc_all_ages.copy()
model = ["bessel", "exponential"]
coordinates = ("x", "y")
coordinate_names = coordinates
variogram_parameters = dict_best_fit_variogram_params
# KrigingParameterInputs.create(**{"correlation_range": correlation_range})
# KrigingAnalysis.create(**{})

##########
transect_df = df_nasc_all_ages.copy()
mesh_df = df_mesh.copy()
##########
tx = transect_df["transect_num"]
uniq_tx = tx.unique()
nt = len(uniq_tx)

from typing import List, Tuple
###################

###################
num_nearest_transects: float = 3
mesh_buffer_distance: float = 2.5
projection: str = "epsg:4326" 
mesh_df = df_mesh.copy()
transect_df = df_nasc_all_ages.copy()
from echopop.nwfsc_feat.projection import wgs84_to_utm
from echopop.nwfsc_feat.spatial import transect_extent
spatial.transect_extent = transect_extent
###################

# Convert mesh DataFrame into a GeoDataframe
mesh_gdf = gpd.GeoDataFrame(
    mesh_df,
    geometry=gpd.points_from_xy(mesh_df["longitude"], mesh_df["latitude"]),
    crs=projection,
)

# Convert the mesh projection to UTM (m)
wgs84_to_utm(mesh_gdf)

# Determine the survey extent by generating the border polygon
survey_polygon = spatial.transect_extent(
    transect_df, projection, num_nearest_transects
)

# Find the mesh coordinates that fall within the buffered polygon
# ---- Convert `grid_buffer` (nmi) to m and add buffer to polygon
survey_polygon_buffered = survey_polygon.buffer(mesh_buffer_distance * 1852)
# ---- Inclusion/union filter mask
within_polygon_mask = mesh_gdf.geometry.within(survey_polygon_buffered)
# ---- Apply mask to the mesh grid
mesh_gdf_masked = mesh_gdf[within_polygon_mask]

# Mask the mesh
mesh_regions = mesh_gdf_masked.drop(columns="geometry")

# Extract the transect regions
# transect_mesh_regions = transect_df.reset_index().loc[interpolated_indices]
# mesh_regions = mesh_df.loc[mesh_indices]
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 6)) # You can set the figure size here
mesh_df.plot.scatter(x="longitude", y="latitude", s=0.2, color="gray", ax=ax, zorder=1)
mesh_regions.plot.scatter(x="longitude", y="latitude", ax=ax, zorder=2)
transect_df.plot.scatter(x="longitude", y="latitude", marker="s", s=1, color="red", ax=ax, zorder=3)
transect_df.loc[transect_df.transect_num == 145].plot.scatter(x="longitude", y="latitude", marker="s", s=1, color="gold", ax=ax, zorder=4)
plt.show()


####################################################################################################
from echopop.survey import Survey
from echopop.utils.validate_dict import KrigingParameterInputs, KrigingAnalysis
from echopop.spatial.transect import correct_transect_intervals
from echopop.acoustics import aggregate_sigma_bs, nasc_to_biomass
# from echopop.biology import (
#     # age1_metric_proportions,
#     # distribute_length_age,
#     # filter_species,
#     # fit_length_weight_relationship,
#     # fit_length_weights,
#     # impute_kriged_values,
#     # # number_proportions,
#     # # partition_transect_age,
#     # quantize_number_counts,
#     # quantize_weights,
#     # reallocate_kriged_age1,
#     # weight_proportions,
# )
from echopop.spatial.krige import kriging
from echopop.spatial.mesh import crop_mesh, mesh_to_transects, stratify_mesh
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import (
    edit_transect_columns,
    save_transect_coordinates,
    summarize_transect_strata,
    transect_spatial_features,
)
from echopop.analysis import (
    acoustics_to_biology,
    apportion_kriged_values,
    krige,
    process_transect_data,
    stratified_summary,
    variogram_analysis,
)
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import edit_transect_columns
from echopop.utils import load as el, load_nasc as eln, message as em
from echopop.utils.load import dataset_integrity
from echopop.spatial.variogram import (
    empirical_variogram,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
)
from echopop.spatial.mesh import griddify_lag_distances
from echopop.spatial.transect import define_western_extent
from echopop.spatial.krige import kriging
from echopop.statistics import stratified_transect_statistic
from echopop.utils.validate_dict import (
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
)
from echopop.spatial.krige import griddify_lag_distances, define_western_extent, adaptive_search_radius, count_within_radius, kriging_interpolation, kriging_lambda, kriging_matrix
from echopop.spatial.krige import kriging_lambda, kriging_matrix
from echopop.nwfsc_feat.spatial import kriging_lambda
survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data(ingest_exports="echoview")
survey.load_survey_data()
survey.transect_analysis(exclude_age1=False)
survey.fit_variogram()
survey.kriging_analysis(extrapolate=True)

self = survey
input_dict, analysis_dict, settings_dict = self.input, self.analysis, self.analysis["settings"]["kriging"]
transect_data, mesh_data, settings_dict = analysis_dict["kriging"]["transect_df"], analysis_dict["kriging"]["mesh_df"], settings_dict

survey.analysis["settings"]["variogram"][""]
self = survey
input_dict, analysis_dict, configuration_dict, settings_dict = self.input, self.analysis["transect"], self.config, self.analysis["settings"]
# Extract the necessary correct strata mean sigma_bs
sigma_bs_strata = analysis_dict["acoustics"]["sigma_bs"]["strata_mean_df"]

# Pull out the length-weight conversion for each stratum
length_weight_strata = analysis_dict["biology"]["weight"]["weight_stratum_df"]

# Get the name of the stratum column
stratum_col = settings_dict["transect"]["stratum_name"]

# Get group-specific columns
age_group_cols = settings_dict["transect"]["age_group_columns"]

# Extract the correct strata dataframe
# ---- Define `strata_df` if KS
if settings_dict["transect"]["stratum"] == "ks":
    strata_df = input_dict["spatial"]["strata_df"].copy()
# Define `inpfc_strata_df` if INPFC
elif settings_dict["transect"]["stratum"] == "inpfc":
    strata_df = input_dict["spatial"]["inpfc_strata_df"].copy()

# Get group-specific column names and create conversion key
name_conversion_key = {age_group_cols["haul_id"]: "haul_num", age_group_cols["nasc_id"]: "nasc"}
# ---- Update if the stratum is not equal to INPFC
if settings_dict["transect"]["stratum"] != "inpfc":
    name_conversion_key.update({age_group_cols["stratum_id"]: stratum_col})

# Rename columns
# ---- Extract NASC data
nasc_data = input_dict["acoustics"]["nasc_df"].copy()
# ---- Change names
nasc_data.rename(columns=name_conversion_key, inplace=True)

# Correct the acoustic survey transect intervals
nasc_interval_df = correct_transect_intervals(nasc_data)

distributions_dict, proportions_dict, TS_L_parameters, settings_dict = (
    input_dict["biology"]["distributions"],
    analysis_dict["biology"]["proportions"],
    configuration_dict["TS_length_regression_parameters"]["pacific_hake"],
    settings_dict
)