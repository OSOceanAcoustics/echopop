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
azimuth_range = 360.
n_lags = 30
lag_resolution = 0.002
sill = 0.91
hole_effect_range = 0.
correlation_range = 0.007
decay_power = 1.5
transect_df = df_nasc_all_ages.copy()
estimates = transect_df["biomass"].to_numpy()
#######

#######
spatial.filter_lag_matrix = filter_lag_matrix
#######
transect_df: pd.DataFrame = df_nasc_all_ages.copy()
variable: str = "biomass_density"
coordinates: List[str] = ["x", "y"]
azimuth_filter: bool = True
azimuth_angle_threshold: float = 180.
lag_resolution = 0.002
n_lags = 30
force_lag_zero: bool = True


##################################################
from lmfit import Minimizer, Parameters
from echopop.spatial.variogram import get_variogram_arguments, variogram

# Set up `lmfit` parameters
variogram_parameters = Parameters()
variogram_parameters.add_many(
    ("nugget", dict_variogram_params["nugget"], True, 0., None),
    ("sill", dict_variogram_params["sill"], True, 0., None),
    ("correlation_range", dict_variogram_params["correlation_range"], True, 0., None),
    ("hole_effect_range", dict_variogram_params["hole_effect_range"], True, 0., None),
    ("decay_power", dict_variogram_params["decay_power"], True, 0., None),
)

# Set up optimization parameters used for fitting the variogram
dict_optimization = {"max_nfev": 500, "ftol": 1e-06, "gtol": 0.0001, "xtol": 1e-06, 
                     "diff_step": 1e-08, "tr_solver": "exact", "x_scale": "jac", 
                     "jac": "3-point"}

model = ["exponential", "bessel"]
lags = lags
lag_counts = lag_counts
gamma = gamma

fit_variogram(lags, lag_counts, gamma, variogram_parameters, model, dict_optimization)

########

####################################################################################################
from echopop.survey import Survey
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

survey = Survey(init_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis()
survey.fit_variogram()
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