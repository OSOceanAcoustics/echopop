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
from echopop.nwfsc_feat import spatial
spatial.get_survey_western_extents = get_survey_western_extents
spatial.adaptive_search_radius = adaptive_search_radius
adaptive_search_radius = spatial.adaptive_search_radius
spatial.western_boundary_search_strategy = western_boundary_search_strategy
spatial.ordinary_kriging_matrix = ordinary_kriging_matrix
spatial.kriging_lambda = kriging_lambda
spatial.parse_stacked_kriging_array = parse_stacked_kriging_array
spatial.kriging_point_estimator = kriging_point_estimator
spatial.krige = krige
spatial.project_kriging_results = project_kriging_results
adaptive_search_strategy = spatial.western_boundary_search_strategy
##########
# Get the survey transects' western extents
transect_western_extents = spatial.get_survey_western_extents(
    transect_df=df_nasc_all_ages,
    coordinate_names=("x", "y"),
    latitude_threshold=51.
)

# Pre-define arguments within a partial function defining the western boundary search strategy
western_strategy = partial(spatial.western_boundary_search_strategy, 
                           western_extent=transect_western_extents,
                           kriging_mesh=kriging_mesh,
                           coordinate_names=("x", "y"))
adaptive_search_strategy = western_strategy
variogram_parameters = {"model": model, **dict_best_fit_variogram_params}
kriging_parameters = {
    "search_radius": search_radius,
    "anisotropy": anisotropy,
    "k_min": k_min,
    "k_max": k_max,
}   


# Apply kriging to interpolate transect values over the mesh grid
kriged_estimates = spatial.krige(
    transect_df=df_nasc_all_ages,
    kriging_mesh=df_mesh,
    coordinate_names=("x", "y"),
    variable="biomass_density",
    kriging_parameters=kriging_parameters,
    variogram_parameters=variogram_parameters,
    adaptive_search_strategy=western_strategy,
)    

# Post-process the kriged estimates

variable = "biomass_density"

###

# Compute the global/survey variance
survey_variance = transect_df[variable].var()
default_mesh_cell_area = 6.25 # nmi

def project_kriging_results(
    kriged_estimates: np.ndarray,
    kriging_mesh: pd.DataFrame,
    transect_df: pd.DataFrame,
    default_mesh_cell_area: float,
    variable: str,
) -> Tuple[pd.DataFrame, float]:
    """
    Compute mesh CV and survey CV from kriged values (assumed to be in absolute units).

    Parameters
    ----------
    kriged_values : np.ndarray, shape (n, 3)
        Kriged values with columns: [estimate, kriging variance, sample variance]
    survey_variance : float
        Variance of the original data (assumed absolute units).

    Returns
    -------
    Tuple[pd.DataFrame, float]
        A tuple containing the output mesh DataFrame with columns including the kriged estimates 
        and variance, sample variance, and cell coefficient of variation (CV). The other value is 
        the overall CV computed for the entire kriging mesh.
    """

    # Create copy
    mesh_results = kriging_mesh.copy()

    # Get the adjusted areas
    if "fraction" in kriging_mesh.columns and default_mesh_cell_area is not None:
        mesh_results["area"] = mesh_results["fraction"] * default_mesh_cell_area

    # Calculate the (global) survey variance
    survey_variance = transect_df[variable].var()

    # Distribute estimates over an area
    survey_estimate = np.nansum(kriged_estimates[:, 0] * mesh_results["area"])

    # Add the kriged array values
    # ---- Kriged estimates
    mesh_results["estimate"] = kriged_estimates[:, 0]
    # ---- Kriged variance
    mesh_results["kriged_variance"] = kriged_estimates[:, 1]
    # ---- Sample variance
    mesh_results["sample_variance"] = kriged_estimates[:, 2]
    # ---- Mesh node CV
    mesh_results["cell_cv"] = (
        mesh_results["area"].mean()
        * np.sqrt(kriged_estimates[:, 1].dot(survey_variance))
        / survey_estimate
        * np.sqrt(len(kriged_estimates[:, 1]))
    )

    # Calculate the global CV
    global_cv = (
        np.sqrt(np.nansum(kriged_estimates[:, 1] * mesh_results["area"]**2) * survey_variance) / 
        survey_estimate
    )

    # Return the DataFrame and global CV estimate
    return mesh_results, global_cv
    



def integrate_grid_cells(
    kriged_estimates: np.ndarray[float],
    kriging_mesh: pd.DataFrame,
    default_mesh_cell_area: Optional[float] = None,
):

    #
    if "area" in kriging_mesh.columns:
        area = kriging_mesh["area"]
    elif "fraction" in kriging_mesh.columns and default_mesh_cell_area is not None:
        area = kriging_mesh["fraction"] * default_mesh_cell_area

    #



# Compute the coefficients of variation (CV)
# ---- Compute the global/survey variance
survey_variance = transect_df[variable].var()
# -------- Drop erroneous negative values along edge
KR = kriged_estimates.copy()
KR[KR[:, 0] < 0, 0] = 0.0
# -------- Distribute biological variable over area
survey_estimate = np.nansum(KR[:, 0] * area)
# ---- Compute the georeferenced CV at each mesh node
mesh_CV1 = (
    area.mean()
    * np.sqrt(KR[:, 1].dot(survey_variance))
    / survey_estimate
    * np.sqrt(len(KR[:, 1]))
)
# ---- Compute the global/survey CV
survey_CV1 = (
    np.sqrt(np.nansum(KR[:, 1] * area**2) * survey_variance) / survey_estimate
)

mesh_CV - mesh_CV1
survey_CV - survey_CV1

KRA = kriged_estimates.copy()
KRA[KRA[:, 0] < 0, 0] = 0.0
KRA[:, 0] *= area 
KRA[:, 1] = (np.sqrt(KRA[:, 1]) * area) ** 2 # scale estimate by area
# KRA[:, 1] *= area**2           # scale variance by area^2

SU_EST = np.nansum(KRA[:, 0])
SU_VAR = (transect_df[variable] * transect_df["area_interval"]).var()

mesh_CV2 = np.sqrt(KRA[:, 1].dot(SU_VAR)) / SU_EST * np.sqrt(len(KRA[:, 1]))
survey_CV2 = np.sqrt(np.nansum(KRA[:, 1]) * SU_VAR) / SU_EST


survey_cv = np.sqrt(np.nansum(variances1 * estimates0**2) * SU_VAR) / survey_estimate1
# -------- Drop erroneous negative values along edge
kriged_estimates[kriged_estimates[:, 0] < 0, 0] = 0.0
# -------- Distribute biological variable over area
SURVEY_estimate1 = np.nansum(kriged_estimates[:, 0] * area)
SU_VAR = (kriged_estimates[:, 0] * area).var()

kriged_values_abs = kriged_estimates.copy()
kriged_values_abs[:, 0] *= area              # scale estimate by area
kriged_values_abs[:, 1] *= area**2           # scale variance by area^2

# ---- Compute the georeferenced CV at each mesh node
mesh_CV1 = (
    np.sqrt(kriged_values_abs[:, 1].dot(SU_VAR))
    / SU_EST
    * np.sqrt(len(kriged_values_abs[:, 1]))
)
# ---- Compute the global/survey CV
survey_CV1= (
    np.sqrt(np.nansum(kriged_estimates[:, 1] * area**2) * survey_variance1) / survey_estimate1
)


kriged_values_abs = kriged_estimates.copy()
kriged_values_abs[:, 0] *= area              # scale estimate by area
kriged_values_abs[:, 1] *= area**2           # scale variance by area^2

estimates0 = kriged_values_abs[:, 0]
variances1 = kriged_values_abs[:, 1]
n = len(estimates0)

survey_estimate1 = np.nansum(estimates0)

mesh_cv = area.mean() * np.sqrt(variances1.dot(survey_variance)) / survey_estimate1 * np.sqrt(n)
survey_cv = np.sqrt(np.nansum(variances1 * estimates0**2) * survey_variance) / survey_estimate1

area.mean()

# -------- Drop erroneous negative values along edge
kriged_estimates[kriged_estimates[:, 0] < 0, 0] = 0.0
# -------- Distribute biological variable over area
SURVEY_estimate1 = np.nansum(kriged_estimates[:, 0] * area)
# ---- Compute the georeferenced CV at each mesh node
mesh_CV1 = (
    area.mean()
    * np.sqrt(kriged_estimates[:, 1].dot(survey_variance))
    / SURVEY_estimate1
    * np.sqrt(len(kriged_estimates[:, 1]))
)
# ---- Compute the global/survey CV
survey_CV1= (
    np.sqrt(np.nansum(kriged_estimates[:, 1] * area**2) * survey_variance1) / survey_estimate1
)




kriged_estimates[:, 2].var()


variogram(local_distance_matrix, {"model": model, **variogram_parameters})
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