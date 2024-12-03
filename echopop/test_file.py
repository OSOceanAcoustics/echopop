import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
import pandas as pd
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
from echopop.utils.load import dataset_integrity
import yaml
from echopop.core import BIODATA_HAUL_MAP, DATA_STRUCTURE, LAYER_NAME_MAP, NAME_CONFIG
from echopop.utils.data_structure_utils import map_imported_datasets
from echopop.utils.validate_df import DATASET_DF_MODEL
from echopop.utils.validate_dict import CONFIG_DATA_MODEL, CONFIG_INIT_MODEL
from echopop.utils.load import map_imported_datasets, read_validated_data, prepare_input_data, preprocess_acoustic_spatial, preprocess_biodata, preprocess_spatial, preprocess_statistics
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

from echopop.survey import Survey

init_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config.yml"
survey_year_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
survey = Survey(init_config_path, survey_year_config_path)
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis(stratum="inpfc")
survey.transect_analysis()



survey.input["spatial"]["strata_df"]
survey.input["acoustics"]["nasc_df"]
survey.analysis["transect"]["acoustics"]
survey.analysis["acoustics"]["adult_transect_df"]



self = survey
species_id: Union[float, list[float]] = 22500
exclude_age1: bool = True
stratum: Literal["inpfc", "ks"] = "inpfc"
verbose: bool = True

# Check dataset integrity
dataset_integrity(self.input, analysis="transect")

# Update settings to reflect the stratum definition
self.analysis["settings"].update(
    {
        "transect": {
            "age_group_columns": {
                "haul_id": "haul_no_age1" if exclude_age1 else "haul_all_ages",
                "nasc_id": "NASC_no_age1" if exclude_age1 else "NASC_all_ages",
                "stratum_id": "stratum_no_age1" if exclude_age1 else "stratum_all_ages",
            },
            "species_id": species_id,
            "stratum": stratum.lower(),
            "stratum_name": "stratum_num" if stratum == "ks" else "stratum_inpfc",
            "unique_strata": (
                np.unique(self.input["spatial"]["strata_df"]["stratum_num"])
                if stratum == "ks"
                else np.unique(self.input["spatial"]["inpfc_strata_df"]["stratum_inpfc"])
            ),
            "exclude_age1": exclude_age1,
        }
    }
)

# Initial data processing of the transect biological and acoustic data
self.analysis["transect"] = process_transect_data(
    self.input, self.analysis["transect"], self.analysis["settings"], self.config
)

# Convert NASC into number density (animals/nmi^2), biomass density (kg/nmi^2), abundance
# (# animals), and biomass (kg) for all fish, sexed (male/female) fish, and unsexed fish
# ---- This further provides the resulting distributions of biomass and abundance over
# ---- length and age for each sex across the entire survey
# biomass_summary, self.analysis["transect"] = acoustics_to_biology(
#     self.input, self.analysis["transect"], self.config, self.analysis["settings"]
# )
input_dict = self.input
analysis_dict = self.analysis["transect"]
configuration_dict = self.config
settings_dict = self.analysis["settings"]

# Convert NASC into number density (animals/nmi^2), biomass density (kg/nmi^2), abundance
# (# animals), and biomass (kg) for all fish, sexed (male/female) fish, and unsexed fish
# strata_adult_proportions, nasc_to_biology = nasc_to_biomass(
#     input_dict, analysis_dict, configuration_dict, settings_dict
# )

# Extract the necessary correct strata mean sigma_bs
sigma_bs_strata = analysis_dict["acoustics"]["sigma_bs"]["strata_mean_df"]

# Pull out the length-weight conversion for each stratum
length_weight_strata = analysis_dict["biology"]["weight"]["weight_stratum_df"]

# Get the name of the stratum column
stratum_col = settings_dict["transect"]["stratum_name"]

# Get group-specific columns
age_group_cols = settings_dict["transect"]["age_group_columns"]

# Extract the correct strata dataframe
# ---- Define `strata_df`
strata_df = input_dict["spatial"]["strata_df"].copy()
# ---- Determine with the default ('ks') needs to be swapped out for 'inpfc'
if settings_dict["transect"]["stratum"] == "inpfc":
    # ---- Get the INPFC strata
    inpfc_df = input_dict["spatial"]["inpfc_strata_df"].copy().set_index(["haul_bin"])
    # ---- Offset starts
    inpfc_df["haul_start"] = inpfc_df["haul_start"] - int(1)
    # ---- Get the `haul_bins`
    haul_bins =  np.unique(inpfc_df.loc[:, "haul_start":"haul_end"].stack().values)
    # NOT CORRECT BINS ! 
    
    # ---- Cut `strata_df` to haul_bins
    strata_df["haul_bin"] = pd.cut(strata_df["haul_num"], haul_bins)
    # ---- Set index
    strata_df.set_index(["haul_bin"], inplace=True)
    # ---- Merge
    strata_df["stratum_inpfc"] = inpfc_df["stratum_inpfc"]
    # ---- Reset index 
    strata_df = strata_df.reset_index().drop(columns=["haul_bin", "stratum_num"])