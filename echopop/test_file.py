from echopop import Survey
import copy
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union
import numpy as np
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
from echopop.graphics import variogram_interactive as egv
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import edit_transect_columns
from echopop.utils import load as el, load_nasc as eln, message as em
from echopop.utils.load import dataset_integrity
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
from echopop.spatial.projection import wgs84_to_utm
from echopop.spatial.transect import transect_bearing, transect_extent
import geopandas as gpd
import geopy.distance
import numpy as np
import pandas as pd
from scipy import interpolate

init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
survey = Survey(init_config, file_config)
survey.load_survey_data(verbose=False)
survey.load_acoustic_data(verbose=False)
survey.transect_analysis(verbose=False)
survey.fit_variogram(verbose=False)
survey.stratified_analysis(verbose=False)
survey.kriging_analysis(verbose=False)

self = survey
cropping_parameters: Dict[str, Any] = {}
kriging_parameters: Dict[str, Any] = {}
coordinate_transform: bool = True
extrapolate: bool = False
best_fit_variogram: bool = False
variable: Literal["biomass"] = "biomass"
variogram_parameters: Optional[Dict[str, Any]] = None
veerbose: bool = True

input_dict = self.input
analysis_dict = self.analysis
settings_dict = self.analysis["settings"]["kriging"]

transect_data = transect_data
mesh_data = mesh_data
cropping_parameters = validated_cropping_methods

transect_data = transect_data.copy()
mesh_data = mesh
cropping_parameters = cropping_parameters