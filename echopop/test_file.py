from echopop.survey import Survey

survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml",
                survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis()
survey.fit_variogram()
##################################
self = survey
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import edit_transect_columns
from echopop.analysis import (
    acoustics_to_biology,
    apportion_kriged_values,
    krige,
    process_transect_data,
    stratified_summary,
    variogram_analysis,
)
import copy
from pathlib import Path
from typing import List, Literal, Optional, Union, TypedDict, Dict, Any
from echopop.utils.validate import (
    VariogramBase, VariogramInitial, VariogramOptimize, MeshCrop, KrigingParameters
)
bearing_tolerance: float = 15.0
coordinate_transform: bool = True
crop_method: Literal["interpolation", "convex_hull"] = "interpolation"
extrapolate: bool = False
best_fit_variogram: bool = True
kriging_parameters: Optional[dict] = None
latitude_resolution: float = 1.25
mesh_buffer_distance: float = 1.25
num_nearest_transects: int = 4
projection: Optional[str] = None
stratum: str = "ks"
variable: str = "biomass_density"
variogram_model: Union[str, List[str]] = ["bessel", "exponential"]
variogram_parameters: Optional[dict] = None
verbose: bool = True

# Create a copy of the existing variogram settings
default_variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()

self.analysis["settings"].update(
    {
        "kriging": {
            "exclude_age1": self.analysis["settings"]["transect"]["exclude_age1"],
            "extrapolate": extrapolate,
            "kriging_parameters": (
                self.input["statistics"]["kriging"]["model_config"]
                if kriging_parameters is None
                else kriging_parameters
            ),
            "projection": (
                self.config["geospatial"]["init"] if projection is None else projection
            ),
            "standardize_coordinates": coordinate_transform,
            "stratum": stratum.lower(),
            "stratum_name": "stratum_num" if stratum == "ks" else "inpfc",
            "variable": variable,
            "variogram_parameters": (
                self.input["statistics"]["variogram"]["model_config"]
                if variogram_parameters is None
                else variogram_parameters
            ),
            "verbose": verbose,
        }
    }
)

variogram_parameters = {
    **self.results["variogram"],
    **self.analysis["settings"]["variogram"],
    **self.results["variogram"]["model_fit"]
}

best_fit_variogram_parameters = {**self.results["variogram"], 
                                 **self.results["variogram"]["model_fit"]}

from echopop.spatial.variogram import initialize_variogram_parameters

self.analysis["settings"]["variogram"]
VariogramEmpirical.create(**self.analysis["settings"]["variogram"])
input_dict = self.input
analysis_dict = self.analysis
settings_dict = self.analysis["settings"]["kriging"]
mesh_cropping_parameters: MeshCrop = MeshCrop({})
kriging_parameters = {}
from echopop.spatial.transect import edit_transect_columns
from echopop.spatial.mesh import crop_mesh

# Validate the mesh cropping method parameters
cropping_params = MeshCrop.create(**mesh_cropping_parameters)

# Validate the empirical variogram parameters
empirical_vario_params = VariogramEmpirical.create(**variogram_parameters)

# Validate the variogram model parameters
model_vario_params = VariogramBase.create(**variogram_parameters)

# Validate the kriging parameters
# ---- Update using the `correlation_range`
kriging_parameters.update({"correlation_range": model_vario_params["correlation_range"]})
# ---- Generate dictionary
kriging_params = KrigingParameters.create(**kriging_parameters)

# Extract kriging mesh data
mesh_data = input_dict["statistics"]["kriging"]["mesh_df"]

# Extract the reference grid (200 m isobath)
isobath_data = input_dict["statistics"]["kriging"]["isobath_200m_df"]

# Define the and prepare the processed and georeferenced transect data
transect_data = edit_transect_columns(analysis_dict["transect"], settings_dict)

# Crop the mesh grid if the kriged data will not be extrapolated
if not settings_dict["extrapolate"]:
    # ---- Compute the cropped mesh
    mesh_full = crop_mesh(transect_data, mesh_data, cropping_params)
    if (settings_dict["verbose"]) & (cropping_params["crop_method"] == "convex_hull"):
        # ---- Print alert
        print(
            f"Kriging mesh cropped to prevent extrapolation beyond the defined "
            f"`mesh_buffer_distance` value ({cropping_params['mesh_buffer_distance']} nmi)."
        )        
else:
    # ---- Else, extract original mesh dataframe
    mesh_df = mesh_data.copy()
    # ---- Extract longitude column name
    mesh_longitude = [col for col in mesh_df.columns if "lon" in col.lower()][0]
    # ---- Latitude
    mesh_latitude = [col for col in mesh_df.columns if "lat" in col.lower()][0]
    # ---- Rename the dataframe
    mesh_full = mesh_df.copy().rename(
        columns={f"{mesh_longitude}": "longitude", f"{mesh_latitude}": "latitude"}
    )

# Standardize the x- and y-coordinates, if necessary
if settings_dict["standardize_coordinates"]:
    # ---- Transform transect data geometry (generate standardized x- and y-coordinates)
    transect_data, d_x, d_y = transform_geometry(transect_data, isobath_data, settings_dict)
    # ---- Transform mesh grid geometry (generate standardized x- and y-coordinates)
    mesh_full, _, _ = transform_geometry(mesh_full, isobath_data, settings_dict, d_x, d_y)
    if settings_dict["verbose"]:
        # ---- Print alert
        print(
            """Longitude and latitude coordinates (WGS84) converted to standardized """
            """coordinates (x and y)."""
        )
else:
    # ---- Else, duplicate the transect longitude and latitude coordinates as 'x' and 'y'
    # -------- x
    transect_data["x"] = transect_data["longitude"]
    # -------- y
    transect_data["y"] = transect_data["latitude"]
    # ---- Duplicate the mesh grid longitude and latitude coordinates as 'x' and 'y'
    # -------- x
    mesh_full["x"] = mesh_full["longitude"]
    # -------- y
    mesh_full["y"] = mesh_full["latitude"]
# --- Append to the analysis attribute
analysis_dict.update({"kriging": {"mesh_df": mesh_full, "transect_df": transect_data}})



def parameterize_variogram(
    default_variogram_parameters: dict,
    variogram_parameters: VariogramBase,
    best_fit_variogram: bool = False, 
    best_fit_variogram_parameters: Optional[VariogramBase] = None,      
):
    """
    Parameterize the variogram that will be used for kriging
    """

    # Use the best-fit parameters
    if best_fit_variogram:
        # ---- Validate each component of the variogram parameters
        # -------- Empirical
        empirical_variogram_params = VariogramEmpirical.create(**variogram_parameters)
        # -------- Base parameters
        kriging_variogram_params = VariogramBase.create(**best_fit_variogram_parameters)
    else:
        # Validate the relevant empirical variogram parameters
        empirical_variogram_params = VariogramEmpirical.create()

        # Initialize and validate the variogram model parameters
        valid_variogram_params = initialize_variogram_parameters(
            variogram_parameters, default_variogram_parameters
        )

        

self.analysis["settings"]["kriging"]

from echopop.utils.validate import VariogramBase, VariogramEmpirical, InitialValues
VariogramBase()
VariogramBase()
VariogramInitial[Union[List[str], Dict[str, InitialValues]]]

KrigingParameters.create(**{"correlation_range": 0.002})
MeshCrop.create(**{"bearing_tolerance": 21, "num_nearest_transects": 4.0})