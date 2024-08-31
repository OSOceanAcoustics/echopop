from typing import (
    Any,
    Callable,
    Dict,
    List,
    Literal,
    Mapping,
    Optional,
    Tuple,
    Type,
    TypedDict,
    Union,
)

import numpy as np
from echopop.survey import Survey
from echopop.spatial.variogram import create_optimization_options

file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
survey = Survey(init_config, file_config)
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis(exclude_age1=True)
self = survey
# survey.fit_variogram()
# self.analysis['settings']['variogram']['optimization']

##################
# ---- Update `VariogramBase`
from echopop.utils.validate import VariogramBase, VariogramOptimize, VariogramInitial, VariogramEmpirical
model = ["bessel", "exponential"]
n_lags = 30
azimuth_range = 360.0
force_lag_zero = True

standardize_coordinates = True

variable = "biomass"
verbose = True
input_variogram_parameters = {"sill": 1.0}
input_optimization_parameters = {"max_fun_evaluations": 400}
input_initial = {
    "nugget": {"min": 0.1, "value": 1.0},
    "sill": {"min": 0.0, "value": 2.0},
    "correlation_range": {"min": 0.0},
    "hole_effect_range": {"min": 0.0, "max": 5.0},
    "decay_power": {"max": np.inf}
}

# Validate "variable" input
if variable not in ["biomass", "abundance"]:
    raise ValueError(
        f"The user input for `variable` ({variable}) is invalid. Only `variable='biomass'` "
        f"and `variable='abundance'` are valid inputs for the `fit_variogram()` method."
    )

# Initialize and validate the variogram parameters
# --------
# Get the default variogram parameters
default_variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
# ---- Update model, n_lags
default_variogram_parameters.update({
    "model": model,
    "n_lags": n_lags
    })
# ---- Update the defaults for the base variogram model parameters
VariogramBase.update_defaults(default_variogram_parameters)
# ---- Initialize and validate the theoretical variogram parameters
variogram_parameters = VariogramBase.create(**input_variogram_parameters)






# variogram_parameters = variogram_parameters
# optimization_parameters = optimization_parameters
transect_dict = self.analysis["transect"]
settings_dict = self.analysis["settings"]["variogram"]
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]
# initial_values = initial_values

# Prepare the transect data
# ---- Create a copy of the transect dictionary
transect_input = copy.deepcopy(transect_dict)
# ---- Edit the transect data
transect_data = edit_transect_columns(transect_input, settings_dict)

# Standardize the transect coordinates, if necessary
if settings_dict["standardize_coordinates"]:
    # ---- Transform geometry
    transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)
    # ---- Print message if verbose
    if settings_dict["verbose"]:
        # ---- Print alert
        print(
            "Longitude and latitude coordinates (WGS84) converted to standardized "
            "coordinates (x and y)."
        )
else:
    # ---- x
    transect_data["x"] = "longitude"
    # ---- y
    transect_data["y"] = "latitude"

from echopop.spatial.variogram import empirical_variogram, optimize_variogram
# Validate the remaining user inputs for the empirical variogram
empirical_variogram_params = VariogramEmpirical.create(**{"azimuth_range": azimuth_range, 
                                                          "force_lag_zero": force_lag_zero, 
                                                          "standardize_coordinates": standardize_coordinates})
# Compute the empirical variogram
lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
    transect_data, {**variogram_parameters, **empirical_variogram_params}, settings_dict
)

optimization_settings = {"parameters": parameters, "config": lmfit_parameters}

# Least-squares fitting
best_fit_variogram, initial_fit, optimized_fit = optimize_variogram(
    lag_counts, lags, gamma_h, variogram_parameters, optimization_settings
)

# Return a dictionary of results
return {
    "best_fit_parameters": best_fit_variogram,
    "initial_fit": {
        "parameters": dict(zip(initial_fit[0], initial_fit[1])),
        "MAD": initial_fit[2],
    },
    "optimized_fit": {
        "parameters": dict(zip(optimized_fit[0], optimized_fit[1])),
        "MAD": optimized_fit[2],
    },
}
