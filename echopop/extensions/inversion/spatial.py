import copy
import warnings
from typing import Any, Dict

import numpy as np
import pandas as pd

from ...spatial.variogram import (
    empirical_variogram,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
)
from ...spatial.projection import transform_geometry
from ...spatial.mesh import crop_mesh
from ...spatial.krige import kriging
from ...utils.validate_dict import (
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
)

def inversion_variogram_analysis(
    transect_df: pd.DataFrame,
    variogram_parameters: dict,
    default_variogram_parameters: dict,
    optimization_parameters: dict,
    initialize_variogram: dict,
    settings_dict: dict, 
):
    """
    Fit optimized variogram to inverted georeferenced data
    """

    # Validate the relevant empirical variogram parameters
    empirical_variogram_params = VariogramEmpirical.create(**settings_dict)

    # Initialize and validate the variogram model parameters
    valid_variogram_params = initialize_variogram_parameters(
        variogram_parameters, default_variogram_parameters
    )

    # Initialize and validate the optimization parameters
    valid_optimization_params = initialize_optimization_config(optimization_parameters)

    # Initialize and validate the initial values/boundary inputs
    valid_initial_values = initialize_initial_optimization_values(
        initialize_variogram, valid_variogram_params
    )

    # Create copy of data
    transect_data = transect_df.copy()

    # Compute the empirical variogram
    lags, gamma_h, lag_counts, _ = empirical_variogram(
        transect_data, {**valid_variogram_params, **empirical_variogram_params}, settings_dict
    )

    # Least-squares fitting
    # ---- Consolidate the optimization dictionaries into a single one
    optimization_settings = {
        "parameters": valid_initial_values,
        "config": valid_optimization_params,
    }
    # ---- Optimize parameters
    best_fit_variogram, initial_fit, optimized_fit = optimize_variogram(
        lag_counts, lags, gamma_h, optimization_settings, **valid_variogram_params
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

def krige_inverted_data(
    transect_df: pd.DataFrame,
    mesh_df: pd.DataFrame,
    analysis_dict: Dict[str, Any],
    settings_dict: Dict[str, Any],
) -> None:
    """
    Krige inverted population estimates
    """
    # Validate cropping method parameters
    validated_cropping_methods = MeshCrop.create(**settings_dict["cropping_parameters"])
    # ---- Update the dictionary
    settings_dict["cropping_parameters"].update({**validated_cropping_methods,
                                                "projection": settings_dict["projection"]})

    # Validate the variogram parameters
    valid_variogram_parameters = VariogramBase.create(**settings_dict["variogram_parameters"])
    # ---- Update the dictionary
    settings_dict["variogram_parameters"].update({**settings_dict["variogram_parameters"],
                                                **valid_variogram_parameters})

    # Validate kriging parameters
    valid_kriging_parameters = {
        **settings_dict["kriging_parameters"],
        **KrigingParameterInputs.create(
            **{**settings_dict["kriging_parameters"], **settings_dict["variogram_parameters"]}
        )
    }
    # ---- Update the dictionary
    settings_dict["kriging_parameters"].update({**valid_kriging_parameters})

    # Crop the mesh grid if the kriged data will not be extrapolated
    if settings_dict["extrapolate"]:
        # ---- Else, extract original mesh dataframe
        mesh_data = mesh_df.copy()
        # ---- Extract longitude column name
        mesh_longitude = [col for col in mesh_data.columns if "lon" in col.lower()][0]
        # ---- Latitude
        mesh_latitude = [col for col in mesh_data.columns if "lat" in col.lower()][0]
        # ---- Rename the dataframe
        mesh_full = mesh_data.copy().rename(
            columns={f"{mesh_longitude}": "longitude", f"{mesh_latitude}": "latitude"}
        )
    else:
        # ---- Compute the cropped mesh
        mesh_full = crop_mesh(
            transect_df, mesh_df, settings_dict["cropping_parameters"]
        )
        # ---- Print message, if verbose
        if (settings_dict["verbose"]) and (
            validated_cropping_methods["crop_method"] == "convex_hull"
        ):
            # ---- Print alert
            print(
                f"Kriging mesh cropped to prevent extrapolation beyond the defined "
                f"`mesh_buffer_distance` value "
                f"({validated_cropping_methods['mesh_buffer_distance']} nmi)."
            )
    # --- Append to the analysis attribute
    analysis_dict.update({"mesh_df": mesh_full})

    # Kriged results
    kriged_results = kriging(
        transect_df, mesh_full, settings_dict
    )

    # Return kriged (interpolated) results
    return kriged_results
    