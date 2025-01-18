"""
Scatterer class used for parameterizing and computing acoustic scattering models
"""

from typing import Any, Dict, Literal

import numpy as np
from .math import length_average, orientation_average

import pandas as pd

from .scattering_models import pcdwba


# class Scatterer:

#     def __init__(
#         self,
#         parameters: Dict[str, Any],
#         scattering_type: Literal[
#             "elastic_shelled",
#             "fluid_like",
#             "resonant",
#         ] = "fluid_like",
#         shape: Literal[
#             "arbitrary", "cylinder", "oblate_spheroid", "prolate_spheroid", "sphere", "umbrella"
#         ] = "arbitrary",
#     ):

#         # Initialize attributes
#         # ---- Model parameters
#         self.parameters = parameters
#         # ---- Scattering-type
#         self.scattering_type = scattering_type
#         # ---- Shape-type
#         self.shape = shape

#         # Build position matrix
#         self.position_vector = construct_position_vector(self.shape, self.parameters)

#     def compute_ts(
#         self, scattering_parameters: Dict[str, Any], position_vector: pd.DataFrame
#     ) -> pd.DataFrame:
#         """
#         Compute theoretical target strength (TS)
#         """

#         # PHASE 1) EXTRACT PARAMETERS
#         # PHASE 2) MAP RELEVANT SCATTERING CALLABLE FUNCTION (e.g. `pcdwba`)
#         # PHASE 3) COMPUTE TS

#         # RETURNS: DataFrame with TS values mapped to specific parameter values
#         pass

#     def predict_Sv(self, scattering_parameters: Dict[str, Any]) -> pd.DataFrame:

#         # PHASE 1) EXTRACT PARAMETERS
#         # PHASE 2) AVERAGE SIGMA_BS (LIKE IN EchoPro_matlab_krill_inversion::DWBAscat_simple1.m)

#         # RETURNS: DataFrame (or np.ndarray?) of predicted Sv values
#         pass


####################################################################################################
# UTILITY FUNCTIONS
####################################################################################################

def compute_Sv(
    number_density: float,
    theta_values: np.ndarray[float],
    theta_mean: float,
    theta_sd: float,
    length_values: np.ndarray[float],
    length_mean: float,
    length_deviation: float,
    form_function: np.ndarray[complex],
    ka: np.ndarray[float],
    ka_center: np.ndarray[float],
) -> Dict[str, Any]: 
    """
    Predict the volumetric backscattering strength from theoretical scattering estimates
    """

    # Orientation-averaged linear scattering coefficient
    orientation_pdf, f_bs_orientation = orientation_average(theta_values, 
                                                            form_function, 
                                                            theta_mean, 
                                                            theta_sd)

    # Length-averaged sigma_bs (normalized to length)
    length_pdf, sigma_bs_length = length_average(length_values, 
                                                 ka, 
                                                 ka_center, 
                                                 f_bs_orientation, 
                                                 length_mean, 
                                                 length_deviation)
    # ---- Convert to sigma_bs (linear backscattering cross-section)
    sigma_bs = sigma_bs_length * (length_mean) ** 2

    # Convert sigma_bs to s_v (linear volumetric backscattering coefficient)
    sv_prediction = sigma_bs.sum(axis=0)

    # Switch to logarithmic domain to compute S_V (volumetric backscattering strength)
    Sv_prediction = 10 * np.log10(number_density * sv_prediction.sum())

    # Return a formatted dictionary
    return {
        "orientation": {
            "estimates": f_bs_orientation,
            "pdf": orientation_pdf,
        },
        "length": {
            "estimates": sigma_bs_length,
            "pdf": length_pdf,
        },
        "predicted_Sv": Sv_prediction,
    }
    

def uniformly_bent_cylinder(
    n_segments: int,
    radius_of_curvature_ratio: float,
    taper_order: float,
) -> pd.DataFrame:
    """
    Generates the normalized position matrix for an uniformly bent cylinder
    """

    # Curvature coefficients (for shape normalization)
    # ---- Gamma
    gamma = 0.5 / radius_of_curvature_ratio
    # ---- Normalization
    norm_ratio = radius_of_curvature_ratio * 2

    # Create normalized horizontal increments along the anterioposterior (z) axis of the body shape
    z = np.linspace(-1, 1, n_segments)

    # Compute the taper vector
    taper = np.sqrt(1-z**taper_order)

    # Bend the cylinder
    # ---- z-axis
    z_curved = np.sin(gamma) * z
    # ---- Dorsoventral axis (x-axis)
    x_curved = 1 - np.sqrt(1-z_curved**2)

    # Normalize the curvature
    # ---- z-axis
    z_norm = z_curved * norm_ratio
    # ---- x-axis
    x_norm = x_curved * norm_ratio

    # Calculate the slope between curved segments
    gamma_tilt = np.arctan2(z_norm, x_norm)

    # Caluculate the orientation angles between curved segments 
    # ---- z-axis differences
    dz = np.diff(z_norm)
    # ---- x-axis differences
    dx = np.diff(x_norm) + np.finfo(float).eps
    # ---- alpha tilt angles
    alpha_tilt = np.append(np.arctan(dz/dx), np.arctan(dz[-1]/dx[-1]))
    # ---- beta tilt angles
    beta_tilt = np.where(alpha_tilt >= 0.0, alpha_tilt - np.pi/2, alpha_tilt + np.pi/2)

    # Compute the along-axis Euclidean distances to construct the position vector
    # ---- Position vector
    r_pos = np.sqrt(x_norm**2 + z_norm**2)
    # ---- Compute the derivative of the position vector derivative
    dr_pos = np.append(np.sqrt(dx[0]*dx[0] + dz[0]*dz[0]), np.sqrt(dx*dx + dz*dz))
    
    # Return the relevant parameters
    return taper, gamma_tilt, beta_tilt, r_pos, dr_pos
