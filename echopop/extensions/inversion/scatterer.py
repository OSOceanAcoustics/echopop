"""
Scatterer class used for parameterizing and computing acoustic scattering models
"""

from typing import Any, Dict, Literal, Union

import numpy as np
import pandas as pd

from .math import length_average, orientation_average, valid_array_row_length
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
    f_bs_orientation = orientation_average(theta_values, form_function, theta_mean, theta_sd)

    # Length-averaged sigma_bs (normalized to length)
    sigma_bs_length = length_average(
        length_values, ka, ka_center, f_bs_orientation, length_mean, length_deviation
    )
    # ---- Convert to sigma_bs (linear backscattering cross-section)
    sigma_bs = sigma_bs_length * (length_mean) ** 2

    # Switch to logarithmic domain to compute S_V (volumetric backscattering strength)
    Sv_prediction = 10 * np.log10(number_density * sigma_bs)

    # Return predicted Sv
    return Sv_prediction


def uniformly_bent_cylinder(
    n_segments: Union[int, np.ndarray[int]],
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
    # ---- If array of values
    if isinstance(n_segments, np.ndarray):
        # ---- Maximum dimension
        max_n = n_segments.max()
        # ---- Generate the linearly spaced indices
        indices = np.arange(max_n)
        # ---- Broadcast a map for valid segment ranges
        valid_mask = indices < n_segments[:, None]
        # ---- Create a placeholder array padded with NaN
        z = np.full((n_segments.size, max_n), np.nan)
        # ---- Fill with the valid values
        z[valid_mask] = np.concatenate([np.linspace(-1.0, 1.0, n) for n in n_segments])
    # ---- If only a single value
    else:
        # ---- Maximum dimension
        max_n = n_segments
        # ---- Create vector
        z = np.linspace(-1.0, 1.0, n_segments)

    # Compute the taper vector
    taper = np.sqrt(1 - z**taper_order)

    # Bend the cylinder
    # ---- z-axis
    z_curved = np.sin(gamma) * z
    # ---- Dorsoventral axis (x-axis)
    x_curved = 1 - np.sqrt(1 - z_curved**2)

    # Normalize the curvature
    # ---- z-axis
    z_norm = z_curved * norm_ratio
    # ---- x-axis
    x_norm = x_curved * norm_ratio

    # Calculate the slope between curved segments
    gamma_tilt = np.arctan2(z_norm, x_norm)

    # Calculate the orientation angles between curved segments
    # ---- z-axis differences
    dz = np.diff(z_norm)
    # ---- x-axis differences
    dx = np.diff(x_norm) + np.finfo(float).eps
    # ---- alpha tilt angles
    alpha_tilt = np.arctan(dz / dx)
    # ---- Get the valid number of values per row
    n_valid = np.apply_along_axis(valid_array_row_length, 1, arr=alpha_tilt)
    # ---- Preallocate array
    new_column = np.array([])
    # ---- Generate values that are appended to the last valid column
    if alpha_tilt.ndim > 1:
        for i in range(alpha_tilt.shape[0]):
            if n_valid[i] == max_n - 1:
                new_column = np.append(new_column, np.arctan(dz[i, -1] / dx[i, -1]))
            else:
                new_column = np.append(new_column, np.nan)
                alpha_tilt[i, n_valid[i]] = np.arctan(dz[i, n_valid[i] - 1] / dx[i, n_valid[i] - 1])
        # ---- Now concatenate the missing column
        alpha_tilt = np.concatenate([alpha_tilt, new_column.reshape(-1, 1)], axis=1)
    else:
        alpha_tilt = np.append(alpha_tilt, np.arctan(dz[-1] / dx[-1]))
    # ---- beta tilt angles
    beta_tilt = np.where(alpha_tilt >= 0.0, alpha_tilt - np.pi / 2, alpha_tilt + np.pi / 2)

    # Compute the along-axis Euclidean distances to construct the position vector
    # ---- Position vector
    r_pos = np.sqrt(x_norm**2 + z_norm**2)
    # ---- Generate values that are appended to the first valid column
    if r_pos.ndim > 1:
        dr_pos = np.concatenate(
            [
                np.sqrt(dx[:, 0] * dx[:, 0] + dz[:, 0] * dz[:, 0]).reshape(-1, 1),
                np.sqrt(dx * dx + dz * dz),
            ],
            axis=1,
        )
    else:
        # ---- Compute the derivative of the position vector derivative
        dr_pos = np.append(np.sqrt(dx[0] * dx[0] + dz[0] * dz[0]), np.sqrt(dx * dx + dz * dz))

    # Return the relevant parameters
    return taper, gamma_tilt, beta_tilt, r_pos, dr_pos
