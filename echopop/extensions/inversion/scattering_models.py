from typing import Any, Dict, List, Union

from numpy.typing import ArrayLike
from pydantic import BaseModel
from scipy.special import j1
import numpy as np

from .math import reflection_coefficient, spherical_hn, valid_array_row_length

def pcdwba(
    taper: float, 
    gamma_tilt: np.ndarray[float], 
    beta_tilt: np.ndarray[float], 
    r_pos: np.ndarray[float], 
    dr_pos: np.ndarray[float],
    length_radius_ratio: float,
    g: Union[np.ndarray[float], float],
    h: Union[np.ndarray[float], float],
    ka: Union[np.ndarray[float], float],
    theta: Union[np.ndarray[float], float],
) -> List[np.ndarray[complex]]:
    """
    Phase-compensated distorted wave Born approximation (DWBA)

    Defined as [1]_,

    ..math::
        f(\phi)=(kaT)^2\frac{J_1(2kaT \cos{\theta})}{2kaT \cos{\theta}}
        e^{i\varepsilon ka(\vec{r}_{pos}/h)\cos{\beta}}

    where :math:`f` is the scattering amplitude, :math:`\phi` is the scatterer orientation angle
    relative to the incident sound wave, :math:`k` is the acoustic wavenumber, :math:`T` is the
    taper coefficient, :math:`a` is the radius at a specific point along the body, :math:`J_1` is
    the cylindrical Bessel function of the first kind, :math:`\theta` is the orientation angle
    relative to the incident sound wave at a specific point along the body, :math:`\varepsilon` is
    the length-to-radius ratio, :math:`\vec{r}_{pos}` is the positional vector, :math:`h` is the
    soundspeed contrast, and :math:`\beta` is the orientation of a specific point along the body.

    References
    ----------
    ..[1] Chu, D., and Ye, Z. (1999). A phase-compensated distorted wave Born approximation
    representation of the bistatic scattering by weakly scattering objects: Application to
    zooplankton

    """

    # Get the array sizes for later broadcasting
    # ---- Number of frequencies/wavenumbers
    n_k = (
        np.apply_along_axis(valid_array_row_length, axis=1, arr=ka) 
        if isinstance(ka, np.ndarray) else len(ka)
    )
    # ---- Number of segments and minimum number of integration points
    n_segments = (
        np.apply_along_axis(valid_array_row_length, axis=1, arr=r_pos) 
        if r_pos.ndim > 1 else len(r_pos)
    )
    # ---- Number of orientation values
    n_theta = len(theta)

    # Compute the reflection coefficient, `C_b`
    C_b = reflection_coefficient(g, h)

    # Pre-allocate output (this will be indexed by center frequency)
    f_bs = []
    
    # Iterate across frequencies
    for i in range(len(n_k)):
        # ---- Reindex 'ka'
        ka_i = ka[i, :n_k[i]]
        # ---- Adjust `ka` to account for body shape tapering 
        ka_tapered = ka_i .reshape(-1, 1) * taper[i, :n_segments[i]] / h
        # ka_tapered = ka.reshape(-1, 1) * taper / h
        # ---- Adjust along-axis tilt angles and slopes to be relative to the incident planar wave
        # -------- Along-axis curvature slopes
        delta_gamma_cos = np.cos(gamma_tilt[i, :n_segments[i]].reshape(-1, 1) - theta)
        # -------- Along-axis tilt angles between segments
        delta_theta_cos = np.abs(np.cos(beta_tilt[i, :n_segments[i]].reshape(-1, 1) - theta))
        # ---- Generate the exponentiated matrix
        M1 = length_radius_ratio * ka_i.reshape(-1, 1) * (r_pos[i, :n_segments[i]] / h)
        # ---- Generate the matrix that accounts for the material properties and position vector 
        # ---- variability
        M2 = h ** 2 * C_b * dr_pos[i, :n_segments[i]] / 4
        # ---- Normalize the `ka` vector
        ka_norm = np.linspace(-2 * ka_i[-1], 2*ka_i[-1], 2*n_segments[i])
        # ---- Pre-compute the cylindrical Bessel function of the first kind of order 1
        J1 = j1(ka_norm)
        # ---- Broadcast `ka_tapered` for multiplication with `delta_gamma_cos`
        ARG = (
            2 * ka_tapered[:, :, np.newaxis] 
            * delta_theta_cos[np.newaxis, :, :] + np.finfo(float).eps
        )
        # -------- Flatten for subsequent interpolation
        ARG_flat = ARG.ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F")
        # ---- Interpolate the values
        J1_interp = np.array([
            np.interp(ARG_flat[:, m], ka_norm, J1) for m in range(n_theta)
        ])
        # -------- Reshape and normalize
        J1_norm = (J1_interp.T / ARG_flat).reshape((n_k[i], n_segments[i], n_theta), order="F")
        # ---- Exponentiate the terms to compute the phase
        phase = np.exp(1j * M1[:, :, np.newaxis] * delta_gamma_cos[np.newaxis, :, :])
        # ---- Combine the adjusted `ka`, interpolated Bessel function output, and phase
        M3 = ka_tapered[:, :, np.newaxis] ** 2 * J1_norm * phase
        # ---- Compute the complex form function (matrix multiplication via Einstein summation)
        f_bs.append([np.einsum("ijk, j->ik", M3, M2) + np.finfo(float).eps])

    # Return the form function
    return f_bs

class validate_pcdwba(BaseModel):
    """
    Pydantic model for validating scattering model parameters specific to the PCDWBA
    """

    # RETURNS: Dict[str, Any]
    pass
