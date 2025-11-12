from typing import Any, Dict, Union

import numpy as np
import pandas as pd
from numba import njit
from scipy.special import j1

from . import utils as ops


@njit(cache=True, fastmath=True)
def fast_phase(M1, delta_gamma_cos):
    r"""
    Compute complex phase exponentials for scattering calculations.

    This helper function calculates the phase component of the scattering form function for the
    PCDWBA using matrix exponentiation.

    Parameters
    ----------
    M1 : numpy.ndarray
        Phase matrix with shape (n_k, n_int) containing wavenumber and
        position-dependent phase terms
    delta_gamma_cos : numpy.ndarray
        Cosine of orientation angles with shape (n_int, n_theta) containing
        the angular dependence relative to the incident wave

    Returns
    -------
    numpy.ndarray
        Complex phase exponentials with shape (n_k, n_int, n_theta)

    Notes
    -----
    This function is JIT-compiled with numba for performance. The phase calculation follows:

    .. math::
        \Theta = \exp\left(i\frac{L}{a} k_fa \left(\frac{\vec{r_0}}{h} \cos \gamma\right) \right)
    """

    # Return
    return np.exp(1j * M1[:, :, np.newaxis] * delta_gamma_cos[np.newaxis, :, :])


@njit(fastmath=True)
def batch_interp(ARG_flat, ka_norm, J1):
    """
    Batched interpolation of cylindrical Bessel functions of the first kind of order 1.

    This helper function performs vectorized interpolation of cylindrical Bessel function of the
    first kind of order 1 values for multiple orientation angles simultaneously.

    Parameters
    ----------
    ARG_flat : numpy.ndarray
        Flattened argument array with shape (n_points, n_theta) containing the arguments for
        Bessel function interpolation
    ka_norm : numpy.ndarray
        Normalized wavenumber array for interpolation grid
    J1 : numpy.ndarray
        Pre-computed cylindrical Bessel function of the first kind of order 1

    Returns
    -------
    numpy.ndarray
        Interpolated Bessel function values with shape (n_points, n_theta)

    Notes
    -----
    This function is JIT-compiled with numba for performance. It performs linear interpolation of
    the Bessel function values for each orientation.
    """
    n_theta = ARG_flat.shape[1]
    n_points = ARG_flat.shape[0]
    out = np.empty((n_points, n_theta))
    for m in range(n_theta):
        out[:, m] = np.interp(ARG_flat[:, m], ka_norm, J1)
    return out


def pcdwba_fbs(
    taper_order: float,
    length_sd_norm: float,
    length_mean: float,
    length_radius_ratio: float,
    radius_of_curvature_ratio: float,
    theta_radians: Union[np.ndarray[float], float],
    k_f: Union[np.ndarray[float], float],
    ka_f: Union[np.ndarray[float], float],
    g: float,
    h: float,
    n_integration: int,
    n_wavelength: int,
):
    r"""
    Compute the linear scattering amplitude form function for the PCDWBA target strength model.

    This function implements the Phase-Compensated Distorted Wave Born Approximation (PCDWBA)[1]_ to
    calculate acoustic backscattering from uniformly bent cylindrical targets such as krill via the
    linear scattering amplitude (:math:`f_{bs}`).

    Parameters
    ----------
    taper_order : float
        Order parameter controlling body tapering
    length_sd_norm : float
        Normalized standard deviation of length distribution
    length_mean : float
        Mean length of targets in meters
    length_radius_ratio : float
        Ratio of length to equivalent cylindrical radius
    radius_of_curvature_ratio : float
        Ratio controlling body curvature
    theta_radians : numpy.ndarray or float
        Orientation angles in radians relative to incident wave
    k_f : numpy.ndarray
        Wavenumber array with shape (n_freq, n_sub_freq) in m^-1
    ka_f : numpy.ndarray
        Dimensionless wavenumber array (k * a) with same shape as k_f
    g : float
        Density contrast ratio (rho_body / rho_water)
    h : float
        Sound speed contrast ratio (c_body / c_water)
    n_integration : int
        Minimum number of integration points along body axis
    n_wavelength : int
        Number of integration points per wavelength

    Returns
    -------
    list
        List of complex form functions for each frequency, with each element containing an array of
        shape (n_freq_sub, n_theta)

    Notes
    -----
    The PCDWBA model accounts for:

    - Body curvature through uniformly bent cylinder geometry
    - Material property contrasts (density and sound speed)
    - Body tapering using power-law profiles
    - Phase compensation for curved geometries

    The form function is computed as:

    .. math::
        f_{\text{bs}}(k_f, \theta_m) =
        \frac{h^2 C_b dr_0}{4}
        \sum\limits_{j=1}^{n_f^\text{int}}
        \left[
            \left(
                \hat{k}_{f}a_j t_j
            \right) ^ 2
        \frac{\text{J}_1(2(\hat{k}_{f}a t)_j \cos(\beta_{jm})}{2(\hat{k}_{f}a t)_j \cos(\beta_{jm})}
        \exp({i\varphi_{fjm}})
        ~\right]

    where :math:`k_f` is the acoustic wavenumber of the surrounding medium (e.g. seawater) at
    frequency :math:`f` (Hz), and :math:`\theta` is the :math:`m^\text{th}` tilt angle (radians) of
    the cylinder subject to a plane wave where broadside incidence is at
    :math:`\theta_m = \frac{\pi}{2}`. The term :math:`\text{J}_1` is the cylindrical Bessel
    function of the first kind of order 1.

    Numerically, the integral is discretized into :math:`n_f^\text{int}` steps; however, the
    user-defined argument (`n_integration`) sets the minimum number of integration points. This
    enables an adaptive discretization rule that adjusts :math:`n_f^\text{int}` for each defined
    frequency. This operates in conjunction with :math:`n_{\lambda}`, which dictates the number of
    integration points per wavelength. Furthermore, a bandwidth surrounding the defined center
    frequencies is based on the variability in body length (:math:`\sigma_L`) bounded by
    :math:`\left[1 - 3.1 * \sigma_L, 1 + 3.1 * \sigma_L \right]`. These new bandwidths are used to
    compute the effective length:

    .. math::
        k_f L_\text{max} = \max_i \Big(k_f L \Big)
        \left( 1 + 3.1 * \sigma_L \right)

    When :math:`k_f L_\text{max} < n_f^\text{int}`, then the value for
    `n_integration` is used. When :math:`k_f L_\text{max} \leq n_f^\text{int}`, then:

    .. math::
        n_f^\text{int} =
        \left\leceil \frac{k_f L_\text{max} n_{\lambda}}{2\pi} \right\rceil

    The :math:`j^\text{th}` element of the position matrix (:math:`\vec{r_0}) is expressed by the
    radius (:math:`a_j`), taper coefficient (:math:`t_j`), and the along-axis tilt angle of the
    curved cylinder[2]_ (:math:`\beta_{jm}) relative to tilt angle :math:`m`. Variability in the
    position matrix cartesian coordinates is expressed as :math:`dr_0`. The cylinder is also
    characterized by its respective material properties:

    .. math::
        C_b =\gamma_k - \gamma_\rho = \frac{1 - g h^2}{g h^2} - \frac{g - 1}{g}

    where :math:`g` and :math:`h` are cylinders' density and sound speed relative to the
    surrounding medium. This includes the cylinder-specific wavenumber (:math:`k_f`) that accounts
    for :math:`h`.

    The phase (:math:`\exp(i\varphi_{fjm})`) is:

    .. math::
        \varphi_{fjm} =
        \left( \frac{L_j}{a_j} \right) k_f a_j \left( \frac{\vec{r_0}_j}{h} \right)
        \cos(\gamma_j - \theta_m)

    where \gamma_j is the slope between integration points relative to the overall tilt angle
    :math:`m`.

    References
    ----------
    .. [1] Chu, D., and Ye, Z. (1999). A phase-compensated distorted wave Born approximation
           representation of the bistatic scattering by weakly scattering objects: Application to
           zooplankton. The Journal of the Acoustical Society of America, 106, 1732-1743.
           doi: 10.1121/1.428036
    .. [2] Stanton, T.K. (1989). Sound scattering by cylinders of finite length. III. Deformed
           cylinders. The Journal of the Acoustical Society of America, 86, 691-705.
           doi: 10.1121/1.398193
    """

    kL_max = np.nanmax(k_f * length_mean, axis=1) * (1 + 3.1 * length_sd_norm)
    n_int = np.where(
        kL_max < n_integration, n_integration, np.ceil(kL_max * n_wavelength / (2 * np.pi))
    ).astype(int)

    taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(
        n_int, radius_of_curvature_ratio, taper_order
    )

    # Get the array sizes for later broadcasting
    # ---- Number of frequencies/wavenumbers
    n_k = np.apply_along_axis(ops.valid_array_row_length, 1, arr=ka_f)
    # ---- Number of orientation values
    n_theta = len(theta_radians)

    # Compute the reflection coefficient, `C_b`
    C_b = ops.reflection_coefficient(g, h)

    # Pre-allocate output (this will be indexed by center frequency)
    f_bs = []

    # Iterate across frequencies
    for i in range(len(n_k)):
        # ---- Reindex 'ka'
        ka_i = ka_f[i, : n_k[i]]
        # ---- Adjust `ka` to account for body shape tapering
        ka_tapered = ka_i.reshape(-1, 1) * taper[i, : n_int[i]] / h
        # ka_tapered = ka.reshape(-1, 1) * taper / h
        # ---- Adjust along-axis tilt angles and slopes to be relative to the incident planar wave
        # -------- Along-axis curvature slopes
        delta_gamma_cos = np.cos(gamma_tilt[i, : n_int[i]].reshape(-1, 1) - theta_radians)
        # -------- Along-axis tilt angles between segments
        delta_theta_cos = np.abs(np.cos(beta_tilt[i, : n_int[i]].reshape(-1, 1) - theta_radians))
        # ---- Generate the exponentiated matrix
        M1 = length_radius_ratio * ka_i.reshape(-1, 1) * (r_pos[i, : n_int[i]] / h)
        # ---- Generate the matrix that accounts for the material properties and position vector
        # ---- variability
        M2 = h**2 * C_b * dr_pos[i, : n_int[i]] / 4
        # ---- Normalize the `ka` vector
        ka_norm = np.linspace(-2 * ka_i[-1], 2 * ka_i[-1], 2 * n_int[i])
        # ---- Pre-compute the cylindrical Bessel function of the first kind of order 1
        J1 = j1(ka_norm)
        # ---- Broadcast `ka_tapered` for multiplication with `delta_gamma_cos`
        ARG = (
            2 * ka_tapered[:, :, np.newaxis] * delta_theta_cos[np.newaxis, :, :]
            + np.finfo(float).eps
        )
        # -------- Flatten for subsequent interpolation
        ARG_flat = ARG.ravel(order="F").reshape((n_k[i] * n_int[i], n_theta), order="F")
        # ---- Interpolate the values
        # J1_interp1 = np.array([np.interp(ARG_flat[:, m], ka_norm, J1) for m in range(n_theta)])
        J1_interp = batch_interp(ARG_flat, ka_norm, J1)
        # -------- Reshape and normalize
        J1_norm = (J1_interp / ARG_flat).reshape((n_k[i], n_int[i], n_theta), order="F")
        # ---- Exponentiate the terms to compute the phase
        phase = fast_phase(M1, delta_gamma_cos)
        # ---- Combine the adjusted `ka`, interpolated Bessel function output, and phase
        M3 = ka_tapered[:, :, np.newaxis] ** 2 * J1_norm * phase
        # ---- Compute the complex form function (matrix multiplication via Einstein summation)
        f_bs.append([np.einsum("ijk, j->ik", M3, M2) + np.finfo(float).eps])

    return f_bs


def pcdwba(
    center_frequencies: np.ndarray[float],
    length_mean: float,
    length_sd_norm: float,
    length_radius_ratio: float,
    taper_order: float,
    radius_of_curvature_ratio: float,
    theta_mean: float,
    theta_sd: float,
    orientation_distribution: Dict[str, Any],
    g: float,
    h: float,
    sound_speed_sw: float,
    frequency_interval: float,
    n_integration: int,
    n_wavelength: int,
    number_density: float,
    length_distribution: Dict[str, Any],
    **kwargs,
):
    """
    Phase-Corrected Distorted Wave Born Approximation (PCDWBA) for acoustic scattering.

    This function implements the PCDWBA[1]_ model for computing acoustic backscattering
    from elongated marine organisms modeled as uniformly bent fluid cylinders with
    tapered ends. The model accounts for organism size and orientation distributions
    through numerical integration.

    Parameters
    ----------
    center_frequencies : np.ndarray[float]
        Array of acoustic frequencies in Hz for scattering calculations.
    length_mean : float
        Mean organism length in meters.
    length_sd_norm : float
        Normalized standard deviation of length distribution (σ_L/μ_L).
    length_radius_ratio : float
        Ratio of organism length to equivalent cylindrical radius.
    taper_order : float
        Parameter controlling the sharpness of body tapering at ends.
    radius_of_curvature_ratio : float
        Ratio describing the radius of an osculating circle relative to the body length, which
        expresses the degree of curvature in a bent cylinder.
    theta_mean : float
        Mean orientation angle in degrees relative to vertical. Broadside incidence is considered
        to be 90°.
    theta_sd : float
        Standard deviation of orientation distribution in degrees.
    orientation_distribution : Dict[str, Any]
        Configuration for orientation averaging:
        - 'family': str, 'gaussian' or 'uniform'
        - 'bins': int, the number of bins each distribution is discretized into
    g : float
        Cylinder density (kg m^-3) relative to seawater.
    h : float
        Cylinder sound speed (m s^-1) relative to seawater.
    sound_speed_sw : float
        Sound speed in seawater in m s^-1.
    frequency_interval : float
        Frequency spacing in Hz for integration over scattering spectrum.
    n_integration : int
        Number of integration segments along organism axis.
    n_wavelength : int
        Number of integration points per acoustic wavelength.
    number_density : float
        Volumetric number density in individuals m^-3.
    length_distribution : Dict[str, Any]
        Configuration for length averaging:
        - 'family': str, 'gaussian' or 'uniform' distribution type
        - 'bins': int, the number of bins each distribution is discretized into
    **kwargs : dict
        Additional arguments passed to subfunctions

    Returns
    -------
    np.ndarray[float]
        Predicted volume backscattering strength (Sv) in dB re 1 m^-1
        for each input frequency. Array has same length as `center_frequencies`.

    Notes
    -----
    The PCDWBA model treats organisms as fluid cylinders with:
    1. **Uniform bending** described by radius_of_curvature_ratio
    2. **Tapered ends** controlled by taper_order parameter
    3. **Weak scattering** approximation (g,h ≈ 1)
    4. **Orientation averaging** over natural swimming orientations
    5. **Length averaging** over population size distributions

    The computational approach involves:
    - Discretization of organism body into integration segments
    - Calculation of scattering from each segment
    - Phase-correct summation accounting for body curvature
    - Statistical averaging over orientation and size distributions

    Model validity requires:
    - Organism length << acoustic wavelength for end-on incidence
    - |g-1|, |h-1| << 1 for weak scattering assumption
    - Sufficient integration points for curved/tapered geometries

    The model is particularly suitable for elongated zooplankton like
    euphausiids (krill) in the geometric scattering regime where
    ka ≈ 1-10 (k=wavenumber, a=characteristic size).

    Examples
    --------
    >>> # Typical krill parameters
    >>> freqs = np.array([38e3, 120e3])
    >>> Sv = pcdwba(
    ...     center_frequencies=freqs,
    ...     length_mean=0.025,  # 25 mm
    ...     length_sd_norm=0.15,
    ...     length_radius_ratio=18.0,
    ...     taper_order=10.0,
    ...     radius_of_curvature_ratio=3.0,
    ...     theta_mean=90.0,  # horizontal
    ...     theta_sd=15.0,
    ...     orientation_distribution={'family': 'gaussian', 'bins': 50},
    ...     g=1.02, h=1.02,
    ...     sound_speed_sw=1500.0,
    ...     frequency_interval=2000.0,
    ...     n_integration=50,
    ...     n_wavelength=10,
    ...     number_density=1000.0,
    ...     length_distribution={'family': 'gaussian', 'bins': 30}
    ... )

    References
    ----------
    .. [1] Chu, D., and Ye, Z. (1999). A phase-compensated distorted wave Born approximation
           representation of the bistatic scattering by weakly scattering objects: Application to
           zooplankton. The Journal of the Acoustical Society of America, 106, 1732-1743.
           doi: 10.1121/1.428036
    """

    # Pre-allocate arrays based on input size
    n_theta = orientation_distribution["bins"]
    n_length = length_distribution["bins"]

    # Generate frequency intervals centered on the central frequencies
    frequencies = ops.generate_frequency_interval(
        center_frequencies, length_sd_norm, frequency_interval
    )

    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_c = ops.wavenumber(center_frequencies, sound_speed_sw)
    # ---- Compute ka (center frequencies)
    ka_c = k_c * length_mean / length_radius_ratio
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k_f = ops.wavenumber(frequencies, sound_speed_sw)
    # -------- Now `ka`
    ka_f = k_f * length_mean / length_radius_ratio

    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        theta_mean - 3.1 * theta_sd,
        theta_mean + 3.1 * theta_sd,
        n_theta,
    )
    # ---- Convert to radians
    theta_radians = theta_values * np.pi / 180.0

    # Compute over vector lengths
    length_values = np.linspace(
        length_mean - 3 * (length_sd_norm * length_mean),
        length_mean + 3 * (length_sd_norm * length_mean),
        n_length,
    )

    f_bs = pcdwba_fbs(
        taper_order,
        length_sd_norm,
        length_mean,
        length_radius_ratio,
        radius_of_curvature_ratio,
        theta_radians,
        k_f,
        ka_f,
        g,
        h,
        n_integration,
        n_wavelength,
    )

    # Orientation averaging
    f_bs_orientation = ops.orientation_average(
        theta_values, f_bs, theta_mean, theta_sd, "uniform"  # orientation_distribution["family"]
    )

    # Length-averaged sigma_bs (normalized to length)
    sigma_bs_length = ops.length_average(
        length_values,
        ka_f,
        ka_c,
        f_bs_orientation,
        length_mean,
        length_mean * length_sd_norm,
        length_distribution["family"],
    )

    # Convert to sigma_bs (linear backscattering cross-section)
    sigma_bs = sigma_bs_length * (length_mean) ** 2

    # Switch to logarithmic domain to compute S_V (volumetric backscattering strength)
    Sv_prediction = 10 * np.log10(number_density * sigma_bs)

    return Sv_prediction


def uniformly_bent_cylinder(
    n_segments: Union[int, np.ndarray[int]],
    radius_of_curvature_ratio: float,
    taper_order: float,
) -> pd.DataFrame:
    """
    Generate geometric parameters for uniformly bent cylinder with tapered ends.

    This function computes the spatial discretization and geometric properties
    of organisms modeled as uniformly bent cylinders with tapering. It calculates
    position vectors, orientation angles, and tapering factors needed for
    acoustic scattering calculations in the PCDWBA model.

    Parameters
    ----------
    n_segments : Union[int, np.ndarray[int]]
        Number of integration segments along the organism axis.
        If int: single organism with n_segments divisions
        If array: multiple organisms with potentially different segment counts
        Must be > 0. Typical range: 20-100 for adequate geometric resolution.
    radius_of_curvature_ratio : float
        Ratio describing body curvature relative to organism length.
        Must be > 0. Larger values indicate straighter bodies:
        - ~1: highly curved (C-shaped)
        - ~3-5: moderately curved (typical krill)
        - >10: nearly straight
    taper_order : float
        Parameter controlling end tapering sharpness.
        Must be > 0. Higher values create more pointed ends:
        - 2: parabolic tapering
        - 10: typical moderate tapering
        - >20: sharp, nearly conical ends

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        taper : np.ndarray
            Tapering factor along body axis. Values in [0,1] where 1 = full
            radius at body center, 0 = pointed ends.
        gamma_tilt : np.ndarray
            Local tilt angles (radians) of curved body segments relative to
            the body coordinate system. Used for phase calculations.
        beta_tilt : np.ndarray
            Local orientation angles (radians) of segments relative to
            incident acoustic wave direction. Critical for scattering.
        r_pos : np.ndarray
            Position vector along curved body axis in normalized coordinates.
            Euclidean distances from body center to each segment.
        dr_pos : np.ndarray
            Differential position vector - incremental distances between
            adjacent segments. Used for integration weighting.

    Notes
    -----
    The uniformly bent cylinder model represents marine organisms with:

    1. **Uniform curvature**: Constant radius of curvature along body axis
    2. **Symmetric tapering**: Gradual radius reduction toward both ends
    3. **Smooth geometry**: Continuous derivatives for stable numerics

    Mathematical formulation:
    - Body axis parameterized as z ∈ [-1, 1] in normalized coordinates
    - Curvature parameter γ = 0.5/radius_of_curvature_ratio
    - Bent coordinates: x = 1 - √(1 - (sin(γz))²), z' = sin(γz)
    - Tapering: taper(z) = √(1 - z^taper_order)

    The local angles are crucial for scattering calculations:
    - gamma_tilt: orientation of body segments in curved geometry
    - beta_tilt: orientation relative to acoustic incidence direction

    These angles determine the effective scattering cross-section and
    phase relationships in the PCDWBA integration.

    The function handles both single organisms (n_segments=int) and
    populations with varying discretization (n_segments=array) for
    efficient batch processing in population-level calculations.

    Examples
    --------
    >>> # Single moderately curved organism with 50 segments
    >>> taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(
    ...     n_segments=50,
    ...     radius_of_curvature_ratio=3.0,  # moderate curvature
    ...     taper_order=10.0  # moderate tapering
    ... )
    >>>
    >>> # Batch processing for multiple organisms
    >>> n_segments_array = np.array([30, 50, 40])  # different discretizations
    >>> results = uniformly_bent_cylinder(
    ...     n_segments=n_segments_array,
    ...     radius_of_curvature_ratio=4.0,
    ...     taper_order=8.0
    ... )

    References
    ----------
    .. [1] Stanton, T.K. (1989). Sound scattering by cylinders of finite length. III. Deformed
           cylinders. The Journal of the Acoustical Society of America, 86, 691-705.
           doi: 10.1121/1.398193
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
    n_valid = np.apply_along_axis(ops.valid_array_row_length, 1, arr=alpha_tilt)
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
