import inspect
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from lmfit import Minimizer, Parameters
from scipy import special

# Set warnings filter
warnings.simplefilter("always")


# Single family models


# ---- Cubic
def cubic(distance_lags, sill, nugget, correlation_range):
    """
    Cubic variogram model with smooth transitions and finite range.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range where correlation becomes zero.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The cubic variogram model is defined as:

    .. math::
        γ(h) = \\begin{cases}
        C_0 + C_1 \\left(7\\left(\\frac{h}{a}\\right)^2 \\right. \\\\
        \\qquad\\qquad \\left. - \\frac{35}{4}\\left(\\frac{h}{a}\\right)^3
        + \\frac{7}{2}\\left(\\frac{h}{a}\\right)^5 \\right. \\\\
        \\qquad\\qquad \\left. - \\frac{3}{4}\\left(\\frac{h}{a}\\right)^7\\right) & h ≤ a \\\\
        C_0 + C_1 & h > a
        \\end{cases}

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the effective range

    The cubic model provides very smooth transitions with continuous first and second derivatives
    at the origin and at the range. It's suitable for highly regular spatial processes requiring
    smooth interpolation.

    References
    ----------
    .. [1] Journel, A.G. & Huijbregts, C.J. (1978). Mining Geostatistics. Academic Press.
    .. [2] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    """

    # Calculate normalized distances
    h_norm = distance_lags / correlation_range

    # Calculate cubic correlation for h <= a
    correlation = np.where(
        h_norm <= 1,
        1 - (7 * h_norm**2 - 35 / 4 * h_norm**3 + 7 / 2 * h_norm**5 - 3 / 4 * h_norm**7),
        0,
    )

    # Calculate variogram
    return nugget + (sill - nugget) * (1 - correlation)


# ---- Exponential
def exponential(distance_lags: np.ndarray, sill: float, nugget: float, correlation_range: float):
    """
    Exponential variogram model for spatial correlation analysis.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range parameter controlling spatial correlation decay.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The exponential variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\frac{h}{a}\\right)\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the correlation range parameter

    The exponential model exhibits monotonic decay and reaches 95% of the sill at distance 3a. It
    represents processes with continuous but non-differentiable spatial correlation.

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-(distance_lags / correlation_range))

    # Compute the exponential semivariogram
    return partial_sill * decay + nugget


# ---- Gaussian
def gaussian(distance_lags: np.ndarray, sill: float, nugget: float, correlation_range: float):
    """
    Gaussian (squared exponential) variogram model for spatial correlation analysis.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range parameter controlling spatial correlation decay.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The Gaussian variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\frac{h^2}{a^2}\\right)\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the correlation range parameter

    The Gaussian model exhibits very smooth spatial transitions with infinitely differentiable
    correlation functions. It approaches the sill more gradually than exponential models and is
    suitable for highly regular spatial processes.

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Rasmussen, C.E. & Williams, C.K.I. (2006). Gaussian Processes for Machine Learning. MIT
           Press.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-(distance_lags**2 / correlation_range**2.0))

    # Compute the Gaussian semivariogram
    return partial_sill * decay + nugget


# ---- J-Bessel
def jbessel(distance_lags: np.ndarray, sill: float, nugget: float, hole_effect_range: float):
    """
    J-Bessel variogram model exhibiting hole-effect (periodic) behavior.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    hole_effect_range : float
        Range parameter controlling the periodicity of oscillations.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The J-Bessel variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\frac{J_1(h/a)}{h/(2a)}\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - J₁ is the first-order Bessel function of the first kind
    - a is the hole-effect range parameter

    This model exhibits damped oscillatory behavior, creating "hole effects" where the variogram
    exceeds the sill before returning. It's suitable for processes with regular spatial patterns or
    cyclical structures.

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - special.j0(hole_effect_range * distance_lags)

    # Compute the J-Bessel semivariogram
    return partial_sill * decay + nugget


# ---- K-Bessel
def kbessel(distance_lags: np.ndarray, sill: float, nugget: float, hole_effect_range: float):
    """
    K-Bessel variogram model with modified Bessel function of the second kind.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range parameter controlling spatial correlation decay.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The K-Bessel variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\frac{2^{1-ν}}{Γ(ν)}
               \\left(\\frac{h}{a}\\right)^ν K_ν\\left(\\frac{h}{a}\\right)\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the correlation range parameter
    - K_ν is the modified Bessel function of the second kind
    - ν is a fixed parameter determined by the specific implementation
    - Γ(ν) is the gamma function

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Stein, M.L. (1999). Interpolation of Spatial Data. Springer.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Avoid the case where `hole_effect_range` = 0.0
    if hole_effect_range == 0.0:
        return np.where(distance_lags == 0.0, 0.0, nugget)

    # Compute the cyclical term
    cycle = np.where(
        distance_lags / hole_effect_range < 1e-4,
        0.0,
        special.kv(1.0, (distance_lags / hole_effect_range)),
    )

    # Calculate the spatial decay term
    decay = np.where(
        distance_lags / hole_effect_range < 1e-4,
        0.0,
        1.0 - (distance_lags / hole_effect_range) * cycle,
    )

    # Compute the J-Bessel semivariogram
    return partial_sill * decay + nugget


# ---- Linear
def linear(distance_lags: np.ndarray, sill: float, nugget: float):
    """
    Linear variogram model for unbounded spatial processes.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Linear slope parameter (variance increases linearly with distance).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The linear variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\cdot h

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the slope parameter (equivalent to sill in this context)

    The linear model represents unbounded spatial processes where variance increases indefinitely
    with distance. It violates the assumption of a finite sill and should be used cautiously,
    typically for trend analysis or detrended residuals.

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Wackernagel, H. (2003). Multivariate Geostatistics. Springer.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Compute the linear semivariogram
    return partial_sill * distance_lags + nugget


# ---- Matern
def matern(distance_lags, sill, nugget, correlation_range, smoothness_parameter):
    """
    Matérn variogram model with flexible smoothness control.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range parameter controlling spatial correlation decay.
    smoothness_parameter : float
        Smoothness parameter (ν) controlling field differentiability (ν > 0).

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The Matérn variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\frac{2^{1-ν}}{Γ(ν)}
               \\left(\\frac{h}{a}\\right)^ν K_ν\\left(\\frac{h}{a}\\right)\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the correlation range parameter
    - ν is the smoothness parameter
    - K_ν is the modified Bessel function of the second kind
    - Γ(ν) is the gamma function

    **Common smoothness parameter interpretations:**
    - ν = 0.5: Exponential variogram (non-differentiable)
    - ν = 1.5: Once differentiable random fields
    - ν = 2.5: Twice differentiable random fields
    - ν → ∞: Approaches Gaussian variogram (infinitely differentiable)

    References
    ----------
    .. [1] Matérn, B. (1986). Spatial Variation. Springer-Verlag.
    .. [2] Rasmussen, C.E. & Williams, C.K.I. (2006). Gaussian Processes for Machine Learning. MIT
           Press.
    """
    # Calculate the argument for the Bessel function
    h_scaled = distance_lags / correlation_range

    # Handle zero distances
    h_scaled = np.where(h_scaled == 0, 1e-10, h_scaled)

    # Calculate Matérn correlation
    bessel_arg = np.sqrt(2 * smoothness_parameter) * h_scaled
    correlation = (
        (2 ** (1 - smoothness_parameter) / special.gamma(smoothness_parameter))
        * (bessel_arg**smoothness_parameter)
        * special.kv(smoothness_parameter, bessel_arg)
    )

    # Handle numerical issues at origin
    correlation = np.where(distance_lags == 0, 1.0, correlation)

    # Calculate variogram
    return nugget + (sill - nugget) * (1 - correlation)


# ---- Nugget
def nugget(distance_lags: np.ndarray, sill: float, nugget: float):
    """
    Pure nugget variogram model representing uncorrelated spatial noise.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (equivalent to nugget in this model).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The pure nugget variogram model is defined as:

    .. math::
        γ(h) = \\begin{cases}
        0 & \\text{if } h = 0 \\\\
        C_0 & \\text{if } h > 0
        \\end{cases}

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect (total variance)

    This model represents completely uncorrelated spatial data where:
    - No spatial correlation exists beyond measurement locations
    - All variance is attributed to measurement error or micro-scale variation
    - The spatial process appears as white noise

    The pure nugget model is useful for:
    - Modeling measurement error components
    - Representing spatially uncorrelated residuals
    - Baseline comparison for other variogram models
    - Data with extremely short correlation ranges

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Journel, A.G. & Huijbregts, C.J. (1978). Mining Geostatistics. Academic Press.
    """

    # Sum together except at lag == 0.0
    return np.where(distance_lags == 0.0, 0.0, sill + nugget)


# --- Pentaspherical
def pentaspherical(distance_lags: np.ndarray, sill: float, nugget: float, correlation_range: float):
    """
    Pentaspherical variogram model with quintic polynomial behavior and finite range.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range where correlation becomes zero.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The pentaspherical variogram model is defined as:

    .. math::
        γ(h) = \\begin{cases}
        C_0 + C_1 \\left(\\frac{15h}{8a} - \\frac{5h^3}{4a^3} \\right. \\\\
        \\qquad\\qquad \\left. + \\frac{3h^5}{8a^5}\\right) & \\text{if } h ≤ a \\\\
        C_0 + C_1 & \\text{if } h > a
        \\end{cases}

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the effective range

    The pentaspherical model exhibits smoother behavior than the spherical model with continuous
    derivatives up to second order. It provides a good balance between flexibility and
    computational efficiency for bounded spatial processes.

    References
    ----------
    .. [1] Journel, A.G. & Huijbregts, C.J. (1978). Mining Geostatistics. Academic Press.
    .. [2] Wackernagel, H. (2003). Multivariate Geostatistics. Springer.
    """

    # Calculate normalized distances
    h_norm = distance_lags / correlation_range

    # Calculate pentaspherical correlation for h <= a
    correlation = np.where(
        h_norm <= 1, 1 - (15 * h_norm / 8 - 5 * h_norm**3 / 4 + 3 * h_norm**5 / 8), 0
    )

    # Calculate variogram
    return nugget + (sill - nugget) * (1 - correlation)


# ---- Power law
def power(distance_lags: np.ndarray, sill: float, nugget: float, power_exponent: float):
    """
    Power law variogram model for fractal spatial processes.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Not used in power law model (set to 0 or ignore).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    power_exponent : float
        Power exponent parameter (0 < β < 2) controlling scaling behavior.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The power law variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\cdot h^β

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the scaling coefficient
    - β is the power exponent (0 < β < 2)

    **Power exponent interpretations:**
    - β = 1: Linear variogram
    - β → 0: Approaches nugget effect only
    - β → 2: Approaches parabolic behavior (Brownian motion)

    This model represents unbounded, self-similar (fractal) spatial processes
    without a finite sill. It's commonly used for modeling phenomena with
    scale-invariant properties.

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Goovaerts, P. (1997). Geostatistics for Natural Resources Evaluation. Oxford University
           Press.
    """
    # Handle zero distances
    h_nonzero = np.where(distance_lags == 0, 0, distance_lags**power_exponent)

    # Calculate variogram (sill parameter acts as scaling coefficient)
    return nugget + sill * h_nonzero


# ---- (Rational) quadratic
def quadratic(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    shape_parameter: float,
):
    """
    Rational quadratic variogram model with polynomial decay behavior.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range parameter controlling spatial correlation decay.
    shape_parameter : float
        Shape parameter (α) controlling decay rate (α > 0).

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The rational quadratic variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\left(1 + \\frac{h^2}{2αa^2}\\right)^{-α}\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the correlation range parameter
    - α is the shape parameter

    The rational quadratic model exhibits polynomial decay and can be viewed as a scale mixture of
    Gaussian processes. It provides intermediate behavior between exponential and Gaussian models.

    References
    ----------
    .. [1] Rasmussen, C.E. & Williams, C.K.I. (2006). Gaussian Processes for Machine Learning. MIT
    Press.
    .. [2] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    """
    # Calculate scaled distances
    h_scaled = distance_lags / correlation_range

    # Calculate correlation
    correlation = (1 + h_scaled**2 / (2 * shape_parameter)) ** (-shape_parameter)

    # Calculate variogram
    return nugget + (sill - nugget) * (1 - correlation)


# ---- Sinc
def sinc(distance_lags: np.ndarray, sill: float, nugget: float, hole_effect_range: float):
    """
    Sinc (cardinal sine) variogram model exhibiting oscillatory hole-effect behavior.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    hole_effect_range : float
        Range parameter controlling the wavelength of oscillations.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The sinc variogram model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\frac{\\sin(\\pi h/a)}{\\pi h/a}\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the hole-effect range parameter
    - The limit as h→0 gives sinc(0) = 1

    The sinc function creates regular oscillatory patterns with:
    - First zero crossing at h = a
    - Subsequent zeros at integer multiples of a
    - Damped oscillations that decay as 1/h

    This model is suitable for processes with:
    - Regular periodic structures
    - Wave-like spatial patterns
    - Alternating zones of positive and negative correlation
    - Spectral characteristics in spatial data

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Wackernagel, H. (2003). Multivariate Geostatistics. Springer.
    """

    # Get machine epsilon
    eps = np.finfo(float).eps

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.sin(hole_effect_range * (distance_lags + eps)) / (distance_lags + eps)

    # Compute the sinc semivariogram
    return partial_sill * decay + nugget


# ---- Spherical
def spherical(distance_lags: np.ndarray, sill: float, nugget: float, correlation_range: float):
    """
    Spherical variogram model with finite effective range.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Effective range where correlation becomes zero.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The spherical variogram model is defined as:

    .. math::
        γ(h) = \\begin{cases}
        C_0 + C_1 \\left(\\frac{3h}{2a} - \\frac{h^3}{2a^3}\\right) & \\text{if } h \\leq a \\\\
        C_0 + C_1 & \\text{if } h > a
        \\end{cases}

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a is the effective range (correlation = 0 beyond this distance)

    The spherical model has a finite range and linear near-origin behavior,
    making it suitable for processes with clear spatial boundaries.

    References
    ----------
    .. [1] Matheron, G. (1963). Principles of geostatistics. Economic Geology, 58(8), 1246-1266.
    .. [2] Journel, A.G. & Huijbregts, C.J. (1978). Mining Geostatistics. Academic Press.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = (3.0 * distance_lags) / (2.0 * correlation_range) - distance_lags**3.0 / (
        2.0 * correlation_range**3.0
    )

    # Compute the spherical semivariogram
    return np.where(
        distance_lags < correlation_range,
        partial_sill * decay + nugget,
        sill + nugget,
    )


# Composite family models (i.e. hole-effects)
# ---- J-Bessel and Gaussian
def bessel_gaussian(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
):
    """
    Composite Bessel-Gaussian variogram model combining periodic and smooth decay.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Range parameter for Gaussian decay component.
    hole_effect_range : float
        Range parameter controlling Bessel oscillation periodicity.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The Bessel-Gaussian composite model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\frac{h^2}{a_1^2}\\right)
               \\cdot \\frac{J_1(h/a_2)}{h/(2a_2)}\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a₁ is the correlation range (Gaussian component)
    - a₂ is the hole-effect range (Bessel component)
    - J₁ is the first-order Bessel function

    This model combines smooth Gaussian decay with oscillatory behavior, creating damped periodic
    patterns suitable for highly regular processes with cyclical spatial structures.

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** 2))

    # Calculate the hole effect
    hole_effect = special.j0(hole_effect_range * distance_lags)

    # Compute the composite J-Bessel and Gaussian semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- J-Bessel and exponential
def bessel_exponential(
    distance_lags: np.ndarray,
    nugget: float,
    sill: float,
    correlation_range: float,
    decay_power: float,
    hole_effect_range: float,
):
    """
    Composite Bessel-exponential variogram model combining periodic and decay behavior.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Range parameter for exponential decay component.
    hole_effect_range : float
        Range parameter controlling Bessel oscillation periodicity.
    decay_power : float
        Decay exponent for the exponential component (0 < α ≤ 2).

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The Bessel-exponential composite model combines periodic and exponential decay:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\left(\\frac{h}{a_1}\\right)^α\\right)
               \\cdot \\frac{J_1(h/a_2)}{h/(2a_2)}\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a₁ is the correlation range (exponential component)
    - a₂ is the hole-effect range (Bessel component)
    - α is the decay power
    - J₁ is the first-order Bessel function

    This model captures both oscillatory patterns and long-range decay, suitable for processes with
    periodic structures that diminish with distance.

    **Decay power interpretations:**
    - α = 1.0: Linear exponential decay with Bessel oscillations
    - α = 2.0: Gaussian decay with Bessel oscillations
    - α ∈ (0,2): Intermediate decay behaviors

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Matérn, B. (1986). Spatial Variation. Springer-Verlag.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** decay_power))

    # Calculate the hole effect
    hole_effect = special.j0(hole_effect_range * distance_lags)

    # Compute the composite J-Bessel and exponential semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- cosine and exponential
def cosine_exponential(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
    enhance_semivariance: bool,
):
    """
    Composite cosine-exponential variogram model with trigonometric modulation.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Range parameter for exponential decay component.
    hole_effect_range : float
        Range parameter controlling cosine oscillation wavelength.
    decay_power : float
        Decay exponent for the exponential component (0 < α ≤ 2).

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The cosine-exponential composite model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\left(\\frac{h}{a_1}\\right)^α\\right)
               \\cdot \\cos\\left(\\frac{\\pi h}{a_2}\\right)\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a₁ is the correlation range (exponential component)
    - a₂ is the hole-effect range (cosine wavelength)
    - α is the decay power

    This model creates sinusoidal modulation of exponential decay, producing regular oscillatory
    patterns. The cosine component creates predictable periodicity, making it suitable for
    processes with known cyclical behavior.

    **Decay power interpretations:**
    - α = 1.0: Linear exponential decay with cosine modulation
    - α = 2.0: Gaussian decay with cosine modulation
    - α ∈ (0,2): Intermediate decay behaviors

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay_modifier = -1.0 if enhance_semivariance is True else 1.0
    decay = decay_modifier * np.exp(-(distance_lags / correlation_range))

    # Calculate the hole effect
    hole_effect = np.cos(hole_effect_range * distance_lags)

    # Compute the composite cosine and exponential semivariogram
    return partial_sill * (1.0 - decay * hole_effect) + nugget


# ---- cosine and Gaussian
def cosine_gaussian(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
):
    """
    Composite cosine-Gaussian variogram model with trigonometric modulation and smooth decay.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance (nugget + partial sill).
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Range parameter for Gaussian decay component.
    hole_effect_range : float
        Range parameter controlling cosine oscillation wavelength.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The cosine-Gaussian composite model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\frac{h^2}{a_1^2}\\right)
               \\cdot \\cos\\left(\\frac{\\pi h}{a_2}\\right)\\right)

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (sill - nugget)
    - a₁ is the correlation range (Gaussian component)
    - a₂ is the hole-effect range (cosine wavelength)

    This model combines:
    - **Gaussian smoothness**: Infinitely differentiable correlation structure
    - **Cosine periodicity**: Regular oscillatory patterns
    - **Smooth transitions**: Gradual approach to the sill

    The Gaussian component provides very smooth spatial correlation while the
    cosine component introduces predictable periodic behavior. This combination
    is ideal for highly regular spatial processes with known cyclical patterns.

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Stein, M.L. (1999). Interpolation of Spatial Data. Springer.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = np.exp(-((distance_lags / correlation_range) ** 2))

    # Calculate the hole effect
    hole_effect = np.cos(hole_effect_range * distance_lags)

    # Compute the composite cosine and Gaussian semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- exponential and linear
def exponential_linear(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
    decay_power: float,
):
    """
    Composite exponential-linear variogram model with bounded and unbounded components.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance for the exponential component.
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Range parameter for exponential decay component.
    hole_effect_range : float
        Range parameter for linear component.
    decay_power : float
        Decay exponent for the exponential component (0 < α ≤ 2).

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The exponential-linear composite model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\left(\\frac{h}{a_1}\\right)^α\\right)\\right)
               + C_2 \\cdot \\frac{h}{a_2}

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (exponential component)
    - a₁ is the correlation range (exponential component)
    - a₂ is the hole-effect range (linear component)
    - α is the decay power
    - C₂ is determined by model implementation

    **Decay power interpretations:**
    - α = 1.0: Linear exponential decay with linear trend
    - α = 2.0: Gaussian decay with linear trend
    - α ∈ (0,2): Intermediate decay behaviors

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Wackernagel, H. (2003). Multivariate Geostatistics. Springer.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** decay_power))

    # Calculate the hole effect
    hole_effect = 1.0 - hole_effect_range * distance_lags**decay_power

    # Compute the composite exponential and linear semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# ---- Gaussian and linear
def gaussian_linear(
    distance_lags: np.ndarray,
    sill: float,
    nugget: float,
    correlation_range: float,
    hole_effect_range: float,
):
    """
    Composite Gaussian-linear variogram model with smooth bounded and unbounded components.

    Parameters
    ----------
    distance_lags : np.ndarray
        Array of spatial lag distances.
    sill : float
        Total variance for the Gaussian component.
    nugget : float
        Nugget effect representing measurement error and micro-scale variation.
    correlation_range : float
        Range parameter for Gaussian decay component.
    hole_effect_range : float
        Range parameter for linear component.

    Returns
    -------
    np.ndarray
        Variogram values at specified lag distances.

    Notes
    -----
    The Gaussian-linear composite model is defined as:

    .. math::
        γ(h) = C_0 + C_1 \\left(1 - \\exp\\left(-\\frac{h^2}{a_1^2}\\right)\\right)
               + C_2 \\cdot \\frac{h}{a_2}

    where:
    - γ(h) is the variogram at lag distance h
    - C₀ is the nugget effect
    - C₁ is the partial sill (Gaussian component)
    - a₁ is the correlation range (Gaussian component)
    - a₂ is the hole-effect range (linear component)
    - C₂ is determined by model implementation

    This model combines smooth Gaussian decay with linear trend components for highly regular
    spatial processes with regional drift.

    References
    ----------
    .. [1] Chilès, J.P. & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
    .. [2] Stein, M.L. (1999). Interpolation of Spatial Data. Springer.
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(-((distance_lags / correlation_range) ** 2))

    # Calculate the hole effect
    hole_effect = 1.0 - hole_effect_range * distance_lags**2

    # Compute the composite Gaussian and linear semivariogram
    return partial_sill * (decay * hole_effect) + nugget


# Variogram function API
VARIOGRAM_MODELS = {
    "single": {
        "cubic": cubic,
        "exponential": exponential,
        "gaussian": gaussian,
        "jbessel": jbessel,
        "kbessel": kbessel,
        "linear": linear,
        "matern": matern,
        "nugget": nugget,
        "pentaspherical": pentaspherical,
        "power": power,
        "quadratic": quadratic,
        "sinc": sinc,
        "spherical": spherical,
    },
    "composite": {
        ("bessel", "exponential"): bessel_exponential,
        ("bessel", "gaussian"): bessel_gaussian,
        ("cosine", "exponential"): cosine_exponential,
        ("cosine", "gaussian"): cosine_gaussian,
        ("exponential", "linear"): exponential_linear,
        ("gaussian", "linear"): gaussian_linear,
    },
}


# Variogram wrapper function
def compute_variogram(
    distance_lags: np.ndarray,
    variogram_parameters: Optional[Dict[str, float]] = None,
    model: Optional[Union[str, List[str]]] = None,
    **kwargs,
):
    """
    Compute the theoretical semivariogram

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: Optional[Dict[str, float]]
        An optional dictionary that can contain values for variogram model parameters (see the
        below table associated with the argument `model`). Alternatively, these parameters can be
        entered directly and are contained within `kwargs`. Possible parameters include:
            - `sill` (Sill): The asymptotic value as lags approach infinity.
            - `nugget` (Nugget): The semivariogram y-intercept that corresponds to variability
            at lag distances shorter than the lag resolution.
            - `correlation_range` (Correlation length scale/range): The ascending rate for the
            semivariogram.
            - `hole_effect_range` (Hole effect range): The (normalized) length scale/range that
            'holes' are observed, which represent 'null' (or very small) points compared to their
            neighboring lags.
            - `decay_power` (Decay term exponent): An exponential term that is used in certain
            generalized exponential (or related) semivariogram models that modulates the ascending
            rate for a semivariogram.
            - `enhance_semivariance` (Semivariance enhancement): A boolean term that determines
            whether the correlation decay in certain cosine-related variogram models are enhanced
            (or not) are further lag distances.
    model: Optional[Union[ str , list ]]
        A string or list of model names. A single name represents a single family model. Two inputs
        represent the desired composite model (e.g. the composite J-Bessel and exponential model).
        Available variogram models and their respective arguments include (alongside
        `distance_lags`):

        +----------------------------+-----------------+--------------------------+
        | :fun:`variogram`           | Input           | Parameters               |
        | model                      |                 |                          |
        +============================+=================+==========================+
        |  :fun:`cubic`              | 'cubic'         | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`exponential`        | 'exponential    | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`gaussian`           | 'gaussian'      | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`jbessel`            | 'jbessel'       | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`kbessel`            | 'kbessel'       | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`linear`             | 'linear'        | - `nugget`               |
        |                            |                 | - `sill`                 |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`matern`             | 'matern'        | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `smoothness_parameter` |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`nugget`             | 'nugget'        | - `nugget`               |
        |                            |                 | - `sill`                 |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`pentaspherical`     | 'pentaspherical'| - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`power`              | 'power'         | - `nugget`               |
        |                            |                 | - `sill`                 |
        |                            |                 | - `power_exponent`       |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`quadratic`          | 'quadratic'     | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `shape_parameter`      |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`sinc`               | 'sinc'          | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`spherical`          | 'spherical'     | - `sill`                 |
        |                            |                 | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`bessel_exponential` | ['bessel',      | - `sill`                 |
        |                            |  'exponential'] | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `decay_power`          |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`bessel_gaussian`    | ['bessel',      | - `sill`                 |
        |                            |  'gaussian']    | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `decay_power`          |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`cosine_exponential` | ['cosine',      | - `sill`                 |
        |                            |  'exponential'] | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        |                            |                 | - `enhance_semivariance` |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`cosine_gaussian`    | ['cosine',      | - `sill`                 |
        |                            |  'gaussian']    | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`exponential_linear` | ['exponential', | - `sill`                 |
        |                            |  'linear']      | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        |                            |                 | - `decay_power`          |
        +----------------------------+-----------------+--------------------------+
        |  :fun:`gaussian_linear`    | ['gaussian',    | - `sill`                 |
        |                            |  'linear']      | - `nugget`               |
        |                            |                 | - `correlation_range`    |
        |                            |                 | - `hole_effect_range`    |
        +----------------------------+-----------------+--------------------------+

    Returns
    ----------
    variogram: np.ndarray
        An array containing the (normalized) semivariance for each lag bin.
    """

    # Determine model source
    if model is not None:
        # ---- Get the variogram arguments and function from `model`
        model_source = model
    elif variogram_parameters is not None:
        # ---- Get the variogram arguments and function from `variogram_parameters`
        model_source = variogram_parameters["model"]
    else:
        raise ValueError("Argument `model` is missing.")

    # Extract function signatures
    variogram_args, variogram_function = get_variogram_arguments(model_source)

    # Evaluate whether required function parameters are present
    if variogram_parameters is not None:
        # ---- Get input arguments
        input_args = variogram_parameters
        # ---- Use `variogram_parameters` as source
        arg_diff = set(list(variogram_args)) - set(input_args) - set(["distance_lags"])
    else:
        # ---- Get input arguments
        input_args = kwargs
        # ---- Use `kwargs` as source
        arg_diff = set(list(variogram_args)) - set(input_args) - set(["distance_lags"])

    # Raise error if any are missing
    if len(arg_diff) > 0:
        raise ValueError(
            f"The following variogram parameters are missing: {', '.join(list(arg_diff))}"
        )

    # Filter out only the variogram parameters required for the model
    required_args = dict((k, input_args[k]) for k in input_args if k in list(variogram_args))

    # Pipe the parameters into the appropriate variogram function
    return variogram_function["model_function"](distance_lags, **required_args)


def get_variogram_arguments(model_name: Union[str, List[str]]):
    """
    Get the variogram function arguments and model function for a given model.

    Parameters
    ----------
    model_name : Union[str, List[str]]
        A string or list of model names. A single name represents a single family model.
        Two inputs represent the desired composite model (e.g. the composite J-Bessel and
        exponential model).

    Returns
    -------
    Tuple[inspect.Signature.parameters, Dict[str, Any]]
        A tuple containing the function signature parameters and a dictionary with the
        model function.
    """

    # Convert to lowercase to match reference model dictionary
    if isinstance(model_name, str):
        model_input = model_name.lower()
    elif isinstance(model_name, list) & len(model_name) == 1:
        model_input = "".join(model_name).lower()
    else:
        model_input = [name.lower() for name in model_name]
        # ---- Alphabetic sort
        model_input.sort()

    # Parse user input from reference model dictionary
    # ---- Check against VARIOGRAM_MODELS API to ensure model exists
    if isinstance(model_input, list) and (tuple(model_input) in VARIOGRAM_MODELS["composite"]):
        # if (len(model_input) > 1) & (tuple(model_input) in VARIOGRAM_MODELS["composite"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["composite"][tuple(model_input)]
    # elif (len([model_input]) == 1) & (model_input in VARIOGRAM_MODELS["single"]):
    elif not isinstance(model_input, list) and (model_input in VARIOGRAM_MODELS["single"]):
        # ---- Parse model function
        model_function = VARIOGRAM_MODELS["single"][model_input]
    else:
        raise LookupError(
            f"The model input ({model_name}) could not be matched to an"
            f" existing variogram method."
        )

    # Check input parameters against the required function arguments
    # ---- Get the function signature
    function_signature = inspect.signature(model_function)
    # ---- Create ordered dictionary of required arguments
    return function_signature.parameters, {"model_function": model_function}


def fit_variogram(
    lags: np.ndarray[float],
    lag_counts: np.ndarray[int],
    gamma: np.ndarray[float],
    model_parameters: Parameters,
    model: Union[str, List[str]] = ["bessel", "exponential"],
    optimizer_kwargs: Dict[str, Any] = {},
) -> Tuple[Dict[str, Any], float, float]:
    """
    Fit theoretical variogram models to empirical semivariogram data using weighted least squares.

    This function performs non-linear optimization to find the best-fitting parameters for
    theoretical variogram models. The optimization uses weighted least squares where weights
    are proportional to the lag counts, giving more influence to lags with more data pairs.

    Parameters
    ----------
    lags : np.ndarray
        Array of lag distances from empirical variogram computation.
    lag_counts : np.ndarray
        Number of data point pairs contributing to each lag estimate. Used as optimization weights.
    gamma : np.ndarray
        Empirical semivariogram values (standardized semivariance) at each lag.
    model : str or List[str]
        Theoretical variogram model specification. Single string for basic models
        ('exponential', 'gaussian', 'spherical', 'jbessel', 'linear'). List of two strings
        for composite models (['bessel', 'exponential'], ['bessel', 'gaussian'],
        ['cosine', 'exponential']).
    model_parameters : lmfit.Parameters
        Parameter object containing initial values, bounds, and constraints for optimization.
        Required parameters depend on the selected model.
    optimizer_kwargs : dict, default={}
        Additional keyword arguments passed to `lmfit.minimize()`. Common options include
        'max_nfev' for maximum function evaluations and solver-specific parameters.

    Returns
    -------
    Tuple[dict, dict, float]
        - Optimized parameter values as dictionary
        - Initial parameter values as dictionary
        - Mean absolute deviation of optimized fit

    Notes
    -----
    The optimization minimizes the weighted objective function:

    .. math::
        \\min_{θ} \\sum_{i=1}^{n} w_i [γ_{emp}(h_i) - γ_{model}(h_i; θ)]^2

    where w_i = lag_counts[i] are the weights and θ represents the model parameters.

    **Available Models:**

    *Single Models:*
    - 'cubic': Smooth cubic polynomial with finite range
    - 'exponential': C₀ + C₁(1 - exp(-h/a))
    - 'gaussian': C₀ + C₁(1 - exp(-h²/a²))
    - 'jbessel': C₀ + C₁(1 - J₁(h/a)/(h/2a)) - exhibits hole effects
    - 'kbessel': K-Bessel function model
    - 'linear': C₀ + C₁·h - unbounded growth
    - 'nugget': Pure nugget effect (no spatial correlation)
    - 'pentaspherical': Quintic polynomial with smooth transitions
    - 'power': Power law model for fractal processes
    - 'quadratic': Rational quadratic with polynomial decay
    - 'sinc': Sinc function with oscillatory behavior
    - 'spherical': Piecewise function with finite range

    *Composite Models:*
    - ['bessel', 'exponential']: Periodic patterns with exponential decay
    - ['bessel', 'gaussian']: Periodic patterns with Gaussian smoothness
    - ['cosine', 'exponential']: Sinusoidal modulation with exponential decay

    The function uses Trust Region Reflective algorithm (default in lmfit) which handles parameter
    bounds robustly and is suitable for the non-linear nature of variogram models.

    References
    ----------
    .. [1] Cressie, N. (1993). Statistics for Spatial Data. Wiley.
    .. [2] Newville, M., et al. (2014). LMFIT: Non-Linear Least-Square Minimization and
           Curve-Fitting for Python. Zenodo.
    """
    # Normalize the lag counts to get the lag weights
    lag_weights = lag_counts / lag_counts.sum()

    # Vertically stack the lags, semivariance, and weights
    data_stack = np.vstack((lags, gamma, lag_weights))

    # Recover the lag resolution
    delta_lag = np.diff(lags).mean()

    # Compute the range
    range = lags.max() + delta_lag

    # Index lag distances that are within the parameterized range
    within_range = np.where(lags <= range)[0]

    # Truncate the data stack
    truncated_stack = data_stack[:, within_range]

    # Create helper cost-function that is weighted using the kriging weights (`w`), lag
    # distances (`x`), and empirical semivariance (`y`)
    def cost_function(parameters, x, y, w, model):
        yr = compute_variogram(x, {**parameters, **{"model": model}})
        return (yr - y) * w

    # Compute the initial fit based on the pre-optimized parameter values
    initial_fit = cost_function(
        model_parameters,
        x=truncated_stack[0],
        y=truncated_stack[1],
        w=truncated_stack[2],
        model=model,
    )

    # Compute the initial mean absolute deviation (MAD)
    mad_initial = np.mean(np.abs(initial_fit))

    # Generate `Minimizer` function class required for bounded optimization
    minimizer = Minimizer(
        cost_function,
        model_parameters,
        fcn_args=(truncated_stack[0], truncated_stack[1], truncated_stack[2], model),
    )

    # Minimize the cost-function to compute the best-fit/optimized variogram parameters
    parameters_optimized = minimizer.minimize(method="least_squares", **optimizer_kwargs)

    # Calculate the optimized MAD
    mad_optimized = np.mean(np.abs(parameters_optimized.residual))

    # Extract the best-fit parameter values
    best_fit_params = parameters_optimized.params.valuesdict()

    # Return the final tuple
    return best_fit_params, mad_initial, mad_optimized
