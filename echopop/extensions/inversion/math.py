"""
Mathematical and numerical utility functions.
"""

from typing import Any, Dict, Literal, Tuple, Union

import numpy as np

from scipy.special import spherical_jn, spherical_yn


def spherical_hn(n, z, derivative=False) -> np.ndarray:
    """
    Spherical Bessel function of the third kind (Hankel function) or its derivative

    Defined as [1]_,

    .. math:: h_n^{(1)}(z)=j_n(z)+in_n(z),

    where :math:`h_n^{(1)}` is the spherical Bessel function of the third kind (or Hankel function
    of the first kind), :math:`j_n` is the spherical Bessel function of the first kind, :math:`n_n`
    is the spherical Bessel function of the second kind (or Neumann function), :math:`n` is the
    order of the function (:math:`n>=0`), :math:`z` is the Bessel function argument value, and
    :math:`i` is an imaginary number.

    Parameters
    ----------
    n: int
        Order of the Bessel function (n >= 0)
    z: Union[float, np.complex]
        Argument of the Bessel function
    derivative: Optional[bool]
        When True, the derivative is computed

    Notes
    -----
    The derivative is computed using the relations [2]_,

    .. math::
        \frac{n}{z} h^{(1)}_n - h^{(1)}_{n+1}(z)

    References
    ----------
    .. [1] https://dlmf.nist.gov/10.47#E5
    .. [2] https://dlmf.nist.gov/10.51#E2

    """
    # == lib/sphhn.m

    # Define internal function
    def _spherical_hn(n, z):
        return spherical_jn(n, z) + 1j * spherical_yn(n, z)

    # Computing derivative
    if derivative:
        return (n / z) * _spherical_hn(n, z) - _spherical_hn(n + 1, z)
    else:
        return _spherical_hn(n, z)

def wavenumber(
    frequency: Union[np.ndarray, float],
    water_sound_speed: float,
) -> np.ndarray[float]:
    """
    Compute the acoustic wavenumber
    """
    
    return 2 * np.pi * frequency / water_sound_speed 

def reflection_coefficient(
    g: Union[np.ndarray, float],
    h: Union[np.ndarray, float],
) -> np.ndarray[float]: 
    """
    Compute the reflection coefficient based on material properties
    """

    return (1-g*h*h)/(g*h*h)-(g-1)/g

def length_average(
    length: np.ndarray[float],
    ka: np.ndarray[float],
    ka_center: np.ndarray[float],
    form_function: np.ndarray[complex],
    length_mean: float,
    length_deviation: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """
    Compute the length-averaged linear backscattering cross-section (:math:`\sigma_{bs}(L)`)
    """

    # Skip weighted averaging if only a single value is present
    if len(form_function) == 1:
        # ---- If complex
        if np.all(np.iscomplex(form_function)):
            return np.sqrt(form_function.real ** 2 + form_function.imag ** 2)  
        # ---- If not complex
        else:
            return form_function            

    # Normalize the length values, if needed
    length_norm = length / length_mean
    # ---- Also normalize the standard deviation
    length_sd_norm = length_deviation / length_mean

    # Weight based on distribution input
    # ---- Gaussian (Normal)    
    if distribution == "gaussian":
        # ---- Get the interval
        length_interval = np.diff(length_norm).mean()
        # ---- Compute the PDF
        PDF = (
            length_interval * np.exp(-0.5*(length_norm-1)**2 / length_sd_norm ** 2) 
            / (np.sqrt(2*np.pi)*length_sd_norm)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(form_function))/len(form_function)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Compute `sigma_bs` (orientation-averaged backscattering cross-section)
    sigma_bs = form_function ** 2

    # Length-weighted averaged sigma_bs
    # ---- Get valid values
    n_vals = np.apply_along_axis(valid_array_row_length, 1, arr=ka)
    # ---- Compute the length-weighted ka 
    ka_weighted = length_norm * ka_center.reshape(-1, 1)
    # ---- Trim values so they fall within the valid/defined bandwidth
    ka_weighted_trim = np.where((ka_weighted >= np.nanmin(ka)) & (ka_weighted <= np.nanmax(ka)), 
                                ka_weighted, 
                                np.nan)
    # ---- Evluate
    sigma_bs_L = np.array([(length_norm ** 2 * PDF 
                  * np.interp(ka_weighted_trim[i], ka[i, :n_vals[i]], form_function[i] ** 2)).sum() 
                  for i in range(len(form_function))])
    # ---- Return the weighted average
    return sigma_bs_L

def orientation_average(
    angle: np.ndarray[float],
    form_function: np.ndarray[complex],
    theta_mean: float,
    theta_sd: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """
    Compute the orientation-averaged linear backscattering cross-section :math:`\sigma_{bs}(\theta)`
    """

    # Skip weighted averaging if only a single value is present
    if len(form_function) == 1:
        # ---- If complex
        if np.all(np.iscomplex(form_function)):
            return np.sqrt(form_function.real ** 2 + form_function.imag ** 2)  
        # ---- If not complex
        else:
            return form_function  

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        # ---- Get interval
        orientation_interval = np.diff(angle).mean()
        # ---- Compute the PDF
        PDF = (
            orientation_interval * np.exp(-0.5*(angle-theta_mean)**2 / theta_sd ** 2) 
            / (np.sqrt(2*np.pi)*theta_sd)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(form_function))/len(form_function)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Return the weighted form function
    # ---- If complex
    return [np.sqrt(np.matmul((f[0].real ** 2 + f[0].imag ** 2), PDF)) 
            if np.all(np.iscomplex(f)) else np.sqrt(np.matmul(f, PDF)) 
            for f in form_function]

def fit_rayleigh_pdf(
    measured: np.ndarray[float],
    density: np.ndarray[float],
    mean: float,
    standard_deviation: float,
    lower_bounds: float,
    upper_bounds: float,
    arg_distribution: Literal["exponential", "gaussian"] = "gaussian",
):
    """
    Fit a single-parameter Rayleigh probability density function to the measured data
    """
    pass


def valid_array_row_length(arr: np.ndarray[float]) -> int:
    """
    Returns the number of valid (i.e. not NaN) length of each row within an array
    """
    
    return np.sum(~np.isnan(arr))