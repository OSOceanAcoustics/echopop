from functools import lru_cache
from typing import Literal, Tuple, Union

import numpy as np
import pandas as pd


def impute_missing_sigma_bs(
    unique_strata: Union[list, np.ndarray], sigma_bs_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Imputes sigma_bs for strata without measurements or values

    Parameters
    ----------
    unique_strata : Union[list, np.ndarray]
        An array comprising all expected stratum numbers.
    sigma_bs_df : pd.DataFrame
        DataFrame containing the mean `sigma_bs` calculated for each stratum.
        Must have stratum as index and 'sigma_bs' column.

    Returns
    -------
    pd.DataFrame
        DataFrame with imputed sigma_bs values for missing strata

    Examples
    --------
    >>> # Create sigma_bs data for some strata
    >>> sigma_bs_df = pd.DataFrame({
    ...     'sigma_bs': [0.001, 0.003, 0.002]
    ... }, index=pd.Index([1, 3, 5], name='stratum'))
    >>>
    >>> # Define all expected strata
    >>> all_strata = [1, 2, 3, 4, 5]
    >>>
    >>> # Impute missing strata (2 and 4)
    >>> result = impute_missing_sigma_bs(all_strata, sigma_bs_df)
    >>> print(result)
           sigma_bs
    stratum
    1       0.001000
    2       0.002000  # interpolated between 1 and 3
    3       0.003000
    4       0.002500  # interpolated between 3 and 5
    5       0.002000

    Notes
    -----
    This function iterates through all stratum layers to impute either the
    nearest neighbor interpolation or mean sigma_bs for strata that are missing values.

    For missing strata, the function:
    1. Finds the nearest stratum below and above the missing stratum
    2. Interpolates the sigma_bs value as the mean of these neighbors
    3. If no neighbors exist on one side, uses the available neighbor value
    """

    # Extract the stratum index name
    stratum_name = sigma_bs_df.index.name

    # Collect present strata
    present_strata = np.unique(sigma_bs_df.index)

    # Invert to retrieve missing strata
    missing_strata = set(unique_strata).difference(set(present_strata))

    # Impute values for missing strata
    if len(missing_strata) > 0:

        # Concatenate the existing data with the missing strata
        sigma_bs_stratum_impute = pd.concat(
            [
                sigma_bs_df,
                pd.DataFrame({stratum_name: list(missing_strata), "sigma_bs": np.nan}).set_index(
                    stratum_name
                ),
            ],
        ).sort_values(stratum_name)

        # Find strata intervals to impute over
        for i in missing_strata:

            # Get the stratum indices below and above the missing stratum
            strata_floor = present_strata[present_strata < i]
            strata_ceil = present_strata[present_strata > i]

            new_stratum_below = (
                np.max(strata_floor) if strata_floor.size > 0 else np.min(strata_ceil)
            )
            new_stratum_above = (
                np.min(strata_ceil) if strata_ceil.size > 0 else np.max(strata_floor)
            )

            # Get the indexed values
            sigma_bs_indexed = sigma_bs_stratum_impute.loc[[new_stratum_below, new_stratum_above]]

            # Impute the mean to the missing value
            sigma_bs_stratum_impute.loc[i] = sigma_bs_indexed.mean()

        # Return the imputed values
        return sigma_bs_stratum_impute
    else:
        return sigma_bs_df


def wavenumber(
    frequency: Union[np.ndarray, float],
    sound_speed_sw: float,
) -> np.ndarray[float]:
    """
    Compute the acoustic wavenumber from frequency and sound speed.

    The acoustic wavenumber relates frequency and wavelength through the
    medium's sound speed, fundamental to acoustic scattering calculations.

    Parameters
    ----------
    frequency : Union[np.ndarray, float]
        Acoustic frequency in Hz. Can be scalar or array for multiple frequencies.
    sound_speed_sw : float
        Sound speed in seawater in m/s, typically ~1500 m/s depending on
        temperature, salinity, and pressure conditions.

    Returns
    -------
    np.ndarray[float]
        Acoustic wavenumber in rad/m. Same shape as input frequency.

    Examples
    --------
    >>> # Single frequency
    >>> k = wavenumber(120e3, 1500.0)  # 120 kHz in 1500 m/s water
    >>> print(f"k = {k:.2f} rad/m")
    k = 502.65 rad/m

    >>> # Multiple frequencies
    >>> freqs = np.array([18e3, 38e3, 120e3])
    >>> k_array = wavenumber(freqs, 1500.0)

    Notes
    -----
    The wavenumber is calculated as:

    .. math::
        k = \\frac{2\\pi f}{c}

    where k is wavenumber (rad/m), f is frequency (Hz), and c is sound speed (m/s).

    References
    ----------
    .. [1] Medwin, H. & Clay, C.S. (1998). Fundamentals of Acoustical
           Oceanography. Academic Press.
    """

    return 2 * np.pi * frequency / sound_speed_sw


def reflection_coefficient(
    g: Union[np.ndarray, float],
    h: Union[np.ndarray, float],
) -> np.ndarray[float]:
    """
    Compute the acoustic reflection coefficient from material properties.

    The reflection coefficient quantifies acoustic impedance mismatch between
    organism tissue and surrounding seawater, crucial for scattering strength
    calculations in biological acoustic models.

    Parameters
    ----------
    g : Union[np.ndarray, float]
        Relative density contrast: ratio of organism density to seawater density.
        Typical values for zooplankton: 1.01-1.10 (slightly denser than water).
    h : Union[np.ndarray, float]
        Relative sound speed contrast: ratio of organism sound speed to
        seawater sound speed. Typical values: 0.99-1.05.

    Returns
    -------
    np.ndarray[float]
        Acoustic reflection coefficient (dimensionless). Same shape as inputs.

    Examples
    --------
    >>> # Typical krill parameters
    >>> g = 1.04  # 4% denser than seawater
    >>> h = 1.02  # 2% faster sound speed than seawater
    >>> R = reflection_coefficient(g, h)
    >>> print(f"Reflection coefficient: {R:.4f}")

    >>> # Array inputs for parameter studies
    >>> g_vals = np.linspace(1.01, 1.08, 5)
    >>> h_vals = np.linspace(1.00, 1.04, 5)
    >>> R_matrix = reflection_coefficient(g_vals[:, None], h_vals)

    Notes
    -----
    The reflection coefficient is calculated using the formula:

    .. math::
        R = \\frac{1 - gh^2}{gh^2} - \\frac{g-1}{g}

    This formulation accounts for both density (g) and compressibility (h)
    contrasts between the scatterer and surrounding medium.
    """

    return (1 - g * h * h) / (g * h * h) - (g - 1) / g


def orientation_average(
    angle: np.ndarray[float],
    form_function: np.ndarray[complex],
    theta_mean: float,
    theta_sd: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """
    Compute orientation-averaged backscattering cross-section.

    This function integrates the complex form function over organism
    orientation angles weighted by a probability distribution, accounting
    for natural variation in organism orientation in the water column.

    Parameters
    ----------
    angle : np.ndarray[float]
        Array of orientation angles in degrees, typically relative to vertical
    form_function : np.ndarray[complex]
        Complex-valued form function f(ka, θ) for each angle and frequency.
        Should have shape compatible with angle array.
    theta_mean : float
        Mean orientation angle in degrees (e.g., 90° for horizontal orientation)
    theta_sd : float
        Standard deviation of orientation distribution in degrees
    distribution : Literal["gaussian", "uniform"], default="gaussian"
        Type of orientation probability distribution:
        - "gaussian": Normal distribution around theta_mean
        - "uniform": Uniform distribution over angle range

    Returns
    -------
    np.ndarray[float]
        Orientation-averaged backscattering cross-section σ_bs in m²

    Examples
    --------
    >>> angles = np.linspace(60, 120, 61)  # ±30° around vertical
    >>> # form_function computed from scattering model
    >>> sigma_bs = orientation_average(angles, form_func, 90.0, 10.0)

    Notes
    -----
    The orientation averaging is performed as:

    .. math::
        \\sigma_{bs} = \\int |f(\\theta)|^2 P(\\theta) d\\theta

    where f(θ) is the complex form function and P(θ) is the orientation
    probability density function (PDF).

    For Gaussian distributions:

    .. math::
        P(\\theta) =
        \\frac{1}{\\sqrt{2\\pi}\\sigma} \\exp\\left(-\\frac{(\\theta-\\mu)^2}{2\\sigma^2}\\right)

    """

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        # ---- Get interval
        orientation_interval = np.diff(angle).mean()
        # ---- Compute the PDF
        PDF = (
            orientation_interval
            * np.exp(-0.5 * (angle - theta_mean) ** 2 / theta_sd**2)
            / (np.sqrt(2 * np.pi) * theta_sd)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(angle)) / len(angle)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Return the weighted form function
    # ---- If complex
    return [
        (
            np.sqrt(np.matmul((f[0].real ** 2 + f[0].imag ** 2), PDF))
            if np.all(np.iscomplex(f))
            else np.sqrt(np.matmul(f, PDF))
        )
        for f in form_function
    ]


def valid_array_row_length(arr: np.ndarray[float]) -> int:
    """
    Returns the number of valid (i.e. not NaN) length of each row within an array
    """

    return np.sum(~np.isnan(arr))


def length_average(
    length_values: np.ndarray[float],
    ka_f: np.ndarray[float],
    ka_c: np.ndarray[float],
    form_function: np.ndarray[complex],
    length_mean: float,
    length_deviation: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """ """

    # Normalize the length values, if needed
    length_norm = length_values / length_mean
    # ---- Also normalize the standard deviation
    length_sd_norm = length_deviation / length_mean

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        # ---- Get the interval
        length_interval = np.diff(length_norm).mean()
        # ---- Compute the PDF
        PDF = (
            length_interval
            * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm**2)
            / (np.sqrt(2 * np.pi) * length_sd_norm)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(length_values)) / len(length_values)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Length-weighted averaged sigma_bs
    # ---- Get valid values
    if ka_f.ndim == 1:
        n_vals = valid_array_row_length(ka_f)
    else:
        n_vals = np.apply_along_axis(valid_array_row_length, 1, arr=ka_f)
    # ---- Compute the length-weighted ka
    ka_weighted = length_norm * ka_c.reshape(-1, 1)
    # ---- Trim values so they fall within the valid/defined bandwidth
    ka_weighted_trim = np.where(
        (ka_weighted >= np.nanmin(ka_f)) & (ka_weighted <= np.nanmax(ka_f)), ka_weighted, np.nan
    )
    # ---- Evaluate
    # -------- One frequency
    if len(form_function) == 1:
        sigma_bs_L = np.array(
            [
                (
                    length_norm**2
                    * PDF
                    * np.interp(ka_weighted_trim[0], ka_f[0], form_function[0] ** 2)
                ).sum()
            ]
        )
    # -------- More than one frequency
    else:
        sigma_bs_L = np.array(
            [
                (
                    length_norm**2
                    * PDF
                    * np.interp(ka_weighted_trim[i], ka_f[i, : n_vals[i]], form_function[i] ** 2)
                ).sum()
                for i in range(len(form_function))
            ]
        )
    # ---- Return the weighted average
    return sigma_bs_L


def _make_freq_key(
    frequency: np.ndarray, length_sd_norm: float, frequency_interval: float, /, ndigits: int = 12
) -> Tuple:
    cf = tuple(round(float(x), ndigits) for x in np.asarray(frequency).ravel())
    return (cf, round(float(length_sd_norm), ndigits), round(float(frequency_interval), ndigits))


@lru_cache(maxsize=128)
def _generate_frequency_interval_cached_key(key: Tuple) -> np.ndarray:
    # key is (cf_tuple, length_sd_norm_rounded, frequency_interval_rounded)
    cf_tuple, length_sd_norm_rounded, frequency_interval_rounded = key
    # rebuild numpy array from tuple
    cf = np.asarray(cf_tuple, dtype=float)
    # perform the original generation
    frequency_lst = [
        np.arange(
            freq * (1 - 3.1 * float(length_sd_norm_rounded)),
            freq * (1 + 3.1 * float(length_sd_norm_rounded)) + float(frequency_interval_rounded),
            float(frequency_interval_rounded),
        )
        for freq in cf
    ]
    # Find the maximum length and pad with NaN to return a regular 2D ndarray
    max_length = max(len(arr) for arr in frequency_lst)
    padded_results = np.full((len(frequency_lst), max_length), np.nan, dtype=float)
    for i, arr in enumerate(frequency_lst):
        padded_results[i, : len(arr)] = arr

    return padded_results


def generate_frequency_interval(
    frequency: np.ndarray, length_sd_norm: float, frequency_interval: float, ndigits: int = 12
) -> np.ndarray:
    """
    Wrapper that builds a stable hashable key and uses an lru_cache-backed worker.
    No module-level globals are created; cache is managed by functools.lru_cache.
    """
    key = _make_freq_key(frequency, length_sd_norm, frequency_interval, ndigits)
    return _generate_frequency_interval_cached_key(key)


# OPTIMIZATION: More efficient parameter extraction
def _extract_parameters_optimized(inverted_data: pd.DataFrame) -> pd.DataFrame:
    """
    Extract parameters more efficiently than using .apply().

    This optimized version avoids the overhead of pandas .apply()
    by using direct iteration and batch DataFrame construction.

    Parameters
    ----------
    inverted_data : pd.DataFrame
        DataFrame with 'parameters' column containing InvParameters objects

    Returns
    -------
    pd.DataFrame
        DataFrame with parameter values as columns
    """
    # Extract all parameter dictionaries at once
    param_dicts = [obj.values for obj in inverted_data["parameters"]]

    # Construct DataFrame in one operation
    return pd.DataFrame(param_dicts, index=inverted_data.index)
