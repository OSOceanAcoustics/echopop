from functools import lru_cache
from typing import Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd


def impute_missing_sigma_bs(
    unique_strata: Union[list, np.ndarray[np.number]], sigma_bs_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Imputes :math:`\\sigma_\\text{bs}` for strata without measurements or values

    Parameters
    ----------
    unique_strata : Union[list, |np.ndarray[np.number]|]
        An array comprising all expected stratum numbers.
    sigma_bs_df : |pd.DataFrame|
        DataFrame containing the mean :math:`\\sigma_\\text{bs}` calculated for each stratum. This
        DataFrame must be indexed by the stratum and must contain the column ``'sigma_bs'``.

    Returns
    -------
    |pd.DataFrame|
        DataFrame with imputed :math:`\\sigma_\\text{bs}` values for missing strata

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
    This function iterates through all stratum layers to impute either the nearest neighbor
    interpolation or mean :math:`\\sigma_\\text{bs}` for strata that are missing values.

    For missing strata, the function:

    1. Finds the nearest stratum below and above the missing stratum
    2. Interpolates the :math:`\\sigma_\\text{bs}` value as the mean of these neighbors
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
    frequency: Union[np.ndarray[float], float],
    sound_speed_sw: float,
) -> np.ndarray[float]:
    """
    Compute the acoustic wavenumber from frequency and sound speed.

    The acoustic wavenumber relates frequency and wavelength through the medium's sound speed,
    fundamental to acoustic scattering calculations.

    Parameters
    ----------
    frequency : Union[|np.ndarray[float]|, float]
        Acoustic frequency in Hz. Can be scalar or array for multiple frequencies.
    sound_speed_sw : float
        Sound speed in seawater in m s :math:`^{-1}`, typically ~1500 m s :math:`^{-1}` depending
        on temperature, salinity, and pressure conditions.

    Returns
    -------
    |np.ndarray[float]|
        Acoustic wavenumber in rad m :math:`^{-1}`. Same shape as input frequency.

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

    where :math:`k` is wavenumber (rad m :math:`^{-1}`), :math:`f` is frequency (Hz), and :math:`c`
    is sound speed (m s :math:`^{-1}`).

    References
    ----------
    .. [1] Medwin, H. & Clay, C.S. (1998). Fundamentals of Acoustical Oceanography. Academic Press.
    """

    return 2 * np.pi * frequency / sound_speed_sw


def reflection_coefficient(
    g: Union[np.ndarray[float], float],
    h: Union[np.ndarray[float], float],
) -> np.ndarray[float]:
    """
    Compute the acoustic reflection coefficient from material properties.

    The reflection coefficient quantifies acoustic impedance mismatch between organism tissue and s
    urrounding seawater, crucial for scattering strength calculations in biological acoustic models.

    Parameters
    ----------
    g : Union[|np.ndarray[float]|, float]
        The density contrast, which is the ratio of organism density to seawater density
        (kg m :math:`^{-3}`). For fluid-like, weak scatterers, values for zooplankton are often
        close to unity (1.00), typically 1.00-1.05.
    h : Union[|np.ndarray[float]|, float]
        The sound speed contrast, which is the ratio of organism sound speed to seawater sound
        speed (m s :math:`^{-1}`). For fluid-like, weak scatterers, values for zooplankton are
        often close to unity (1.00), typically 1.00-1.05.

    Returns
    -------
    |np.ndarray[float]|
        Acoustic reflection coefficient (:math:`\\mathscr{R}_{12}`, dimensionless) with the same
        shape as the inputs.

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
    Let the material properties of the surrounding seawater be denoted as :math:`(\\rho_1, c_1)`
    and those of the scatterer :math:`(\\rho_2, c_2)`. The density contrast :math:`g` and sound
    speed contrast :math:`h` are defined as:

    .. math::
        g = \\frac{\\rho_2}{\\rho_1}, \\quad h = \\frac{c_2}{c_1}

    These are then used to compute the reflection coefficient :math:`\\mathscr{R}_{12}`:

    .. math::
        \\mathscr{R}_{12} = \\frac{1 - gh^2}{gh^2} - \\frac{g-1}{g}

    where :math:`\\mathscr{R}_{12}` is dimensionless. This quantity constitutes a simplification of
    the boundary conditions at the interface between the surrounding seawater and scatterer. It
    encapsulates how acoustic waves are partially reflected an transmitted due to this impedance
    mismatch.
    """

    # Convert to arrays for consistent operations
    g_arr = np.asarray(g)
    h_arr = np.asarray(h)

    # Validate shapes for broadcasting
    if g_arr.shape != h_arr.shape:
        if g_arr.size == 1 or h_arr.size == 1:
            pass
        else:
            raise ValueError(
                "'g' and 'h' must be scalars, or vectors of the same shape, or one scalar and one "
                "vector."
            )

    return (1 - g_arr * h_arr * h_arr) / (g_arr * h_arr * h_arr) - (g_arr - 1) / g_arr


def orientation_average(
    theta: np.ndarray[float],
    form_function: Union[np.ndarray[complex], np.ndarray[float]],
    theta_mean: Optional[float] = None,
    theta_sd: Optional[float] = None,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
    output_type: Literal["sigma_bs", "f_bs"] = "sigma_bs",
    convert_type: bool = True,
) -> list[np.ndarray[float]]:
    """
    Compute the orientation-averaged form function (:math:`\\mathscr{f}(f, \\theta)`) for
    complex linear scattering coefficient :math:`f_\\text{bs}` or backscattering cross-section
    :math:`\\sigma_\\text{bs}`.

    This function integrates the form function, :math:`\\mathscr{f}(f, \\theta)`, over the
    scatterer's orientation distribution. This effectively weights :math:`\\mathscr{f}(f, \\theta)`
    using a probability density function (PDF) to account for natural variation in scatterer
    orientation in the water column. The result is the orientation-averaged form function
    :math:`\\bar{\\mathscr{f}}(f)`.

    Parameters
    ----------
    theta : |np.ndarray[float]|
        Array of orientation angles (:math:`\\theta`, °) in degrees, typically relative to the
        incident acoustic wave where broadside incidence is considered to be 90°.
    form_function : |np.ndarray[complex]| or |np.ndarray[float]|
        The scattering form function that is either the complex linear scattering coefficient
        (:math:`f_\\text{bs}(f, \\theta)`) or backscattering cross-section
        (:math:`\\sigma_\\text{bs}(f, \\theta)`).
        This must be an a complex or float array with the same shape as the array for ``theta``.
    theta_mean : Optional[float]
        Mean orientation angle (:math:`\\bar{\\theta}`, °) in degrees.
    theta_sd : Optional[float]
        Standard deviation of orientation values (:math:`\\sigma_{\\theta}`, °).
    distribution : Literal["gaussian", "uniform"], default="gaussian"
        Orientation probability distribution function (PDF) used for averaging. There are two
        available options:

        - ``"gaussian"``: Normal/Gaussian distribution around :math:`\\bar{\\theta}`
          (``theta_mean``). This is the default distribution, and requires both ``theta_mean`` and
          ``theta_sd`` to be specified.
        - ``"uniform"``: Uniform distribution over all orientation values in ``theta``.

    output_type : Literal["sigma_bs", "f_bs"], default="sigma_bs"
        Specifies the desired output type after length averaging:

        - ``"sigma_bs"``: Outputs the length-averaged backscattering cross-section
          (:math:`\\sigma_\\text{bs}(f)`, m²). When ``form_function`` is complex, it is converted to
          :math:`\\sigma_\\text{bs}` before averaging before averaging when ``convert_type`` is
          True. This is the default option.
        - ``"f_bs"``: Outputs the length-averaged linear scattering coefficient
          (:math:`f_\\text{bs}(f)`, m).

    convert_type : bool, default=True
        If ``True`` and ``output_type`` is set to ``"sigma_bs"``, the function converts the input
        ``form_function`` from complex linear scattering coefficient (:math:`f_\\text{bs}`) to
        backscattering cross-section (:math:`\\sigma_\\text{bs}`) before averaging. If ``False``,
        the function assumes the input ``form_function`` is already in the desired output type.

    Returns
    -------
    List[|np.ndarray[float]|]
        Orientation-averaged backscattering cross-section :math:`\\sigma_\\text{bs}` (m²). When
        ``form_function`` represents the linear scattering coefficient (i.e., a complex array), the
        output remains :math:`f_\\text{bs}` after averaging.

    Examples
    --------
    >>> angles = np.linspace(60, 120, 61)  # ±30° around vertical
    >>> # form_function computed from scattering model
    >>> sigma_bs = orientation_average(angles, form_func, 90.0, 10.0)

    Notes
    -----
    The orientation averaging when :math:`\\mathscr{f}(f, \\theta)` is real is defined as:

    .. math::
        \\bar{\\mathscr{f}}(f) = \\int \\mathscr{f}(f, \\theta) P(\\theta) d\\theta

    where :math:`\\mathscr{f}(f, \\theta)` is the form function at frequency :math:`f` and
    orientation :math:`\\theta`, and :math:`P(\\theta)` is the probability density function (PDF).
    Since :math:`\\mathscr{f}(f, \\theta)` can also be the linear scattering coefficient, the
    averaged form function calculation is modified to account for this, defined as:

    .. math::
        \\bar{\\mathscr{f}}(f) = \\sqrt{
            \\int (|\\mathscr{f}(f, \\theta)|^2) P(\\theta) d\\theta
        }

    which corresponds to the real-valuted root mean square (RMS). When using a Gaussian
    distribution, the PDF is defined as:

    .. math::
        P(\\theta) = \\frac{1}{\\sqrt{2\\pi}\\sigma_{\\theta}}
        \\exp\\left(-\\frac{(\\theta-\\bar{\\theta})^2}{2\\sigma_{\\theta}^2}\\right)

    When using a uniform distribution, the PDF is constant over the range of :math:`\\theta`:

    .. math::
        P(\\theta) = \\frac{1}{n_\\theta}

    where :math:`n_\\theta` is the number of discrete orientations in :math:`\\theta`.

    For computational purposes, this integral is approximated using a discrete sum over
    :math:`\\theta`:

    .. math::
        \\begin{align*}
            \\bar{\\mathscr{f}}(f) &\\approx
                \\sum\\limits_{\\theta} \\mathscr{f}(f, \\theta) P(\\theta) \\Delta\\theta
                \\quad \\text{when }
                \\bar{\\mathscr{f}}(f) \\rightarrow \\bar{\\mathscr{f}}(f) \\\\
            \\bar{\\mathscr{f}}(f) &\\approx \\sqrt{
                \\sum\\limits_{\\theta} |\\mathscr{f}(f, \\theta)|^2 P(\\theta) \\Delta\\theta
            }
            \\quad \\text{when }
            \\mathscr{f}(f, L) \\in \\mathbb{C} \\text{ and }
            \\bar{\\mathscr{f}}(f) \\rightarrow \\bar{\\mathscr{f}}(f) \\\\
            \\bar{\\mathscr{f}}(f) &\\approx
                \\sum\\limits_{\\theta} |\\mathscr{f}(f, \\theta)|^2 P(\\theta) \\Delta\\theta
                \\quad \\text{when }
                \\bar{f}_\\text{bs}(f) \\rightarrow \\bar{\\sigma}_\\text{bs}(f) \\\\
        \\end{align*}

    where :math:`\\Delta\\theta` is the interval between discrete orientation angles in
    :math:`\\theta`. In the case of the bottom equation, the squared magnitude is used when
    converting the mean linear scattering coefficient :math:`\\bar{f}_\\text{bs}(f)` to the mean
    backscattering cross-section :math:`\\bar{\\sigma}_\\text{bs}(f)`.
    """

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        if theta_mean is None or theta_sd is None:
            raise ValueError(
                "Both 'theta_mean' and 'theta_sd' must be specified for Gaussian distribution."
            )
        # ---- Get interval
        orientation_interval = np.diff(theta).mean()
        # ---- Compute the PDF
        PDF = (
            orientation_interval
            * np.exp(-0.5 * (theta - theta_mean) ** 2 / theta_sd**2)
            / (np.sqrt(2 * np.pi) * theta_sd)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(theta)) / len(theta)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Return the weighted form function
    result = []
    for f in form_function:
        # ---- Output: f_bs
        if output_type == "f_bs":
            val = np.sqrt(np.matmul(np.abs(f[0]) ** 2, PDF))
        # ---- Output: sigma_bs
        else:
            if convert_type:
                # ---- Apply conversion
                val = np.matmul(np.abs(f[0]) ** 2, PDF)
            else:
                val = np.matmul(f[0], PDF)
        result.append(val)
    return result


def valid_array_row_length(arr: np.ndarray[float]) -> int:
    """
    Returns the number of valid (i.e. not NaN) length of each row within an array
    """

    return np.sum(~np.isnan(arr))


def length_average(
    length_values: np.ndarray[float],
    ka_f: np.ndarray[float],
    ka_c: np.ndarray[float],
    form_function: Union[np.ndarray[complex], np.ndarray[float]],
    length_mean: float,
    length_deviation: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
    output_type: Literal["sigma_bs", "f_bs"] = "sigma_bs",
    convert_type: bool = True,
) -> list[np.ndarray[float]]:
    """
    Compute length-averaged backscattering cross-section (:math:`\\sigma_\\text{bs}(f, L)`) or
    linear scattinering coefficient :math:`f_\\text{bs}` (m) or backscattering cross-section
    :math:`\\sigma_\\text{bs}` (m²).

    This function integrates the form function, :math:`\\mathscr{f}(f, L)`, over the scatterer's
    length distribution. This effectively weights :math:`\\mathscr{f}(f, L)` using a probability
    density function (PDF) to account for natural variation in scatterer length in the water
    column. The result is the length-averaged form function :math:`\\bar{\\mathscr{f}}(f)`.

    Parameters
    ----------
    length_values : |np.ndarray[float]|
        Array of scatterer lengths (:math:`L`, m).
    ka_f : |np.ndarray[float]|
        Array of dimensionless wavenumber-length products for each frequency (ka).
    ka_c : |np.ndarray[float]|
        Array of reference ka values for each length, typically computed from central frequency.
    form_function : |np.ndarray[float]| or |np.ndarray[float]|
        The scattering form function, either the complex linear scattering coefficient
        (:math:`f_\\text{bs}(f, L)`) or backscattering cross-section
        (:math:`\\sigma_\\text{bs}(f, L)`).
    length_mean : float
        Mean scatterer length (:math:`\\bar{L}`, m).
    length_deviation : float
        Standard deviation of scatterer lengths (:math:`\\sigma_L`, m).
    distribution : Literal["gaussian", "uniform"], default="gaussian"
        Length probability distribution function (PDF) used for averaging. There are two available
        options:

        - ``"gaussian"``: Normal/Gaussian distribution around :math:`\\bar{L}`
          (``length_mean``). This is the default distribution, and requires both ``length_mean`` and
          ``length_sd`` to be specified.
        - ``"uniform"``: Uniform distribution over all orientation values in ``length_values``.

    output_type : Literal["sigma_bs", "f_bs"], default="sigma_bs"
        Specifies the desired output type after length averaging:

        - ``"sigma_bs"``: Outputs the length-averaged backscattering cross-section
          (:math:`\\sigma_\\text{bs}(f)`, m²). When ``form_function`` is complex, it is converted to
          :math:`\\sigma_\\text{bs}` before averaging when ``convert_type`` is True. This is the
          default option.
        - ``"f_bs"``: Outputs the length-averaged linear scattering coefficient
          (:math:`f_\\text{bs}(f)`, m).

    convert_type : bool, default=True
        If ``True`` and ``output_type`` is set to ``"sigma_bs"``, the function converts the input
        ``form_function`` from complex linear scattering coefficient (:math:`f_\\text{bs}`) to
        backscattering cross-section (:math:`\\sigma_\\text{bs}`) before averaging. If ``False``,
        the function assumes the input ``form_function`` is already in the desired output type.

    Returns
    -------
    List[|np.ndarray[float]|]
        Length-averaged backscattering cross-section :math:`\\sigma_\\text{bs}` (m²). When
        ``form_function`` is complex, it is converted to :math:`\\sigma_\\text{bs}` before
        averaging. The output quantity yields a single length-averaged :math:`\\sigma_\\text{bs}`
        value for each frequency.

    Examples
    --------
    >>> lengths = np.linspace(0.01, 0.03, 10)  # 1 to 3 cm
    >>> sigma_bs_L = length_average(lengths, ka_f, ka_c, form_func, 0.02, 0.005)

    Notes
    -----
    The length averaging when :math:`\\mathscr{f}(f, L)` is real is defined as:

    .. math::
        \\bar{\\mathscr{f}}^*(f) = \\int \\mathscr{f}(f, L) P(L^*) dL^*

    where :math:`\\mathscr{f}(f, L)` is the form function at frequency :math:`f` and normalized
    length :math:`L^*`, and :math:`P(L^*)` is the normalized probability density function (PDF)
    such that:

    .. math::
        P(L^*) = \\frac{P(L)}{\\bar{L}}

    Since :math:`\\mathscr{f}(f, L)` can also be complex, the form function can be modified to
    calculate the squared magnitude:

    .. math::
        \\bar{\\mathscr{f}}^*(f) = \\int (|\\mathscr{f}(f, L)|^2) P(L^*) dL^*

    For Gaussian distributions, the PDF is defined as:

    .. math::
        P(L^*) = \\frac{1}{\\sqrt{2\\pi}\\sigma_{L^*}}
        \\exp\\left(-\\frac{(L^*-\\bar{L^*})^2}{2\\sigma_{L^*}^2}\\right)

    For uniform distributions, the PDF is constant over the range of :math:`L^*`:

    .. math::
        P(L^*) = \\frac{1}{n_{L^*}}

    where :math:`n_{L^*}` is the number of discrete lengths in :math:`L^*`.

    For computational purposes, this integral is approximated using a discrete sum over :math:`L^*`:

    .. math::
        \\begin{align*}
            \\bar{\\mathscr{f}}^*(f) &\\approx \\sum_{L^*} |\\mathscr{f}(f, L)| P(L^*) \\Delta L^*
            \\quad \\text{when }
            \\bar{\\mathscr{f}}^*(f) \\rightarrow \\bar{\\mathscr{f}}^*(f) \\\\
            \\bar{\\mathscr{f}}^*(f) &\\approx \\sum_{L^*} |\\mathscr{f}(f, L)|^2 P(L^*) \\Delta L^*
            \\quad \\text{when }
            \\bar{f}_\\text{bs}(f) \\rightarrow \\bar{\\sigma}_\\text{bs}^*(f)
        \\end{align*}

    where :math:`\\Delta L^*` is the interval between discrete length values in :math:`L^*`. In the
    case of the bottom equation, the squared magnitude is used when converting the mean linear
    scattering coefficient :math:`\\bar{f}_\\text{bs}(f)` to the mean backscattering
    cross-section :math:`\\bar{\\sigma}_\\text{bs}^*(f)`. Because
    :math:`\\bar{\\sigma}_\\text{bs}^*(f)` is normalized to length, it is in units of m, not m².
    Consequently, an additional scaling is required to convert to the true mean backscattering
    cross-section :math:`\\bar{\\sigma}_\\text{bs}(f)` in m²:

    .. math::
        \\bar{\\sigma}_\\text{bs}(f) = \\bar{\\sigma}_\\text{bs}^*(f) \\times \\bar{L}^2
    """

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

    # Helper function to copmute the weighted averages
    def _length_weighted_sum(form, ka_wt, ka_freq):
        # ---- Interpolate to get the form function at the weighted ka values
        form_interp = np.interp(ka_wt, ka_freq, form)
        # ---- Define the weights
        weights = length_norm**2 * PDF
        # ---- Compute the weighted sum
        if output_type == "f_bs":
            return (weights * form_interp).sum()
        else:
            if convert_type:
                val = (weights * np.abs(form_interp) ** 2).sum()
                val *= length_mean**2
                return val
            else:
                return (weights * form_interp).sum()

    # Evaluate
    result = []
    for i in range(len(form_function)):
        form_func = form_function[i]
        # ---- Get weighted sum
        form_wt = _length_weighted_sum(form_func, ka_weighted_trim[i], ka_f[i, : n_vals[i]])
        result.append(form_wt)
    # ---- Return the weighted average
    return result


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
