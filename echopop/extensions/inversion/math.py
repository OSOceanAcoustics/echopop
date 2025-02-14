"""
Mathematical and numerical utility functions.
"""

from typing import Any, Dict, Literal, Tuple, Union
import pandas as pd
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

    return (1 - g * h * h) / (g * h * h) - (g - 1) / g


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
    # if len(form_function) == 1:
    #     # ---- If complex
    #     if np.all(np.iscomplex(form_function)):
    #         return np.sqrt(form_function.real**2 + form_function.imag**2)
    #     # ---- If not complex
    #     else:
    #         return form_function

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
            length_interval
            * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm**2)
            / (np.sqrt(2 * np.pi) * length_sd_norm)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(form_function)) / len(form_function)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Length-weighted averaged sigma_bs
    # ---- Get valid values
    n_vals = np.apply_along_axis(valid_array_row_length, 1, arr=ka)
    # ---- Compute the length-weighted ka
    ka_weighted = length_norm * ka_center.reshape(-1, 1)
    # ---- Trim values so they fall within the valid/defined bandwidth
    ka_weighted_trim = np.where(
        (ka_weighted >= np.nanmin(ka)) & (ka_weighted <= np.nanmax(ka)), ka_weighted, np.nan
    )
    # ---- Evaluate
    # -------- One frequency
    if len(form_function) == 1:
        sigma_bs_L = np.array(
            [
                (
                    length_norm**2 
                    * PDF 
                    * np.interp(ka_weighted_trim[0], ka[0], form_function[0] ** 2)
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
                    * np.interp(ka_weighted_trim[i], ka[i, : n_vals[i]], form_function[i] ** 2)
                ).sum()
                for i in range(len(form_function))
            ]
        )
    # ---- Return the weighted average
    return sigma_bs_L

def inverse_normalize_series(series: pd.Series, ranges_dict: Dict[str, Any]):
    """
    Apply inverse min-max normalization to a `pandas.Series`
    """

    # Create copy
    series_copy = series.copy() 
    
    # Get the columns that will be inverse transformed
    transform_cols = list(ranges_dict.keys())
    
    # Convert `ranges_dict` to a `pd.DataFrame`
    ranges_df = pd.DataFrame(ranges_dict).T

    # Apply inverse transformation
    series_copy[transform_cols] = (
        series_copy[transform_cols] * (ranges_df["high"] - ranges_df["low"]) + ranges_df["low"]
    )

    # Return the inverse min-max normalized series
    return series_copy

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

    # # Skip weighted averaging if only a single value is present
    # if len(form_function) == 1:
    #     # ---- If complex
    #     if np.all(np.iscomplex(form_function)):
    #         return [np.sqrt(f.real**2 + f.imag**2) for f in form_function[0]]
    #     # ---- If not complex
    #     else:
    #         return form_function

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
        PDF = np.ones(len(form_function)) / len(form_function)
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

def transect_mask(coordinates_df: pd.DataFrame,
                  mask: pd.Series,
                  filter_dict: Dict[str, Any]) -> pd.Series:
    """
    Generate transect mask
    """
    # Transect-based filtering
    if "transect_min" in filter_dict:
        mask &= coordinates_df.index.get_level_values(0) >= filter_dict["transect_min"]
    if "transect_max" in filter_dict:
        mask &= coordinates_df.index.get_level_values(0) <= filter_dict["transect_max"]
    if "transect" in filter_dict:
        mask &= coordinates_df.index.get_level_values(0).isin(filter_dict["transect"])

    return mask

def coordinate_mask(coordinates_df: pd.DataFrame,
                    mask: pd.Series,
                    filter_dict: Dict[str, Any]) -> pd.Series:
    """
    Generate coordinate mask
    """
    # Longitude-based filtering
    if "longitude_min" in filter_dict:
        mask &= coordinates_df["longitude"] >= filter_dict["longitude_min"]
    if "longitude_max" in filter_dict:
        mask &= coordinates_df["longitude"] <= filter_dict["longitude_max"]

    # Latitude-based filtering
    if "latitude_min" in filter_dict:
        mask &= coordinates_df["latitude"] >= filter_dict["latitude_min"]
    if "latitude_max" in filter_dict:
        mask &= coordinates_df["latitude"] <= filter_dict["latitude_max"]

    return mask

def create_composite_mask(coordinates_df: pd.DataFrame, 
                          filter_dict: Dict[str, Any]) -> pd.Series:
    """
    Create composite transect and coordinate mask
    """

    # Start with all True mask
    mask = pd.Series(True, index=coordinates_df.index)

    # Transect-based filtering
    mask &= transect_mask(coordinates_df, mask, filter_dict)

    # Coordinate-based filtering
    mask &= coordinate_mask(coordinates_df, mask, filter_dict)

    # Return the composite mask
    return mask    

def reduce_dataset(metadata_dfs: Dict[str, pd.DataFrame], 
                   measurements_df: pd.DataFrame,
                   filter_dict: Dict[str, Any]) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Mask the dataset to only include specified transects, longitudes, and latitudes
    """
    # Get coordinates
    coordinates_df = metadata_dfs["coordinates_df"]

    # Create composite mask
    df_mask = create_composite_mask(coordinates_df, filter_dict)

    # Apply mask to metadata attributes
    updated_metadata = {name: df.loc[df_mask] for name, df in metadata_dfs.items()}

    # Apply mask to dataset
    # ---- For MultiIndex
    if isinstance(measurements_df.index, pd.MultiIndex):
        updated_measurements = measurements_df[df_mask.reindex_like(measurements_df)]
    # ---- For single index
    else:
        # ---- Initialize transect mask
        mask_transect = pd.Series(True, index=measurements_df.index)
        # ---- Transect-based filtering
        mask_transect &= transect_mask(measurements_df, mask_transect, filter_dict)
        # ---- Apply mask
        updated_measurements = measurements_df[mask_transect.reindex_like(measurements_df)]

    # Return both
    return updated_metadata, updated_measurements

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


def generate_frequency_interval(
    frequency: np.ndarray[float], length_sd_norm: float, frequency_interval: float
) -> np.ndarray[float]:
    """
    Generate frequency interval 2D array centered on an array input of center frequencies.
    """

    frequency_lst = [
        np.arange(
            freq * (1 - 3.1 * length_sd_norm),
            freq * (1 + 3.1 * length_sd_norm) + frequency_interval,
            frequency_interval,
        )
        for freq in frequency
    ]

    # Find the maximum length of the generated arrays
    max_length = max(len(arr) for arr in frequency_lst)

    # Create a padded 2D array with NaN for shorter arrays
    padded_results = np.full((len(frequency_lst), max_length), np.nan)
    for i, arr in enumerate(frequency_lst):
        padded_results[i, : len(arr)] = arr

    return padded_results
