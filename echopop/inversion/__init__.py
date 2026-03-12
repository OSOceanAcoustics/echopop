"""
Acoustic inversion subpackage for estimating fish population characteristics from backscatter.

Exposes the core inversion classes and helpers: ``InversionBase``, ``InversionLengthTS``,
``InversionMatrix``, ``InvParameters``, and the ``ts_length_regression`` utility.
"""

from .inversion_base import InversionBase, InvParameters
from .inversion_length_TS import InversionLengthTS, ts_length_regression
from .inversion_matrix import InversionMatrix
from .post_inversion import estimate_population
from .utils import (
    generate_frequency_interval,
    impute_missing_sigma_bs,
    length_average,
    orientation_average,
    reflection_coefficient,
    wavenumber,
)

__all__ = [
    "InvParameters",
    "InversionBase",
    "InversionLengthTS",
    "InversionMatrix",
    "estimate_population",
    "generate_frequency_interval",
    "impute_missing_sigma_bs",
    "length_average",
    "orientation_average",
    "reflection_coefficient",
    "wavenumber",
    "ts_length_regression",
]
