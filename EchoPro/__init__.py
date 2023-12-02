from .computation.semivariogram import SemiVariogram
from .survey import Survey

__all__ = ["Survey", "SemiVariogram"]

from _echopro_version import version as __version__  # noqa