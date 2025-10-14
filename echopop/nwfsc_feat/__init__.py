"""
NWFSC FEAT-specific utilities and workflow code.

This sub-package contains modules, functions, and classes used exclusively for the
Northwest Fisheries Science Center (NWFSC) Fisheries Ecosystem Assessment Team (FEAT) workflow.
These are not intended for general use outside of FEAT applications.
"""

from . import parameters
from .functions import (
    get_survey_western_extents, 
    transect_ends_crop, 
    western_boundary_search_strategy
)
from .reporter import Reporter

__all__ = [
    # Generic functions
    "get_survey_western_extents",
    "transect_ends_crop",
    "western_boundary_search_strategy",
    # Report generation class,
    "Reporter",
    # Submodules
    "functions",
    "parameters",
    "reporter",
]
