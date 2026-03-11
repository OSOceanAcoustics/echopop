"""
Utility functions for Echopop.
"""

from .base import (
    apply_filters,
    binify,
    binned_distribution,
    create_grouped_table,
    create_pivot_table,
    group_interpolator_creator,
    round_half_up
)

from .feat_functions import (
    get_survey_western_extents,
    transect_ends_crop,
    western_boundary_search_strategy
)

from . import feat_parameters

__all__ = [
    # Generic utilities
    "apply_filters",
    "binify",
    "binned_distribution",
    "create_grouped_table",
    "create_pivot_table",
    "group_interpolator_creator",
    "round_half_up",
    # FEAT-specific functions
    "get_survey_western_extents",
    "transect_ends_crop",
    "western_boundary_search_strategy",

    # Submodules
    "base",
    "feat_functions",
    "feat_parameters"
]