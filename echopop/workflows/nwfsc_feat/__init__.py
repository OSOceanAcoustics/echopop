"""
NWFSC FEAT-specific utilities and workflow code.

This sub-package contains modules, functions, and classes used exclusively for the
Northwest Fisheries Science Center (NWFSC) Fisheries Ecosystem Assessment Team (FEAT) workflow.
These are not intended for general use outside of FEAT applications.
"""

from . import parameters
from .apportionment import (
    distribute_population_estimates,
    distribute_unaged_from_aged,
    mesh_biomass_to_nasc,
    reallocate_excluded_estimates,
    remove_group_from_estimates,
    sum_population_tables,
)
from .biology import compute_abundance, compute_biomass, length_binned_weights
from .functions import (
    get_survey_western_extents,
    transect_ends_crop,
    western_boundary_search_strategy,
)
from .reporter import Reporter
from .workflows.year_specific import cli_utils

__all__ = [
    # Generic functions
    "cli_utils",
    "get_survey_western_extents",
    "transect_ends_crop",
    "western_boundary_search_strategy",
    # Apportionment functions
    "distribute_population_estimates",
    "distribute_unaged_from_aged",
    "mesh_biomass_to_nasc",
    "reallocate_excluded_estimates",
    "remove_group_from_estimates",
    "sum_population_tables",
    # Biology functions
    "compute_abundance",
    "compute_biomass",
    "length_binned_weights",
    # Report generation class,
    "Reporter",
    # Submodules
    "apportionment",
    "biology",
    "functions",
    "parameters",
    "reporter",
]
