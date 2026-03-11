"""
NWFSC FEAT-specific utilities and workflow code.

This sub-package contains modules, functions, and classes used exclusively for the
Northwest Fisheries Science Center (NWFSC) Fisheries Ecosystem Assessment Team (FEAT) workflow.
These are not intended for general use outside of FEAT applications.
"""

from ...survey.biology import compute_abundance, length_binned_weights

from ...survey.apportionment import distribute_unaged_from_aged, mesh_biomass_to_nasc, reallocate_excluded_estimates, remove_group_from_estimates, sum_population_tables

from . import feat_parameters
from ...survey.apportionment import (
    distribute_population_estimates,
)
from ...survey.biology import compute_biomass
from .comparisons import (
    compute_dataset_differences,
    load_all_geodata_reports,
    plot_dataset_differences,
)
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
    # Comparison functions
    "compute_dataset_differences",
    "load_all_geodata_reports",
    "plot_dataset_differences",
    # Report generation class,
    "Reporter",
    # Submodules
    "apportionment",
    "biology",
    "comparison",
    "functions",
    "feat_parameters",
    "reporter",
]
