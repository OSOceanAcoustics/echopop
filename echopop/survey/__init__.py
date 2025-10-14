"""
Survey analysis and methodology functions.

This package contains functions for acoustic survey design, analysis, and population estimation
including stratified sampling, bootstrap statistics, biological proportions, and population
apportionment.
"""

from .apportionment import (
    distribute_population_estimates,
    distribute_unaged_from_aged,
    mesh_biomass_to_nasc,
    reallocate_excluded_estimates,
    remove_group_from_estimates,
    sum_population_tables,
)
from .biology import (
    compute_abundance,
    compute_biomass,
    fit_length_weight_regression,
    length_binned_weights,
)
from .proportions import (
    binned_weights,
    compute_binned_counts,
    get_nasc_proportions_slice,
    get_number_proportions_slice,
    get_weight_proportions_slice,
    number_proportions,
    scale_weight_proportions,
    scale_weights_by_stratum,
    stratum_averaged_weight,
    weight_proportions,
)
from .stratified import JollyHampton
from .transect import compute_interval_distance

__all__ = [
    # Main class
    "JollyHampton",
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
    "fit_length_weight_regression",
    "length_binned_weights",
    # Proportion functions
    "binned_weights",
    "compute_binned_counts",
    "get_nasc_proportions_slice",
    "get_number_proportions_slice",
    "get_weight_proportions_slice",
    "number_proportions",
    "scale_weight_proportions",
    "scale_weights_by_stratum",
    "stratum_averaged_weight",
    "weight_proportions",
    # Statistics functions
    "confidence_interval",
    # Transect functions
    "compute_interval_distance",
    # Submodules
    "apportionment",
    "biology",
    "proportions",
    "statistics",
    "stratified",
    "transect",
]
