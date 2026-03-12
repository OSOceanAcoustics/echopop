"""
Survey analysis and methodology functions.

This package contains functions for acoustic survey design, analysis, and population estimation
including stratified sampling, bootstrap statistics, biological proportions, and population
apportionment.
"""

from .apportionment import (
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
    quantize_length_data,
)
from .proportions import (
    binned_weights,
    compute_binned_counts,
    fitted_weight_proportions,
    get_nasc_proportions_slice,
    get_number_proportions_slice,
    get_weight_proportions_slice,
    number_proportions,
    stratum_averaged_weight,
    weight_proportions,
)
from .stratified import JollyHampton
from .transect import compute_interval_distance

__all__ = [
    # Class [Stratified]
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
    "quantize_length_data",
    # Proportion functions
    "binned_weights",
    "compute_binned_counts",
    "fitted_weight_proportions",
    "get_nasc_proportions_slice",
    "get_number_proportions_slice",
    "get_weight_proportions_slice",
    "number_proportions",
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
