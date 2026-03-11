"""
Report generation and comparisons.

This package contains functions for generating FEAT-specific reports and visualizing differences 
between EchoPro (MATLAB) and Echopop (Python) outputs.
"""

from .compare import (
    compute_dataset_differences,
    load_all_geodata_reports,
    plot_dataset_differences,
    plot_geodata,
    plot_haul_count_comparisons,
    plot_population_table_comparisons,
    read_aged_geodata,
    read_geodata,
    read_pivot_table_report,
)
from .reporter import Reporter

__all__ = [
    # Report generation class,
    "Reporter",
    # Comparison functions
    "compute_dataset_differences",
    "load_all_geodata_reports",
    "plot_dataset_differences",
    "plot_geodata",
    "plot_haul_count_comparisons",
    "plot_population_table_comparisons",
    "read_aged_geodata",
    "read_geodata",
    "read_pivot_table_report", 
    # Submodules
    "compare"
]