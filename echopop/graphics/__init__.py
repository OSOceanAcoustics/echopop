"""
Visualization sub-package for echopop.

Collects plotting functions for survey transects, kriged mesh outputs, age-length distributions,
variogram diagnostics, and interactive GUI widgets.
"""

from .age_length_heatmap import plot_age_length_heatmap
from .diagnostics import Diagnostics
from .kriged_mesh import plot_kriged_mesh
from .transect_map import plot_transect_map
from .variogram_gui import VariogramGUI

__all__ = [
    "plot_age_length_heatmap",
    "plot_kriged_mesh",
    "plot_transect_map",
    "Diagnostics",
    "VariogramGUI",
]
