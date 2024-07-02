import ipywidgets as widgets
import matplotlib.pyplot as plt
from typing import Dict
from IPython.display import display, clear_output

import pandas as pd
import numpy as np

from ..spatial.variogram import (
    get_variogram_arguments, empirical_variogram
)

def click_theoretical_variogram_button(button, variogram_settings, plot_output_widget):

    # Plot/update the plot widget
    with plot_output_widget:
        # Clear 
        clear_output(wait=True)
        # Figure and axis objects
        fig, ax = plot_empirical_variogram()
        # Add the theoretical variogram line to the plot
        plot_theoretical_variogram(fig, ax, variogram_settings)
        # Edit canvas settings
        fig.canvas.toolbar_visible = False
        fig.canvas.header_visible = False
        fig.canvas.footer_visible = False
        fig.canvas.resizable = False
        # Show the plot
        plt.show(fig)


def click_empirical_variogram_button(button, transect_data, general_settings, settings_dict):

    # Compute the distance lags
    distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
    # ---- Compute the maximum range
    max_range = general_settings["n_lags"] * general_settings["lag_resolution"]

    # Update the general settings for computing the empirical variogram
    general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags, 
                             "range": max_range})
    
    # Compute the empirical variogram
    lags, gamma_h, lag_counts, _ = empirical_variogram(
        transect_data, general_settings, settings_dict
    )

    # Return the outputs
    return lags, gamma_h, lag_counts, general_settings