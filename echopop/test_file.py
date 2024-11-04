from echopop import Survey
import inspect
from typing import Any, Dict, Optional, Tuple, Union
import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union
import numpy as np
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import verde as vd
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.ticker import FixedLocator
from echopop.graphics.validate_plot import PlotModel
from echopop.graphics import plotting as egp, variogram_interactive as egv
from echopop.graphics.plotting import (
    get_survey_bounds, 
    get_coastline, 
    apply_aspect_ratio, 
    prune_args, 
    add_alongtransect_data, 
    add_colorbar, 
    add_transect_lines,
    format_axes,
    format_heatmap_ticks,
    add_heatmap_grid,
)

init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
survey = Survey(init_config, file_config)
survey.load_survey_data()
survey.load_acoustic_data()
survey.transect_analysis()
survey.fit_variogram()
survey.kriging_analysis(variogram_parameters={"n_lags": 30}, variable="biomass_density")
survey.plot(kind="mesh", 
            variable="biomass", 
            plot_parameters={"log_base": 2,  
                             "vmax": 1e7},
            plot_type="pcolormesh")


self = survey
kind: Literal["age_length_distribution", "mesh", "transect"] = "mesh"
variable = "biomass"
plot_parameters: Dict[str, Any] = {"log_base": 10, 
                                   "vmin": 1e3, 
                                   "vmax": 1e8,
                                   }
plot_type: Optional[Literal["heatmap", "hexbin", "scatter", "pcolormesh"]] = "hexbin"

# Get associated plotting function information
plot_info = egp.PLOT_MAP(self, kind)

# Initialize 'parameters' dictionary
parameters = plot_parameters.copy()

# Proceed with plotting
# ---- Type: spatial
if plot_info["type"] == "spatial":
    # ---- Get the geospatial configuration
    geo_config = self.config["geospatial"]
    # ---- Get the coastline, if it exists, from `plot_parameters`
    geo_input = plot_parameters.get("geo_config", None)
    # ---- Advance to see if 'coastline' exists
    if geo_input is not None and "coastline" not in geo_input:
        # ---- Get the default coastline
        egp.get_coastline(geo_config)
    # ---- Update the parameterization
    parameters.update(dict(geo_config=geo_config))

# Prepare the data range for validation testing
parameters["data_range"] = (
    plot_parameters.get("vmin", None),
    plot_parameters.get("vmax", None),
)

# Add the primary arguments into the dictionary
parameters.update(dict(kind=kind, plot_type=plot_type, variable=variable))

# Prepare plotting parameters
validated_parameters = egp.validate_plot_args(**parameters)

############
dataset = plot_info.get("data")
dataset["kriged_variance"]
z = dataset["sample_cv"]
10 ** np.round(np.log10(z.max()), 1)
kwargs = validated_parameters

# Get the dataset variable name
variable = kwargs.get("variable", None)

# Get sex
sex = kwargs.get("sex", None)

# Get the correct dataset for plotting
# ---- Abundance
if variable == "abundance":
    dataset = data_dict["abundance"]["aged_abundance_df"].copy()
    cmap = kwargs.get("cmap", "viridis")
    colorbar_label = kwargs.get("colorbar_label", 
                                (f"Abundance ({sex} fish, # animals)").capitalize())
# ---- Biomass
elif variable == "biomass":
    dataset = data_dict["biomass"]["aged_biomass_df"].copy()
    cmap = kwargs.get("cmap", "plasma")
    colorbar_label = kwargs.get("colorbar_label", 
                                (f"Biomass ({sex} fish, kg)").capitalize())

# Prepare the dataset
# ---- Stack sum
dataset_stk = dataset.stack(future_stack=True).sum(axis=1).reset_index(name=variable)
# ---- Convert the length and age columns from intervals into numerics
# -------- Length
dataset_stk["length"] = dataset_stk["length_bin"].apply(lambda x: x.mid).astype(float)
# -------- Age
dataset_stk["age"] = dataset_stk["age_bin"].apply(lambda x: x.mid).astype(float)
# ---- Subset the data based on the sex
dataset_sub = dataset_stk.loc[dataset_stk["sex"] == sex, :]

# Convert into a heatmap
heatmap_data = dataset_sub.pivot_table(
    values=variable, columns="age", index="length", aggfunc="sum"
)

# Get limits
# ---- vmin
vmin = kwargs.get("vmin", 0.0)
# ---- vmax
vmax = kwargs.get("vmax", 10 ** np.round(np.log10(heatmap_data.stack().max())))

# Get axis limits
if not kwargs.get("axis_limits", None):
    # ---- Default to a rectangle
    axis_limits = dict(
        x=dict(xmin=dataset_sub["age"].min(), xmax=dataset_sub["age"].max()),
        y=dict(
            xmin=dataset_sub["length"].min() * 0.25, xmax=dataset_sub["length"].max() * 0.25
        ),
    )

# Prepare default parameters
# ---- x
kwargs.update(dict(xlabel=kwargs.get("xlabel", "Age (years)")))
# ---- y
kwargs.update(dict(ylabel=kwargs.get("ylabel", "Fork length (cm)")))

# Initialize figure
# ---- Update the 'figsize' if it doesn't yet exist
if not kwargs.get("figsize", None):
    kwargs.update(dict(figsize=kwargs.get("figsize", apply_aspect_ratio(5.5, axis_limits))))
# Initialize figure
fig, ax = plt.subplots(**prune_args(plt.subplots, **kwargs))
# ---- Format the x- and y-axis tick labels and spacing
plot_data = format_heatmap_ticks(ax, heatmap_data)
# ---- Normalize the colormapping
colormap_norm = add_colorbar(
    ax,
    **dict(kwargs, cmap=cmap, colorbar_label=colorbar_label, vmin=vmin, vmax=vmax, x_pad=0.025)
)
# ---- Plot the heatmap
ax.imshow(
    plot_data["heatmap_array"],
    norm=colormap_norm,
    cmap=cmap,
    interpolation=None,
    extent=plot_data["extent"],
    aspect="auto",
)
# ---- Add a heatmap grid, if requested
if kwargs.get("grid_heatmap"):
    add_heatmap_grid(ax=ax, **plot_data)
# ---- Format the figure area axes
format_axes(
    ax, axis_limits=axis_limits, xlabel=kwargs.get("xlabel"), ylabel=kwargs.get("ylabel")
)
# ---- Tighten the layout and display
plt.tight_layout()
plt.show()