import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from echopop import Survey
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import AutoLocator, LogFormatterSciNotation, LogLocator,  ScalarFormatter
import geopandas as gpd
import pandas as pd
from typing import Any, Dict, Optional, Tuple, Union
import numpy as np

def scale_sizes(values, min_value, max_value, min_size=25, max_size=250):
    """
    Scale point size
    """
    # Censor values if needed
    sizes = values.copy()
    sizes.loc[sizes < min_value] = min_value
    sizes.loc[sizes > max_value] = max_value

    return ((sizes - min_value) / (max_value - min_value)) * (max_size - min_size) + min_size

def get_survey_bounds(dataframe: pd.DataFrame, geo_config: Dict[str, Any]) -> gpd.GeoDataFrame:
    """
    Extract the survey data boundary coordinates along the x- and y-axes
    """

    # Get projection
    projection = geo_config["init"]

    # Convert to a GeoDataFrame
    geodataframe = gpd.GeoDataFrame(
        dataframe,
        geometry=gpd.points_from_xy(dataframe["longitude"], dataframe["latitude"]),
        crs=projection,    
    )
    
    # Get total boundaries
    return geodataframe.total_bounds

def add_colorbar(axes: GeoAxes, 
                 colormap: str,
                 label: str,
                 limits: Tuple[float], 
                 log: Optional[float] = None):
    """
    Add colorbar to the plot
    """

    # Apply log-transformation, if defined
    if log:
        # ---- Define locator
        locator = LogLocator(base=log)
        # ---- Define formatter
        formatter = LogFormatterSciNotation(base=log)
        # ---- Define scale
        scale = LogNorm(vmin=np.min(limits), vmax=np.max(limits))
    else: 
        # ---- Define locator
        locator = AutoLocator()
        # ---- Define formatter
        formatter = ScalarFormatter()
        # ---- Return 'None' for scale
        scale = Normalize(vmin=np.min(limits), vmax=np.max(limits))

    # Validate the colormap existence and apply the scaling
    try:        
        # ---- Create the colormap
        scaled_colormap = plt.cm.ScalarMappable(cmap=colormap, norm=scale)
    except ValueError as e:
        # ---- Drop traceback
        e.__traceback__ = None
        # ---- Raise Error
        raise(e)

    # Create the colorbar
    cbar = plt.colorbar(scaled_colormap, ax=axes, ticks=locator, format=formatter, shrink=0.5, 
                        fraction=0.075, pad=0.025)
    # ---- Render with label
    cbar.set_label(label)

    # Return scale
    return  scale

def format_axes(axes: GeoAxes, axis_labels: Dict[str, Any], axis_limits: Dict[str, float]):
    """
    Format plotting axis labels and ticks
    """
    
    # Set axis tick marks/spacing
    # ---- x
    axes.set_xticks(axes.get_xticks())
    # ---- y
    axes.set_yticks(axes.get_yticks())

    # Set labels
    # ---- x
    axes.set_xlabel(axis_labels["x"])
    # ---- y
    axes.set_ylabel(axis_labels["y"])

    # Set the axis limits
    axes.set_extent([np.min(axis_limits["x"]), np.max(axis_limits["x"]),
                     np.min(axis_limits["y"]), np.max(axis_limits["y"])]) 

    # Remove margin padding
    axes.margins(0, 0)

def add_transect_lines(dataset: pd.DataFrame, plot_order: int):
    """
    Add transect lines from the dataset
    """

    # Add transect lines layer
    (
        # ---- Group by each transect
        dataset.groupby(["transect_num"])
        # ---- Plot each transect line
        .apply(lambda transect: plt.plot(transect["longitude"], 
                                         transect["latitude"],
                                         c="#848884",
                                         linewidth=0.5,
                                         zorder=plot_order),
                include_groups=False)
    )

def apply_aspect_ratio(axis_limits: np.ndarray, figure_width: float) -> Tuple[float]:
    """
    Apply the defined aspect ratio on the figure window size
    """

    # Compute the aspect ratio of the plotting data
    phi = (
        (np.max(axis_limits["y"]) - np.min(axis_limits["y"])) 
        / (np.max(axis_limits["x"]) - np.min(axis_limits["x"]))
    )

    # Adjust the figure width and height
    return (figure_width, figure_width * phi)

def get_coastline(container_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Get coastline datafor plotting
    """

    # Get the `cartopy` coastline, if it doesn't already exist
    if (
        "coastline" not in container_dict 
        and not isinstance(container_dict["coastline"], cfeature.NaturalEarthFeature)
    ):
        # ---- Get the coastline information
        coastlines = cfeature.NaturalEarthFeature(
            category="physical",
            name="land", 
            scale="10m",
            facecolor="#C5C9C7",
            edgecolor="#929591",
            linewidth=0.5,
            zorder=1,
        )
        # ---- Update the dictionary
        container_dict.update({
            "coastline": coastlines
        })

def add_alongtransect_data(ax: GeoAxes, dataset: pd.DataFrame, variable: str, 
                           data_range: Tuple[float], colormap: str,
                           colormap_norm: Union[LogNorm, Normalize], plot_order: int):
    """
    Plot transect estimates
    """

    # Sort the data based on coordinates
    data_geo = (
        dataset.copy().sort_values(by=["transect_num", "longitude", "latitude"])
        .reset_index(drop=True)
    )
    
    # 'Hurdle' the data (separate presence/absence and magnitude processes)
    data_hurdle = data_geo.loc[data_geo[variable] > 0.0].copy()
    # ---- Sort by the variable magnitudes
    data_hurdle = data_hurdle.sort_values(by=[variable]).reset_index(drop=True)
    # ---- Normalize the point sizes
    data_hurdle.loc[:, "normalized_var"] = scale_sizes(data_hurdle.loc[:, variable],
                                                       np.min(data_range),
                                                       np.max(data_range),
                                                       1, 75)

    # Add the scatter plot
    ax.scatter(
        x=data_hurdle["longitude"],
        y=data_hurdle["latitude"],
        c=data_hurdle[variable],
        norm=colormap_norm,
        cmap=colormap,
        s=data_hurdle["normalized_var"],
        edgecolor="black",
        linewidth=0.2,
        zorder=plot_order,
    )

def plot_transect(dataset: pd.DataFrame, variable: str, figure_width: float, 
                  geo_config: Dict[str, Any],
                  colormap: Optional[str] = None, 
                  data_range: Optional[Tuple[float]] = None,
                  log_base: Optional[float] = None, 
                  axis_limits: Optional[Dict[str, Tuple[float]]] = None):

    # Determine axis limits
    if not axis_limits:
        # ---- Get survey bounds
        survey_bounds = get_survey_bounds(dataset, geo_config)
        # ---- Additional buffering
        axis_limits = dict(
            x=(survey_bounds[0] * 1.005, survey_bounds[2] * 0.995),
            y=(survey_bounds[1] * 0.995, survey_bounds[3] * 1.005)
        )

    # Compute the figure dimensions based on the data aspect ratio
    figure_dimensions = apply_aspect_ratio(axis_limits, figure_width)

    # Add prepend label string, if required
    if "_male" in variable:
        label_prepend = "Male "
    elif "_female" in variable:
        label_prepend = "Female "
    else:
        label_prepend = ""

    # Prepare default parameters
    axis_labels = dict(x="Longitude (\u00B0E)", y="Latitude (\u00B0N)")

    # Prepare parameters for plotting
    # ---- Abundance
    if "abundance" in variable:
        colormap = colormap if colormap else "viridis"
        data_range = data_range if data_range else (1e0, 1e7)
        label = (label_prepend + "abundance\n$#~animals").capitalize()     
    # ---- Biomass density
    elif "biomass" in variable:
        colormap = colormap if colormap else "plasma"
        data_range = data_range if data_range else (1e0, 1e7)
        label = (label_prepend + "biomass\n$\\mathregular{kg}$").capitalize()
    # ---- Biomass density
    elif "biomass_density" in variable:
        colormap = colormap if colormap else "inferno"
        data_range = data_range if data_range else (1e0, 1e7)
        label = (label_prepend + "biomass density\n$\\mathregular{kg~nmi^{-2}}$").capitalize()
    # ---- NASC
    elif "nasc" in variable:
        colormap = colormap if colormap else "cividis"
        data_range = data_range if data_range else (1e0, 1e4)
        label = "NASC\n$\\mathregular{m^{2}~nmi^{-2}}$"
    # ---- Number density
    elif "number_density" in variable:
        colormap = colormap if colormap else "magma"
        data_range = data_range if data_range else (1e0, 1e7)
        label = (label_prepend + "number density\n$\\mathregular{animals~nmi^{-2}}$").capitalize()

    # Initialize figure
    plt.figure(figsize=figure_dimensions)
    # ---- Definte GeoAxes
    ax = plt.axes(projection=ccrs.PlateCarree())
    # ---- Add coastline
    ax.add_feature(geo_config["coastline"])
    # ---- Add transect lines
    add_transect_lines(dataset, plot_order=1)
    # ---- Normalize the colormapping
    colormap_norm = add_colorbar(ax, colormap, label, data_range, log_base)
    # ---- Add transect data
    add_alongtransect_data(ax, dataset, variable, data_range, colormap, colormap_norm, plot_order=2)
    # ---- Format the figure area axes
    format_axes(ax, axis_labels, axis_limits)
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()

def plot_mesh(dataset: pd.DataFrame, variable: str, figure_width: float, geo_config: Dict[str, Any], 
              colormap: Optional[str] = None, data_range: Optional[Tuple[float]] = None,
              log_base: Optional[float] = None, 
              axis_limits: Optional[Dict[str, Tuple[float]]] = None):

    # Determine axis limits
    if not axis_limits:
        # ---- Get survey bounds
        survey_bounds = get_survey_bounds(dataset, geo_config)
        # ---- Additional buffering
        axis_limits = dict(
            x=(survey_bounds[0] * 1.005, survey_bounds[2] * 0.995),
            y=(survey_bounds[1] * 0.995, survey_bounds[3] * 1.005)
        )

    # Compute the figure dimensions based on the data aspect ratio
    figure_dimensions = apply_aspect_ratio(axis_limits, figure_width)

    # Prepare default parameters
    axis_labels = dict(x="Longitude (\u00B0E)", y="Latitude (\u00B0N)")

    # Prepare units
    units = "kg"
    kriged_variable = "biomass"

    # Adjust values, if required
    dataset.loc[dataset["biomass"] < 0.0, "biomass"] = 0.0
    dataset.loc[dataset["kriged_mean"] < 0.0, "kriged_mean"] = 0.0
    if log_base and dataset[variable].min() == 0.0:
        # ---- Adjust offset
        dataset.loc[:, variable] = dataset.loc[:, variable] + 1e-1

    # Prepare parameters for plotting
    if variable == "biomass":
        colormap = colormap if colormap else "plasma"
        data_range = data_range if data_range else (1, 1e6)
        label = "Kriged biomass\n$\\mathregular{kg}$"
    elif variable == "kriged_mean": 
        colormap = colormap if colormap else "inferno"
        data_range = data_range if data_range else (1, 1e5)
        label = "Kriged " + kriged_variable + f" density\n{units} " + "nmi$^{-2}$"
    elif variable == "kriged_variance":
        colormap = colormap if colormap else "hot"
        data_range = data_range if data_range else (0, 5)
        label = f"Kriged {kriged_variable} density variance" + f"\n({units} " + "nmi$^{-2})^{2}$"
    elif variable == "sample_cv":
        colormap = colormap if colormap else "magma"
        data_range = data_range if data_range else (0, np.ceil(np.max(dataset[variable])/0.1)*0.1)
        label = "Kriged $CV$"
    elif variable == "sample_variance":
        colormap = colormap if colormap else "cividis"
        data_range = data_range if data_range else (0, 1e1)
        label = f"Sample {kriged_variable} variance" + f"\n{units}" + "$^{-2}$"
        
    # Initialize figure
    plt.figure(figsize=figure_dimensions)
    # ---- Definte GeoAxes
    ax = plt.axes(projection=ccrs.PlateCarree())
    # ---- Add coastline
    ax.add_feature(geo_config["coastline"])
    # ---- Normalize the colormapping
    colormap_norm = add_colorbar(ax, colormap, label, data_range, log_base)
    # ---- Add meshed data
    ax.scatter(
        dataset["longitude"], dataset["latitude"], c=dataset[variable], s=2, marker="H", 
        cmap=colormap, norm=colormap_norm
    )
    # ---- Format the figure area axes
    format_axes(ax, axis_labels, axis_limits)
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()

def PLOT_MAP(cls: Survey, kind: str):

    # Reference dictionary
    PLOT_MAP = {
        "mesh": {
            "data": lambda v: v.results["kriging"]["mesh_results_df"].copy(),
            "function": plot_mesh,
            "type": "spatial",
        },
        "transect": {
            "data": lambda v: v.analysis["transect"]["acoustics"]["adult_transect_df"].copy(),
            "function": plot_transect,
            "type": "spatial",
        },
    }

    # Get the associated plotting parameterization/information
    plotting_info = PLOT_MAP.get(kind)
    # ---- Fill the data key
    plotting_info.update({"data": plotting_info["data"](cls)})

    # Return the dictionary
    return plotting_info

    