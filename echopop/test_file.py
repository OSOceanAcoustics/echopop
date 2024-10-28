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
import echopop.graphics.plotting as egp
# plt.close('all')
# matplotlib.use("TkAgg")

init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
survey = Survey(init_config, file_config)
survey.load_survey_data()
survey.load_acoustic_data()
survey.transect_analysis()
survey.fit_variogram()
survey.kriging_analysis(variogram_parameters={"n_lags": 30}, variable="biomass_density")
data = survey.analysis["transect"]["acoustics"]["adult_transect_df"].copy()
geo_config = survey.config["geospatial"]

class survey_copy(Survey):

    def 

survey2 = survey_copy(init_config, file_config)
survey2.load_survey_data()
survey2.load_acoustic_data()
survey2.transect_analysis()

def scale_sizes(values, min_value, max_value, min_size=25, max_size=250):

    # Censor values if needed
    sizes = values.copy()
    sizes.loc[sizes < min_value] = min_value
    sizes.loc[sizes > max_value] = max_value

    return ((sizes - min_value) / (max_value - min_value)) * (max_size - min_size) + min_size

def get_survey_bounds(dataframe: pd.DataFrame, geo_config: Dict[str, Any]) -> gpd.GeoDataFrame:

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



def hurdle_data(dataframe: pd.DataFrame, variable: str):

    # Sort DataFrame
    data_sort = (
        dataframe.copy()
        .sort_values(by=["transect_num", "longitude", "latitude"]).reset_index(drop=True)
    )

    # Extract the non-zero data (ignore the presence/absence process)
    data_subset = data_sort.copy().loc[data_sort.copy()[variable] > 0.0]


    
ax.set_extent([total_bounds[0] * 1.005, 
               total_bounds[2] * 0.995, 
               total_bounds[1] * 0.995, 
               total_bounds[3] * 1.005])


axis_labels = dict(x="Longitude (\u00B0E)", y="Latitude (\u00B0N)")
axis_limits = dict(x=(total_bounds[0] * 1.005, total_bounds[2] * 0.995), y=(total_bounds[1] * 0.995, total_bounds[3] * 1.005))

colormap = mcolors.LogNorm(vmin=1e0, vmax=1e4)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=colormap)
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, fraction=0.075, pad=0.025)
# Compute the aspect ratio
phi = (total_bounds[3] - total_bounds[1]) / (total_bounds[2] - total_bounds[0])



data = data.sort_values(by=["transect_num", "longitude", "latitude"]).reset_index(drop=True)
data_subset = data.copy().loc[data.copy().nasc > 1]
data_subset["normalized_size"] = scale_sizes(data_subset["nasc"], 1e0, 1e4, 1, 75)

data_gdf = gpd.GeoDataFrame(
    data,
    geometry=gpd.points_from_xy(data["longitude"], data["latitude"]),
    crs="epsg:4326",    
)
total_bounds = data_gdf.total_bounds
alpha = (total_bounds[3] - total_bounds[1]) / (total_bounds[2] - total_bounds[0])
width = 9
height = width * alpha
label = "NASC\n$\\mathregular{m^{2}~nmi^{-2}}$"
limits = (1e0, 1e4)
log = 10
colormap = "viridis"

coastlines = cfeature.NaturalEarthFeature(
    category='physical',
    name='land',  # Use 'land' to fill the land portion
    scale='10m',
    facecolor="#C5C9C7",
    edgecolor="#929591",
    linewidth=0.5,
    zorder=1,
)
colormap = mcolors.LogNorm(vmin=1e0, vmax=1e4)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=colormap)
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, fraction=0.075, pad=0.025)

plt.figure(figsize=(width, height*0.75))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(coastlines)
(
    data.groupby('transect_num')
    .apply(lambda group: plt.plot(group["longitude"], group["latitude"], 
                                  c="#848884", linewidth=0.5, 
                                  zorder=2), 
           include_groups=False)
)
data_sort = data_subset.sort_values(by=["nasc"]).reset_index(drop=True)
colormap_norm = add_colorbar(ax, "viridis", label, limits, log)
ax.scatter(
    x=data_sort["longitude"],
    y=data_sort["latitude"],
    c=data_sort["nasc"],
    norm=colormap_norm,
    cmap="viridis",
    s=data_sort["normalized_size"],
    edgecolor="black",
    linewidth=0.2,
    zorder=3
)

format_axes(ax, axis_labels, axis_limits)
plt.tight_layout()
plt.show()

data_range = (1e0, 1e4)
variable = "nasc"
figure_width = 6
colormap = "viridis"

def get_coastline(container_dict: Dict[str, Any]) -> Dict[str, Any]:

    # Get the `cartopy` coastline
    coastlines = cfeature.NaturalEarthFeature(
        category="physical",
        name="land",  # Use 'land' to fill the land portion
        scale="10m",
        facecolor="#C5C9C7",
        edgecolor="#929591",
        linewidth=0.5,
        zorder=1,
    )

    # Add to the appropriate dictionary
    container_dict.update({
        "coastline": coastlines
    })

get_coastline(geo_config)
transect_data = survey.analysis["transect"]["acoustics"]["adult_transect_df"].copy()
egp.plot_transect(transect_data, variable="nasc", figure_width=6.0, geo_config=geo_config, data_range=(1e1, 2e4))
def add_transect_lines(dataset: pd.DataFrame, plot_order: int):

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


def add_alongtransect_data(ax: GeoAxes, dataset: pd.DataFrame, variable: str, 
                           data_range: Tuple[float], colormap: str,
                           colormap_norm: Union[LogNorm, Normalize], plot_order: int):

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

def apply_aspect_ratio(axis_limits: np.ndarray, figure_width: float) -> Tuple[float]:

    # Compute the aspect ratio of the plotting data
    phi = (
        (np.max(axis_limits["y"]) - np.min(axis_limits["y"])) 
        / (np.max(axis_limits["x"]) - np.min(axis_limits["x"]))
    )

    # Adjust the figure width and height
    return (figure_width, figure_width * phi)

def plot_transect(dataset: pd.DataFrame, variable: str, figure_width: float, 
                  colormap: str, geo_config: Dict[str, Any],
                  data_range: Optional[Tuple[float]] = None,
                  log_base: Optional[float] = None, 
                  axis_limits: Optional[Dict[str, Tuple[float]]] = None):

    # Determine axis limits
    if not axis_limits:
        # ---- Get survey bounds
        survey_bounds = get_survey_bounds(dataset, geo_config)
        # ---- Additional buffering
        axis_limits = dict(
            x=(survey_bounds[0], survey_bounds[2]),
            y=(survey_bounds[1], survey_bounds[3])
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

    # Prepare default parameter dictionary
    defaults = {
        "axis_labels": dict(x="Longitude (\u00B0E)", y="Latitude (\u00B0N)"),
    }

    # Prepare parameters for plotting
    # ---- Abundance
    if "abundance" in variable:
        defaults.update({
            "data_range": data_range if data_range else (1e0, 1e7),
            "label": (
                label_prepend + "abundance\n$\\mathregular{#~animals}$"
            ).capitalize()            
        })
    # ---- Biomass density
    elif "biomass" in variable:
        defaults.update({
            "data_range": data_range if data_range else (1e0, 1e7),
            "label": (
                label_prepend + "biomass\n$\\mathregular{kg}$"
            ).capitalize()            
        })
    # ---- Biomass density
    elif "biomass_density" in variable:
        defaults.update({
            "data_range": data_range if data_range else (1e0, 1e7),
            "label": (
                label_prepend + "biomass density\n$\\mathregular{kg~nmi^{-2}}$"
            ).capitalize()            
        })
    # ---- NASC
    elif "nasc" in variable:
        defaults.update({
            "data_range": data_range if data_range else (1e0, 1e4),
            "label": "NASC\n$\\mathregular{m^{2}~nmi^{-2}}$"            
        })
    # ---- Number density
    elif "number_density" in variable:
        defaults.update({
            "data_range": data_range if data_range else (1e0, 1e7),
            "label": (
                label_prepend + "number density\n$\\mathregular{animals~nmi^{-2}}$"
            ).capitalize()            
        })

    # Initialize figure
    plt.figure(figsize=figure_dimensions)
    # ---- Definte GeoAxes
    ax = plt.axes(projection=ccrs.PlateCarree())
    # ---- Add coastline
    ax.add_feature(geo_config["coastline"])
    # ---- Add transect lines
    add_transect_lines(data, plot_order=1)
    # ---- Normalize the colormapping
    colormap_norm = add_colorbar(ax, colormap, label, data_range, log_base)
    # ---- Add transect data
    add_alongtransect_data(ax, dataset, variable, data_range, colormap, colormap_norm, plot_order=2)
    # ---- Format the figure area axes
    format_axes(ax, axis_labels, axis_limits)
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()
    

"NASC\n$\\mathregular{m^{2}~nmi^{-2}}$"         
(label_prepend + "number density\n$\\mathregular{animals~nmi^{-2}}$").capitalize()
label = f"{{{label_prepend}}} number density\n$\\mathregular{animals~nmi^{-2}}$"
data.number_density.max()
1e6
1e7
re.match(r"biomass_density$+_", "biomass_density")
r"_male" in "biomass_density_female"

import re
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
from pydantic import BaseModel, Field, RootModel, ValidationError, field_validator, model_validator

class PlotModel(BaseModel):
    """
    Base Pydantic model for plotting parameter inputs
    """

    # Validator method
    @classmethod
    def judge(cls, **kwargs):
        """
        Validator method
        """
        try:
            return cls(**kwargs)
        except ValidationError as e:
            e.__traceback__ = None
            raise e

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method

        Notes
        ----------
        This is for `pytest` testing.
        """

        return cls.judge(**kwargs).model_dump(exclude_none=True)

class PlotParameters(PlotModel):
    axis_label: Dict[str, str]
    axis_limits: Optional[Tuple[float]] = Field(default=None)
    colormap: str = Field(default="viridis")
    data_range: Tuple[float]
    figure_width: float = Field(default=6.0)
    label: str
    log_base: Optional[float] = Field(default=None)

    @field_validator("axis_limits", "data_range", mode="before")
    def validate_a(cls, v):
        # Check length
        if isinstance(v, tuple) and len(v) != 2:
            raise ValueError(
                "Tuple input must have a length of 2."
            )

class TransectEstimateParameters(PlotModel):

    @classmethod
    def _VALID_ARGS(cls)
    
    # Factory method
    @classmethod
    def parameterize(cls, **kwargs):
        """
        Factory creation method
        """

        return cls.judge(**kwargs).model_dump(exclude_none=True)

cls = 

class DD(PlotModel):
    a: Tuple[float] = Field(default=(0.0, 5.0))

    @field_validator("a", mode="before")
    def validate_a(cls, v):
        # Check length
        if isinstance(v, tuple) and len(v) != 2:
            raise ValueError(
                "Tuple input must have a length of 2."
            )

cls = DD
cls.model_fields
v = (2.0)
DD.judge(**dict(a=(2.0,)))

class TransectEstimateParameters(PlotModel):

    @classmethod
    def _VALID_ARGS(cls)
    
    # Factory method
    @classmethod
    def parameterize(cls, **kwargs):
        """
        Factory creation method
        """

        return cls.judge(**kwargs).model_dump(exclude_none=True)





colormap = mcolors.LogNorm(vmin=1e0, vmax=1e4)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=colormap)
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, fraction=0.075, pad=0.025)
cbar.set_label("NASC\n$\\mathregular{m^{2}~nmi^{-2}}$")
ax.set_xticks(ax.get_xticks())
ax.set_yticks(ax.get_yticks())
ax.set_xlabel("Longitude (\u00B0E)")
ax.set_ylabel("Latitude (\u00B0N)")
ax.margins(0, 0)
ax.set_extent([total_bounds[0] * 1.005, 
               total_bounds[2] * 0.995, 
               total_bounds[1] * 0.995, 
               total_bounds[3] * 1.005])
ax.margins(0, 0)
plt.tight_layout()
plt.show()


# Create a figure and axes using the PlateCarree projection
ax = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(1, 1, projection=ccrs.PlateCarree())
colormap = mcolors.LogNorm(vmin=1e0, vmax=1e4)
# colormap = plt.colormaps.get_cmap(norm).resampled(256)
# custom_cmap = ListedColormap(colormap)
projection = ccrs.PlateCarree()
# Add coastlines to the plot, filling the land portion with gray
coastlines = cfeature.NaturalEarthFeature(
    category='physical',
    name='land',  # Use 'land' to fill the land portion
    scale='10m',
    facecolor='gray',
    edgecolor='black',
    linewidth=0.5,
    zorder=1,
    transform="epsg:4326"
)
ax.add_feature(coastlines)
# data_grp = data.groupby(["transect_num"])
# data_grp.plot(x=data_grp["longitude"], y=data_grp["latitude"], s=2)
# pivot_df = data.pivot(index="longitude", columns="nasc", values="latitude")
data.groupby('transect_num').apply(lambda group: plt.plot(group["longitude"], group["latitude"], c="black", zorder=2), include_groups=False)

ax.scatter(
    x=data_subset["longitude"],
    y=data_subset["latitude"],
    c=data_subset["nasc"],
    norm=colormap,
    cmap="viridis",
    s=data_subset["normalized_size"],
    zorder=3
)

# Set the extent of the plot
ax.set_extent([-135, -120, 33, 58])
ax.set_aspect('equal')
# Show the plot
plt.show()



data = survey.results["kriging"]["mesh_results_df"].copy()
figure_dimensions = apply_aspect_ratio(axis_limits, figure_width)

# Initialize figure
plt.figure(figsize=figure_dimensions)
# ---- Definte GeoAxes
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
ax.add_feature(geo_config["coastline"])
# ---- Add mesh
ax.pcolormesh(data["longitude"], data["latitude"], data["biomass"], vmin=np.min(data_range), vmax=np.max(data_range))
# ---- Format the figure area axes
format_axes(ax, axis_labels, axis_limits)
plt.tight_layout()
plt.show()



import seaborn as sns
import matplotlib.tri as tri
# data_pvt = data.pivot(index="latitude", columns="longitude", values="biomass")
x = data["longitude"]
y = data["latitude"]
z = data["biomass"]
triang = tri.Triangulation(x, y)
X,Y = np.meshgrid(x,y)
Z=z.values.reshape(len(y),len(x))

lon, lat = np.meshgrid(data["longitude"], data["latitude"])                                                      
plt.figure()
plt.tripcolor(triang, z, shading='gouraud')
plt.colorbar()
plt.title('Heatmap of Unevenly Spaced Data')
plt.show()

survey.analysis["kriging"]
data = survey.results["kriging"]["mesh_results_df"]

axis_labels = dict(x="Longitude (\u00B0E)", y="Latitude (\u00B0N)")
axis_limits = None
# Determine axis limits
if not axis_limits:
    # ---- Get survey bounds
    survey_bounds = egp.get_survey_bounds(data, geo_config)
    # ---- Additional buffering
    axis_limits = dict(
        x=(survey_bounds[0], survey_bounds[2]),
        y=(survey_bounds[1], survey_bounds[3])
    )

egp.get_coastline(geo_config)
import echopop.core as ecc

from typing import Literal

dataset = survey.results["kriging"]["mesh_results_df"]
dataset.kriged_variance.max()
variable: Literal["kriged_mean", "kriged_variance", "sample_variance", "sample_cv", "biomass"] = "kriged_mean"
figure_width = 6.0
colormap: Optional[str] = None
data_range: Optional[Tuple[float]] = None
log_base: Optional[float] = None
axis_limits: Optional[Dict[str, Tuple[float]]] = None

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
figure_dimensions = egp.apply_aspect_ratio(axis_limits, figure_width)

# Prepare default parameters
axis_labels = dict(x="Longitude (\u00B0E)", y="Latitude (\u00B0N)")

# Prepare units
units = "kg"
kriged_variable = "biomass"

# Prepare parameters for plotting
log = None
data_range = None
colormap = None
variable = "biomass"
# ---- Abundance
data.loc[data["biomass"] < 0.0, "biomass"] = 0.0
data.loc[data["kriged_mean"] < 0.0, "kriged_mean"] = 0.0
if log and data[variable].min() == 0.0:
    # ---- Adjust offset
    data.loc[:, variable] = data.loc[:, variable] + 1e-1

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
    data_range = data_range if data_range else (0, np.ceil(np.max(data[variable])/0.1)*0.1)
    label = "Kriged $CV$"
elif variable == "sample_variance":
    colormap = colormap if colormap else "cividis"
    data_range = data_range if data_range else (0, 1e1)
    label = f"Sample {kriged_variable} variance" + f"\n{units}" + "$^{-2}$"
    
# Initialize figure
plt.figure(figsize=(6, 8.5))
# ---- Definte GeoAxes
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
ax.add_feature(geo_config["coastline"])
colormap_norm = egp.add_colorbar(ax, colormap, label=label, limits=data_range, log=log)

ax.scatter(
    data["longitude"], data["latitude"], c=data[variable], s=2, marker="H", 
    cmap=colormap, norm=colormap_norm
)

egp.format_axes(ax, axis_labels, axis_limits)

plt.tight_layout()
plt.show()

distribution_data = survey.analysis["transect"]["biology"]["distributions"]["binned_aged_counts_df"].copy()
df = survey.analysis["transect"]["biology"]["population"]["tables"]["abundance"]["aged_abundance_df"].copy()
distribution_data = df.stack(future_stack=True).sum(axis=1).reset_index(name="abundance")
distribution_data["length"] = distribution_data["length_bin"].apply(lambda x: x.mid).astype(float)
distribution_data["age"] = distribution_data["age_bin"].apply(lambda x: x.mid).astype(float)
distribution_sub = distribution_data.loc[distribution_data["sex"] == "all", :]
heatmap_data = distribution_sub.pivot_table(values="abundance", columns="age", index="length", aggfunc="sum")
length_labels = heatmap_data.index.values
age_labels = heatmap_data.columns.values
heatmap_array = heatmap_data.values * 1e-3

from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FixedLocator, FixedFormatter

delta_age = np.ceil(np.diff(age_labels).mean())
delta_length = np.ceil(np.diff(length_labels).mean())

centers = [age_labels.min(), age_labels.max(), length_labels.max(), length_labels.min()]
dx, = np.diff(centers[:2])/(heatmap_array.shape[1]-1)
dy, = -np.diff(centers[2:])/(heatmap_array.shape[0]-1)
extent = [centers[0]-dx/2, centers[1]+dx/2, centers[2]+dy/2, centers[3]-dy/2]

age_axis = np.arange(start=age_labels.min(), stop=age_labels.max(), step=delta_age * 2, dtype=int)
length_axis = np.arange(start=length_labels.min(), stop=length_labels.max(), step=delta_length * 4)

fig, ax = plt.subplots()
im = ax.imshow(heatmap_array, cmap="viridis", vmin=1, vmax=3.5e3, interpolation=None, extent=extent, aspect="auto")
cbar = ax.figure.colorbar(im, ax=ax)
plt.xticks(np.arange(centers[0], centers[1]+dx, dx))
plt.yticks(np.arange(centers[3], centers[2]+dy, dy))
y_locator = FixedLocator(length_axis)
ax.yaxis.set_major_locator(y_locator)
x_locator = FixedLocator(age_axis)
ax.xaxis.set_major_locator(x_locator)
ax.hlines((length_labels - delta_length / 2).tolist(), #transform=ax.get_xaxis_transform(), 
          xmin=age_labels.min() - delta_age / 2, xmax=age_labels.max() + delta_age / 2, colors="black")
ax.vlines((age_labels - delta_age / 2).tolist(), #transform=ax.get_xaxis_transform(), 
          ymin=length_labels.min() - delta_length / 2, ymax=length_labels.max() + delta_length / 2, colors="black")
plt.show()


fig, ax = plt.subplots()
im = ax.imshow(heatmap_array, cmap="viridis", vmin=1, vmax=3.5e3, aspect="auto", zorder=1, extent=extent)

# Set axis tick marks/spacing
# ---- x
# delta_age = np.ceil(np.diff(age_labels).mean())
# xticks = ax.get_xticks()
# xticks_adj = xticks - xticks.min() + age_labels.min()
# ax.set_xticks(xticks_adj)
# # ---- y
# delta_length = np.ceil(np.diff(length_labels).mean())
# yticks = ax.get_yticks()
# yticks_adj = yticks - yticks.min() + length_labels.min()
# ax.set_yticks(ax.get_yticks())

# Automatically reduce the number of ticks if there are too many labels

age_axis = np.arange(start=age_labels.min(), stop=age_labels.max(), step=delta_age * 2, dtype=int)
ax.set_xticks(age_axis.astype(int).tolist())


length_axis = np.arange(start=length_labels.min(), stop=length_labels.max(), step=delta_length * 4)
ax.set_yticks(length_axis.astype(int).tolist())
# xticks = int(age_labels) - 1
# ax.set_xticks(age_labels)


# ax.set_yticklabels([length_labels[i] for i in y_ticks], fontsize=8)

age_ind = np.arange(age_labels.min(), age_labels.max(), 1, dtype=int)
length_ind = np.arange(length_labels.min(), length_labels.max(), 1, dtype=int)
# x_locator = FixedLocator(age_axis)
# y_locator = FixedLocator([.85, 1.15, 1.28, 1.9])
# ax.set_xticklabels([age_labels[i] for i in age_axis], fontsize=8)
# ax.xaxis.set_major_locator(x_locator)

# ax.set_xticklabels(age_axis, fontsize=12)
# ax.vlines(age_ind + 0.5, #transform=ax.get_xaxis_transform(), 
#           ymin=length_labels.min() - delta_length, ymax=length_labels.max(), colors="white")
# ax.hlines(length_ind + 0.5, #transform=ax.get_xaxis_transform(), 
#           xmin=age_labels.min() - delta_age, xmax=age_labels.max(), colors="white")
# ax.hlines(y=np.concatenate([[0.0], length_labels.values]) - 1.0, #transform=ax.get_yaxis_transform(), 
#           xmin=-1.0, xmax=39.0, colors="white")
# ax.vlines([2.0, 3.0, 4.0], 
#           ymin=length_labels.min(), 
#           ymax=length_labels.max(), 
#           transform=ax.get_yaxis_transform())
cbar = ax.figure.colorbar(im, ax=ax)

# ax.grid(b=True, which='major', color='white', linestyle='-')
# ax.grid(which="minor", color="black", linestyle="-", linewidth=2)
# ax.tick_params(which="minor", size=0)
# Set the x-axis limit to 80
plt.show()