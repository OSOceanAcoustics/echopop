from echopop import Survey
import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union
import numpy as np
from echopop.graphics import plotting as egp, variogram_interactive as egv
import holoviews as hv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib
import matplotlib.colors as colors
from scipy.interpolate import griddata
import matplotlib.gridspec as gridspec
import verde as vd
import pyproj
from pyproj import Transformer, CRS
from echopop.spatial.projection import utm_string_generator

init_config = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config.yml"
file_config = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
survey = Survey(init_config, file_config)
survey.load_survey_data()
survey.load_acoustic_data()
survey.transect_analysis()
survey.fit_variogram()
survey.kriging_analysis(variogram_parameters={"n_lags": 30}, variable="biomass_density")

self = survey
kind = "mesh"
variable = "biomass"
figure_width = 5.5
plot_parameters: Dict[str, Any] = {}
axis_limits: Optional[Dict[str, Tuple[float]]] = None
colormap = "viridis"
data_range = None
log_base = None
geo_config = validated_parameters.get("geo_config")
dataset = plot_info.get("data")

# Determine axis limits
if not axis_limits:
    # ---- Get survey bounds
    survey_bounds = get_survey_bounds(dataset, geo_config)
    # ---- Additional buffering
    axis_limits = dict(
        x=(survey_bounds[0] * 1.005, survey_bounds[2] * 0.995),
        y=(survey_bounds[1] * 0.995, survey_bounds[3] * 1.005),
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
    data_range = (
        data_range if data_range else (0, np.ceil(np.max(dataset[variable]) / 0.1) * 0.1)
    )
    label = "Kriged $CV$"
elif variable == "sample_variance":
    colormap = colormap if colormap else "cividis"
    data_range = data_range if data_range else (0, 1e1)
    label = f"Sample {kriged_variable} variance" + f"\n{units}" + "$^{-2}$"


x = dataset['longitude'].values
y = dataset['latitude'].values
z = dataset[variable].values

x_unique = np.unique(x)
y_unique = np.unique(y)

xg = np.unique(np.round(x, 1))
yg = np.unique(np.round(y, 1))

x_unique.size
norm = mcolors.SymLogNorm(linthresh=1, linscale=0.5, vmin=1, vmax=np.max(data_range))

plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
pc = plt.hexbin(x, y, C=z, gridsize=(yg.size - 1, xg.size - 1), reduce_C_function=np.sum, cmap="viridis", norm=norm, edgecolor="black", linewidth=0.05)
ax.add_feature(geo_config["coastline"])
cb = plt.colorbar(ax=pc.axes, norm=mcolors.SymLogNorm(linthresh=1, linscale=1, vmin=1, vmax=np.max(data_range)))
plt.tight_layout()
plt.show()



utm_string_generator(data.longitude.mean(), data.latitude.mean())

proj_epsg = CRS.from_string("EPSG:4326")
projection_epsg = pyproj.Proj(proj_epsg)

proj_utm = CRS.from_string("EPSG:32609")
projection_utm = pyproj.Proj(proj_utm)

data = dataset
proj_coords = projection_epsg(data.longitude.values, data.latitude.values)
# proj_coords = projection_utm(data.longitude.values, data.latitude.values)

projection = pyproj.Proj(proj="merc", lat_ts=proj_coords[1].mean())
proj_coords = projection(proj_coords[0], proj_coords[1])

# projection = pyproj.Proj(proj="merc", lat_ts=data.latitude.mean())
# projection = pyproj.Proj(proj=proj_utm)
# proj_coords = projection(data.longitude.values, data.latitude.values)
# proj_coords = crs_trans.transform(data.longitude.values, data.latitude.values)
# proj_coords[1].max()

spacing = 5 / 60

# chain = vd.Chain(
#     [
#         ("reduce", vd.BlockReduce("mean", spacing * 111e3)),
#         # ("mean", vd.Spline()),
#         ("trend", vd.Trend(degree=1)),
#     ]
# )
chain = vd.Chain(
    [
        ("reduce", vd.BlockReduce("mean", spacing * 111e3)),
        ("spline", vd.KNeighbors(k=5, reduction=np.median)),
        
    ]
)

chain.fit(proj_coords, data[variable])

# reducer = vd.BlockReduce(np.median, spacing=spacing * 111e3)
# filter_coords, filter_bathy = reducer.filter(proj_coords, data[variable])
# spline = vd.Spline().fit(filter_coords, filter_bathy)
region = vd.get_region((data.longitude, data.latitude))
# grid = spline.grid(region=region, spacing=spacing * 111e3, data_names=f"{variable}")
# grid = vd.distance_mask(proj_coords, maxdist=10e3, grid=grid)

# chain.fit(proj_coords, data[variable])
grid = chain.grid(
    region=region,
    spacing=spacing,
    projection=projection,
    data_names=f"{variable}",
)
grid = vd.distance_mask(
    data_coordinates=(data.longitude, data.latitude),
    maxdist=spacing * 111e3,
    grid=grid,
    projection=projection,
)

from matplotlib.colors import SymLogNorm, Normalize

def add_colorbar(
    axes: GeoAxes,
    colormap: str,
    label: str,
    limits: Tuple[float, float],
    log: Optional[float] = None,
    x_pad: Optional[float] = None,
):
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
        scale = SymLogNorm(linthresh=min=np.min(limits), vmax=np.max(limits))
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
        raise (e)

    # Add pad, if needed
    x_pad = x_pad if x_pad else 0.0

    # Create the colorbar
    cbar = plt.colorbar(
        scaled_colormap,
        ax=axes,
        ticks=locator,
        format=formatter,
        shrink=0.5,
        fraction=0.075,
        pad=0.025 + x_pad,
    )
    # ---- Render with label
    cbar.set_label(label)

    # Return scale
    return scale

from typing import Literal
from matplotlib.colors import LogNorm, Normalize, SymLogNorm

cmap = "plasma"
vmin, vmax = (0, 1e6)
figsize = (7, 7)
kwargs = dict(cmap=cmap)
            #   , add_colorbar=True, norm=SymLogNorm(linthresh=1))
type: Literal["hexbins", "pcolor", "scatter"] = "hexbins"
norm = SymLogNorm(linthresh = 1, vmin=vmin, vmax=vmax)

x = dataset['longitude'].values
y = dataset['latitude'].values
z = dataset[variable].values

xg = np.unique(np.round(x, 1))
yg = np.unique(np.round(y, 1))

import inspect

data_args = dict(x=x, y=y)
plot_args = {}

class MeshPlot(SpatialPlot):

    type: Literal["hexbin", "pcolormesh", "scatter"]
    variable: Literal["biomass", "kriged_mean", "kriged_variance", "sample_cv", "sample_variance"]

class HexbinPlot(BaseModel):

    cmap: str
    extent: Tuple[float, float, float, float]
    gridsize: Union[int, Tuple[int, int]]
    vmin: float
    vmax: float
    x: np.ndarray[float]
    y: np.ndarray[float]

    @classmethod
    def _DEFAULT_ARG_FACTORY(cls, type):
        
        # Create reference default parameters
        
        return {
            "age_length_distribution": BiologicalHeatmapPlot,
            "mesh": MeshPlot,
            "transect": TransectPlot,
        }
    

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
        data_range = (
            data_range if data_range else (0, np.ceil(np.max(dataset[variable]) / 0.1) * 0.1)
        )
        label = "Kriged $CV$"
    elif variable == "sample_variance":
        colormap = colormap if colormap else "cividis"
        data_range = data_range if data_range else (0, 1e1)
        label = f"Sample {kriged_variable} variance" + f"\n{units}" + "$^{-2}$"


if type == "hexbin":
    # ---- Define plotting function
    plot_func = plt.hexbin
    # ---- Get function arguments
    func_args = dict(inspect.signature(plot_func).parameters)
    # ---- Update the data arguments
    data_args.update(dict(gridsize = (yg.size - 1, xg.size - 1)))
    # ----


    plot_args = dict(
        x=x, y=y, gridsize=(yg.size - 1, xg.size - 1), reduce_C_function=np.sum, cmap=cmap, 
        norm=norm, edgecolor="black", linewidth=0.05
    )
    plot_args.update({param: kwargs.get(param, None) for param in kwargs if param in plt.hexbin.__code__.co_varnames})   
    plot_func = plt.hexbin

dd = inspect.signature(plot_func)
dict(dd.parameters)
inspect.getfullargspec(plot_func)
inspect.getcallargs(plt.hexbin)
plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
# plot_func(**plot_args)
plot_func(x, y, C=z, gridsize=(yg.size - 1, xg.size - 1), reduce_C_function=np.sum, cmap="viridis", 
          norm=norm, edgecolor="black", linewidth=0.05)
plt.show()

pt =  plt.hexbin(x, y, C=z, gridsize=(yg.size - 1, xg.size - 1), reduce_C_function=np.sum, cmap="viridis", norm=norm, edgecolor="black", linewidth=0.05)

plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
pt
plt.show()


x = dataset['longitude'].values
y = dataset['latitude'].values
z = dataset[variable].values

xg = np.unique(np.round(x, 1))
yg = np.unique(np.round(y, 1))

x_unique.size
norm = mcolors.SymLogNorm(linthresh=1, linscale=0.5, vmin=1, vmax=np.max(data_range))

plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
pc = plt.hexbin(x, y, C=z, gridsize=(yg.size - 1, xg.size - 1), reduce_C_function=np.sum, cmap="viridis", norm=norm, edgecolor="black", linewidth=0.05)
ax.add_feature(geo_config["coastline"])
cb = plt.colorbar(ax=pc.axes, norm=mcolors.SymLogNorm(linthresh=1, linscale=1, vmin=1, vmax=np.max(data_range)))
plt.tight_layout()
plt.show()




# faults.
reducer = vd.BlockReduce(np.median, spacing=spacing * 111e3)
filter_coords, filter_bathy = reducer.filter(proj_coords, data[variable])
grd = vd.KNeighbors(k=10, reduction=np.mean)
grd.fit(proj_coords, data[variable])
grid = grd.grid(
    region=region,
    spacing=spacing,
    projection=projection,
    dims=["latitude", "longitude"],
    data_names=f"{variable}",
)
grid = vd.distance_mask(
    data_coordinates=(data.longitude, data.latitude),
    maxdist=spacing * 111e3,
    grid=grid,
    projection=projection,
)

plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
grid[variable].plot.pcolormesh(**kwargs)
ax.add_feature(geo_config["coastline"])
plt.show()


ax = plt.axes(projection=ccrs.PlateCarree())
# ax.add_feature(geo_config["coastline"])
grid[variable].plot.pcolormesh(**kwargs)
plt.show()


plt.__dir__()
plt.pcolormesh.__code__.co_varnames
    
plt.figure()

plt.figure(figsize=(7, 6))
plt.title("Gridded bathymetry in Cartesian coordinates")
pc = grid[variable].plot.pcolormesh(cmap="viridis", vmax=1e6, add_colorbar=False)
plt.colorbar(pc).set_label("bathymetry (m)")
# plt.plot(filter_coords[0], filter_coords[1], ".k", markersize=0.5)
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.show()



import pygmt

plt.figure(figsize=(7, 6))
plt.title("Gridded bathymetry in Cartesian coordinates")
pc = grid[variable].plot.pcolormesh(cmap="viridis", vmax=0, add_colorbar=False)
plt.colorbar(pc).set_label("bathymetry (m)")
plt.plot(filter_coords[0], filter_coords[1], ".k", markersize=0.5)
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.show()

plt.figure(figsize=(7, 6))
ax = plt.axes(projection=ccrs.Mercator())
ax.set_title("Geographic grid of bathymetry")
pc = grid_geo[variable].plot.pcolormesh(
    ax=ax, transform=ccrs.PlateCarree(), vmax=0, zorder=-1, add_colorbar=False
)
plt.colorbar(pc).set_label("meters")
plt.show()

reducer = vd.BlockReduce("mean", spacing=spacing * 111e3)
filter_coords, filter_bathy = reducer.filter(proj_coords, data[variable])

spline = vd.BlockMean().fit(filter_coords, filter_bathy)


spline = vd.Spline().fit(filter_coords, filter_bathy)
grid = spline.grid(spacing=spacing * 11e3, data_names=f"{variable}")



region = vd.get_region((dataset.longitude, dataset.latitude))
grid_geo = spline.grid(
    region=region,
    spacing=spacing,
    projection=projection,
    dims=["latitude", "longitude"],
    data_names="bathymetry",
)

plt.figure(figsize=(7, 6))
plt.title("Gridded bathymetry in Cartesian coordinates")
pc = grid.bathymetry.plot.pcolormesh(cmap="viridis", vmin=0, vmax=1e5, add_colorbar=True)
plt.plot(filter_coords[0], filter_coords[1], ".k", markersize=0.5)
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.show()

plt.figure(figsize=(7, 6))
plt.title("Gridded bathymetry in Cartesian coordinates")
pc = grid[variable].plot.pcolormesh(cmap="viridis", vmax=1e5, add_colorbar=True)
plt.colorbar(pc).set_label("bathymetry (m)")
# plt.plot(filter_coords[0], filter_coords[1], ".k", markersize=0.5)
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.show()


x = dataset["longitude"].values
y = dataset["latitude"].values
z = dataset[variable].values
dataset_xr = dataset.filter(["longitude", "latitude", variable]).set_index(["longitude", "latitude"]).to_xarray()


lon2d, lat2d = np.meshgrid(x, y)
proj = ccrs.PlateCarree()
extent = [
        np.min(axis_limits["x"]),
        np.max(axis_limits["x"]),
        np.min(axis_limits["y"]),
        np.max(axis_limits["y"]),
]

ax = plt.axes(projection=proj)
ax.set_extent(extent, crs=proj)
dataset_xr[variable].plot.pcolormesh(ax=ax, transform=proj, robust=True)
plt.show()

dataset[variable].plot.pcolormesh(ax=ax, transform=proj, robust=True)




numcols, numrows = 240, 240
xi = np.linspace(x.min(), x.max(), numcols)
yi = np.linspace(y.min(), y.max(), numrows)
xi, yi = np.meshgrid(xi, yi)

# interpolate, there are better methods, especially if you have many datapoints
zi = griddata((x, y), z, (xi, yi), method='cubic')

fig, ax = plt.subplots(figsize=(12, 12))
ax.pcolormesh(xi, yi, zi, 500, cmap='magma', zorder = 2)
plt.show()

extent = [
        np.min(axis_limits["x"]),
        np.max(axis_limits["x"]),
        np.min(axis_limits["y"]),
        np.max(axis_limits["y"]),
]
z[np.isnan(z)] = 0.0

dataset_xr = dataset.filter(["longitude", "latitude", variable]).set_index(["longitude", "latitude"]).to_xarray()
dataset_xr[variable].plot(norm=colors.SymLogNorm(linthresh=1.0, vmin=1, vmax=1e5))
plt.show()

albo = ccrs.AlbersEqualArea(central_latitude=0,
 false_easting=0, 
 false_northing=0, 
 central_longitude=x.mean(), 
 standard_parallels=(20, 50) )

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection=albo)
ax.set_extent(extent, crs=albo)
ax.coastlines( zorder=9)

pcm = ax.imshow(z,
                extent=extent,
                cmap="viridis",
                origin="lower",
                norm=colors.SymLogNorm(linthresh=1.0, vmin=1, vmax=1e5))
cbar = fig.colorbar(pcm, ax=ax,  shrink=0.5, pad=0.1, extend='min')
ax.gridlines()
plt.show()

xy = np.concatenate([x.ravel()[:, None], y.ravel()[:, None]], axis=1)
grid_xy = np.meshgrid(x, y)
grid_xy = np.vstack((grid_xy[0].ravel(), grid_xy[1].ravel())).transpose()

import xarray as xr

import matplotlib.tri as mtri
dataset_xr = dataset.filter(["longitude", "latitude", variable]).set_index(["longitude", "latitude"]).to_xarray()

z = dataset_xr[variable].values
x = dataset_xr['longitude'].values
y = dataset_xr['latitude'].values

fig = plt.figure(figsize=(6,6))    
ax = plt.axes(projection=ccrs.PlateCarree())
ax.pcolormesh(lon, lat, topo_data ,transform=ccrs.PlateCarree())
ax.set_extent([-134.8, -119.95, 34.25, 55.1], ccrs.PlateCarree())
ax.coastlines()
plt.show()

fig = plt.figure(figsize=(6,6))    
ax = plt.axes(projection=ccrs.PlateCarree())
pm = plt.pcolormesh(x, y, z)
ax.set_extent([-134.8, -119.95, 34.25, 55.1], ccrs.PlateCarree())
plt.show()

# Create a mesh grid from the unique values of longitude and latitude
x_unique = np.unique(x)
y_unique = np.unique(y)
norm = mcolors.SymLogNorm(linthresh=1, linscale=0.5, vmin=1, vmax=np.max(data_range))

import wradlib as wrl
import matplotlib.pyplot as plt
import warnings
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_dem(ax):
    filename = wrl.util.get_wradlib_data_file("geo/bangladesh.tif")
    ds = wrl.io.open_raster(filename)
    # pixel_spacing is in output units (lonlat)
    ds = wrl.georef.reproject_raster_dataset(ds, spacing=0.005)
    rastervalues, rastercoords, proj = wrl.georef.extract_raster_dataset(ds)
    # specify kwargs for plotting, using terrain colormap and LogNorm
    dem = ax.pcolormesh(
        rastercoords[..., 0],
        rastercoords[..., 1],
        rastervalues,
        cmap=plt.cm.terrain,
        norm=LogNorm(vmin=1, vmax=3000),
    )
    # make some space on the right for colorbar axis
    div1 = make_axes_locatable(ax)
    cax1 = div1.append_axes("right", size="5%", pad=0.1)
    # add colorbar and title
    # we use LogLocator for colorbar
    cb = plt.gcf().colorbar(dem, cax=cax1, ticks=ticker.LogLocator(subs=range(10)))
    cb.set_label("terrain height [m]")

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, aspect="equal")
plot_dem(ax)

gridded = wrl.comp.togrid(
    xy,
    grid_xy,
    128000.0,
    np.array([x.mean(), y.mean()]),
    z.ravel(),
    wrl.ipol.Nearest,
)

plt.figure(figsize=figure_dimensions)
pc = plt.hexbin(x, y, C=z, gridsize=200, reduce_C_function=np.sum, cmap="viridis", norm=norm)
cb = plt.colorbar(ax=pc.axes, norm=mcolors.SymLogNorm(linthresh=1, linscale=1, vmin=1, vmax=np.max(data_range)))
plt.tight_layout()
plt.show()

from scipy import interpolate
# Create interpolator:
N = 10000
a = np.array([-10, -10, 0])
b = np.array([15, 15, 0])
x0 = 3*np.random.randn(N, 3) + a
x1 = 5*np.random.randn(N, 3) + b
xx = np.vstack([x0, x1])
v0 = np.exp(-0.01*np.linalg.norm(x0-a, axis=1)**2)
v1 = np.exp(-0.01*np.linalg.norm(x1-b, axis=1)**2)
v = np.hstack([v0, v1])

ndpol = interpolate.LinearNDInterpolator(np.vstack([x, y])[:,:2], z)
X, Y = np.meshgrid(x, y, sparse=True)
V = ndpol(list(zip(X.ravel(),Y.ravel()))).reshape(X.shape)






ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
ax.add_feature(geo_config["coastline"])
topo_plot = ax.scatter(c=topo_data, vmin=0, vmax=1000, transform=data_proj)
ax.set_extent(
    [
        np.min(axis_limits["x"]),
        np.max(axis_limits["x"]),
        np.min(axis_limits["y"]),
        np.max(axis_limits["y"]),
    ]
)
axpos = ax.get_position()
plt.tight_layout()
plt.show()





# Create a 2D grid for the biomass values
Z = dataset.pivot_table(index='latitude', columns='longitude', values=variable)
Z = Z.values
X, Y = np.meshgrid(x_unique, y_unique)
norm = mcolors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=np.min(data_range), vmax=np.max(data_range))





nrows, ncols = y_unique.size, x_unique.size
grid = z.reshape((nrows, ncols))

map_proj = ccrs.PlateCarree()
data_proj = ccrs.PlateCarree()
geodetic_proj = ccrs.PlateCarree()

fig = plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
ax.add_feature(geo_config["coastline"])
topo_plot = ax.contourf(lon,\
                        lat,\
                        topo_data,\
                        levels=[0,200,400,600,800,1000],\
                        extend='max',\
                        transform=data_proj)
ax.set_extent(
    [
        np.min(axis_limits["x"]),
        np.max(axis_limits["x"]),
        np.min(axis_limits["y"]),
        np.max(axis_limits["y"]),
    ]
)
axpos = ax.get_position()
cbar_ax = fig.add_axes([axpos.x1+0,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(topo_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
plt.tight_layout()
plt.show()

fig = plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
ax.add_feature(geo_config["coastline"])
topo_plot = ax.scatter(c=topo_data, vmin=0, vmax=1000, transform=data_proj)
ax.set_extent(
    [
        np.min(axis_limits["x"]),
        np.max(axis_limits["x"]),
        np.min(axis_limits["y"]),
        np.max(axis_limits["y"]),
    ]
)
axpos = ax.get_position()
plt.tight_layout()
plt.show()

# scatter = ax.scatter(dataset["longitude"], dataset["latitude"], c=dataset[variable], cmap='viridis', norm=norm, transform=ccrs.PlateCarree())
pmesh = ax.pcolormesh(X, Y, Z, shading='auto', cmap='viridis', norm=norm)
plt.colorbar(pmesh, ax=ax, orientation='vertical', shrink=0.7)
# cbar.set_label('SymLog-scaled color values')
plt.tight_layout()
plt.show()


plt.figure(figsize=figure_dimensions)
ax.imshow(lon, lat, topo_data, origin='lower')
plt.tight_layout()
plt.show()

ax = plt.axes(projection=ccrs.PlateCarree())
plt.contourf(x, y, z, 60, transform=ccrs.PlateCarree())
plt.tight_layout()
plt.show()

plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(geo_config["coastline"])
topo_plot, aa = ax.pcolormesh(
    [x, y], z, vmin=0.0, vmax=1e5, transform=ccrs.PlateCarree()
)

plt.tight_layout()
plt.show()

# make data
X, Y = np.meshgrid(np.linspace(-3, 3, 16), np.linspace(-3, 3, 16))
Z = (1 - X/2 + X**5 + Y**3) * np.exp(-X**2 - Y**2)

# plot
fig, ax = plt.subplots()

ax.imshow(Z, origin='lower')

plt.show()


X, Y = np.meshgrid(dataset["longitude"], dataset["latitude"])

plt.figure(figsize=figure_dimensions)
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
ax.add_feature(geo_config["coastline"])
norm = mcolors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=np.min(data_range), vmax=np.max(data_range))
# scatter = ax.scatter(dataset["longitude"], dataset["latitude"], c=dataset[variable], cmap='viridis', norm=norm, transform=ccrs.PlateCarree())
pmesh = ax.pcolormesh(X, Y, dataset[variable].values, shading='auto', cmap='viridis', norm=norm)
plt.colorbar(pmesh, ax=ax, orientation='vertical', shrink=0.7)
# cbar.set_label('SymLog-scaled color values')
plt.tight_layout()
plt.show()

ax.scatter(
    dataset["longitude"],
    dataset["latitude"],
    c=dataset[variable],
    s=2,
    marker="H",
    cmap="viridis",
    norm=norm,
)
plt.colorbar(ax=ax)
plt.tight_layout()
plt.show()


c = ax.imshow(dataset["longitude"], dataset["latitude"], , norm=norm, cmap='viridis')


# Initialize figure
plt.figure(figsize=figure_dimensions)
# ---- Define GeoAxes
ax = plt.axes(projection=ccrs.PlateCarree())
# ---- Add coastline
ax.add_feature(geo_config["coastline"])
# ---- Normalize the colormapping
colormap_norm = add_colorbar(ax, colormap, label, data_range, log_base)
# ---- Add meshed data
ax.scatter(
    dataset["longitude"],
    dataset["latitude"],
    c=dataset[variable],
    s=2,
    marker="H",
    cmap=colormap,
    norm=colormap_norm,
)
# ---- Format the figure area axes
format_axes(ax, axis_labels, axis_limits)
# ---- Tighten the layout and display
plt.tight_layout()
plt.show()

import matplotlib.tri as tri
# Sample irregular grid data
x = np.array([0, 1, 2, 3, 4])
y = np.array([0, 1, 2, 3])  # y does not have to match the shape of x
data = np.array([[1, 2, 3, 4],
                 [5, 6, 7, 8],
                 [9, 10, 11, 12],
                 [13, 14, 15, 16]])

# Create a mesh grid for pcolormesh
X, Y = np.meshgrid(x, y)

# Use SymLogNorm for color normalization
norm = mcolors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=data.min(), vmax=data.max())

# Create the plot
plt.figure(figsize=(10, 5))
c = plt.pcolormesh(X, Y, data, shading='flat', cmap='viridis', norm=norm)

# Add colorbar
plt.colorbar(c, label='SymLog-scaled Data values')

# Label axes
plt.xlabel('X-axis')
plt.ylabel('Y-axis')

plt.title('Irregular Mesh Data with pcolormesh and SymLogNorm')
plt.show()

kind="mesh"
variable="biomass"
log_base=10
survey.plot(kind="age_length_distribution", variable="biomass", log_base=10)

dataset = survey.analysis["transect"]["biology"]["population"]["tables"]["biomass"]["aged_biomass_df"]
data_dict = survey.analysis["transect"]["biology"]["population"]["tables"]
variable = "biomass"
sex = "all"
dataset_stk = dataset.stack(future_stack=True).sum(axis=1).reset_index(name=variable)
# ---- Convert the length and age columns from intervals into numerics
# -------- Length
dataset_stk["length"] = dataset_stk["length_bin"].apply(lambda x: x.mid).astype(float)
# -------- Age
dataset_stk["age"] = dataset_stk["age_bin"].apply(lambda x: x.mid).astype(float)
# ---- Subset the data based on the sex
dataset_sub = dataset_stk.loc[dataset_stk["sex"] == sex, :]

dataset_sub["biomass"].max()
5e8