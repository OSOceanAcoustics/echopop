
import geoviews as gv
import holoviews as hv
import hvplot.pandas
import cartopy.crs as ccrs
from echopop import Survey
# setting bokeh as backend
# hv.extension('bokeh')
from geoviews import opts
import geoviews.tile_sources as gts
import numpy as np
# going to use show() to open plot in browser
from bokeh.plotting import show
init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
survey = Survey(init_config, file_config)
survey.load_survey_data()
survey.load_acoustic_data()
survey.transect_analysis()
survey.input["biology"]["specimen_df"]["weight"].mean()

def scale_sizes(values, min_value, max_value, min_size=25, max_size=250):

    # Censor values if needed
    sizes = values.copy()
    sizes.loc[sizes < min_value] = min_value
    sizes.loc[sizes > max_value] = max_value

    return ((sizes - min_value) / (max_value - min_value)) * (max_size - min_size) + min_size

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import geopandas as gpd

data = survey.analysis["transect"]["acoustics"]["adult_transect_df"].copy()
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
ax.scatter(
    x=data_sort["longitude"],
    y=data_sort["latitude"],
    c=data_sort["nasc"],
    norm=colormap,
    cmap="viridis",
    s=data_sort["normalized_size"],
    edgecolor="black",
    linewidth=0.2,
    zorder=3
)
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



# data = survey.analysis["transect"]["acoustics"]["adult_transect_df"].copy()
# data["normalized_size"] = scale_sizes(data["nasc"], 1e1, 1e3, 5, 100)
# data_sorted = data.sort_values('normalized_size')
# plot = data_sorted.hvplot(
#     kind="scatter", 
#     x="longitude", 
#     y="latitude",
#     color="nasc",
#     size="normalized_size",
#     cmap="viridis",
#     legend=True
# )
# # land = gts.OSM.opts(xlim=(-134.0, -120), ylim=(30, 58))
# land = gv.feature.land.opts(
#     scale='50m', 
#     projection=ccrs.PlateCarree(), 
#     alpha=1, 
#     fill_color="gray", 
#     line_color="black",
#     # xlim=(-134.0, -120),
#     # ylim=(30, 58)
# ).redim.range(longitude=(-134.0, -120), latitude=(30, 58))
# # Overlay the plots
# combined_plot = land * plot
# # Set the plot size (width and height) and display the plot
# # combined_plot.opts(
# #     opts.Overlay(width=800, height=600),  # Set the size of the window
# #     opts.Points(size=5)  # Adjust the point size if needed
# # )
# show(
#     hv.render(
#         combined_plot.opts(
#             opts.Overlay(width=800, height=600)
#         )
#     )
# )

# import matplotlib.pyplot as plt
# import cartopy
# import cartopy.feature as cfeature
# import cartopy.crs as ccrs
# import numpy as np

# plt.figure(figsize=(8, 8))
# ax=plt.axes(projection=ccrs.PlateCarree())
# ax.add_feature(cartopy.feature.LAND.with_scale("10m"), facecolor="gray")
# ax.coastlines("10m", color="black", linewidth=1.0)
# plt.show()