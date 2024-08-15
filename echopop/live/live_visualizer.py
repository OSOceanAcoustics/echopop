from echopop.live.sql_methods import SQL
from shapely import wkt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import geopandas as gpd
from typing import Union, Optional
from pathlib import Path

def plot_livesurvey_grid(grid_db: Union[Path, pd.DataFrame],
                         projection: str,
                         coast_db: Optional[Union[Path, pd.DataFrame]] = None):

    # Extract grid data from database if needed
    if isinstance(grid_db, Path):
        # ---- SELECT
        grid_data = SQL(grid_db, "select", table_name="grid_df")
    elif not isinstance(grid_db, pd.DataFrame):
        raise TypeError(
            "Grid data input (`grid_data`) must either be a `Path` or `pandas.DataFrame` object."
        )
    else:
        grid_data = grid_db
    
    # Extract coast data from database if needed
    if isinstance(coast_db, Path):
        # ---- SELECT
        coast_data = SQL(coast_db, "select", table_name="coastline_df")  
    elif coast_data is None:
        # ---- SELECT from `grid_data`
        coast_data = SQL(grid_db, "select", table_name="coastline_df")  
    elif not isinstance(coast_db, pd.DataFrame):
        raise TypeError(
            "Coast data input (`coast_data`) must either be a `Path` or `pandas.DataFrame` object, "
            "or exist within the SQL database as a table (`'coastline_df'`) within the `grid_data` "
            "input (i.e. `grid_data.db`)."
        )      
    else:
        coast_data = coast_db  
    
    # Format columns if needed (well-known-text to Polygon)
    # ---- `grid_data`
    if isinstance(grid_data["geometry"][0], str):
        grid_data["geometry"] = grid_data["geometry"].apply(wkt.loads)
    # ---- `coastline_data`
    if isinstance(coast_data["geometry"][0], str):
        coast_data["geometry"] = coast_data["geometry"].apply(wkt.loads)    
    
    # Generate GeoDataFrames
    # ---- `grid`
    grid_gdf = gpd.GeoDataFrame(grid_data, geometry="geometry", crs=projection)
    # ---- `coast`
    coast_gdf = gpd.GeoDataFrame(coast_data, geometry="geometry", crs=projection)

    # Get appropriate plot axis-limits
    axis_limits = grid_gdf.total_bounds

    # Variable label dictionary map
    VARIABLE_MAP = {
        "number_density_mean": {
            "name": "Mean number density",
            "units": "fish $\\mathregular{nmi^{-2}}$",
            "colormap": "viridis",
        }, 
        "biomass_density_mean": {
            "name": "Mean biomass density",
            "units": "kg $\\mathregular{nmi^{-2}}$",
            "colormap": "plasma",
        },     
        "biomass": {
            "name": "Biomass",
            "units": "kg",
            "colormap": "cividis",
        },
        "abundance": {
            "name": "Abundance",
            "units": "$\\it{N}$",
            "colormap": "inferno",
        }
    }

    # Create a figure and a 2x2 grid of subplots
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    # List of variables to plot
    variables = list(VARIABLE_MAP.keys())

    # Iterate through and plot all subplots
    for ax, var in zip(axes.flat, variables):
        # ---- Get the colormap
        colormap = plt.colormaps.get_cmap(VARIABLE_MAP[var]["colormap"]).resampled(256)
        # ---- Invert
        newcolors = colormap (np.linspace(0, 1, 256))[::-1]
        # ---- Define `white`
        white = np.array([1, 1, 1, 1])
        # ---- Replace "start" color
        newcolors[0, :] = white
        # ---- Create the new custom colormap
        custom_cmap = ListedColormap(newcolors)
        # ---- Normalize colorscale
        norm=plt.Normalize(vmin=grid_gdf[var].min(), vmax=grid_gdf[var].max())
        # ---- Plot the polygons with color fills based on the variable (non-zero)
        grid_gdf.plot(column=var, ax=ax, edgecolor="gainsboro", legend=False, cmap=custom_cmap,
                      norm=norm,
                      markersize=0, linewidth=0.5)        
        # ---- Add coastline data layer
        coast_gdf.plot(ax=ax, linewidth=1.2, color='gray', edgecolor="black")
        # ---- Set axis limits
        ax.set_xlim(axis_limits[0]*1.005, axis_limits[2]*1.01)
        ax.set_ylim(axis_limits[1]*0.98, axis_limits[3]*1.005)
        # ---- Trim down the margins
        ax.margins(0,0)
        # ---- Set adjustable aspect ratio
        # ax.set_aspect('equal', adjustable='box')
        # ---- Set the title and labels
        var_info = VARIABLE_MAP[var]
        ax.set_title(f"{var_info['name']}")
        # ---- Set axis labels
        plt.xlabel(u'Longitude (\u00B0E)')
        plt.ylabel(u'Latitude (\u00B0N)')
        # ---- Add colorbar
        sm = plt.cm.ScalarMappable(cmap=custom_cmap, 
                                   norm=plt.Normalize(vmin=grid_gdf[var].min(), 
                                                      vmax=grid_gdf[var].max()))
        sm._A = []  # fake up the array of the scalar mappable
        cbar = fig.colorbar(sm, ax=ax, shrink=0.5)
        cbar.set_label(f"{var_info['units']}")
        # ---- Add scalebar
        scalebar_length = 250  # Length of scale bar in km
        scalebar_length_in_degrees = scalebar_length / 111  # Assuming 1 degree = 111 km
        # ---- Transform scale bar coordinates to axis units
        # scalebar_x = axis_limits[0]*1.005 + (axis_limits[2]*1.01 - axis_limits[0]*1.005) * 0.1
        # scalebar_y = axis_limits[1]*0.98 + (axis_limits[3]*1.005 - axis_limits[1]*0.98) * 0.1
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        x_scale = (x1 - x0) * 0.1
        y_scale = (y1 - y0) * 0.1
        # scalebar_y_offset = (axis_limits[3]*1.005 - axis_limits[1]*0.98) * 0.05
        # ---- Plot scalebar        
        # ax.plot([scalebar_x, scalebar_x + scalebar_length / 100], 
        #         [scalebar_y, scalebar_y], color='black', lw=2)
        ax.plot([x0 + x_scale, x0 + x_scale + scalebar_length_in_degrees], 
                [y0 + y_scale, y0 + y_scale], color='black', lw=2)
        # ---- Add scale text
        ax.text(x0 + x_scale + scalebar_length_in_degrees / 2, y0 + y_scale - (y1 - y0) * 0.025, 
                f'{scalebar_length} km', ha='center', va='top', color='black')

        # ax.text(scalebar_x + (scalebar_length / 200), 
        #         scalebar_y - scalebar_y_offset, 
        #         f'{scalebar_length} km', ha='center', va='bottom', color='black')

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()

def plot_livesurvey_track(survey_data_db: Union[Path, pd.DataFrame],
                          projection: str,
                          coast_db: Optional[Union[Path, pd.DataFrame]] = None):

    # Extract grid data from database if needed
    if isinstance(survey_data_db, Path):
        # ---- SELECT
        survey_data = SQL(survey_data_db, "select", table_name="survey_data_df")
    elif not isinstance(survey_data_db, pd.DataFrame):
        raise TypeError(
            "Grid data input (`grid_data`) must either be a `Path` or `pandas.DataFrame` object."
        )
    else:
        survey_data = survey_data_db
    
    # Extract coast data from database if needed
    if isinstance(coast_db, Path):
        # ---- SELECT
        coast_data = SQL(coast_db, "select", table_name="coastline_df")   
    elif not isinstance(coast_db, pd.DataFrame):
        raise TypeError(
            "Coast data input (`coast_data`) must either be a `Path` or `pandas.DataFrame` object."
        )      
    else:
        coast_data = coast_db
    
    # Format columns if needed (well-known-text to Polygon)
    # ---- `coastline_data`
    if isinstance(coast_data["geometry"][0], str):
        coast_data["geometry"] = coast_data["geometry"].apply(wkt.loads)    
    
    # Generate GeoDataFrames
    # ---- `grid`
    survey_gdf = gpd.GeoDataFrame(survey_data, 
                                  geometry=gpd.points_from_xy(survey_data["longitude"], 
                                                              survey_data["latitude"]),
                                                              crs=projection)
    # ---- `coast`
    coast_gdf = gpd.GeoDataFrame(coast_data, geometry="geometry", crs=projection)

    # Get appropriate plot axis-limits
    axis_limits = survey_gdf.total_bounds

    # Variable label dictionary map
    VARIABLE_MAP = {
        "number_density": {
            "name": "Mean number density",
            "units": "fish $\\mathregular{nmi^{-2}}$",
            "colormap": "inferno",
            "minimum": 0.0,
            "cbar_reverse": True,
            "size": [25, 250]
        }, 
        "biomass_density": {
            "name": "Mean biomass density",
            "units": "kg $\\mathregular{nmi^{-2}}$",
            "colormap": "plasma",
            "minimum": 0.0,
            "cbar_reverse": True,
            "size": [25, 250]
        },     
        "nasc": {
            "name": "Nautical area scattering coefficient",
            "units": "$\\mathregular{m^{2}~nmi^{-2}}$",
            "colormap": "viridis",
            "minimum": 0.0,
            "cbar_reverse": False,
            "size": [25, 250]
        },
        "max_Sv": {
            "name": "Max $\\mathregular{S_V}$",
            "units": "dB re. 1 $\\mathregular{m^-1}$",
            "colormap": "viridis",
            "minimum": -999,
            "cbar_reverse": True,
            "color_threshold": {
                "minimum": -80.0,
                "maximum": -36.0
            },
            "size": [5, 200]
        },
        # "mean_Sv": {
        #     "name": "$Mean \\mathregular{S_V}$",
        #     "units": "dB re. 1 $\\mathregular{m^-1}$",
        #     "colormap": "viridis",
        #     "minimum": -999,
        #     "cbar_reverse": True,
        #     "color_threshold": {
        #         "minimum": -80.0,
        #         "maximum": -36.0
        #     }
        # },
    }

    # List of variables to plot
    variables = list(VARIABLE_MAP.keys())

    def scale_sizes(values, min_value, max_value, min_size=25, max_size=250):

        # Censor values if needed
        sizes = values.copy()
        sizes.loc[sizes < min_value] = min_value
        sizes.loc[sizes > max_value] = max_value

        return (
            ((sizes - min_value) / (max_value - min_value))
            * (max_size - min_size) + min_size
        )    
    
    # Create a figure and a 2x2 grid of subplots
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    # Iterate through and plot all subplots
    for ax, var in zip(axes.flat, variables):
        # ---- Get the colormap
        colormap = plt.colormaps.get_cmap(VARIABLE_MAP[var]["colormap"]).resampled(256)
        # ---- Invert
        if VARIABLE_MAP[var]["cbar_reverse"]:
            newcolors = colormap(np.linspace(0, 1, 256))[::-1]
        # ---- Create the new custom colormap
        custom_cmap = ListedColormap(newcolors)
        # ---- Plot cruisetrack
        # survey_gdf.plot(ax=ax, color="dimgray", linewidth=0.25, linestyle="-")
        ax.plot(survey_gdf.geometry.x, survey_gdf.geometry.y, color="dimgray", 
                linewidth=0.25, linestyle="-")
        # ---- Drop "empty" values
        sub_gdf = survey_gdf[survey_gdf[var] > VARIABLE_MAP[var]["minimum"]]
        # ---- Assign color range
        if "color_threshold" in VARIABLE_MAP[var].keys():
            min_value = VARIABLE_MAP[var]["color_threshold"]["minimum"]
            max_value = VARIABLE_MAP[var]["color_threshold"]["maximum"]
        else:
            min_value = sub_gdf[var].min()
            max_value = sub_gdf[var].max()
        # ---- Normalize colorscale
        norm=plt.Normalize(vmin=min_value, vmax=max_value)
        # ---- Plot the points with color fills based on the variable (non-zero)
        ax.scatter(
            [geom.x for geom in sub_gdf.geometry],
            [geom.y for geom in sub_gdf.geometry],
            c=sub_gdf[var],
            s=scale_sizes(values=sub_gdf[var], 
                          min_value=min_value, 
                          max_value=max_value,
                          min_size=VARIABLE_MAP[var]["size"][0],
                          max_size=VARIABLE_MAP[var]["size"][1]),
            cmap=custom_cmap,
            norm=norm,
            edgecolor="black",
            linewidths=0.5
        )    
        # ---- Add coastline data layer
        coast_gdf.plot(ax=ax, linewidth=1.2, color='gray', edgecolor="black")
        # ---- Set axis limits
        ax.set_xlim(axis_limits[0]*1.005, axis_limits[2]*0.995)
        ax.set_ylim(axis_limits[1]*0.98, axis_limits[3]*1.005)
        # ---- Trim down the margins
        ax.margins(0,0)
        # ---- Set adjustable aspect ratio
        # ax.set_aspect('equal', adjustable='box')
        # ---- Set the title and labels
        var_info = VARIABLE_MAP[var]
        ax.set_title(f"{var_info['name']}")
        # ---- Set axis labels
        plt.xlabel(u'Longitude (\u00B0E)')
        plt.ylabel(u'Latitude (\u00B0N)')
        # ---- Add colorbar
        sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
        sm._A = []  # fake up the array of the scalar mappable
        cbar = fig.colorbar(sm, ax=ax, shrink=0.5)
        cbar.set_label(f"{var_info['units']}")
        # ---- Add scalebar
        scalebar_length = 250  # Length of scale bar in km
        scalebar_length_in_degrees = scalebar_length / 111  # Assuming 1 degree = 111 km
        # ---- Transform scale bar coordinates to axis units
        # scalebar_x = axis_limits[0]*1.005 + (axis_limits[2]*1.01 - axis_limits[0]*1.005) * 0.1
        # scalebar_y = axis_limits[1]*0.98 + (axis_limits[3]*1.005 - axis_limits[1]*0.98) * 0.1
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        x_scale = (x1 - x0) * 0.1
        y_scale = (y1 - y0) * 0.1
        # scalebar_y_offset = (axis_limits[3]*1.005 - axis_limits[1]*0.98) * 0.05
        # ---- Plot scalebar        
        # ax.plot([scalebar_x, scalebar_x + scalebar_length / 100], 
        #         [scalebar_y, scalebar_y], color='black', lw=2)
        ax.plot([x0 + x_scale, x0 + x_scale + scalebar_length_in_degrees], 
                [y0 + y_scale, y0 + y_scale], color='black', lw=2)
        # ---- Add scale text
        ax.text(x0 + x_scale + scalebar_length_in_degrees / 2, y0 + y_scale - (y1 - y0) * 0.025, 
                f'{scalebar_length} km', ha='center', va='top', color='black')

        # ax.text(scalebar_x + (scalebar_length / 200), 
        #         scalebar_y - scalebar_y_offset, 
        #         f'{scalebar_length} km', ha='center', va='bottom', color='black')

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()
