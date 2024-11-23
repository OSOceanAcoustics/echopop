"""
Diagnostic data visualizations
"""

from typing import Tuple

import bokeh as bok
import geopandas as gpd
import geoviews as gv
import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn
from bokeh.models import HoverTool
from shapely.geometry import LineString

from ..survey import Survey

####################################################################################################
# Utility functions
# --------------------------------------------------------------------------------------------------


def get_discrete_color_palette(num_colors: int) -> Tuple[str, ...]:
    """
    Import appropriate discrete colormap
    """

    # Bookend values
    num_colors = 3 if num_colors < 3 else 8 if num_colors > 8 else num_colors

    # Generate color palette name
    palette_name = f"Colorblind{num_colors}"

    # Attempt to retrieve palette
    try:
        return getattr(bok.palettes, palette_name)
    except ImportError:
        raise ValueError(f"Colorblind `bokeh` palette for {num_colors} colors not found.")


####################################################################################################
# Diagnostic extension class
# --------------------------------------------------------------------------------------------------


class DiagnosticPlot:
    """
    Echopop diagnostic plots used for visually scrutinizing inputs and results

    This class includes visualizations used by the FEAT team

    Parameters
    ----------
    survey: echopop.survey.Survey
        An `echopop.survey.Survey` class object containing all of the relevant results required for
        producing various `*.xlsx`-formatted reports

    Attributes
    ----------
    data: echopop.survey.Survey
        The `echopop.survey.Survey`-class object containing all of the data used to generate the
        requested diagnostic visualizations
    tile: geoviews.elements.geo.WMTS
        A GeoViews tile element comprising an ESRI ocean basemap for plotting
    """

    # Class-level flag to track renderer initialization
    _is_bokeh_initialized = False

    # Class-level flag to track renderer initialization
    _is_theme_initialized = False

    def __init__(self, survey: Survey):

        # Evaluate whether `bokeh` has been globally initialized
        if not DiagnosticPlot._is_bokeh_initialized:
            # ---- Initialize the extension
            hv.extension("bokeh")
            # ---- Mark the class flag to prevent redundant initializations
            DiagnosticPlot._is_bokeh_initialized = True

        # Evaluate whether `bokeh.theme` has been globally initialized
        if not DiagnosticPlot._is_theme_initialized:
            # ---- Initialize the theme
            hv.renderer("bokeh").theme = bok.themes.Theme(json=self._json_theme)

        # Assign the data
        self.data = survey

        # Store the ESRI baselayer map to avoid redundant calls
        self.tile = gv.tile_sources.EsriOceanBase

    def mesh_cropping_results(self) -> pn.Column:

        # Get the datasets
        DATAINPUT = {
            "Latitude/longitude": {
                "cropped": (
                    self.data.analysis["kriging"]["mesh_df"].copy()
                    if "kriging" in self.data.analysis
                    else {}
                ),
                "full": self.data.input["statistics"]["kriging"]["mesh_df"].copy(),
            },
            "Isobath-referenced": {
                "cropped": self.data.analysis["kriging"]["mesh_df"].copy(),
            },
        }

        # Get the isobath data
        isobath_df = self.data.input["statistics"]["kriging"]["isobath_200m_df"]

        # Amend the column names, if needed
        [
            v.rename(
                columns={
                    col: (
                        "longitude"
                        if "longitude" in col
                        else "latitude" if "latitude" in col else col
                    )
                    for col in v.columns
                },
                inplace=True,
            )
            for k, v in DATAINPUT["Latitude/longitude"].items()
        ]

        # Initialize widget options
        # ---- Variable options
        var_options = ["Latitude/longitude", "Isobath-referenced"]
        # ---- Get column names
        var_columns = {
            "Latitude/longitude": ["longitude", "latitude"],
            "Isobath-referenced": ["x", "y"],
        }
        # ---- Tile options
        var_tiles = {
            "Latitude/longitude": self.tile,
            "Isobath-referenced": None,
        }
        # ---- Axis options
        var_axes = {
            "Latitude/longitude": {
                "xlab": "Longitude (\u00B0E)",
                "ylab": "Latitude (\u00B0N)",
            },
            "Isobath-referenced": {
                "xlab": "Normalized x",
                "ylab": "Normalized y",
            },
        }

        # Initialize `panel` extension
        pn.extension()

        # Create dropdown selector
        var_dropdown = pn.widgets.Select(name="Variable", options=var_options)

        # Create plotting function
        @pn.depends(selection=var_dropdown)
        def update_plot(selection):

            # Parse dataset
            data = DATAINPUT[selection]
            # ---- Get tile source
            tile = var_tiles[selection]
            # ----  Get x-axis label
            xlab = var_axes[selection]["xlab"]
            # ---- Get y-axis label
            ylab = var_axes[selection]["ylab"]
            # ---- Get coordinates
            coord = var_columns[selection]

            # Initialize plotting layers
            plot_layers = []

            # Plot the data
            # ---- Base tile layer
            if tile:
                plot_layers.append(tile)
            # ---- Full mesh
            if "full" in data:
                plot_layers.append(
                    gv.Points(data["full"], kdims=coord, label="Full mesh").opts(
                        size=8, marker="s", color="gray"
                    )
                )
            # ---- Cropped mesh
            plot_layers.append(
                gv.Points(data["cropped"], kdims=coord, label="Cropped mesh").opts(
                    size=8, marker="s", color="black"
                )
            )
            # ---- Add isobath, if applicable
            if selection == "Latitude/longitude":
                plot_layers.append(
                    gv.Path(
                        isobath_df, kdims=["longitude", "latitude"], label="200 m isobath"
                    ).opts(
                        line_dash="dashed", line_width=2, color="red", alpha=0.75, show_legend=True
                    )
                )
            # ---- Layer the data
            return hv.Overlay(plot_layers).opts(height=800, width=800, xlabel=xlab, ylabel=ylab)

        # Return the stacked Column object
        return pn.Column(var_dropdown, update_plot)

    def mesh_regions(self) -> hv.core.overlay.Overlay:

        # Get the full dataset
        DATAINPUT = self.data.analysis["kriging"]["transect_mesh_regions_df"].copy()
        # ---- Get color code for each region
        colormap = get_discrete_color_palette(len(DATAINPUT.mesh_region.unique()))
        # ---- Apply the colormap
        DATAINPUT["region_color"] = DATAINPUT["mesh_region"].apply(lambda x: colormap[x - 1])

        # Format tooltip
        TOOLTIPS = """
            <div style="position: relative;">
                <style>
                    div.bk-tooltip-content > div > div:not(:first-child) {{
                        display:none !important;
                    }}
                </style>
                <div>
                    <b>Transect: </b> <span>@transect_num</span> <br>
                </div>
            </div>
            """
        # ---- Create hovertool
        hover_tool = HoverTool(tooltips=TOOLTIPS, mode="mouse")

        # Create hovertool data
        HOVERPATH = gv.Points(
            DATAINPUT, kdims=["longitude", "latitude"], vdims=["transect_num"]
        ).opts(alpha=0.0, tools=[hover_tool])

        # Initialize plot layers
        plot_layers = []

        # Iterate and add colored regions
        for region in DATAINPUT["mesh_region"].unique():
            # ---- Format region label
            region_label = f"Region {region}"
            # ---- Subset the data
            subset = DATAINPUT.loc[DATAINPUT["mesh_region"] == region]
            # ---- Pull the color label
            id_color = subset["region_color"].iloc[0]

            #
            plot_layers.append(
                gv.Points(
                    subset,
                    kdims=["longitude", "latitude"],
                    vdims=["transect_num"],
                    label=region_label,
                ).opts(
                    color=id_color,
                    show_legend=True,
                )
            )

        # Stack the plotting layers for visualization
        # ---- Stack
        plot_stacked = (self.tile * HOVERPATH * gv.Overlay(plot_layers)).opts(
            width=800, height=800, xlabel="Longitude (\u00B0E)", ylabel="Latitude (\u00B0N)"
        )
        # ---- Return
        return plot_stacked

    def nasc_map(self) -> hv.core.overlay.Overlay:

        # Get the full dataset
        DATAINPUT = self.data.analysis["transect"]["acoustics"]["adult_transect_df"].copy()
        # ---- Convert to a GeoDataFrame
        gdf = gpd.GeoDataFrame(
            DATAINPUT, geometry=gpd.points_from_xy(DATAINPUT.longitude, DATAINPUT.latitude)
        )

        # Create `shapely` LineString objects
        lines_gdf = (
            gdf.groupby("transect_num")
            .apply(
                lambda x: LineString(x.geometry.tolist()),
                include_groups=False,
            )
            .reset_index(name="geometry")
        )
        # ---- Convert back to GeoDataFrame
        transect_gdf = gpd.GeoDataFrame(lines_gdf)

        # Aggregate the data
        DATACOPY = DATAINPUT.loc[DATAINPUT.nasc > 0].copy()
        # ---- Apply min-max scaling
        DATACOPY["nasc_scale"] = (DATACOPY.nasc - DATACOPY.nasc.min()) / (
            DATACOPY.nasc.max() - DATACOPY.nasc.min()
        )
        # ---- Set sizing
        DATACOPY["nasc_size"] = DATACOPY["nasc_scale"] * 25

        # Create TOOLTIP
        TOOLTIPS = """
            <div style="position: relative;">
                <style>
                    div.bk-tooltip-content > div > div:not(:first-child) {{
                        display:none !important;
                    }}
                </style>
                <div>
                    <b>Transect: </b> <span>@transect_num</span> <br>
                    <b>NASC: </b> <span>@nasc</span> m^2 nmi^-2 <br>
                </div>
            </div>
            """
        # ---- Create hovertool
        hover_tool = HoverTool(tooltips=TOOLTIPS, mode="mouse")

        # Create transect layer
        transect_gv = gv.Path(transect_gdf, label="Transect").opts(color="black", show_legend=True)

        # Create point layer
        points_gv = gv.Points(
            DATACOPY,
            kdims=["longitude", "latitude"],
            vdims=["transect_num", "nasc", "nasc_scale", "nasc_size"],
        ).opts(
            # ---- COLOR
            colorbar=True,
            colorbar_position="right",
            cmap="reds",
            color="nasc_scale",
            clim=(0, 1),
            colorbar_opts={
                "title": "Min-max scaled NASC [0, 1]",
                "title_text_font_style": "bold",
                "title_text_color": "black",
                "major_label_text_color": "black",
            },
            # AESTHETICS
            line_color="black",
            line_width=0.25,
            size="nasc_size",
            tools=[hover_tool],
        )

        # Stack all of the layers to complete the plot
        full_plot = (transect_gv * self.tile * points_gv).opts(
            height=800, width=800, xlabel="Longitude (\u00B0E)", ylabel="Latitude (\u00B0N)"
        )
        # ---- Return plot
        return full_plot

    def stratified_results(self) -> pn.Column:

        # Get the full datasets
        # ---- Get stratification information
        strata_info = self.data.input["spatial"]["inpfc_strata_df"].copy()
        # ---- Initialize
        DATAINPUT = {"Biomass": {}, "Biomass density": {}}
        # ---- Kriging
        if "kriging" in self.data.results:
            strat_krig = self.data.results["stratified"]["kriging"].copy()
            # ---- Add 'biomass'
            DATAINPUT["Biomass"].update(
                {
                    "kriging": pd.DataFrame(
                        {
                            "stratum": strata_info.stratum_inpfc,
                            "total": strat_krig["estimate"]["strata"]["total"] * 1e-6,
                            "lower": [k[0] * 1e-6 for k in strat_krig["ci"]["strata"]["total"]],
                            "upper": [k[1] * 1e-6 for k in strat_krig["ci"]["strata"]["total"]],
                        }
                    )
                }
            )
            # ---- Add 'biomass density'
            DATAINPUT["Biomass density"].update(
                {
                    "kriging": pd.DataFrame(
                        {
                            "stratum": strata_info.stratum_inpfc,
                            "density": strat_krig["estimate"]["strata"]["density"] * 1e-6,
                            "lower": [k[0] * 1e-6 for k in strat_krig["ci"]["strata"]["density"]],
                            "upper": [k[1] * 1e-6 for k in strat_krig["ci"]["strata"]["density"]],
                        }
                    )
                }
            )
        # ---- Transect
        if "transect" in self.data.results:
            strat_krig = self.data.results["stratified"]["transect"].copy()
            # ---- Add 'biomass'
            DATAINPUT["Biomass"].update(
                {
                    "transect": pd.DataFrame(
                        {
                            "stratum": strata_info.stratum_inpfc,
                            "total": strat_krig["estimate"]["strata"]["total"] * 1e-6,
                            "lower": [k[0] * 1e-6 for k in strat_krig["ci"]["strata"]["total"]],
                            "upper": [k[1] * 1e-6 for k in strat_krig["ci"]["strata"]["total"]],
                        }
                    )
                }
            )
            # ---- Add 'biomass density'
            DATAINPUT["Biomass density"].update(
                {
                    "transect": pd.DataFrame(
                        {
                            "stratum": strata_info.stratum_inpfc,
                            "density": strat_krig["estimate"]["strata"]["density"] * 1e-6,
                            "lower": [k[0] * 1e-6 for k in strat_krig["ci"]["strata"]["density"]],
                            "upper": [k[1] * 1e-6 for k in strat_krig["ci"]["strata"]["density"]],
                        }
                    )
                }
            )

        # Initialize `panel` extension
        pn.extension()

        # Initialize widgets
        # ---- Variable options
        var_options = ["Biomass", "Biomass density"]
        # ---- Get column names
        var_columns = {"Biomass": "total", "Biomass density": "density"}
        # ---- y-axis names
        y_axis_label = {
            "Biomass": "Biomass (kmt) [mean ± 95% CI]",
            "Biomass density": "Biomass density (kmt nmi⁻²) [mean ± 95% CI]",
        }
        # ---- Create dropdown selector
        var_dropdown = pn.widgets.Select(name="Variable", options=var_options)

        # Create plotting function
        @pn.depends(selection=var_dropdown)
        def update_plot(selection):

            # Parse dataset
            data = DATAINPUT[selection]
            # ---- Get data column
            value = var_columns[selection]

            # Initialize plotting layers
            plot_layers = []

            # Plot the data
            # ---- Kriging
            if "kriging" in data:
                # ---- Add errorbar
                kriging_pt_err = hv.Segments(
                    data["kriging"], kdims=["stratum", "lower", "stratum", "upper"], label="Kriging"
                ).opts(color="red")
                # ---- Point estimates
                kriging_pts = hv.Points(
                    data["kriging"], kdims=["stratum", value], label="Kriging"
                ).opts(
                    size=12,
                    color="red",
                    ylim=(0.0, np.inf),
                    xticks=data["kriging"]["stratum"].tolist(),
                )
                # ---- Add to layer
                plot_layers.append(kriging_pt_err * kriging_pts)
                # ---- Transect
                if "transect" in data:
                    # ---- Add errorbar
                    transect_pt_err = hv.Segments(
                        data["transect"],
                        kdims=["stratum", "lower", "stratum", "upper"],
                        label="Transect",
                    ).opts(color="black")
                    # ---- Point estimates
                    transect_pts = hv.Points(
                        data["transect"],
                        kdims=["stratum", var_columns[selection]],
                        label="Transect",
                    ).opts(
                        size=12,
                        color="black",
                    )
                    # ---- Add to layer
                    plot_layers.append(transect_pt_err * transect_pts)
                # ---- Layer the data
                return hv.Overlay(plot_layers).opts(
                    height=600, width=800, xlabel="INPFC stratum", ylabel=y_axis_label[selection]
                )

        # Return the stacked Column object
        return pn.Column(var_dropdown, update_plot)

    # Reference attributes
    _json_theme = {
        "attrs": {
            "Axis": {
                "major_label_text_font_size": "12.5pt",
                "major_label_text_color": "black",
                "axis_label_text_font_size": "16pt",
                "axis_label_text_font_style": "bold",
                "axis_label_text_color": "black",
            },
            "BaseColorBar": {
                "title_text_font_size": "16pt",
                "major_label_text_font_size": "11pt",
            },
        },
    }
