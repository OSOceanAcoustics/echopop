from typing import Any, Dict, Optional, Union

import geopandas as gpd
import geoviews as gv
import holoviews as hv
import hvplot.pandas  # noqa: F401; import is required even though it is not directly called
import pandas as pd
import shapely as sg
from bokeh.models import HoverTool
from bokeh.themes import Theme

from .transect_map import get_transect_lines
from .utils import dataframe_to_geodataframe


class Diagnostics:
    """
    Diagnostics plotting utilities for survey and mesh data.

    Parameters
    ----------
    json_theme : dict, optional
        A Bokeh theme dictionary to apply to plots.
    """

    # Class-level flag to track renderer initialization
    _is_bokeh_initialized = False

    # Class-level flag to track renderer initialization
    _is_theme_initialized = False

    # Internal default JSON theme
    _json_theme = {
        "attrs": {
            "Title": {
                "align": "left",
                "text_font_size": "20px",
                "text_color": "black",
            },
            "Axis": {
                "axis_label_text_font_style": "bold",
                "axis_label_text_color": "black",
                "axis_label_text_font_size": "18px",
                "major_label_text_font_size": "15px",
                "major_label_text_color": "black",
            },
            "ColorBar": {
                "title_text_font_style": "normal",
                "title_text_font_size": "18px",
                "title_text_color": "black",
                "major_label_text_font_size": "16px",
                "major_label_text_color": "black",
            },
            "Legend": {
                "title_text_font_style": "bold",
                "title_text_font_size": "16px",
                "title_text_color": "black",
            },
        },
    }

    def __init__(self, json_theme: Optional[Dict[str, Any]] = None):
        """
        Initialize Diagnostics with optional Bokeh theme.

        Parameters
        ----------
        json_theme : dict, optional
            A Bokeh theme dictionary to apply to plots.
        """
        # Store the theme
        if json_theme is not None:
            self._json_theme = json_theme

        # Evaluate whether `bokeh` has been globally initialized
        if not Diagnostics._is_bokeh_initialized:
            # ---- Initialize the extension
            hv.extension("bokeh")
            # ---- Mark the class flag to prevent redundant initializations
            Diagnostics._is_bokeh_initialized = True

        # Evaluate whether `bokeh.theme` has been globally initialized
        if not Diagnostics._is_theme_initialized:
            # ---- Initialize the theme
            hv.renderer("bokeh").theme = Theme(json=self._json_theme)
            # ---- Mark the class flag
            Diagnostics._initialized = True

    def _format_transect_labels(self, transect_data: gpd.GeoDataFrame) -> gv.Labels:
        """
        Format and add transect numbers to the centroid of each transect line.

        Parameters
        ----------
        transect_data : geopandas.GeoDataFrame
            GeoDataFrame containing transect line data with 'transect_num' and 'mesh_region'.

        Returns
        -------
        geoviews.Labels
            Labels object with transect numbers at line centroids.
        """

        # Convert the point data for each line to LINESTRING objects
        transect_lines = get_transect_lines(transect_data)

        # Merge the two datasets
        transect_lines = transect_lines.merge(
            transect_data.filter(["transect_num", "mesh_region"]).drop_duplicates()
        )

        # Create copy for centroids
        centroids = transect_lines.copy()

        # Get the centroids
        centroids["centroid"] = centroids.to_crs("+proj=cea").geometry.centroid.to_crs(
            centroids.crs
        )

        # Extract the coordinates
        # ---- X/Longitude
        centroids["x"] = centroids["centroid"].x
        # ---- Y/Latitude
        centroids["y"] = centroids["centroid"].y

        # Format the labels
        # ---- Format as integers
        if (centroids["transect_num"] % 1 == 0).all():
            centroids["labels"] = centroids["transect_num"].astype(int).astype(str)
        # ---- Format as floats
        else:
            centroids["labels"] = centroids["transect_num"].astype(float).astype(str)

        # Return the geoviews.Labels object
        return gv.Labels(centroids, ["x", "y"], "labels").opts(
            text_font_size="8pt", text_color="black"
        )

    def plot_transect_mesh_regions(
        self,
        transect_data: Union[pd.DataFrame, gpd.GeoDataFrame],
        mesh_data: Union[pd.DataFrame, gpd.GeoDataFrame],
        projection: str = "epsg:4326",
    ) -> hv.Overlay:
        """
        Plot the mesh region assignment to transects.

        Parameters
        ----------
        transect_data : pandas.DataFrame or geopandas.GeoDataFrame
            Transect data with 'longitude', 'latitude', 'transect_num', and 'mesh_region'.
        mesh_data : pandas.DataFrame or geopandas.GeoDataFrame
            Mesh data with 'longitude' and 'latitude'.
        projection : str, default "epsg:4326"
            CRS projection string.

        Returns
        -------
        holoviews.Overlay
            Overlay of mesh points, transect lines, and transect labels.
        """

        # Type-checking
        # ---- Transect data
        if isinstance(transect_data, pd.DataFrame):
            transect_data = dataframe_to_geodataframe(
                transect_data, projection, ("longitude", "latitude")
            )
        # ---- Mesh data    if isinstance(transect_data, pd.DataFrame):
        if isinstance(mesh_data, pd.DataFrame):
            mesh_data = dataframe_to_geodataframe(mesh_data, projection, ("longitude", "latitude"))

        # Convert the transect data into LINESTRINGs
        transect_lines = get_transect_lines(transect_data)

        # Merge the two datasets
        transect_lines = transect_lines.merge(
            transect_data.filter(["transect_num", "mesh_region"]).drop_duplicates()
        )

        # Helper function for setting the legend title
        def set_legend_title(plot, element):
            if plot.state.legend:
                plot.state.legend.title = "Mesh region"

        # Compute the centroids of each transect and generate text labels
        transect_labels_layer = self._format_transect_labels(transect_data)

        # Plot the transect lines layer
        transect_lines_layer = transect_lines.hvplot.paths(
            geo=True, by="mesh_region", line_width=5
        ).opts(show_legend=True, hooks=[set_legend_title])

        # Plot the mesh point layer
        mesh_points_layer = mesh_data.hvplot.pandas.points(geo=True, label="Mesh points").opts(
            size=2, color="#A39D9D"
        )

        # Combine the layers
        layers = (
            gv.tile_sources.CartoLight
            * mesh_points_layer
            * transect_lines_layer
            * transect_labels_layer
        )

        # Return the layers with defined options
        return layers.opts(
            width=700,
            height=700,
            xlabel="Longitude (\u00b0E)",
            ylabel="Latitude (\u00b0N)",
            title="Transect-mesh region mapping",
            active_tools=["pan", "wheel_zoom"],
        )

    def plot_mesh_cropping(
        self,
        mesh_data: Union[pd.DataFrame, gpd.GeoDataFrame],
        cropped_mesh_data: Union[pd.DataFrame, gpd.GeoDataFrame],
        isobath_data: Union[pd.DataFrame, gpd.GeoDataFrame],
        projection: str = "epsg:4326",
    ) -> hv.Overlay:
        """
        Plot the results of the crop meshing.

        Parameters
        ----------
        mesh_data : pandas.DataFrame or geopandas.GeoDataFrame
            Original mesh data.
        cropped_mesh_data : pandas.DataFrame or geopandas.GeoDataFrame
            Cropped mesh data.
        isobath_data : pandas.DataFrame or geopandas.GeoDataFrame
            Isobath line data.
        projection : str, default "epsg:4326"
            CRS projection string.

        Returns
        -------
        holoviews.Overlay
            Overlay of mesh points, cropped mesh points, and isobath line.
        """
        # Type-checking
        # ---- Mesh data
        if isinstance(mesh_data, pd.DataFrame):
            mesh_data = dataframe_to_geodataframe(mesh_data, projection, ("longitude", "latitude"))
        # ---- Cropped mesh data
        if isinstance(cropped_mesh_data, pd.DataFrame):
            cropped_mesh_data = dataframe_to_geodataframe(
                cropped_mesh_data, projection, ("longitude", "latitude")
            )
        # ---- Isobath data
        if isinstance(isobath_data, pd.DataFrame):
            isobath_data = dataframe_to_geodataframe(
                isobath_data, projection, ("longitude", "latitude")
            )

        # Convert isobath data into a LINESTRING
        isobath_line = sg.LineString(isobath_data["geometry"].tolist())

        # Plot the mesh point layer
        mesh_points_layer = mesh_data.hvplot.points(geo=True, label="Mesh points").opts(
            size=2, color="#A39D9D"
        )

        # Plot the cropped mesh point layer
        cropped_mesh_points_layer = cropped_mesh_data.hvplot.points(
            geo=True, label="Cropped mesh points"
        ).opts(size=4, color="black")

        # Isobath layer
        isobath_line_layer = (
            gpd.GeoDataFrame({"geometry": isobath_line}, index=[0])
            .hvplot.paths(geo=True, label="200m isobath")
            .opts(line_width=2, color="red")
        )

        # Combine the layers
        layers = (
            gv.tile_sources.CartoLight
            * mesh_points_layer
            * cropped_mesh_points_layer
            * isobath_line_layer
        )

        # Return the layers with defined options
        return layers.opts(
            width=700,
            height=700,
            xlabel="Longitude (\u00b0E)",
            ylabel="Latitude (\u00b0N)",
            title="Mesh cropping",
            active_tools=["pan", "wheel_zoom"],
        )

    def plot_nasc_map(
        self,
        transect_data: Union[pd.DataFrame, gpd.GeoDataFrame],
        projection: str = "epsg:4326",
    ) -> None:
        """
        Plot NASC (Nautical Area Scattering Coefficient) values for transects.

        Parameters
        ----------
        transect_data : pandas.DataFrame or geopandas.GeoDataFrame
            Transect data with 'longitude', 'latitude', 'transect_num', and 'nasc'.
        projection : str, default "epsg:4326"
            CRS projection string.

        Returns
        -------
        holoviews.Overlay
            Overlay of transect lines and NASC points.
        """
        # Type-checking
        # ---- Transect data
        if isinstance(transect_data, pd.DataFrame):
            transect_data = dataframe_to_geodataframe(
                transect_data, projection, ("longitude", "latitude")
            )

        # Scale NASC values [minmax normalization]
        transect_data["nasc_scaled"] = (transect_data["nasc"] - transect_data["nasc"].min()) / (
            transect_data["nasc"].max() - transect_data["nasc"].min()
        )

        # Size the points accordingly
        transect_data["nasc_size"] = transect_data["nasc_scaled"] * 25 + 2

        # Convert the point data for each line to LINESTRING objects
        transect_lines = get_transect_lines(transect_data)

        # Transect lines layer
        transect_lines_layer = transect_lines.hvplot.paths(
            geo=True, line_width=2, label="Transect lines"
        ).opts(color="black", legend_position="top_left")

        # Create TOOLTIP
        TOOLTIPS = """
            <div style="position: relative;">
                <div>
                    <b>Transect: </b> <span>@transect_num</span> <br>
                    <b>NASC: </b> <span>@nasc</span> m² nmi⁻² <br>
                </div>
            </div>
            """
        # ---- Create hovertool
        hover_tool = HoverTool(
            tooltips=TOOLTIPS, mode="mouse", point_policy="snap_to_data", line_policy="nearest"
        )

        # NASC points -- zeros
        zero_nasc_points_layer = gv.Points(
            transect_data[transect_data["nasc"] == 0.0], vdims=["nasc_scaled"]
        ).opts(
            colorbar=True,
            colorbar_position="right",
            cmap="reds",
            color="nasc_scaled",
            colorbar_opts={
                "title": "Min-max scaled NASC [0, 1]",
                "title_text_font_style": "bold",
                "title_text_color": "black",
                "major_label_text_color": "black",
            },
            line_color="black",
            line_width=0.25,
            size=0.5,
        )

        # NASC points -- nonzeros
        nasc_points_layer = gv.Points(
            transect_data[transect_data["nasc"] > 0.0],
            vdims=["transect_num", "nasc", "nasc_scaled", "nasc_size"],
        ).opts(
            colorbar=True,
            colorbar_position="right",
            cmap="reds",
            color="nasc_scaled",
            colorbar_opts={
                "title": "Min-max scaled NASC [0, 1]",
                "title_text_font_style": "bold",
                "title_text_color": "black",
                "major_label_text_color": "black",
            },
            line_color="black",
            line_width=0.25,
            size="nasc_size",
            tools=[hover_tool],
        )

        # Combine the layers
        layers = (
            gv.tile_sources.CartoLight
            * transect_lines_layer
            * zero_nasc_points_layer
            * nasc_points_layer
        )

        # Return the layers with defined options
        return layers.opts(
            width=700,
            height=700,
            xlabel="Longitude (\u00b0E)",
            ylabel="Latitude (\u00b0N)",
            title="'Extreme' NASC mapping",
            active_tools=["pan", "wheel_zoom"],
        )

    def plot_stratified_results(
        self,
        stratum_results: Dict[str, pd.DataFrame],
    ) -> None:
        """
        Plot stratified biomass results with error bars and relative bias.

        Parameters
        ----------
        stratum_results : dict of str to pandas.DataFrame
            Dictionary mapping aggregation type to DataFrame with index as stratum and columns for
            'biomass' and 'cv'.

        Returns
        -------
        holoviews.Overlay
            Overlay of error bars, points, and CV hover layer.
        """
        # Ensure unique_strata is sorted for consistent plotting
        unique_strata = sorted(
            {str(item) for v in stratum_results.values() for item in v.index.tolist()}
        )

        # Create x-axis mapping for tick spacing
        stratum_to_x = {str(s): i for i, s in enumerate(unique_strata)}
        # ---- Define x-axis tick spacing
        xticks = [(v, k.capitalize() if k == "survey" else k) for k, v in stratum_to_x.items()]

        # Define the horizontal offsets
        n_types = len(stratum_results)
        x_offset = 1 / (n_types + 1)  # or 1/n_types for two types
        # ---- Map the offsets to each key
        agg_offsets = {
            agg_type: (i - (n_types - 1) / 2) * x_offset
            for i, agg_type in enumerate(stratum_results.keys())
        }

        # Iterate through the biomass values
        # ---- Initialize
        all_data = []
        # ---- Populate list
        for agg_type, df in stratum_results.items():
            temp_data = pd.DataFrame(
                {
                    "stratum": df.index,
                    "mean": df["biomass"]["mean"] * 1e-9,
                    "lower": (df["biomass"]["mean"] - df["biomass"]["low"]) * 1e-9,
                    "upper": (df["biomass"]["high"] - df["biomass"]["mean"]) * 1e-9,
                    "rbias": df["biomass"]["bias"] / df["biomass"]["mean"] * 1e2,
                    "agg_type": agg_type,
                }
            )
            # ---- Format CV
            cv_data = pd.DataFrame(
                {
                    "stratum": "cv",
                    "mean": df["cv"]["mean"],
                    "lower": df["cv"]["low"],
                    "upper": df["cv"]["high"],
                    "rbias": df["cv"]["bias"] / df["cv"]["mean"] * 1e2,
                    "agg_type": agg_type,
                }
            ).dropna()
            # ---- Add CV
            temp_data = pd.concat([temp_data, cv_data], axis=0)
            # ---- Make 'stratum' into strings
            temp_data["stratum"] = temp_data["stratum"].astype(str)
            # ---- Assign x-axis mapping
            temp_data["x_plot"] = temp_data["stratum"].map(stratum_to_x) + temp_data[
                "agg_type"
            ].map(agg_offsets)
            # ---- Add to list
            all_data.append(temp_data)

        # Concatenate into new DataFrame
        plot_data = pd.concat(all_data, ignore_index=True)

        # Format error bar layer
        errorbars_layer = (
            plot_data[plot_data["stratum"] != "cv"]
            .hvplot.errorbars(
                x="x_plot",
                y="mean",
                yerr1="lower",
                yerr2="upper",
                by="agg_type",
                label="Confidence intervals",
                color="black",
                line_width=2,
            )
            .opts(
                show_legend=True,
            )
        )

        # Format and define the stratum-specific layers
        stratum_points_layer = hv.Points(
            plot_data[~plot_data["stratum"].isin(["survey", "cv"])],
            kdims=["x_plot", "mean"],
            vdims=["rbias", "agg_type"],
        ).opts(
            color="agg_type",
            cmap="Category10",
            size=8,
            line_color="black",
            line_width=0.25,
            tools=["hover"],
            hover_tooltips=[("Relative bias", "@rbias{0.00}%")],
        )

        # Format and define the survey-specific layer
        # ---- Extract data
        survey_data = plot_data[plot_data["stratum"] == "survey"]
        # ---- Define layer
        survey_points_layer = hv.Points(
            survey_data, kdims=["x_plot", "mean"], vdims=["agg_type"]
        ).opts(
            color="agg_type",
            cmap="Category10",
            size=8,
            line_color="black",
            line_width=0.25,
        )

        # Spoof CV layer for hover tooltip
        cv_data = plot_data[plot_data["stratum"] == "cv"].copy()
        # ---- Add x- and y-coordinates
        cv_data["x"] = survey_data["x_plot"].to_numpy()
        cv_data["y"] = survey_data["mean"].to_numpy()
        # ---- Create tool
        hover_tool = HoverTool(
            tooltips="""
            <div>
            <span style="font-style: italic;">CV</span> [mean ± <i>CI</i>]:
            @mean{0.000} [@lower{0.000}, @upper{0.000}]
            </div>
            """
        )
        # ---- Define layer
        cv_points_layer = hv.Points(
            cv_data, kdims=["x", "y"], vdims=["agg_type", "lower", "mean", "upper"]
        ).opts(
            color="agg_type",
            cmap="Category10",
            size=8,
            line_color="black",
            line_width=0.25,
            tools=[hover_tool],
            alpha=0.0,
        )

        # Combine layers
        layers = errorbars_layer * stratum_points_layer * survey_points_layer * cv_points_layer

        # Return the layers with defined options
        return layers.opts(
            width=700,
            height=700,
            xticks=xticks,
            xlabel="Stratum",
            ylabel="Biomass (kmt) [mean ± CI]",
            title="Stratified biomass estimates and variability",
            active_tools=["pan", "wheel_zoom"],
            legend_position="top_left",
        )
