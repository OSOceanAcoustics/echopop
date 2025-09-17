from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple, Union

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely.geometry as sg

from . import utils as gutils


def get_transect_lines(
    data: gpd.GeoDataFrame,
):
    """
    Convert point data to LINESTRING transects grouped by 'transect_num'.

    Parameters
    ----------
    data : geopandas.GeoDataFrame
        Input data with 'geometry' and 'transect_num' columns.

    Returns
    -------
    lines_gdf : geopandas.GeoDataFrame
        GeoDataFrame with one LINESTRING per transect.

    Examples
    --------
    >>> lines = get_transect_lines(gdf)
    >>> lines.plot()

    Notes
    -----
    Each transect is constructed by sorting points within each 'transect_num' group and connecting
    them into a :class:`shapely.geometry.LineString` (`LineString docs
    <https://shapely.readthedocs.io/en/stable/manual.html#linestring>`_).
    """

    # Input validation and type-checking
    if not isinstance(data, gpd.GeoDataFrame):
        raise TypeError("Data must be a geopandas.GeoDataFrame.")
    if "transect_num" not in data.columns:
        raise KeyError("The column 'transect_num' is required for creating transect LINESTRINGs.")
    if "geometry" not in data.columns:
        raise KeyError("Column 'geometry' is required for creating transect LINESTRINGs.")
    if len(data) == 0:
        raise ValueError("Input data is empty.")

    # Sort the lines
    data_sorted = (
        data.groupby("transect_num")
        .apply(
            lambda g: g.sort_values("geometry", key=lambda s: s.apply(lambda p: (p.x, p.y))),
            include_groups=False,
        )
        .reset_index()
    )

    # Convert to lines
    lines = data_sorted.groupby("transect_num")["geometry"].apply(
        lambda pts: sg.LineString(pts.tolist()), include_groups=False
    )

    # Convert back to a full GeoDataFrame
    lines_gdf = gpd.GeoDataFrame(lines, geometry="geometry", crs=data.crs).reset_index()

    # Return
    return lines_gdf


def plot_transect_map(
    data: Union[pd.DataFrame, gpd.GeoDataFrame],
    variable: str,
    projection: str = "EPSG:4326",
    coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
    scatter_kwargs: Optional[Dict[str, Any]] = None,
    transect_kwargs: Optional[Dict[str, Any]] = None,
    coast_kwargs: Optional[Dict[str, Any]] = None,
    axis_kwargs: Optional[Dict[str, Any]] = None,
    plot_kwargs: Optional[Dict[str, Any]] = None,
    colorbar_kwargs: Optional[Dict[str, Any]] = None,
    savepath: Optional[Path] = None,
    savefig_kwargs: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Plot survey transects and variable values on a map.

    Parameters
    ----------
    data : pandas.DataFrame or geopandas.GeoDataFrame
        Input data with coordinates and variable to plot.
    variable : str
        Name of the column to plot.
    projection : str, default='EPSG:4326'
        CRS for the plot. This should be either a projected or geodetic coordinate system definition
        compatible with :class:`geopandas.GeoDataFrame`
        (see https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html).
    coordinate_names : tuple[str], default=('longitude', 'latitude')
        Names of the coordinate columns. This is a tuple with an expected order of 'x' and then 'y'.
    scatter_kwargs : dict, optional
        Additional keyword arguments passed directly to :meth:`matplotlib.axes.Axes.scatter`
        ([docs](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.scatter.html)).
        For example, you can control marker size, color, alpha, etc.
        Example: `scatter_kwargs={'s': 20, 'c': 'red', 'alpha': 0.7}`
    transect_kwargs : dict, optional
        Additional keyword arguments passed to :meth:`geopandas.GeoDataFrame.plot` for transect
        lines. For example, you can control line color, width, style, etc.
        Example: `transect_kwargs={'color': 'black', 'linewidth': 1.5}`
    coast_kwargs : dict, optional
        Additional keyword arguments passed to the coastline plotting function (e.g., Cartopy's
        `ax.add_feature`). These control the appearance of the coastline overlay.
        Example: `coast_kwargs={'edgecolor': 'black', 'linewidth': 0.5}`
    axis_kwargs : dict, optional
        Additional keyword arguments passed to axis formatting functions (e.g., setting
        axis limits, labels, or grid). These are merged with any defaults and passed to
        the axis/axes object.
        Example: `axis_kwargs={'xlabel': 'Longitude', 'ylabel': 'Latitude', 'xlim': (-130, -120)}`
    plot_kwargs : dict, optional
        Additional keyword arguments passed to the main plotting function. These can
        override values in `scatter_kwargs` or `transect_kwargs`.
        Example: `plot_kwargs={'alpha': 0.8}`
    colorbar_kwargs : dict, optional
        Additional keyword arguments passed to :func:`matplotlib.pyplot.colorbar`
        ([docs](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html))
        for customizing the colorbar. These control label, orientation, ticks, etc.
        Example: `colorbar_kwargs={'label': 'Biomass (kg)', 'orientation': 'vertical'}`
    save_path : Path, optional
        Filepath for saving the figure.
    savefig_kwargs : dict, optional
        Keyword arguments used by `matplotlib.pyplot.savefig` for saving the figure to the
        associated save filepath.

    Returns
    -------
    None

    Examples
    --------
    >>> plot_transect_map(df, 'biomass',
    ...                  scatter_kwargs={'s': 10, 'c': 'blue'},
    ...                  transect_kwargs={'color': 'black', 'linewidth': 1.5},
    ...                  coast_kwargs={'edgecolor': 'gray'},
    ...                  colorbar_kwargs={'label': 'Biomass (kg)'})

    Notes
    -----
    All keyword argument dictionaries are optional and are passed directly to the underlying
    plotting functions. If a keyword is present in both a specific kwargs dict and `plot_kwargs`,
    the value in `plot_kwargs` takes precedence.
    """

    # Create copies
    scatter_kwargs = {} if scatter_kwargs is None else scatter_kwargs.copy()
    transect_kwargs = {} if transect_kwargs is None else transect_kwargs.copy()
    coast_kwargs = {} if coast_kwargs is None else coast_kwargs.copy()
    axis_kwargs = {} if axis_kwargs is None else axis_kwargs.copy()
    plot_kwargs = {} if plot_kwargs is None else plot_kwargs.copy()
    colorbar_kwargs = {} if colorbar_kwargs is None else colorbar_kwargs.copy()
    savefig_kwargs = {} if savefig_kwargs is None else savefig_kwargs.copy()

    # Input validation and type-checking
    if not isinstance(data, (pd.DataFrame, gpd.GeoDataFrame)):
        raise TypeError("Data must be a pandas.DataFrame or geopandas.GeoDataFrame.")
    if variable not in data.columns:
        raise KeyError(f"Input column '{variable}' missing from the input dataset.")
    if len(data) == 0:
        raise ValueError("Input data is empty.")
    if not all(name in data.columns for name in coordinate_names):
        raise KeyError(f"Coordinate columns {coordinate_names} not found in data.")

    # Convert to GeoDataFrame, if necessary
    if isinstance(data, pd.DataFrame):
        data = gutils.dataframe_to_geodataframe(data, projection, coordinate_names)

    # Unpack coordinate names
    x, y = coordinate_names

    # Get the coastline
    _, coast_clipped, _ = gutils.get_coastline(
        gdf=data, projection=projection, coast_kwargs=coast_kwargs
    )

    # Subset the dataset to only include non-zero values
    data_nonzero = data.loc[data[variable] > 0.0]

    # Get the `vmax`
    vmax = scatter_kwargs.pop("vmax", 10 ** np.round(np.log10(data_nonzero.loc[:, variable].max())))

    # Get the `vmin`
    vmin = scatter_kwargs.pop("vmin", 0.0)

    # Get the overall data boundaries
    # ---- Compute the total survey extent
    x0, y0, x1, y1 = [
        axis_kwargs.pop(pt, data.total_bounds[idx])
        for pt, idx in zip(["x0", "y0", "x1", "y1"], [0, 1, 2, 3])
    ]

    # Check for axis labels
    ylabel = axis_kwargs.pop("ylabel", x)
    xlabel = axis_kwargs.pop("xlabel", y)

    # Check for colorbar label
    colorbar_label = colorbar_kwargs.pop("label", variable)

    # Get `cmap`
    cmap = colorbar_kwargs.pop("cmap", "viridis")

    # Get point scaling
    scaling = scatter_kwargs.pop("s", gutils.scale_sizes)
    if isinstance(scaling, Callable):
        s = scaling(data_nonzero.loc[:, variable], vmin, vmax, 0.1, 100)

    # Get figure size
    figsize = plot_kwargs.pop("figsize", gutils.apply_aspect_ratio(5.5, x0, x1, y0, y1))

    # Get the transect lines
    lines_gdf = get_transect_lines(data)

    # Initialize figure
    gutils.call_with_pruned(plt.figure, {"figsize": figsize, **plot_kwargs})

    # Initialize axes
    ax = gutils.call_with_pruned(plt.axes, {"extent": [x0, x1, y0, y1], **axis_kwargs})
    # ---- Set axes
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Add the coastline
    coast_clipped.plot(ax=ax, **{**{"edgecolor": "black", "facecolor": "#C3C7C3"}, **coast_kwargs})

    # Plot the transect lines
    lines_gdf.plot(ax=ax, zorder=2, **{**{"color": "#696B69", "linewidth": 0.5}, **transect_kwargs})

    # Plot the scatter points
    SCATTER_PLOT = ax.scatter(
        zorder=2,
        x=data_nonzero.geometry.x,
        y=data_nonzero.geometry.y,
        c=data_nonzero[variable],
        s=s,
        **scatter_kwargs,
    )

    # Define the colorbar
    # ---- Check for a mappable object
    mappable = colorbar_kwargs.pop("mappable", SCATTER_PLOT)
    # ---- Generate the colorbar
    plt.colorbar(mappable=mappable, ax=ax, label=colorbar_label, cmap=cmap, **colorbar_kwargs)

    # Remove margin padding
    ax.margins(0, 0)

    # Print plot with tight margins/padding
    plt.tight_layout()

    # Save?
    if savepath is not None:
        # ---- Filetyping
        if isinstance(savepath, Path):
            plt.savefig(savepath, **savefig_kwargs)
        else:
            raise TypeError(
                f"The filepath for 'savepath' must be type `pathlib.Path`, not "
                f"'{type(savepath).__name__}'."
            )

    # Show
    plt.show()
