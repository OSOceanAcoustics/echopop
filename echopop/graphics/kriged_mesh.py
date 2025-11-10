from pathlib import Path
from typing import Any, Dict, Literal, Optional, Tuple, Union

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import verde as vd
import xarray as xr

from . import utils as gutils


def interpolation_mesh(
    x: pd.Series,
    y: pd.Series,
    z: pd.Series,
    crs: pyproj.CRS,
    pseudocolormesh_kwargs: Optional[Dict[str, Any]] = None,
) -> xr.DataArray:
    """
    Interpolate scattered data onto a regular mesh using Verde.

    Parameters
    ----------
    x : pd.Series
        X-coordinates of the data points.
    y : pd.Series
        Y-coordinates of the data points.
    z : pd.Series
        Values at each (x, y) location.
    crs : pyproj.CRS
        Coordinate reference system for the input data.
    pseudocolormesh_kwargs : dict, optional
        Additional keyword arguments for mesh creation (e.g., 'spacing').

    Returns
    -------
    grid : xarray.DataArray
        Masked grid of interpolated values.

    Examples
    --------
    >>> grid = interpolation_mesh(df['x'], df['y'], df['z'], pyproj.CRS("EPSG:4326"))
    >>> grid.plot.pcolormesh()

    Notes
    -----
    This function uses :mod:`verde` (`verde documentation
    <https://www.fatiando.org/verde/latest/>`_)to perform spatial interpolation and returns an
    :class:`xarray.DataArray` (`xarray.DataArray docs
    <https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html>`_).
    """

    # Assume no user-specified kwargs when not supplied
    if pseudocolormesh_kwargs is None:
        pseudocolormesh_kwargs = {}

    # Input validation and type-checking
    if not isinstance(x, pd.Series) or not isinstance(y, pd.Series) or not isinstance(z, pd.Series):
        raise TypeError("x, y, and z must all be pandas.Series objects.") from None
    if len(x) != len(y) or len(x) != len(z):
        raise ValueError("x, y, and z must all be the same length.") from None
    if len(x) == 0:
        raise ValueError("Input series x, y, and z must not be empty.") from None
    if x.isnull().any() or y.isnull().any() or z.isnull().any():
        raise ValueError("Input series x, y, and z must not contain NaN values.") from None

    # Get the geo_config settings for transformation
    # ---- Create a projection object
    projection_init = pyproj.Proj(crs)
    # ---- Transform the coordinates
    coords_init = projection_init(x, y)

    # Swap to mercator
    projection = pyproj.Proj(proj="merc", lat_ts=coords_init[1].mean())
    # ---- Retransform
    coords = projection(coords_init[0], coords_init[1])

    # Get resolution/spacing
    spacing = pseudocolormesh_kwargs.pop("spacing", 2.5 / 60.0)

    # Define processing chain
    # ---- Chain
    chain = vd.Chain(
        [
            ("reduce", vd.BlockReduce("mean", spacing * 111e3)),
            ("spline", vd.KNeighbors(k=5, reduction=np.mean)),
        ]
    )
    # ---- Fit the chain
    chain.fit(coords, z)

    # Create grid object
    # ---- Define region
    region = vd.get_region((x, y))
    # ---- Create grid
    coord_grid = chain.grid(
        region=region,
        spacing=spacing,
        projection=projection,
    )
    # ---- Mask the grid
    coord_mask = vd.distance_mask(
        data_coordinates=(x, y),
        maxdist=spacing * 111e3,
        grid=coord_grid,
        projection=projection,
    )

    # Return the masked coordinate values
    return coord_mask["scalars"]


def plot_kriged_mesh(
    data: Union[pd.DataFrame, gpd.GeoDataFrame],
    variable: str,
    projection: str = "EPSG:4326",
    coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
    plot_type: Literal["hexbin", "pcolormesh", "scatter"] = "hexbin",
    scatter_kwargs: Optional[Dict[str, Any]] = None,
    hexbin_kwargs: Optional[Dict[str, Any]] = None,
    pseudocolormesh_kwargs: Optional[Dict[str, Any]] = None,
    coast_kwargs: Optional[Dict[str, Any]] = None,
    axis_kwargs: Optional[Dict[str, Any]] = None,
    plot_kwargs: Optional[Dict[str, Any]] = None,
    colorbar_kwargs: Optional[Dict[str, Any]] = None,
    savepath: Optional[Path] = None,
    savefig_kwargs: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Plot a kriged mesh or survey data using various plot types.

    Parameters
    ----------
    data : pandas.DataFrame or geopandas.GeoDataFrame
        Input data with coordinates and the variable to plot.
    variable : str
        Name of the column to plot.
    projection : str, default='EPSG:4326'
        CRS for the plot. This should be either a projected or geodetic coordinate system definition
        compatible with with :class:`geopandas.GeoDataFrame`. See the `GeoPandas GeoDataFrame
        documentation
        <https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html>`_ for more
        details.
    coordinate_names : tuple[str], default=('longitude', 'latitude')
        Names of the coordinate columns. This is a tuple with an expected order of 'x' and then 'y'.
    plot_type : {'hexbin', 'pcolormesh', 'scatter'}, default='hexbin'
        Type of plot to produce. Options are:

        - 'hexbin': Creates a hexagonal binned plot using :func:`matplotlib.pyplot.hexbin`
          ([docs](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hexbin.html)), which
          visualizes a summary statistic defined by the user (via the `reduce_C_function` argument,
          which defaults to :func:`numpy.mean`) over a two-dimensional grid of hexagons. The
          plotted hexagonal bins can be configured using keyword arguments supplied to the
          `hexbin_kwargs` argument.
        - 'pcolormesh': Interpolates the variable onto a regular grid using
          :func:`interpolation_mesh`, which internally uses the :mod:`verde` library (`verde
          documentation <https://www.fatiando.org/verde/latest/>`_) for spatial interpolation that
          returns an :class:`xarray.DataArray`. The resulting grid is then displayed as a
          pseudocolor mesh using :meth:`xarray.DataArrray.plot.pcolormesh`. The plotted
          pseudoclor mesh can be configured using keyword arguments supplied to the
          `pseudocolormesh_kwargs` argument.
        - 'scatter': Plots each data point individually using :meth:`matplotlib.axes.Axes.scatter`
          (`matplotlib.axes.Axes.scatter docs
          <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.scatter.html>`_) coloring
          points based on the variable value. The plotted points can be configured using keyword
          arguments supplied to the `scatter_kwargs` argument.

    scatter_kwargs : dict, optional
        Additional keyword arguments passed directly to `matplotlib.pyplot.scatter` if
        `plot_type='scatter'`. For example, you can control marker size, color, alpha, etc.
        Example: `scatter_kwargs={'s': 20, 'c': 'red', 'alpha': 0.7}`
    hexbin_kwargs : dict, optional
        Additional keyword arguments passed directly to `matplotlib.pyplot.hexbin` if
        `plot_type='hexbin'`. For example, you can control grid size, colormap, etc.
        Example: `hexbin_kwargs={'gridsize': 50, 'cmap': 'viridis'}`
    pseudocolormesh_kwargs : dict, optional
        Additional keyword arguments passed to the mesh interpolation and to
        `xarray.DataArray.plot.pcolormesh` if `plot_type='pcolormesh'`. For example, you can
        specify mesh spacing, colormap, shading, etc.
        Example: `pseudocolormesh_kwargs={'spacing': 0.01, 'cmap': 'plasma'}`
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
        Additional keyword arguments passed to the main plotting function, depending on
        `plot_type`. These are merged with the specific kwargs above and can override them.
        Example: `plot_kwargs={'alpha': 0.8}`
    colorbar_kwargs : dict, optional
        Additional keyword arguments passed to `matplotlib.pyplot.colorbar` for customizing
        the colorbar. These control label, orientation, ticks, etc.
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
    >>> plot_kriged_mesh(df, 'biomass', plot_type='hexbin',
    ...                  hexbin_kwargs={'gridsize': 40, 'cmap': 'viridis'},
    ...                  coast_kwargs={'edgecolor': 'gray'},
    ...                  colorbar_kwargs={'label': 'Biomass (kg)'})

    Notes
    -----
    All keyword argument dictionaries are optional. If provided, they are passed directly to the
    underlying matplotlib or Cartopy plotting functions. If the same keyword is present in both
    a specific kwargs dict (e.g., `scatter_kwargs`) and `plot_kwargs`, the value in `plot_kwargs`
    takes precedence.
    """

    # Create copies
    scatter_kwargs = {} if scatter_kwargs is None else scatter_kwargs.copy()
    hexbin_kwargs = {} if hexbin_kwargs is None else hexbin_kwargs.copy()
    pseudocolormesh_kwargs = {} if pseudocolormesh_kwargs is None else pseudocolormesh_kwargs.copy()
    coast_kwargs = {} if coast_kwargs is None else coast_kwargs.copy()
    axis_kwargs = {} if axis_kwargs is None else axis_kwargs.copy()
    plot_kwargs = {} if plot_kwargs is None else plot_kwargs.copy()
    colorbar_kwargs = {} if colorbar_kwargs is None else colorbar_kwargs.copy()
    savefig_kwargs = {} if savefig_kwargs is None else savefig_kwargs.copy()

    # Input validation and type-checking
    if not isinstance(data, (pd.DataFrame, gpd.GeoDataFrame)):
        raise TypeError(
            "Data input must be a pandas.DataFrame or geopandas.GeoDataFrame."
        ) from None
    if variable not in data.columns:
        raise KeyError(f"The variable '{variable}' not found in data columns.") from None
    if len(data) == 0:
        raise ValueError("Input data is empty.") from None
    if not all(name in data.columns for name in coordinate_names):
        raise KeyError(f"Coordinate columns {coordinate_names} not found in data.") from None

    # Check that plotting type is valid
    if plot_type not in ["hexbin", "scatter", "pcolormesh"]:
        raise ValueError(
            f"Defined `plot_type` '{plot_type}' is invalid. Options are limited to: 'hexbin' "
            f"(default), 'scatter', and 'pcolormesh'"
        ) from None
    if variable not in data.columns:
        raise KeyError(f"Input column '{variable}' missing from the input dataset.") from None

    # Convert to GeoDataFrame, if necessary
    if isinstance(data, pd.DataFrame):
        data = gutils.dataframe_to_geodataframe(data, projection, coordinate_names)

    # Unpack coordinate names
    x, y = coordinate_names

    # Get the coastline
    _, coast_clipped, _ = gutils.get_coastline(
        gdf=data, projection=projection, coast_kwargs=coast_kwargs
    )

    # Get the `vmax`
    vmax = scatter_kwargs.pop("vmax", 10 ** np.round(np.log10(data.loc[:, variable].max())))

    # Get the `vmin`
    vmin = scatter_kwargs.pop("vmin", 0.0)

    # Get the overall data boundaries
    # ---- Compute the total survey extent
    x0, y0, x1, y1 = [
        axis_kwargs.pop(pt, (data.total_bounds * np.array([1.005, 0.995, 0.995, 1.005]))[idx])
        for pt, idx in zip(["x0", "y0", "x1", "y1"], [0, 1, 2, 3])
    ]

    # Check for axis labels
    ylabel = axis_kwargs.pop("ylabel", x)
    xlabel = axis_kwargs.pop("xlabel", y)

    # Check for colorbar label
    colorbar_label = colorbar_kwargs.pop("label", variable)

    # Get `cmap`
    cmap = colorbar_kwargs.pop("cmap", "viridis")

    # Get figure size
    figsize = plot_kwargs.pop("figsize", gutils.apply_aspect_ratio(5.5, x0, x1, y0, y1))

    # Initialize figure
    gutils.call_with_pruned(plt.figure, {"figsize": figsize, **plot_kwargs})

    # Initialize axes
    ax = gutils.call_with_pruned(plt.axes, axis_kwargs)
    # ---- Set axes
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ---- Set limits
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

    # Add the coastline
    coast_clipped.plot(ax=ax, **{**{"edgecolor": "black", "facecolor": "#C3C7C3"}, **coast_kwargs})

    # Add plotting layer
    if plot_type == "hexbin":
        # ---- Get the unique x- and y-coordinates
        xu = np.round(data.geometry.x, 1).unique()
        yu = np.round(data.geometry.y, 1).unique()
        # ---- Check against user `grid_size`
        gridsize = hexbin_kwargs.pop("gridsize", (xu.size, yu.size))
        # ---- Produce plot
        PLOT = plt.hexbin(
            x=data.geometry.x,
            y=data.geometry.y,
            C=data[variable],
            cmap=cmap,
            gridsize=gridsize,
            vmin=vmin,
            vmax=vmax,
            **hexbin_kwargs,
        )
    elif plot_type == "scatter":
        # ---- Produce plot
        PLOT = ax.scatter(
            x=data.geometry.x,
            y=data.geometry.y,
            c=data[variable],
            vmin=vmin,
            vmax=vmax,
            **scatter_kwargs,
        )
    elif plot_type == "pcolormesh":
        # ---- Prepare interpolation mesh
        grid_z = interpolation_mesh(
            data[x], data[y], data[variable], projection, pseudocolormesh_kwargs
        )
        # ---- Plot
        PLOT = grid_z.plot.pcolormesh(
            add_colorbar=False, cmap=cmap, vmin=vmin, vmax=vmax, **pseudocolormesh_kwargs
        )

    # Define the colorbar
    # ---- Check for a mappable object
    mappable = colorbar_kwargs.pop("mappable", PLOT)
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
