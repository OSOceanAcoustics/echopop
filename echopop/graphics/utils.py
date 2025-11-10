import inspect
import warnings
from typing import Any, Dict, Optional, Tuple

import cartopy.feature as cfeature
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry as sg
from cartopy.mpl.geoaxes import GeoAxes


def call_with_pruned(func, kwargs, /, **overrides):
    """
    Call a function with only the keyword arguments it accepts.

    Parameters
    ----------
    func : callable
        The function to call.
    kwargs : dict
        Dictionary of keyword arguments. Only those accepted by `func` will be used.
    overrides : dict
        Additional keyword arguments to override or add to `kwargs`. These take precedence.

    Returns
    -------
    result
        The result of calling `func` with the accepted keyword arguments.

    Examples
    --------
    >>> def foo(a, b=1): return a + b
    >>> call_with_pruned(foo, {'a': 2, 'b': 3, 'c': 4})
    5

    Notes
    -----
    This utility is useful when you have a dictionary of parameters but only want to pass those
    that are valid for a given function. Any keys in `overrides` will override those in `kwargs`.
    """

    # Get call signature
    sig = inspect.signature(func)

    # Prune
    pruned = {k: kwargs.pop(k) for k in list(kwargs) if k in sig.parameters}

    # Warn about pruned keys
    pruned_keys = set(kwargs) - set(pruned)
    if pruned_keys:
        warnings.warn(
            f"The following keys were not accepted by '{func.__name__}' and were pruned: "
            f"{sorted(pruned_keys)}",
            UserWarning,
            stacklevel=2,
        )

    # Return the pruned keyword arguments
    return func(**{**pruned, **overrides})


def get_coastline(
    gdf, resolution="10m", projection="epsg:4326", coast_kwargs: Optional[Dict[str, Any]] = None
):
    """
    Get and clip coastline data from cartopy to the bounds of a GeoDataFrame.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Boundary GeoDataFrame. Used to determine the region for coastline clipping.
    resolution : str, default='10m'
        Cartopy resolution ('10m', '50m', or '110m').
    projection : str, default="epsg:4326"
        CRS projection string for the returned coastline.
    coast_kwargs : dict, optional
        Additional keyword arguments passed to :class:`cartopy.feature.NaturalEarthFeature`
        ([docs](https://scitools.org.uk/cartopy/docs/latest/reference/generated/
        cartopy.feature.NaturalEarthFeature.html)) for coastline appearance (e.g., `edgecolor`,
        `facecolor`, `linewidth`). These control how the coastline is rendered when added to a plot.

    Returns
    -------
    tuple
        (full_coast, clipped_coast_original, boundary_box_unbuffered_gdf)

    Examples
    --------
    >>> coast, coast_clip, bbox = get_coastline(gdf, coast_kwargs={'edgecolor': 'black'})

    Notes
    -----
    The `coast_kwargs` dictionary is passed directly to Cartopy's feature creation and can be
    used to customize the appearance of the coastline overlay.
    """

    # Update kwargs if none are supplied
    if coast_kwargs is None:
        coast_kwargs = {}

    # Input validation and type-checking
    if not hasattr(gdf, "geometry"):
        raise TypeError("The input dataset for `gdf` must be a `GeoDataFrame`.")
    if gdf.empty:
        raise ValueError("Input GeoDataFrame is empty.")

    # Get original boundaries
    xmin0, ymin0, xmax0, ymax0 = gdf.total_bounds

    # Create boundary boxes
    boundary_box_unbuffered = sg.box(xmin0, ymin0, xmax0, ymax0)
    boundary_box_unbuffered_gdf = gpd.GeoDataFrame(
        geometry=[boundary_box_unbuffered], crs=projection
    )

    # Get coastline from cartopy
    coast_feature = call_with_pruned(
        cfeature.NaturalEarthFeature,
        {
            **{
                "category": "physical",
                "name": "land",
                "scale": resolution,
                "edgecolor": "black",
                "facecolor": "none",
            },
            **coast_kwargs,
        },
    )
    coast_geometries = list(coast_feature.geometries())
    full_coast = gpd.GeoDataFrame(geometry=coast_geometries, crs=projection)

    # Clip the coastline
    clipped_coast_original = gpd.clip(
        full_coast, sg.box(xmin0 - 1, ymin0 - 1, xmax0 + 1, ymax0 + 1)
    )

    return full_coast, clipped_coast_original, boundary_box_unbuffered_gdf


def format_geoaxes(
    axes: GeoAxes,
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    xlabel: str,
    ylabel: str,
    **kwargs,
):
    """
    Format axis labels, ticks, and extent for a GeoAxes.

    Parameters
    ----------
    axes : cartopy.mpl.geoaxes.GeoAxes
        The axes to format.
    x0, x1, y0, y1 : float
        Axis limits.
    xlabel, ylabel : str
        Axis labels.
    kwargs : dict
        Additional keyword arguments passed to axis formatting.

    Returns
    -------
    None

    Examples
    --------
    >>> format_geoaxes(ax, 0, 10, 0, 5, 'Longitude', 'Latitude')
    """

    # Input validation and type-checking
    if not hasattr(axes, "set_extent"):
        raise TypeError("Axes must be a GeoAxes instance.")
    if x1 == x0 or y1 == y0:
        raise ValueError("Axis limits must define a non-zero range.")

    # Set axis tick marks/spacing
    # ---- x
    axes.set_xticks(axes.get_xticks())
    # ---- y
    axes.set_yticks(axes.get_yticks())

    # Set labels
    # ---- x
    axes.set_xlabel(xlabel)
    # ---- y
    axes.set_ylabel(ylabel)

    # Set the axis limits
    axes.set_extent([x0, x1, y0, y1])

    # Remove margin padding
    axes.margins(0, 0)


def scale_sizes(
    values: pd.Series,
    min_value: float,
    max_value: float,
    min_size: float = 25,
    max_size: float = 250,
):
    """
    Scale values to a range of point sizes for plotting.

    Parameters
    ----------
    values : pandas.Series
        Values to scale.
    min_value : float
        Minimum value for scaling.
    max_value : float
        Maximum value for scaling.
    min_size : float, optional
        Minimum point size (default 25).
    max_size : float, optional
        Maximum point size (default 250).

    Returns
    -------
    sizes : numpy.ndarray
        Scaled point sizes.

    Examples
    --------
    >>> scale_sizes(np.array([1, 2, 3]), 1, 3)
    array([ 25., 137.5, 250.])

    Notes
    -----
    The function linearly rescales `values` to the range [`min_size`, `max_size`]. If all values
    are equal, all points will be assigned `min_size`.
    """

    # Create copy
    sizes = values.copy()

    # Censor values if needed
    sizes.loc[sizes < min_value] = min_value
    sizes.loc[sizes > max_value] = max_value

    # If the values are equal
    if max_value == min_value:
        return np.full_like(values, min_size, dtype=float)

    return ((sizes - min_value) / (max_value - min_value)) * (max_size - min_size) + min_size


def apply_aspect_ratio(
    figure_width: float, x0: float, x1: float, y0: float, y1: float
) -> Tuple[float, float]:
    """
    Compute figure size to maintain aspect ratio.

    Parameters
    ----------
    figure_width : float
        Desired figure width.
    x0, x1, y0, y1 : float
        Data bounds.

    Returns
    -------
    (width, height) : tuple of float
        Figure size.

    Examples
    --------
    >>> apply_aspect_ratio(6, 0, 10, 0, 5)
    (6, 3)
    """

    # Validate inputs
    if x1 == x0:
        raise ValueError("Variables `x1` and `x0` must not be equal for aspect ratio calculation.")

    # Compute the aspect ratio of the plotting data
    phi = abs((y1 - y0) / (x1 - x0))

    # Adjust the figure width and height
    return (figure_width, figure_width * phi)


def dataframe_to_geodataframe(
    data: pd.DataFrame, projection: str, coordinate_names: Tuple[str, str], **kwargs
) -> gpd.GeoDataFrame:
    """
    Convert a DataFrame to a GeoDataFrame using coordinate columns.

    Parameters
    ----------
    data : pandas.DataFrame
        Input DataFrame.
    projection : str
        CRS projection string.
    coordinate_names : tuple of str
        Names of the coordinate columns.
    kwargs : dict
        Additional keyword arguments for :class:`geopandas.GeoDataFrame`.

    Returns
    -------
    data_gdf : geopandas.GeoDataFrame
        GeoDataFrame with geometry column.

    Examples
    --------
    >>> gdf = dataframe_to_geodataframe(df, "EPSG:4326", ("lon", "lat"))
    """

    if not all(name in data.columns for name in coordinate_names):
        raise KeyError(f"Coordinate columns {coordinate_names} not found in data.")

    # Unpack coordinate names
    x, y = coordinate_names

    # Convert to GeoDataFrame
    data_gdf = gpd.GeoDataFrame(
        data=data, geometry=gpd.points_from_xy(x=data[x], y=data[y]), crs=projection, **kwargs
    )

    return data_gdf
