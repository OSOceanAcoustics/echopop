from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator

from .. import utils
from . import utils as gutils


def add_heatmap_grid(
    ax: plt.Axes,
    age_labels: np.ndarray,
    delta_age: float,
    delta_length: float,
    length_labels: np.ndarray,
) -> None:
    """
    Add grid lines to a heatmap plot for age-length data.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to which the grid lines will be added.
    age_labels : numpy.ndarray
        Array of age bin centers.
    delta_age : float
        Width of each age bin.
    delta_length : float
        Width of each length bin.
    length_labels : numpy.ndarray
        Array of length bin centers.

    Returns
    -------
    None

    Examples
    --------
    >>> add_heatmap_grid(ax, age_labels, delta_age, delta_length, length_labels)
    """

    # Input validation
    if not hasattr(ax, "hlines") or not hasattr(ax, "vlines"):
        raise TypeError("Axes object `ax` must be a matplotlib.axes.Axes instance.")
    if len(age_labels) == 0 or len(length_labels) == 0:
        raise ValueError("Age and length labels must not be empty.")
    if delta_age <= 0 or delta_length <= 0:
        raise ValueError("The Delta-age and Delta-length incremeents must be positive.")

    # Create linear offsets for grid centering
    # ---- Age
    age_offset = delta_age / 2
    # ---- Length
    length_offset = delta_length / 2

    # Create list increments for each axis
    # ---- Age
    age_increments = (age_labels - age_offset).tolist()
    # ---- Length
    length_increments = (length_labels - length_offset).tolist()

    # Add to the plot
    # ---- Vertical lines (age, x)
    ax.vlines(
        age_increments,
        ymin=length_labels.min() - length_offset,
        ymax=length_labels.max() + length_offset,
        colors="black",
    )
    # ---- Horizontal lines (length, y)
    ax.hlines(
        length_increments,
        xmin=age_labels.min() - age_offset,
        xmax=age_labels.max() + age_offset,
        colors="black",
    )


def format_heatmap_mapping(
    ax: plt.Axes, data: pd.DataFrame
) -> Tuple[List[float], float, float, pd.Series, pd.Series]:
    """
    Format axis ticks, labels, and extent for an age-length heatmap.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to format.
    data : pandas.DataFrame
        DataFrame with MultiIndex (length_bin) and MultiColumn (age_bin).

    Returns
    -------
    extent : list of float
        The extent of the heatmap for imshow.
    delta_age : float
        Width of each age bin.
    delta_length : float
        Width of each length bin.
    age_labels : pandas.Series
        Centers of age bins.
    length_labels : pandas.Series
        Centers of length bins.

    Examples
    --------
    >>> extent, delta_age, delta_length, age_labels, length_labels = format_heatmap_mapping(ax, df)
    """

    # Input validation and type-checking
    if not hasattr(ax, "set_xticks") or not hasattr(ax, "set_yticks"):
        raise TypeError("Axes object `ax` must be a matplotlib.axes.Axes instance.")
    if not isinstance(data, pd.DataFrame):
        raise TypeError("Data must be a pandas.DataFrame.")
    # if not hasattr(data.index, "mid") or not hasattr(data.columns, "mid"):
    #     raise TypeError("Data index and columns must be pandas.IntervalIndex.")

    # Extract the values for the indices and columns
    # ---- Length
    length_labels = pd.Series(data.index).apply(lambda x: x.mid).astype(float)
    # ---- Age
    age_labels = pd.Series(data.columns).apply(lambda x: x.mid).astype(float)
    # ---- Get heatmap values
    heatmap_array = data.values

    # Get the centers of all bins
    # ---- Centers list
    centers = [age_labels.min(), age_labels.max(), length_labels.max(), length_labels.min()]
    # ---- Compute change in x
    (dx,) = np.diff(centers[:2]) / (heatmap_array.shape[1] - 1)
    # ---- Compute change in y
    (dy,) = -np.diff(centers[2:]) / (heatmap_array.shape[0] - 1)
    # ---- Compute the extent
    extent = [centers[0] - dx / 2, centers[1] + dx / 2, centers[2] + dy / 2, centers[3] - dy / 2]

    # Get the bin spacing for each axis
    # ---- Age (x)
    delta_age = np.ceil(np.diff(age_labels).mean())
    # ---- Length (y)
    delta_length = np.ceil(np.diff(length_labels).mean())

    # Create the formatted tick spacings for each axis
    # ---- Age (x)
    age_axis = np.arange(
        start=age_labels.min(), stop=age_labels.max(), step=delta_age * 2, dtype=int
    )
    # ---- Length (y)
    length_axis = np.arange(
        start=length_labels.min(), stop=length_labels.max(), step=delta_length * 4
    )

    # Assign x-axis major tick labels
    x_locator = FixedLocator(age_axis)
    ax.xaxis.set_major_locator(x_locator)

    # Assign y-axis major tick labels
    y_locator = FixedLocator(length_axis)
    ax.yaxis.set_major_locator(y_locator)

    # Return the tuple
    return extent, delta_age, delta_length, age_labels, length_labels


def plot_age_length_heatmap(
    data: pd.DataFrame,
    include_filter: Dict[str, Any] = {},
    exclude_filter: Dict[str, Any] = {},
    replace_value=None,
    axis_kwargs: Optional[Dict[str, Any]] = None,
    plot_kwargs: Optional[Dict[str, Any]] = None,
    colorbar_kwargs: Optional[Dict[str, Any]] = None,
    imshow_kwargs: Optional[Dict[str, Any]] = None,
    savepath: Optional[Path] = None,
    savefig_kwargs: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Plot an age-length heatmap from a DataFrame.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame indexed by 'length_bin' and with columns 'age_bin'.
    include_filter : dict, optional
        Dictionary of filters to include specific data. Passed to :func:`utils.apply_filters`.
    exclude_filter : dict, optional
        Dictionary of filters to exclude specific data. Passed to :func:`utils.apply_filters`.
    replace_value : any, optional
        Value to use for missing or filtered data.
    axis_kwargs : dict, optional
        Additional keyword arguments for axis formatting (e.g., labels).
        Example: `axis_kwargs={'xlabel': 'Age', 'ylabel': 'Length'}`
    plot_kwargs : dict, optional
        Additional keyword arguments for :func:`matplotlib.pyplot.subplots`.
    colorbar_kwargs : dict, optional
        Additional keyword arguments for :func:`matplotlib.pyplot.colorbar`
        ([docs](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html)).
    imshow_kwargs : dict, optional
        Additional keyword arguments for :meth:`matplotlib.axes.Axes.imshow`
        ([docs](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.imshow.html)).
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
    >>> plot_age_length_heatmap(df, axis_kwargs={'xlabel': 'Age', 'ylabel': 'Length'})

    Notes
    -----
    All keyword argument dictionaries are optional and are passed directly to the underlying
    plotting functions. If a keyword is present in both a specific kwargs dict and `plot_kwargs`,
    the value in `plot_kwargs` takes precedence.
    """
    # Create copies
    axis_kwargs = {} if axis_kwargs is None else axis_kwargs.copy()
    plot_kwargs = {} if plot_kwargs is None else plot_kwargs.copy()
    colorbar_kwargs = {} if colorbar_kwargs is None else colorbar_kwargs.copy()
    imshow_kwargs = {} if imshow_kwargs is None else imshow_kwargs.copy()
    savefig_kwargs = {} if savefig_kwargs is None else savefig_kwargs.copy()

    # Input validation and type-checking
    if not isinstance(data, pd.DataFrame):
        raise TypeError("Data must be a pandas.DataFrame.")
    if data.empty:
        raise ValueError("Input data is empty.")
    # if not hasattr(data.index, "mid") or not hasattr(data.columns, "mid"):
    #     raise TypeError("Data index and columns must be pandas.IntervalIndex.")

    # Index check
    if "length_bin" not in data.index.names:
        raise IndexError(
            f"The input DataFrame is expected to be indexed by 'length_bin'. Please "
            f"make sure to set the DataFrame index before plotting. Current index/indices: "
            f"{', '.join(data.index.names)}."
        ) from None

    # Column check
    if "age_bin" not in data.columns.names:
        raise KeyError(
            f"The input DataFrame is expected to be have the column 'age_bin'. The DataFrame has "
            f"the current column(s): {', '.join(data.columns.names)}."
        ) from None

    # Filter the dataset
    data_subset = utils.apply_filters(
        data,
        include_filter=include_filter,
        exclude_filter=exclude_filter,
        replace_value=replace_value,
    )

    # Stack and isolate just the age-length table
    if len(data_subset.columns.names) > 1:
        age_length_df = (
            data_subset.stack(future_stack=True).unstack("length_bin").sum().unstack("age_bin")
        )
    else:
        age_length_df = data_subset

    # Get the `vmax`
    if "norm" in imshow_kwargs and "vmax" not in imshow_kwargs:
        vmax = None
    else:
        vmax = imshow_kwargs.pop("vmax", 10 ** np.round(np.log10(age_length_df.values.max())))

    # Get the `vmin`
    if "norm" in imshow_kwargs and "vmin" not in imshow_kwargs:
        vmin = None
    else:
        vmin = imshow_kwargs.pop("vmin", 0.0)

    # Check for axis labels
    xlabel = axis_kwargs.pop("xlabel", "Age (years)")
    ylabel = axis_kwargs.pop("ylabel", "Fork length (cm)")

    # Check for colorbar label
    colorbar_label = colorbar_kwargs.pop("label", "value")

    # Get `cmap`
    cmap = imshow_kwargs.pop("cmap", "viridis")

    # Initialize figure
    fig, ax = gutils.call_with_pruned(plt.subplots, {**plot_kwargs, **axis_kwargs})
    # ---- Set axes
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Map the heatmap grid
    extent, delta_age, delta_length, age_labels, length_labels = format_heatmap_mapping(
        ax, age_length_df
    )

    # Produce the plot
    PLOT = ax.imshow(
        X=age_length_df,
        aspect="auto",
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        **imshow_kwargs,
    )

    # Define the colorbar
    # ---- Check for a mappable object
    mappable = colorbar_kwargs.pop("mappable", PLOT)
    # ---- Generate the colorbar
    # -------- Match 'cmap' if missing
    if "cmap" not in colorbar_kwargs:
        colorbar_kwargs["cmap"] = cmap
    plt.colorbar(mappable=mappable, ax=ax, label=colorbar_label, **colorbar_kwargs)

    # Add the grid
    add_heatmap_grid(ax, age_labels, delta_age, delta_length, length_labels)

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
