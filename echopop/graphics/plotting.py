import inspect
import warnings
from typing import Any, Callable, Dict, Optional, Tuple, Union

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import verde as vd
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.ticker import FixedLocator

from .validate_plot import PlotModel


def get_survey_bounds(dataframe: pd.DataFrame, geo_config: Dict[str, Any]) -> gpd.GeoDataFrame:
    """
    Extract the survey data boundary coordinates along the x- and y-axes
    """

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


def add_colorbar(
    axes: GeoAxes,
    cmap: str,
    colorbar_label: str,
    vmin: float,
    vmax: float,
    log_base: Optional[float] = None,
    norm: Optional[Normalize] = None,
    x_pad: Optional[float] = None,
    **kwargs,
):
    """
    Add colorbar to the plot
    """

    # Adjust locator/formatter/scale
    if not norm:
        # ---- Apply log-transformation, if defined
        if log_base:
            # ---- Define scale
            norm = SymLogNorm(linthresh=1.0, vmin=vmin, vmax=vmax, base=log_base)
        else:
            # ---- Return 'None' for scale
            norm = Normalize(vmin=vmin, vmax=vmax)

    # Add pad, if needed
    x_pad = x_pad if x_pad else 0.0

    # Validate the colormap existence and apply the scaling
    try:
        # ---- Create the colormap
        scaled_colormap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    except ValueError as e:
        # ---- Drop traceback
        e.__traceback__ = None
        # ---- Raise Error
        raise (e)

    # Add pad, if needed
    x_pad = x_pad if x_pad else 0.0

    # Create the colorbar
    # ---- Prune any additional kwargs
    colorbar_pruned = {k: kwargs.pop(k) for k in prune_args(plt.colorbar, **kwargs)}
    # ---- Create colorbar
    cbar = plt.colorbar(
        scaled_colormap,
        ax=axes,
        shrink=0.5,
        fraction=0.075,
        pad=0.025 + x_pad,
        **colorbar_pruned,
    )
    # ---- Render with label
    cbar.set_label(colorbar_label)

    # Return scale
    return norm


def format_axes(
    axes: Union[Axes, GeoAxes],
    axis_limits: Dict[str, Any],
    xlabel: str,
    ylabel: str,
    **kwargs,
):
    """
    Format plotting axis labels and ticks
    """

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

    # Get the array of values for each axis
    # ---- x
    xvals = np.array(list(axis_limits["x"].values()))
    # ---- y
    yvals = np.array(list(axis_limits["y"].values()))

    # Set the axis limits
    if isinstance(axes, GeoAxes):
        axes.set_extent(
            [
                xvals.min(),
                xvals.max(),
                yvals.min(),
                yvals.max(),
            ]
        )

    # Remove margin padding
    axes.margins(0, 0)


def add_transect_lines(dataset: pd.DataFrame, plot_order: int):
    """
    Add transect lines from the dataset
    """

    # Add transect lines layer
    (
        # ---- Group by each transect
        dataset.groupby(["transect_num"])
        # ---- Plot each transect line
        .apply(
            lambda transect: plt.plot(
                transect["longitude"],
                transect["latitude"],
                c="#848884",
                linewidth=0.5,
                zorder=plot_order,
            ),
            include_groups=False,
        )
    )


def apply_aspect_ratio(figure_width: float, axis_limits: Dict[str, Any]) -> Tuple[float, float]:
    """
    Apply the defined aspect ratio on the figure window size
    """

    # Get the array of values for each axis
    # ---- x
    xvals = np.array(list(axis_limits["x"].values()))
    # ---- y
    yvals = np.array(list(axis_limits["y"].values()))

    # Compute the aspect ratio of the plotting data
    phi = (yvals.max() - yvals.min()) / (xvals.max() - xvals.min())

    # Adjust the figure width and height
    return (figure_width, figure_width * phi)


def add_alongtransect_data(
    ax: GeoAxes,
    dataset: pd.DataFrame,
    variable: str,
    vmin: float,
    vmax: float,
    cmap: str,
    norm: Union[LogNorm, Normalize],
    plot_order: int,
    s: Union[float, Callable],
    **kwargs,
):
    """
    Plot transect estimates
    """

    # Sort the data based on coordinates
    data_geo = (
        dataset.copy()
        .sort_values(by=["transect_num", "longitude", "latitude"])
        .reset_index(drop=True)
    )

    # 'Hurdle' the data (separate presence/absence and magnitude processes)
    data_hurdle = data_geo.loc[data_geo[variable] > 0.0].copy()
    # ---- Sort by the variable magnitudes
    data_hurdle = data_hurdle.sort_values(by=[variable]).reset_index(drop=True)

    # Add the scatter plot
    # ---- Format, if Callable
    if isinstance(s, Callable):
        s = s(data_hurdle.loc[:, variable], vmin, vmax)
    # ---- Prune
    scatter_pruned = {k: kwargs.pop(k) for k in prune_args(plt.scatter, **kwargs)}
    # ---- Add scatter points
    ax.scatter(
        x=data_hurdle["longitude"],
        y=data_hurdle["latitude"],
        c=data_hurdle[variable],
        norm=norm,
        cmap=cmap,
        s=s,
        zorder=plot_order,
        **scatter_pruned,
    )

    # Return pruned
    return scatter_pruned


def plot_transect(
    dataset: pd.DataFrame,
    variable: str,
    geo_config: Dict[str, Any],
    **kwargs,
):

    # Get the dataset variable name
    z = dataset[variable].values

    # Get axis limits
    if not kwargs.pop("axis_limits"):
        # ---- Get survey bounds
        survey_bounds = get_survey_bounds(dataset, geo_config)
        # ---- Additional buffering
        axis_limits = dict(
            x=dict(xmin=survey_bounds[0] * 1.005, xmax=survey_bounds[2] * 0.995),
            y=dict(ymin=survey_bounds[1] * 0.995, ymax=survey_bounds[3] * 1.005),
        )
    else:
        axis_limits = kwargs.pop("axis_limits")

    # Get `vmin` and `vmax`
    # ---- `vmax`
    vmax = kwargs.pop("vmax")
    # ---- Format `vmax`, if needed
    if isinstance(vmax, Callable):
        vmax = vmax(z)
    # ---- `vmin`
    vmin = kwargs.pop("vmin")

    # Prepare axis labels and other parameters
    # ---- x
    xlabel = kwargs.pop("xlabel") if "xlabel" in kwargs else "Longitude (\u00B0E)"
    # ---- y
    ylabel = kwargs.pop("ylabel") if "ylabel" in kwargs else "Latitude (\u00B0N)"
    # ---- colorbar
    colorbar_label = kwargs.pop("colorbar_label")
    # ---- cmap
    cmap = kwargs.pop("cmap")

    # Initialize figure
    # ---- Update the 'figsize' if it doesn't yet exist
    fig_size = (
        kwargs.pop("figsize") if "figsize" in kwargs else apply_aspect_ratio(5.5, axis_limits)
    )
    # ---- Prune the kwargs
    figure_pruned = {k: kwargs.pop(k) for k in prune_args(plt.figure, **kwargs)}
    # ---- Prepare figure
    plt.figure(**{**{"figsize": fig_size}, **figure_pruned})

    # Initialize GeoAxes
    # ---- Prune the kwargs
    geoaxes_pruned = {k: kwargs.pop(k) for k in prune_args(plt.axes, **kwargs)}
    # ---- Define GeoAxes
    ax = plt.axes(projection=geo_config["plot_projection"], **geoaxes_pruned)

    # Add coastline
    ax.add_feature(geo_config["coastline"])

    # Add transect lines
    add_transect_lines(dataset, plot_order=1)

    # Normalize the colormapping
    # ---- Get `log_base` if available
    log_base = kwargs.pop("log_base")
    # ---- Prune the kwargs
    colormap_norm_pruned = {k: kwargs.pop(k) for k in prune_args(plt.figure, **kwargs)}
    # ---- Create colormap
    colormap_norm = add_colorbar(
        ax,
        cmap=cmap,
        colorbar_label=colorbar_label,
        vmin=vmin,
        vmax=vmax,
        log_base=log_base,
        **colormap_norm_pruned,
    )

    # Add transect data
    # ---- Get size value (or function)
    s = kwargs.pop("s")
    # ---- Plot
    pruned_kwargs = add_alongtransect_data(
        ax,
        dataset,
        variable=variable,
        cmap=cmap,
        norm=colormap_norm,
        plot_order=2,
        vmin=vmin,
        vmax=vmax,
        s=s,
        **kwargs,
    )
    # ---- Update **kwargs
    kwargs = {k: v for k, v in kwargs.items() if k not in pruned_kwargs}

    # Format the figure area axes
    format_axes(ax, axis_limits=axis_limits, xlabel=xlabel, ylabel=ylabel)
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()

    # Check for unused parameters
    plot_unused_param_check(kwargs)


def plot_mesh(
    dataset: pd.DataFrame,
    variable: str,
    geo_config: Dict[str, Any],
    **kwargs,
):

    # Get the dataset variable name
    variable_col = (
        "kriged_mean"
        if variable == "biomass_density"
        else (
            "sample_cv"
            if variable == "kriged_cv"
            else "sample_variance" if variable == "local_variance" else variable
        )
    )

    # Get the x-axis values
    x = dataset["longitude"].values
    # ---- Get the y-axis values
    y = dataset["latitude"].values
    # ---- Get the z-axis values
    z = dataset[variable_col].values

    # Get axis limits
    if not kwargs.pop("axis_limits"):
        # ---- Get survey bounds
        survey_bounds = get_survey_bounds(dataset, geo_config)
        # ---- Additional buffering
        axis_limits = dict(
            x=dict(xmin=survey_bounds[0] * 1.005, xmax=survey_bounds[2] * 0.995),
            y=dict(ymin=survey_bounds[1] * 0.995, ymax=survey_bounds[3] * 1.005),
        )
    else:
        axis_limits = kwargs.pop("axis_limits")

    # Get `vmin` and `vmax`
    # ---- `vmax`
    vmax = kwargs.pop("vmax")
    # ---- Format `vmax`, if needed
    if isinstance(vmax, Callable):
        vmax = vmax(z)
    # ---- `vmin`
    vmin = kwargs.pop("vmin")

    # Prepare axis labels and other parameters
    # ---- x
    xlabel = kwargs.pop("xlabel") if "xlabel" in kwargs else "Longitude (\u00B0E)"
    # ---- y
    ylabel = kwargs.pop("ylabel") if "ylabel" in kwargs else "Latitude (\u00B0N)"
    # ---- colorbar
    colorbar_label = kwargs.pop("colorbar_label")
    # ---- cmap
    cmap = kwargs.pop("cmap")

    # Initialize figure
    # ---- Update the 'figsize' if it doesn't yet exist
    fig_size = (
        kwargs.pop("figsize") if "figsize" in kwargs else apply_aspect_ratio(5.5, axis_limits)
    )
    # ---- Prune the kwargs
    figure_pruned = {k: kwargs.pop(k) for k in prune_args(plt.figure, **kwargs)}
    # ---- Prepare figure
    plt.figure(**{**{"figsize": fig_size}, **figure_pruned})

    # Initialize GeoAxes
    # ---- Prune the kwargs
    geoaxes_pruned = {k: kwargs.pop(k) for k in prune_args(plt.axes, **kwargs)}
    # ---- Define GeoAxes
    ax = plt.axes(projection=geo_config["plot_projection"], **geoaxes_pruned)

    # Add coastline
    ax.add_feature(geo_config["coastline"])

    # Normalize the colormapping
    # ---- Get `log_base` if available
    log_base = kwargs.pop("log_base")
    # ---- Prune the kwargs
    colormap_norm_pruned = {k: kwargs.pop(k) for k in prune_args(plt.figure, **kwargs)}
    # ---- Create colormap
    colormap_norm = add_colorbar(
        ax,
        cmap=cmap,
        colorbar_label=colorbar_label,
        vmin=vmin,
        vmax=vmax,
        log_base=log_base,
        **colormap_norm_pruned,
    )

    # Plot
    # ---- Scatter
    if kwargs.get("plot_type") == "scatter":
        # ---- Prune
        scatter_pruned = {k: kwargs.pop(k) for k in prune_args(plt.scatter, **kwargs)}
        # ---- Plot
        ax.scatter(x=x, y=y, c=z, cmap=cmap, norm=colormap_norm, **scatter_pruned)
    # ---- Hexbin
    elif kwargs.get("plot_type") == "hexbin":
        # ---- Get unique x
        xg = np.unique(np.round(x, 1))
        # ---- Get unique y
        yg = np.unique(np.round(y, 1))
        # ---- Plot
        plt.hexbin(
            x=x,
            y=y,
            C=z,
            gridsize=(yg.size, xg.size),
            linewidth=kwargs.pop("linewidth"),
            cmap=cmap,
            norm=colormap_norm,
            reduce_C_function=kwargs.pop("reduce_C_function"),
        )
    # ---- Pseudocolor mesh
    elif kwargs.get("plot_type") == "pcolormesh":
        # ---- Create gridded (interpolated) dataset
        grid_z = spatial_mesh(x, y, z, geo_config=geo_config)
        # ---- Plot
        grid_z.plot.pcolormesh(
            norm=colormap_norm, cmap=cmap, add_colorbar=kwargs.pop("add_colorbar")
        )

    # Format the figure area axes
    format_axes(ax, axis_limits=axis_limits, xlabel=xlabel, ylabel=ylabel)
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()

    # Check for unused parameters
    plot_unused_param_check(kwargs)


def format_heatmap_ticks(
    ax: Axes,
    heatmap_data: pd.DataFrame,
) -> Dict[str, Any]:

    # Extract the values for the indices and columns
    # ---- Length
    length_labels = heatmap_data.index.values
    # ---- Age
    age_labels = heatmap_data.columns.values
    # ---- Get heatmap values
    heatmap_array = heatmap_data.values

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

    # Return plotting parameter/data dictionary
    return dict(
        extent=extent,
        length_labels=length_labels,
        age_labels=age_labels,
        heatmap_array=heatmap_array,
        delta_age=delta_age,
        delta_length=delta_length,
    )


def spatial_mesh(
    x: np.typing.NDArray[float],
    y: np.typing.NDArray[float],
    z: np.typing.NDArray[float],
    geo_config: Dict[str, Any],
):

    # Get the geo_config settings for transformation
    # ---- Starting
    init = geo_config.get("init")
    # ---- Create a projection object
    projection_init = pyproj.Proj(init)
    # ---- Transform the coordinates
    coords_init = projection_init(x, y)

    # Swap to mercator
    projection = pyproj.Proj(proj="merc", lat_ts=coords_init[1].mean())
    # ---- Retransform
    coords = projection(coords_init[0], coords_init[1])

    # Define processing chain
    # ---- Spacing
    spacing = 2.5 / 60
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


def add_heatmap_grid(
    ax: Axes,
    age_labels: np.ndarray,
    delta_age: float,
    delta_length: float,
    length_labels: np.ndarray,
    **kwargs,
):

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


def plot_age_length_distribution(
    data_dict: Dict[str, Any],
    variable: str,
    sex: str,
    grid_heatmap: bool,
    **kwargs,
):

    # Get the correct dataset for plotting
    dataset = data_dict[variable][f"aged_{variable}_df"].copy()

    # Prepare the dataset
    # ---- Stack sum
    dataset_stk = dataset.stack(future_stack=True).sum(axis=1).reset_index(name=variable)
    # ---- Convert the length and age columns from intervals into numerics
    # -------- Length
    dataset_stk["length"] = dataset_stk["length_bin"].apply(lambda x: x.mid).astype(float)
    # -------- Age
    dataset_stk["age"] = dataset_stk["age_bin"].apply(lambda x: x.mid).astype(float)
    # ---- Subset the data based on the sex
    dataset_sub = dataset_stk.loc[dataset_stk["sex"] == sex, :]

    # Convert into a heatmap
    heatmap_data = dataset_sub.pivot_table(
        values=variable, columns="age", index="length", aggfunc="sum"
    )

    # Get axis limits
    if not kwargs.pop("axis_limits"):
        axis_limits = dict(
            x=dict(xmin=dataset_sub["age"].min(), xmax=dataset_sub["age"].max()),
            y=dict(
                ymin=dataset_sub["length"].min() * 0.25, ymax=dataset_sub["length"].max() * 0.25
            ),
        )
    else:
        axis_limits = kwargs.pop("axis_limits")

    # Get `vmin` and `vmax`
    # ---- `vmax`
    vmax = kwargs.pop("vmax")
    # ---- Format `vmax`, if needed
    if isinstance(vmax, Callable):
        vmax = vmax(heatmap_data.stack())
    # ---- `vmin`
    vmin = kwargs.pop("vmin")

    # Prepare axis labels and other parameters
    # ---- x
    xlabel = kwargs.pop("xlabel") if "xlabel" in kwargs else "Age (years)"
    # ---- y
    ylabel = kwargs.pop("ylabel") if "ylabel" in kwargs else "Fork length (cm)"
    # ---- colorbar
    colorbar_label = kwargs.pop("colorbar_label")
    # -------- Reformat
    colorbar_label = colorbar_label.replace("\n", f" ({sex} fish, ") + ")"
    # ---- cmap
    cmap = kwargs.pop("cmap")

    # Initialize figure
    # ---- Update the 'figsize' if it doesn't yet exist
    fig_size = (
        kwargs.pop("figsize") if "figsize" in kwargs else apply_aspect_ratio(5.5, axis_limits)
    )
    # ---- Prune the kwargs
    figure_pruned = {k: kwargs.pop(k) for k in prune_args(plt.subplots, **kwargs)}
    # ---- Prepare figure
    fig, ax = plt.subplots(**{**{"figsize": fig_size}, **figure_pruned})

    # Format the heatmap tick spacing
    plot_data = format_heatmap_ticks(ax, heatmap_data)

    # Normalize the colormapping
    # ---- Get `log_base` if available
    log_base = kwargs.pop("log_base")
    # ---- Prune the kwargs
    colormap_norm_pruned = {k: kwargs.pop(k) for k in prune_args(plt.figure, **kwargs)}
    # ---- Create colormap
    colormap_norm = add_colorbar(
        ax,
        cmap=cmap,
        colorbar_label=colorbar_label,
        vmin=vmin,
        vmax=vmax,
        log_base=log_base,
        **colormap_norm_pruned,
    )

    # Plot the heatmap
    # ---- Prune the kwargs
    imshow_pruned = {k: kwargs.pop(k) for k in prune_args(plt.imshow, **kwargs)}
    # ---- Plot
    ax.imshow(
        plot_data["heatmap_array"],
        norm=colormap_norm,
        cmap=cmap,
        extent=plot_data["extent"],
        **imshow_pruned,
    )

    # Overlay a grid
    if grid_heatmap:
        add_heatmap_grid(ax=ax, **plot_data)

    # Format the figure area axes
    format_axes(ax, axis_limits=axis_limits, xlabel=xlabel, ylabel=ylabel)
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()

    # Check for unused parameters
    plot_unused_param_check(kwargs)


def PLOT_MAP(cls, kind: str):

    # Reference dictionary
    PLOT_MAP = {
        "age_length_distribution": {
            "data": lambda v: v.analysis["transect"]["biology"]["population"]["tables"],
            "function": plot_age_length_distribution,
            "type": "biological",
        },
        "mesh": {
            "data": lambda v: v.results["kriging"]["mesh_results_df"].copy(),
            "function": plot_mesh,
            "type": "spatial",
        },
        "transect": {
            "data": lambda v: v.analysis["transect"]["acoustics"]["adult_transect_df"].copy(),
            "function": plot_transect,
            "type": "spatial",
        },
    }

    # Get the associated plotting parameterization/information
    plotting_info = PLOT_MAP.get(kind)
    # ---- Fill the data key
    plotting_info.update({"data": plotting_info["data"](cls)})

    # Return the dictionary
    return plotting_info


def prune_args(func, **kwargs):
    return {k: v for k, v in kwargs.items() if k in dict(inspect.signature(func).parameters)}


def validate_plot_args(**kwargs):
    """
    Validate the plotting arguments.
    """

    return PlotModel.create(**kwargs)


def plot_unused_param_check(kwargs):
    """
    Check for unused parameters
    """

    if set(kwargs) != set(["kind", "plot_type"]):
        # ---- Get differing parameters
        mismatch = list(set(kwargs) - set(["kind", "plot_type"]))
        # ---- Join as a list
        mismatch_str = ", ".join(f"'{param}'" for param in mismatch)
        # ---- Print unused parameters (as a warning)
        warnings.warn(f"The following plotting parameters were NOT used: {mismatch_str}!")
