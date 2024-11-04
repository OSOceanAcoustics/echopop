import inspect
from typing import Any, Dict, Optional, Tuple, Union

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


def scale_sizes(values, min_value, max_value, min_size=25, max_size=250):
    """
    Scale point size
    """
    # Censor values if needed
    sizes = values.copy()
    sizes.loc[sizes < min_value] = min_value
    sizes.loc[sizes > max_value] = max_value

    return ((sizes - min_value) / (max_value - min_value)) * (max_size - min_size) + min_size


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
    cbar = plt.colorbar(
        scaled_colormap,
        ax=axes,
        shrink=0.5,
        fraction=0.075,
        pad=0.025 + x_pad,
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
    # ---- Normalize the point sizes
    data_hurdle.loc[:, "normalized_var"] = scale_sizes(
        data_hurdle.loc[:, variable], vmin, vmax, 1, 75
    )

    # Add the scatter plot
    ax.scatter(
        x=data_hurdle["longitude"],
        y=data_hurdle["latitude"],
        c=data_hurdle[variable],
        norm=norm,
        cmap=cmap,
        s=data_hurdle["normalized_var"],
        edgecolor="black",
        linewidth=0.2,
        zorder=plot_order,
    )


def plot_transect(
    dataset: pd.DataFrame,
    **kwargs,
):

    # Get the dataset variable name
    variable = kwargs.get("variable", None)
    # ---- Get the z-axis values
    z = dataset[variable].values

    # Get axis limits
    if not kwargs.get("axis_limits", None):
        # ---- Get survey bounds
        survey_bounds = get_survey_bounds(dataset, kwargs.get("geo_config", None))
        # ---- Additional buffering
        axis_limits = dict(
            x=dict(xmin=survey_bounds[0] * 1.005, xmax=survey_bounds[2] * 0.995),
            y=dict(ymin=survey_bounds[1] * 0.995, ymax=survey_bounds[3] * 1.005),
        )
    else:
        axis_limits = kwargs.get("axis_limits")

    # Add prepend label string, if required
    if "_male" in variable:
        label_prepend = "Male "
    elif "_female" in variable:
        label_prepend = "Female "
    else:
        label_prepend = ""

    # Prepare parameters for plotting
    # ---- Abundance
    if "abundance" in variable:
        cmap = kwargs.get("cmap", "viridis")
        colorbar_label = kwargs.get(
            "colorbar_label", (label_prepend + "abundance\n$#~animals$").capitalize()
        )
    # ---- Biomass density
    elif "biomass_density" in variable:
        cmap = kwargs.get("cmap", "inferno")
        colorbar_label = kwargs.get(
            "colorbar_label",
            (label_prepend + "biomass density\n$\\mathregular{kg~nmi^{-2}}$").capitalize(),
        )
    # ---- Biomass density
    elif "biomass" in variable:
        cmap = kwargs.get("cmap", "plasma")
        colorbar_label = kwargs.get(
            "colorbar_label", (label_prepend + "biomass\n$\\mathregular{kg}$").capitalize()
        )
    # ---- NASC
    elif "nasc" in variable:
        cmap = kwargs.get("cmap", "cividis")
        colorbar_label = kwargs.get("colorbar_label", "NASC\n$\\mathregular{m^{2}~nmi^{-2}}$")
    # ---- Number density
    elif "number_density" in variable:
        cmap = kwargs.get("cmap", "magma")
        colorbar_label = kwargs.get(
            "colorbar_label",
            (label_prepend + "number density\n$\\mathregular{animals~nmi^{-2}}$").capitalize(),
        )

    # Get vmin and vmax
    # ---- vmin
    vmin = kwargs.get("vmin", 0.0)
    # ---- vmax
    vmax = kwargs.get("vmax", 10 ** np.round(np.log10(z.max())))

    # Prepare default parameters
    # ---- x
    kwargs.update(dict(xlabel=kwargs.get("xlabel", "Longitude (\u00B0E)")))
    # ---- y
    kwargs.update(dict(ylabel=kwargs.get("ylabel", "Latitude (\u00B0N)")))

    # Initialize figure
    # ---- Update the 'figsize' if it doesn't yet exist
    if not kwargs.get("figsize", None):
        kwargs.update(dict(figsize=kwargs.get("figsize", apply_aspect_ratio(5.5, axis_limits))))
    # ---- Prepare figure
    plt.figure(**prune_args(plt.figure, **kwargs))
    # ---- Define GeoAxes
    ax = plt.axes(
        projection=kwargs.get("geo_config")["plot_projection"], **prune_args(plt.axes, **kwargs)
    )
    # ---- Add coastline
    ax.add_feature(kwargs.get("geo_config")["coastline"])
    # ---- Add transect lines
    add_transect_lines(dataset, plot_order=1)
    # ---- Normalize the colormapping
    colormap_norm = add_colorbar(
        ax,
        **dict(kwargs, cmap=cmap, colorbar_label=colorbar_label, vmin=vmin, vmax=vmax),
    )
    # ---- Add transect data
    add_alongtransect_data(
        ax,
        dataset,
        **dict(kwargs, cmap=cmap, norm=colormap_norm, plot_order=2, vmin=vmin, vmax=vmax),
    )
    # ---- Format the figure area axes
    format_axes(
        ax, axis_limits=axis_limits, xlabel=kwargs.get("xlabel"), ylabel=kwargs.get("ylabel")
    )
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()


def plot_mesh(
    dataset: pd.DataFrame,
    **kwargs,
):

    # Prepare units
    units = "kg"
    kriged_variable = "biomass"

    # Adjust values, if required
    dataset.loc[dataset["biomass"] < 0.0, "biomass"] = 0.0
    dataset.loc[dataset["kriged_mean"] < 0.0, "kriged_mean"] = 0.0

    # Get the dataset variable name
    variable = kwargs.get("variable", None)
    # ---- Get the x-axis values
    x = dataset["longitude"].values
    # ---- Get the y-axis values
    y = dataset["latitude"].values
    # ---- Get the z-axis values
    z = dataset[variable].values

    # Get axis limits
    if not kwargs.get("axis_limits", None):
        # ---- Get survey bounds
        survey_bounds = get_survey_bounds(dataset, kwargs.get("geo_config", None))
        # ---- Additional buffering
        axis_limits = dict(
            x=dict(xmin=survey_bounds[0] * 1.005, xmax=survey_bounds[2] * 0.995),
            y=dict(ymin=survey_bounds[1] * 0.995, ymax=survey_bounds[3] * 1.005),
        )
    else:
        axis_limits = kwargs.get("axis_limits")

    # Adjust the plotting properties
    if variable == "biomass":
        cmap = kwargs.get("cmap", "plasma")
        reduce_C_function = kwargs.get("reduce_C_function", np.sum)
        colorbar_label = kwargs.get("colorbar_label", "Kriged biomass\n$\\mathregular{kg}$")
        vmax = kwargs.get("vmax", 10 ** np.round(np.log10(z.max())))
    elif variable == "kriged_mean":
        cmap = kwargs.get("cmap", "inferno")
        reduce_C_function = kwargs.get("reduce_C_function", np.mean)
        colorbar_label = kwargs.get(
            "colorbar_label", "Kriged " + kriged_variable + f" density\n{units} " + "nmi$^{-2}$"
        )
        vmax = kwargs.get("vmax", 10 ** np.round(np.log10(z.max())))
    elif variable == "kriged_variance":
        cmap = kwargs.get("cmap", "hot")
        reduce_C_function = kwargs.get("reduce_C_function", np.mean)
        colorbar_label = kwargs.get(
            "colorbar_label",
            f"Kriged {kriged_variable} density variance" + f"\n({units} " + "nmi$^{-2})^{2}$",
        )
        vmax = kwargs.get("vmax", 10 ** np.round(np.log10(z.max()), 1))
    elif variable == "sample_cv":
        cmap = kwargs.get("cmap", "magma")
        reduce_C_function = kwargs.get("reduce_C_function", np.mean)
        colorbar_label = kwargs.get("colorbar_label", "Kriged $CV$")
        vmax = kwargs.get("vmax", np.ceil(z.max() / 0.1) * 0.1)
    elif variable == "sample_variance":
        cmap = kwargs.get("cmap", "cividis")
        reduce_C_function = kwargs.get("reduce_C_function", np.mean)
        colorbar_label = kwargs.get(
            "colorbar_label", f"Sample {kriged_variable} variance" + f"\n{units}" + "$^{-2}$"
        )
        vmax = kwargs.get("vmax", 10 ** np.round(np.log10(z.max())))

    # Get vmin
    vmin = kwargs.get("vmin", 0.0)

    # Prepare default parameters
    # ---- x
    kwargs.update(dict(xlabel=kwargs.get("xlabel", "Longitude (\u00B0E)")))
    # ---- y
    kwargs.update(dict(ylabel=kwargs.get("ylabel", "Latitude (\u00B0N)")))

    # Get the plot type
    plot_type = kwargs.get("plot_type", None)

    # Initialize figure
    # ---- Update the 'figsize' if it doesn't yet exist
    if not kwargs.get("figsize", None):
        kwargs.update(dict(figsize=kwargs.get("figsize", apply_aspect_ratio(5.5, axis_limits))))
    # ---- Prepare figure
    plt.figure(**prune_args(plt.figure, **kwargs))
    # ---- Define GeoAxes
    ax = plt.axes(
        projection=kwargs.get("geo_config")["plot_projection"], **prune_args(plt.axes, **kwargs)
    )
    # ---- Add coastline
    ax.add_feature(kwargs.get("geo_config")["coastline"])
    # ---- Normalize the colormapping
    colormap_norm = add_colorbar(
        ax, **dict(kwargs, cmap=cmap, colorbar_label=colorbar_label, vmin=vmin, vmax=vmax)
    )
    # ---- Scatter
    if plot_type == "scatter":
        ax.scatter(
            x=x,
            y=y,
            c=z,
            s=2,
            marker="s",
            cmap=cmap,
            norm=colormap_norm,
        )
    # ---- Hexbin
    elif plot_type == "hexbin":
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
            linewidth=0.05,
            cmap=cmap,
            norm=colormap_norm,
            reduce_C_function=reduce_C_function,
        )
    # ---- Pseudocolor mesh
    elif plot_type == "pcolormesh":
        # ---- Create gridded (interpolated) dataset
        grid_z = spatial_mesh(x, y, z, geo_config=kwargs.get("geo_config"))
        # ---- Plot
        grid_z.plot.pcolormesh(norm=colormap_norm, cmap=cmap, add_colorbar=False)
    # ---- Format the figure area axes
    format_axes(
        ax, axis_limits=axis_limits, xlabel=kwargs.get("xlabel"), ylabel=kwargs.get("ylabel")
    )
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()


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
    **kwargs,
):

    # Get the dataset variable name
    variable = kwargs.get("variable", None)

    # Get sex
    sex = kwargs.get("sex", None)

    # Get the correct dataset for plotting
    # ---- Abundance
    if variable == "abundance":
        dataset = data_dict["abundance"]["aged_abundance_df"].copy()
        cmap = kwargs.get("cmap", "viridis")
        colorbar_label = kwargs.get(
            "colorbar_label", (f"Abundance ({sex} fish, # animals)").capitalize()
        )
    # ---- Biomass
    elif variable == "biomass":
        dataset = data_dict["biomass"]["aged_biomass_df"].copy()
        cmap = kwargs.get("cmap", "plasma")
        colorbar_label = kwargs.get("colorbar_label", (f"Biomass ({sex} fish, kg)").capitalize())

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

    # Get limits
    # ---- vmin
    vmin = kwargs.get("vmin", 0.0)
    # ---- vmax
    vmax = kwargs.get("vmax", 10 ** np.round(np.log10(heatmap_data.stack().max())))

    # Get axis limits
    if "axis_limits" not in kwargs or kwargs.get("axis_limits") is None:
        axis_limits = dict(
            x=dict(xmin=dataset_sub["age"].min(), xmax=dataset_sub["age"].max()),
            y=dict(
                ymin=dataset_sub["length"].min() * 0.25, ymax=dataset_sub["length"].max() * 0.25
            ),
        )
    else:
        axis_limits = kwargs.get("axis_limits")

    # Prepare default parameters
    # ---- x
    kwargs.update(dict(xlabel=kwargs.get("xlabel", "Age (years)")))
    # ---- y
    kwargs.update(dict(ylabel=kwargs.get("ylabel", "Fork length (cm)")))

    # Initialize figure
    # ---- Update the 'figsize' if it doesn't yet exist
    if not kwargs.get("figsize", None):
        kwargs.update(dict(figsize=kwargs.get("figsize", apply_aspect_ratio(5.5, axis_limits))))
    # Initialize figure
    fig, ax = plt.subplots(**prune_args(plt.subplots, **kwargs))
    # ---- Format the x- and y-axis tick labels and spacing
    plot_data = format_heatmap_ticks(ax, heatmap_data)
    # ---- Normalize the colormapping
    colormap_norm = add_colorbar(
        ax,
        **dict(kwargs, cmap=cmap, colorbar_label=colorbar_label, vmin=vmin, vmax=vmax, x_pad=0.025),
    )
    # ---- Plot the heatmap
    ax.imshow(
        plot_data["heatmap_array"],
        norm=colormap_norm,
        cmap=cmap,
        interpolation=None,
        extent=plot_data["extent"],
        aspect="auto",
    )
    # ---- Add a heatmap grid, if requested
    if kwargs.get("grid_heatmap"):
        add_heatmap_grid(ax=ax, **plot_data)
    # ---- Format the figure area axes
    format_axes(
        ax, axis_limits=axis_limits, xlabel=kwargs.get("xlabel"), ylabel=kwargs.get("ylabel")
    )
    # ---- Tighten the layout and display
    plt.tight_layout()
    plt.show()


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
