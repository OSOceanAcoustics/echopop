import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from echopop.workflows.nwfsc_feat import graphics as egra

# ==================================================================================================
# ==================================================================================================
# DEFINE DEMO DATA ROOT DIRECTORY
# -------------------------------

# Estabalish workflow directory
WORKFLOW_DIR = Path(os.getcwd()) / "echopop/workflow"
# ---- Validate existence
WORKFLOW_DIR.exists()

# Demo folder
DEMO_DIR = WORKFLOW_DIR / "demo"
# ---- Validate existence
DEMO_DIR.exists()

# Assign sub-folder for files
FILES_DIR = DEMO_DIR / "files"

# Assign sub-folder for output figures
FIGURES_DIR = DEMO_DIR / "figures"

# ==================================================================================================
# Load the pickled dataframes
# ---------------------------
try:
    # NASC - transect data
    df_nasc_noage1_prt = pd.read_pickle(FILES_DIR / "df_nasc_no_age1_prt.pkl")
    # Mesh - kriged data
    df_kriged_results = pd.read_pickle(FILES_DIR / "df_kriged_results.pkl")
    # Abundance table - kriged data
    df_kriged_abundance_table = pd.read_pickle(FILES_DIR / "df_kriged_abundance_table.pkl")
    # Biomass table - kriged data
    df_kriged_biomass_table = pd.read_pickle(FILES_DIR / "df_kriged_biomass_table.pkl")
    # Verbose validation upon success
    print("Pickled demo DataFrames successfully 'unpickled'.")
except Exception as e:
    raise e from None

# ==================================================================================================
# Plot and store the kriged mesh plot
# -----------------------------------
# Base kriged mesh [hexbins]
# --------------------------
# """
# -----------------
# plot_kriged_mesh(
#     ...,
#     plot_type='hexbin',
#     hexbin_kwargs,
#     coast_kwargs,
#     colorbar_kwargs,
#     axis_kwargs,
#     plot_kwargs,
#     savefig_kwargs
#     )

# - hexbin_kwargs: dict
#     Passed to matplotlib.pyplot.hexbin (if plot_type="hexbin").
#     Common options:
#         - gridsize: int or (int, int), default=100
#             Number of hexagons in the x-direction, or (nx, ny).
#         - mincnt: int, default=None
#             Minimum number of points in a hexagon for it to be colored.
#         - linewidths: float, default=None
#             Line width of hexagon edges.
#         - edgecolors: color, default=None
#             Color of hexagon edges.
#         - reduce_C_function: callable, default=np.mean
#             Function to aggregate values in each hexagon.

# - coast_kwargs: dict
#     Passed to GeoAxes.add_geometries or GeoDataFrame.plot for coastlines.
#     Common options:
#         - edgecolor: color, default='black'
#             Color of coastline edges.
#         - facecolor: color, default='none'
#             Fill color for land.
#         - linewidth: float, default=1.0
#             Width of coastline lines.
#         - alpha: float, default=1.0
#             Transparency.

# - colorbar_kwargs: dict
#     Passed to matplotlib.pyplot.colorbar.
#     Common options:
#         - cmap: str or Colormap, default=None
#             Colormap.
#         - label: str, default=None
#             Label for the colorbar.
#         - orientation: {'vertical', 'horizontal'}, default='vertical'
#             Orientation of the colorbar.
#         - shrink: float, default=1.0
#             Fraction by which to shrink the colorbar.
#         - aspect: float, default=20
#             Aspect ratio of the colorbar.

# - axis_kwargs: dict
#     Passed to matplotlib Axes set methods.
#     Common options:
#         - xlabel: str, default=None
#             X-axis label.
#         - ylabel: str, default=None
#             Y-axis label.
#         - title: str, default=None
#             Plot title.
#         - xlim: tuple, default=None
#             X-axis limits.
#         - ylim: tuple, default=None
#             Y-axis limits.

# - plot_kwargs: dict
#     Passed to matplotlib.pyplot.figure or plt.subplots.
#     Common options:
#         - figsize: tuple, default=(8, 6)
#             Figure size in inches.

# - savefig_kwargs: dict
#     Pass to `matplotlib.pyplot.savefig` if argument 'savepath' is supplied
#     Common options:
#         - dpi: float
#             Figure resolution in dots per inch.
# -----------------
# """

# Base usage
egra.plot_kriged_mesh(
    data=df_kriged_results,
    variable="biomass_density",
    savepath=FIGURES_DIR / "kriged_mesh_hexbin_base.png",
)

# 'Extended' usage
egra.plot_kriged_mesh(
    data=df_kriged_results,
    variable="biomass_density",
    savepath=FIGURES_DIR / "kriged_mesh_hexbin_plus.png",
    hexbin_kwargs={"gridsize": (30, 30)},
    coast_kwargs={"edgecolor": "black", "linewidth": 0.5},
    colorbar_kwargs={"cmap": "plasma", "label": "Biomass density ($\\mathregular{kg~nmi^{-2}}$)"},
    axis_kwargs={"xlabel": "Longitude (\u00b0E)", "ylabel": "Latitude (\u00b0N)"},
    plot_kwargs={"figsize": (8, 6)},
    savefig_kwargs={"dpi": 300},
)

# ---------------------
# Kriged mesh [scatter]
# ---------------------
# """
# -----------------
# plot_kriged_mesh(
#     ...,
#     plot_type='scatter',
#     scatter_kwargs
#     )

# - scatter_kwargs: dict
#     Passed to matplotlib.pyplot.scatter for point plotting.
#     Common options:
#         - s: float or array, default=20
#             Marker size.
#         - c: color or array, default=None
#             Marker color.
#         - alpha: float, default=1.0
#             Transparency.
#         - marker: str, default='o'
#             Marker style.

# Note: the other keyword argument dictionaries otherwise remain applicable.
# -----------------
# """

egra.plot_kriged_mesh(
    data=df_kriged_results,
    variable="biomass_density",
    plot_type="scatter",
    savepath=FIGURES_DIR / "kriged_mesh_scatter.png",
    scatter_kwargs={"s": 0.1, "alpha": 0.5},
    coast_kwargs={"edgecolor": "black", "linewidth": 0.5},
    colorbar_kwargs={"cmap": "plasma", "label": "Biomass density ($\\mathregular{kg~nmi^{-2}}$)"},
    axis_kwargs={"xlabel": "Longitude (\u00b0E)", "ylabel": "Latitude (\u00b0N)"},
    plot_kwargs={"figsize": (8, 6)},
    savefig_kwargs={"dpi": 300},
)

# -----------------------------
# Kriged mesh [pseudocolormesh]
# -----------------------------
# """
# -----------------
# plot_kriged_mesh(
#     ...,
#     plot_type='pcolormesh',
#     pseudocolormesh_kwargs
#     )

# - pcolormesh_kwargs: dict
#     Passed to xarray.DataArray.plot.pcolormesh for point plotting.
#     Common options:
#         - spacing : float, optional
#             The grid spacing (in coordinate units) for the interpolation mesh. This is used
#         - shading : {'flat', 'nearest', 'gouraud', 'auto'}, optional
#             Shading algorithm for the mesh. Default is 'auto'.
#         - alpha : float, optional
#             Transparency of the mesh (0.0 transparent through 1.0 opaque).
#         - add_labels : bool, optional
#             Whether to add axis labels automatically. Default is True.
#         - rasterized : bool, optional
#             Rasterize the mesh for vector output (e.g., PDF/SVG). Default is False.
#         - antialiased : bool, optional
#             Whether to use antialiasing. Default is True.
#         - robust : bool, optional
#             If True, use the 2nd and 98th percentiles for color limits instead of the min/max.
#         - extend : {'neither', 'both', 'min', 'max'}, optional
#             Direction(s) to extend the colorbar beyond its normal range.
#         - levels : int or array-like, optional
#             Number or positions of discrete color levels (for contour-like plots).
#         - vmin, vmax : float, optional
#             Lower and upper bounds for color scaling. If not set, determined from data.

# Note: the other keyword argument dictionaries otherwise remain applicable.
# -----------------
# """

egra.plot_kriged_mesh(
    data=df_kriged_results,
    variable="biomass_density",
    plot_type="pcolormesh",
    savepath=FIGURES_DIR / "kriged_mesh_pcolormesh.png",
    pseudocolormesh_kwargs={
        "shading": "auto",
        "alpha": 0.5,
        "add_labels": False,
        "rasterized": True,
        "levels": 10,
        "robust": True,
        "extend": "both",
    },
    coast_kwargs={"edgecolor": "black", "linewidth": 0.5},
    colorbar_kwargs={"cmap": "plasma", "label": "Biomass density ($\\mathregular{kg~nmi^{-2}}$)"},
    axis_kwargs={"xlabel": "Longitude (\u00b0E)", "ylabel": "Latitude (\u00b0N)"},
    plot_kwargs={"figsize": (8, 6)},
    savefig_kwargs={"dpi": 300},
)

# ==================================================================================================
# Plot and store the transect map
# -------------------------------

# """
# plot_transect_map(
#     ...,
#     transect_kwargs
# )

# transect_kwargs: dict
#     Passed to geopandas.GeoDataFrame.plot when plotting the transect lines.
#     Common options include:
#         - color : str
#             Line color.
#         - linewidth : float
#             Line width.

# Note: the other keyword argument dictionaries otherwise remain applicable.
# """

# Base
egra.plot_transect_map(
    data=df_nasc_noage1_prt,
    variable="biomass",
    savepath=FIGURES_DIR / "transect_map_base.png",
)

# 'Extended' usage
egra.plot_transect_map(
    data=df_nasc_noage1_prt,
    variable="nasc",
    savepath=FIGURES_DIR / "transect_map_plus.png",
    scatter_kwargs={"alpha": 0.25},
    coast_kwargs={"edgecolor": "black", "linewidth": 0.5},
    colorbar_kwargs={"cmap": "inferno", "label": "Biomass (kg)"},
    axis_kwargs={"xlabel": "Longitude (\u00b0E)", "ylabel": "Latitude (\u00b0N)"},
    plot_kwargs={"figsize": (8, 6)},
    savefig_kwargs={"dpi": 300},
)

# ==================================================================================================
# Plot and store the age-length table heatmaps
# --------------------------------------------
# """
# plot_age_length_heatmap(
#     ...,
#     imshow_kwargs
# )

# imshow_kwargs: dict
#     Passed to matplotlib.axes.Axes.imshow when plotting the heatmap.
#     Common options include:
#         - vmin, vmax : float, optional
#             Lower and upper bounds for color scaling. If not set, determined from data.

# Note: the other keyword argument dictionaries otherwise remain applicable.
# """

# Base
egra.plot_age_length_heatmap(
    data=df_kriged_abundance_table,
    savepath=FIGURES_DIR / "kriged_abundance_age_length_heatmap.png",
)

# 'Extended' usage
# ---- Futz around with the color mapping, e.g applying a log-transform
# ---- This can be added to `colorbar_kwargs`
from matplotlib.colors import SymLogNorm  # noqa: E402

# ~ Define scale
norm = SymLogNorm(linthresh=1.0, vmin=0.0, vmax=df_kriged_biomass_table.max().max())

# ~ Create colormap
scaled_colormap = plt.cm.ScalarMappable(cmap="inferno", norm=norm)

egra.plot_age_length_heatmap(
    data=df_kriged_biomass_table,
    savepath=FIGURES_DIR / "kriged_biomass_age_length_heatmap_plus.png",
    colorbar_kwargs={"mappable": scaled_colormap, "label": "Biomass (kg)"},
    imshow_kwargs={"cmap": "inferno", "norm": norm},
    axis_kwargs={
        "xlabel": r"$\mathregular{\ell}$ (cm)",
        "ylabel": r"$\mathregular{\alpha}$ (years)",
    },
    savefig_kwargs={"dpi": 300},
)

# ------------------------------
# Using the `*_filter` arguments
# ------------------------------
# """
# plot_age_length_heatmap(
#     ...,
#     include_filter,
#     exclude_filter,
#     replace_value=None
# )

# - include_filter: dict
#     Used to **include only rows** where the column matches the given value(s).
#     Example: {"age_bin": [2, 3]} will only plot data for age bins 2 and 3.

# - exclude_filter: dict
#     Used to **exclude rows** where the column matches the given value(s).
#     Example: {"age_bin": [2, 3]} will REMOVE data for age bins 2 and 3.

# - replace_value: float, default=None
#     While `replace_value=None`, values filtered out via `exclude_filter` will be completely
#     excised
#     from the dataset. However, when `replace_value` is set to a float, all values that would
#     otherwise be filtered out are instead set to the designated replacement value.
#     Example: `{"age_bin": [1]}` with `'replace_value'=0.0` would replace all values in
#     `age_bin==1.0` to `0.0` instead of removing them entirely.

#     **Pattern:** These filters are a recurrent pattern in Echopop. They allow you to flexibly
#     subset your data for plotting, making it easy to focus on specific ages, lengths, or other
#     groupings. You can use lists, single values, or callables for advanced filtering.
# """

# Drop age-1
egra.plot_age_length_heatmap(
    data=df_kriged_biomass_table,
    savepath=FIGURES_DIR / "kriged_biomass_age_length_heatmap_noage1.png",
    exclude_filter={"age_bin": [0, 1]},
    replace_value=0.0,
    colorbar_kwargs={"mappable": scaled_colormap, "label": "Age-2+ biomass (kg)"},
    imshow_kwargs={"cmap": "inferno", "norm": norm},
    axis_kwargs={
        "xlabel": r"$\mathregular{\ell}$ (cm)",
        "ylabel": r"$\mathregular{\alpha}$ (years)",
    },
    savefig_kwargs={"dpi": 300},
)

# Make sex-specific
# ---- Male
egra.plot_age_length_heatmap(
    data=df_kriged_biomass_table,
    savepath=FIGURES_DIR / "kriged_biomass_age_length_heatmap_male_noage1.png",
    include_filter={"sex": "male"},
    exclude_filter={"age_bin": [1]},
    replace_value=0.0,
    colorbar_kwargs={"mappable": scaled_colormap, "label": "Age-2+ male-specific biomass (kg)"},
    imshow_kwargs={"cmap": "inferno", "norm": norm},
    axis_kwargs={
        "xlabel": r"$\mathregular{\ell}$ (cm)",
        "ylabel": r"$\mathregular{\alpha}$ (years)",
    },
    savefig_kwargs={"dpi": 300},
)
# ---- Female
egra.plot_age_length_heatmap(
    data=df_kriged_biomass_table,
    savepath=FIGURES_DIR / "kriged_biomass_age_length_heatmap_female_noage1.png",
    include_filter={"sex": "female"},
    exclude_filter={"age_bin": [1]},
    replace_value=0.0,
    colorbar_kwargs={"mappable": scaled_colormap, "label": "Age-2+ female-specific biomass (kg)"},
    imshow_kwargs={"cmap": "inferno", "norm": norm},
    axis_kwargs={
        "xlabel": r"$\mathregular{\ell}$ (cm)",
        "ylabel": r"$\mathregular{\alpha}$ (years)",
    },
    savefig_kwargs={"dpi": 300},
)
