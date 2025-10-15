import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import xarray as xr

from echopop.graphics.age_length_heatmap import add_heatmap_grid, format_heatmap_mapping
from echopop.graphics.kriged_mesh import interpolation_mesh
from echopop.graphics.transect_map import get_transect_lines


def test_interpolation_mesh_valid():
    x = pd.Series(np.linspace(0, 1, 10))
    y = pd.Series(np.linspace(0, 1, 10))
    z = pd.Series(np.random.rand(10))
    arr = interpolation_mesh(x, y, z, pyproj.CRS("EPSG:4326"))
    assert isinstance(arr, xr.DataArray)
    assert hasattr(arr, "plot")


def test_get_transect_lines_valid():
    gdf = gpd.GeoDataFrame(
        {"transect_num": [1, 1, 2, 2], "geometry": gpd.points_from_xy([0, 1, 2, 3], [0, 1, 2, 3])}
    )
    lines = get_transect_lines(gdf)
    assert isinstance(lines, gpd.GeoDataFrame)
    assert isinstance(lines.geometry, gpd.GeoSeries)
    assert all(lines.geometry.geom_type == "LineString")


def test_add_heatmap_grid_and_format():
    fig, ax = plt.subplots()
    age_labels = np.array([1, 2, 3])
    length_labels = np.array([10, 20, 30])
    add_heatmap_grid(ax, age_labels, 1, 10, length_labels)
    plt.close()
    # No assertion needed; just check for errors


def test_format_heatmap_mapping_valid():
    fig, ax = plt.subplots()
    idx = pd.IntervalIndex.from_breaks([0, 1, 2, 3])
    cols = pd.IntervalIndex.from_breaks([10, 20, 30, 40])
    df = pd.DataFrame(np.random.rand(3, 3), index=idx, columns=cols)
    extent, da, dl, ages, lengths = format_heatmap_mapping(ax, df)
    plt.close()
    assert isinstance(extent, list)
