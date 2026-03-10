import geopandas as gpd
import numpy as np
import pandas as pd
import pytest

from echopop.graphics import utils as gutils


def test_scale_sizes_basic():
    """Test scale_sizes with basic ascending values."""
    arr = pd.Series(np.array([1, 2, 3]))
    result = gutils.scale_sizes(arr, 1, 3)
    assert np.allclose(result, [25, 137.5, 250])


def test_scale_sizes_equal():
    """Test scale_sizes when all values are equal."""
    arr = pd.Series(np.array([5, 5, 5]))
    result = gutils.scale_sizes(arr, 5, 5)
    assert np.all(result == 25)


def test_call_with_pruned():
    """Test call_with_pruned filters extra kwargs."""

    def foo(a, b=1):
        return a + b

    assert gutils.call_with_pruned(foo, {"a": 2, "b": 3, "c": 4}) == 5


def test_apply_aspect_ratio():
    """Test apply_aspect_ratio returns correct width and height."""
    w, h = gutils.apply_aspect_ratio(6, 0, 10, 0, 5)
    assert w == 6
    assert h == 3


def test_apply_aspect_ratio_zero_division():
    """Test apply_aspect_ratio raises ValueError when range is zero."""
    with pytest.raises(ValueError):
        gutils.apply_aspect_ratio(6, 0, 0, 0, 5)


def test_dataframe_to_geodataframe_success():
    """Test dataframe_to_geodataframe returns a GeoDataFrame."""
    df = pd.DataFrame({"lon": [0, 1], "lat": [0, 1]})
    gdf = gutils.dataframe_to_geodataframe(df, "EPSG:4326", ("lon", "lat"))
    assert isinstance(gdf, gpd.GeoDataFrame)


def test_dataframe_to_geodataframe_missing():
    """Test dataframe_to_geodataframe raises KeyError for missing columns."""
    df = pd.DataFrame({"x": [0, 1], "y": [0, 1]})
    with pytest.raises(KeyError):
        gutils.dataframe_to_geodataframe(df, "EPSG:4326", ("lon", "lat"))


def test_get_coastline_empty():
    """Test get_coastline raises ValueError for empty GeoDataFrame."""
    gdf = gpd.GeoDataFrame({"geometry": []})
    with pytest.raises(ValueError):
        gutils.get_coastline(gdf)


def test_format_geoxes_type_error():
    """Test format_geoaxes raises TypeError for invalid axes type."""

    class Dummy:
        pass

    with pytest.raises(TypeError):
        gutils.format_geoaxes(Dummy(), 0, 1, 0, 1, "x", "y")
