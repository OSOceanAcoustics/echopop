import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point


@pytest.fixture
def sample_transect_data_cropping():
    """Sample transect data for testing cropping functions."""
    return pd.DataFrame(
        {
            "transect_num": [
                1,
                1,
                1,
                2,
                2,
                2,
                3,
                3,
                3,
                121,
                121,
                121,
                122,
                122,
                122,
                129,
                129,
                129,
                130,
                130,
                130,
            ],
            "longitude": [
                -125.0,
                -125.1,
                -125.2,
                -125.5,
                -125.6,
                -125.7,
                -125.0,
                -125.1,
                -125.2,
                -126.0,
                -126.1,
                -126.2,
                -126.5,
                -126.6,
                -126.7,
                -127.0,
                -127.1,
                -127.2,
                -127.5,
                -127.6,
                -127.7,
            ],
            "latitude": [
                48.0,
                48.1,
                48.2,
                48.0,
                48.1,
                48.2,
                49.0,
                49.1,
                49.2,
                50.0,
                50.1,
                50.2,
                50.0,
                50.1,
                50.2,
                51.0,
                51.1,
                51.2,
                51.0,
                51.1,
                51.2,
            ],
            "interval": [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
            "max_depth": [
                100,
                90,
                80,
                110,
                95,
                85,
                105,
                88,
                92,
                120,
                100,
                90,
                115,
                98,
                88,
                125,
                105,
                95,
                130,
                110,
                100,
            ],
            "layer_depth_min": [
                15,
                25,
                35,
                15,
                25,
                35,
                15,
                25,
                35,
                15,
                25,
                35,
                15,
                25,
                35,
                15,
                25,
                35,
                15,
                25,
                35,
            ],
            "layer_depth_max": [
                45,
                55,
                65,
                45,
                55,
                65,
                45,
                55,
                65,
                45,
                55,
                65,
                45,
                55,
                65,
                45,
                55,
                65,
                45,
                55,
                65,
            ],
            "biomass_density": [
                10.5,
                12.3,
                8.7,
                15.2,
                11.1,
                9.8,
                13.4,
                16.7,
                14.2,
                18.1,
                12.5,
                10.3,
                14.8,
                11.9,
                13.6,
                16.2,
                14.5,
                12.8,
                17.3,
                15.1,
                13.9,
            ],
        }
    )


@pytest.fixture
def sample_mesh_data_cropping():
    """Sample mesh data for testing cropping functions."""
    # Create a regular grid of mesh points
    lons = np.linspace(-128.0, -124.0, 20)
    lats = np.linspace(47.0, 52.0, 25)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    return pd.DataFrame(
        {
            "et_id": range(len(lon_grid.flatten())),
            "longitude": lon_grid.flatten(),
            "latitude": lat_grid.flatten(),
            "area (km^2)": np.random.uniform(0.5, 2.0, len(lon_grid.flatten())),
            "x": np.random.uniform(-0.5, 0.5, len(lon_grid.flatten())),
            "y": np.random.uniform(-0.5, 0.5, len(lon_grid.flatten())),
            "transect_num": np.random.choice(
                [1, 2, 3, 121, 122, 129, 130], len(lon_grid.flatten())
            ),
            "area": 10.0,
        }
    )


@pytest.fixture
def sample_region_function():
    """Sample region function for testing transect_ends_crop."""

    def transect_mesh_region_sample(region):
        """Sample implementation of transect mesh region function."""
        if region == 1:
            return 1, 3, [1.1, 2.1, 3.1], [1.4, 2.4, 3.4]
        elif region == 2:
            return 121, 122, [121.6, 122.6], [121.9, 122.9]
        else:  # region == 3
            return 129, 130, [129.1, 130.1], [129.4, 130.4]

    return transect_mesh_region_sample


@pytest.fixture
def sample_transect_geodataframe():
    """Sample transect GeoDataFrame for testing."""
    data = pd.DataFrame(
        {
            "transect_num": [1, 1, 1, 2, 2, 2, 3, 3, 3],
            "longitude": [-125.0, -125.1, -125.2, -125.0, -125.1, -125.2, -125.0, -125.1, -125.2],
            "latitude": [48.0, 48.1, 48.2, 48.5, 48.6, 48.7, 49.0, 49.1, 49.2],
            "biomass_density": [10.5, 12.3, 8.7, 15.2, 11.1, 9.8, 13.4, 16.7, 14.2],
        }
    )

    return gpd.GeoDataFrame(
        data, geometry=gpd.points_from_xy(data["longitude"], data["latitude"]), crs="epsg:4326"
    )


@pytest.fixture
def sample_utm_coordinates():
    """Sample UTM coordinates for testing."""
    return pd.DataFrame(
        {
            "longitude": [-125.0, -125.1, -125.2],
            "latitude": [48.0, 48.1, 48.2],
            "x": [500000, 495000, 490000],  # Approximate UTM coordinates
            "y": [5300000, 5311000, 5322000],
        }
    )


@pytest.fixture
def sample_projection_string():
    """Sample projection string for testing."""
    return "epsg:4326"


@pytest.fixture
def sample_utm_projection_string():
    """Sample UTM projection string for testing."""
    return "epsg:32610"  # UTM Zone 10N


@pytest.fixture
def sample_mesh_region_assignments():
    """Sample mesh region assignments for testing."""
    return pd.DataFrame(
        {
            "transect_num": [1, 2, 3, 121, 122, 129, 130],
            "mesh_region": [1, 1, 1, 2, 2, 3, 3],
            "transect_lower_bound": [1.1, 2.1, 3.1, 121.6, 122.6, 129.1, 130.1],
            "transect_upper_bound": [1.4, 2.4, 3.4, 121.9, 122.9, 129.4, 130.4],
        }
    )


@pytest.fixture
def sample_spatial_grouped():
    """Sample spatial grouped data for testing transect_coordinate_centroid."""
    points = [Point(-125.0, 48.0), Point(-125.1, 48.1), Point(-125.2, 48.2)]
    return gpd.GeoSeries(points, crs="epsg:4326")


@pytest.fixture
def sample_coordinates_for_utm():
    """Sample coordinates for UTM string generation."""
    return {"longitude": -125.0, "latitude": 48.0, "expected_utm": "32610"}


@pytest.fixture
def sample_southern_hemisphere_coordinates():
    """Sample southern hemisphere coordinates for UTM testing."""
    return {"longitude": -125.0, "latitude": -48.0, "expected_utm": "32710"}


@pytest.fixture
def sample_geodataframe():
    """Sample GeoDataFrame for testing."""
    data = pd.DataFrame(
        {
            "longitude": [-125.0, -125.1, -125.2],
            "latitude": [48.0, 48.1, 48.2],
            "name": ["Point A", "Point B", "Point C"],
        }
    )

    return gpd.GeoDataFrame(
        data, geometry=gpd.points_from_xy(data["longitude"], data["latitude"]), crs="epsg:4326"
    )


@pytest.fixture
def sample_geoseries():
    """Sample GeoSeries for testing."""
    points = [Point(-125.0, 48.0), Point(-125.1, 48.1), Point(-125.2, 48.2)]
    return gpd.GeoSeries(points, crs="epsg:4326")
