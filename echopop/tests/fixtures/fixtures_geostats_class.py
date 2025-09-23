import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_mesh_data_with_lonlat():
    """Create sample mesh data with longitude/latitude columns for geostatistics tests."""
    np.random.seed(999)

    # Create a simple grid of mesh points
    lon_range = np.linspace(-125.0, -123.0, 10)
    lat_range = np.linspace(46.0, 48.0, 10)

    mesh_points = []
    for i, lon in enumerate(lon_range):
        for j, lat in enumerate(lat_range):
            mesh_points.append(
                {"longitude": lon, "latitude": lat, "transect_num": i % 5}  # Some transect numbers
            )

    return pd.DataFrame(mesh_points)


@pytest.fixture
def sample_coordinates_df_with_biomass():
    """Create sample coordinate data with biomass_density column for geostatistics tests."""
    # Use same pattern as working variogram tests
    return pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6, -123.4, -123.2, -123.0, -122.8, -122.6],
            "latitude": [46.0, 46.2, 46.4, 46.6, 46.8, 47.0, 47.2, 47.4],
            "biomass_density": [10.5, 15.2, 12.8, 18.1, 14.3, 16.7, 11.9, 13.4],
            "transect_num": [1, 1, 2, 2, 3, 3, 4, 4],
        }
    )


@pytest.fixture
def sample_mesh_data_latlon():
    """Create sample mesh data with x/y columns for testing custom coordinate names."""
    np.random.seed(999)

    # Create sample transect data
    n_points = 1000
    x_coords = np.random.uniform(-122.5, -123.5, n_points)
    y_coords = np.random.uniform(32.5, 48.5, n_points)

    # Collect into DataFrame
    df = pd.DataFrame(
        {
            "longitude": x_coords,
            "latitude": y_coords,
        }
    )

    return df


@pytest.fixture
def sample_mesh_data_xy():
    """Create sample mesh data with x/y columns for testing custom coordinate names."""
    np.random.seed(999)

    # Create a simple grid of mesh points
    x_range = np.linspace(-125.0, -123.0, 10)
    y_range = np.linspace(46.0, 48.0, 10)

    mesh_points = []
    for i, x in enumerate(x_range):
        for j, y in enumerate(y_range):
            mesh_points.append({"x": x, "y": y, "transect_num": i % 5})

    return pd.DataFrame(mesh_points)


@pytest.fixture
def sample_coordinates_df_latlon():
    """Create sample coordinate data with x/y columns for testing custom coordinate names."""
    np.random.seed(999)

    # Create sample transect data
    n_points = 500
    x_coords = np.random.uniform(-124.5, -123.5, n_points)
    y_coords = np.random.uniform(35.5, 47.5, n_points)
    biomass_values = np.random.normal(10.0, 3.0, n_points)

    # Collect into DataFrame
    df = pd.DataFrame(
        {
            "longitude": x_coords,
            "latitude": y_coords,
            "biomass_density": biomass_values,
            "transect_num": np.random.randint(1, 6, n_points),
        }
    )

    return df.sort_values(["transect_num", "longitude"])


@pytest.fixture
def sample_coordinates_df_xy():
    """Create sample coordinate data with x/y columns for testing custom coordinate names."""
    np.random.seed(999)

    # Create sample transect data
    n_points = 50
    x_coords = np.random.uniform(-124.5, -123.5, n_points)
    y_coords = np.random.uniform(46.5, 47.5, n_points)
    biomass_values = np.random.normal(10.0, 3.0, n_points)

    return pd.DataFrame(
        {
            "x": x_coords,
            "y": y_coords,
            "biomass_density": biomass_values,
            "transect_num": np.random.randint(1, 6, n_points),
        }
    )
