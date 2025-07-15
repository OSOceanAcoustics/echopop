import pandas as pd
from shapely.geometry import Point

from echopop.nwfsc_feat.FEAT import transect_mesh_region_2019
from echopop.nwfsc_feat.mesh import hull_crop, transect_ends_crop
from echopop.nwfsc_feat.projection import utm_string_generator, wgs84_to_utm
from echopop.nwfsc_feat.spatial import transect_coordinate_centroid, transect_extent


# ==================================================================================================
# Test transect_mesh_region_2019
# -------------------------------
def test_transect_mesh_region_2019_region_1():
    """Test transect_mesh_region_2019 for region 1."""
    start, end, lower, upper = transect_mesh_region_2019(1)

    # Check start and end values
    assert start == 1
    assert end == 119

    # Check boundary lists have correct length
    assert len(lower) == 119
    assert len(upper) == 119

    # Check boundary values follow expected pattern
    assert lower[0] == 1.1  # First transect lower bound
    assert upper[0] == 1.4  # First transect upper bound
    assert lower[-1] == 119.1  # Last transect lower bound
    assert upper[-1] == 119.4  # Last transect upper bound


def test_transect_mesh_region_2019_region_2():
    """Test transect_mesh_region_2019 for region 2."""
    start, end, lower, upper = transect_mesh_region_2019(2)

    # Check start and end values
    assert start == 121
    assert end == 127

    # Check boundary lists have correct length
    assert len(lower) == 7
    assert len(upper) == 7

    # Check boundary values follow expected pattern
    assert lower[0] == 121.6  # First transect lower bound
    assert upper[0] == 121.9  # First transect upper bound


def test_transect_mesh_region_2019_region_3():
    """Test transect_mesh_region_2019 for region 3."""
    start, end, lower, upper = transect_mesh_region_2019(3)

    # Check start and end values
    assert start == 129
    assert end == 145

    # Check boundary lists have correct length
    assert len(lower) == 17
    assert len(upper) == 17

    # Check boundary values follow expected pattern
    assert lower[0] == 129.1  # First transect lower bound
    assert upper[0] == 129.4  # First transect upper bound


def test_transect_mesh_region_2019_invalid_region():
    """Test transect_mesh_region_2019 with invalid region."""
    start, end, lower, upper = transect_mesh_region_2019(99)

    # Should default to region 3 behavior
    assert start == 129
    assert end == 145


# ==================================================================================================
# Test utm_string_generator
# --------------------------
def test_utm_string_generator_northern_hemisphere():
    """Test utm_string_generator for northern hemisphere."""
    # Test coordinates in Pacific Northwest (UTM Zone 10N)
    utm_code = utm_string_generator(-125.0, 48.0)
    assert utm_code == "32610"

    # Test coordinates in different zone
    utm_code = utm_string_generator(-120.0, 45.0)
    assert utm_code == "32611"


def test_utm_string_generator_southern_hemisphere():
    """Test utm_string_generator for southern hemisphere."""
    # Test coordinates in southern hemisphere
    utm_code = utm_string_generator(-125.0, -48.0)
    assert utm_code == "32710"


def test_utm_string_generator_single_digit_zone():
    """Test utm_string_generator for single digit UTM zones."""
    # Test a location that should be in UTM zone 1
    utm_code = utm_string_generator(-177.0, 60.0)
    assert utm_code == "32601"  # Should be zero-padded


def test_utm_string_generator_edge_cases():
    """Test utm_string_generator for edge cases."""
    # Test longitude at zone boundary
    utm_code = utm_string_generator(-126.0, 50.0)
    assert utm_code == "32610"

    # Test at equator
    utm_code = utm_string_generator(-120.0, 0.0)
    assert utm_code == "32611"


# ==================================================================================================
# Test wgs84_to_utm
# ------------------
def test_wgs84_to_utm_basic(sample_geodataframe):
    """Test basic wgs84_to_utm functionality."""
    original_crs = sample_geodataframe.crs

    # Apply transformation
    wgs84_to_utm(sample_geodataframe)

    # Check that CRS has changed
    assert sample_geodataframe.crs != original_crs
    assert "32610" in str(sample_geodataframe.crs)  # Should be UTM zone 10N for Pacific Northwest


def test_wgs84_to_utm_preserves_data(sample_geodataframe):
    """Test that wgs84_to_utm preserves all data."""
    original_shape = sample_geodataframe.shape
    original_columns = sample_geodataframe.columns.tolist()

    # Apply transformation
    wgs84_to_utm(sample_geodataframe)

    # Check that data is preserved
    assert sample_geodataframe.shape == original_shape
    assert sample_geodataframe.columns.tolist() == original_columns


# ==================================================================================================
# Test transect_coordinate_centroid
# ----------------------------------
def test_transect_coordinate_centroid_basic(sample_geoseries):
    """Test basic transect_coordinate_centroid functionality."""
    centroid = transect_coordinate_centroid(sample_geoseries)

    # Check that result is a Point
    assert isinstance(centroid, Point)

    # Check that centroid coordinates are reasonable
    assert isinstance(centroid.x, float)
    assert isinstance(centroid.y, float)


def test_transect_coordinate_centroid_single_point():
    """Test transect_coordinate_centroid with single point."""
    import geopandas as gpd

    single_point = gpd.GeoSeries([Point(1.0, 2.0)])

    centroid = transect_coordinate_centroid(single_point)

    # Centroid of single point should be the point itself
    assert centroid.x == 1.0
    assert centroid.y == 2.0


def test_transect_coordinate_centroid_collinear_points():
    """Test transect_coordinate_centroid with collinear points."""
    import geopandas as gpd

    collinear_points = gpd.GeoSeries([Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0)])

    centroid = transect_coordinate_centroid(collinear_points)

    # Centroid should be at middle point
    assert abs(centroid.x - 1.0) < 1e-10
    assert abs(centroid.y - 0.0) < 1e-10


# ==================================================================================================
# Test transect_extent
# ---------------------
def test_transect_extent_basic(sample_transect_data_cropping):
    """Test basic transect_extent functionality."""
    extent = transect_extent(sample_transect_data_cropping, "epsg:4326", 2)

    # Check that result is a geometry object
    assert hasattr(extent, "bounds")
    assert hasattr(extent, "area")

    # Check that extent has positive area
    assert extent.area > 0


def test_transect_extent_single_transect():
    """Test transect_extent with single transect."""
    single_transect = pd.DataFrame(
        {
            "longitude": [-125.0, -125.1, -125.2],
            "latitude": [48.0, 48.1, 48.2],
            "transect_num": [1, 1, 1],
        }
    )

    extent = transect_extent(single_transect, "epsg:4326", 1)

    # Should still create a valid extent
    assert hasattr(extent, "bounds")
    assert extent.area > 0


def test_transect_extent_num_nearest_transects(sample_transect_data_cropping):
    """Test transect_extent with different num_nearest_transects values."""
    extent1 = transect_extent(sample_transect_data_cropping, "epsg:4326", 1)
    extent2 = transect_extent(sample_transect_data_cropping, "epsg:4326", 2)

    # Both should create valid extents
    assert hasattr(extent1, "bounds")
    assert hasattr(extent2, "bounds")
    assert extent1.area > 0
    assert extent2.area > 0


# ==================================================================================================
# Test hull_crop
# ---------------
def test_hull_crop_basic(sample_transect_data_cropping, sample_mesh_data_cropping):
    """Test basic hull_crop functionality."""
    cropped_mesh = hull_crop(sample_transect_data_cropping, sample_mesh_data_cropping)

    # Check that result is a DataFrame
    assert isinstance(cropped_mesh, pd.DataFrame)

    # Check that mesh has been cropped (should be smaller)
    assert len(cropped_mesh) <= len(sample_mesh_data_cropping)

    # Check that required columns are preserved
    assert "longitude" in cropped_mesh.columns
    assert "latitude" in cropped_mesh.columns


def test_hull_crop_parameters(sample_transect_data_cropping, sample_mesh_data_cropping):
    """Test hull_crop with different parameters."""
    cropped_mesh1 = hull_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, num_nearest_transects=2
    )
    cropped_mesh2 = hull_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, num_nearest_transects=3
    )

    # Both should return valid DataFrames
    assert isinstance(cropped_mesh1, pd.DataFrame)
    assert isinstance(cropped_mesh2, pd.DataFrame)

    # Check that both preserve required columns
    for df in [cropped_mesh1, cropped_mesh2]:
        assert "longitude" in df.columns
        assert "latitude" in df.columns


def test_hull_crop_buffer_distance(sample_transect_data_cropping, sample_mesh_data_cropping):
    """Test hull_crop with different buffer distances."""
    cropped_mesh1 = hull_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, mesh_buffer_distance=1.0
    )
    cropped_mesh2 = hull_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, mesh_buffer_distance=5.0
    )

    # Larger buffer should include more mesh points
    assert len(cropped_mesh2) >= len(cropped_mesh1)


# ==================================================================================================
# Test transect_ends_crop
# ------------------------
def test_transect_ends_crop_basic(sample_transect_data_cropping, sample_mesh_data_cropping):
    """Test basic transect_ends_crop functionality."""
    cropped_mesh, annotated_transects = transect_ends_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, 0.05, transect_mesh_region_2019
    )

    # Check that both outputs are DataFrames
    assert isinstance(cropped_mesh, pd.DataFrame)
    assert isinstance(annotated_transects, pd.DataFrame)

    # Check that mesh has been cropped
    assert len(cropped_mesh) <= len(sample_mesh_data_cropping)

    # Check that transect data has been annotated
    assert "mesh_region" in annotated_transects.columns
    assert "transect_lower_bound" in annotated_transects.columns
    assert "transect_upper_bound" in annotated_transects.columns


def test_transect_ends_crop_latitude_resolution(
    sample_transect_data_cropping, sample_mesh_data_cropping
):
    """Test transect_ends_crop with different latitude resolutions."""
    cropped_mesh1, _ = transect_ends_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, 0.01, transect_mesh_region_2019
    )
    cropped_mesh2, _ = transect_ends_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, 0.1, transect_mesh_region_2019
    )

    # Both should return valid DataFrames
    assert isinstance(cropped_mesh1, pd.DataFrame)
    assert isinstance(cropped_mesh2, pd.DataFrame)

    # Check that results have different sizes (different resolution affects cropping)
    assert len(cropped_mesh1) != len(cropped_mesh2)


def test_transect_ends_crop_region_annotation(
    sample_transect_data_cropping, sample_mesh_data_cropping
):
    """Test that transect_ends_crop correctly annotates regions."""
    _, annotated_transects = transect_ends_crop(
        sample_transect_data_cropping, sample_mesh_data_cropping, 0.05, transect_mesh_region_2019
    )

    # Check that regions are properly assigned
    assert "mesh_region" in annotated_transects.columns
    unique_regions = annotated_transects["mesh_region"].unique()
    assert all(region in [1, 2, 3] for region in unique_regions)

    # Check that boundary values are assigned
    assert "transect_lower_bound" in annotated_transects.columns
    assert "transect_upper_bound" in annotated_transects.columns
    assert not annotated_transects["transect_lower_bound"].isna().all()
    assert not annotated_transects["transect_upper_bound"].isna().all()
