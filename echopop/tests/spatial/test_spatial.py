import numpy as np

from echopop.nwfsc_feat.spatial import (
    filter_lag_matrix,
    lag_distance_matrix,
    standardize_coordinates,
)


# ==================================================================================================
# Test standardize_coordinates
# ----------------------------
def test_standardize_coordinates_basic(sample_coordinates_df):
    """Test basic coordinate standardization."""
    result_df, delta_lon, delta_lat = standardize_coordinates(
        sample_coordinates_df, x_offset=0.0, y_offset=0.0
    )

    # Check that x and y columns are created
    assert "x" in result_df.columns, "x column should be created"
    assert "y" in result_df.columns, "y column should be created"

    # Check that delta values are returned
    assert delta_lon is not None, "delta_longitude should be returned"
    assert delta_lat is not None, "delta_latitude should be returned"

    # Check that values are numeric
    assert all(
        isinstance(x, (int, float, np.number)) for x in result_df["x"]
    ), "x values should be numeric"
    assert all(
        isinstance(y, (int, float, np.number)) for y in result_df["y"]
    ), "y values should be numeric"


def test_standardize_coordinates_with_offsets(sample_coordinates_df):
    """Test coordinate standardization with offsets."""
    result_df, delta_lon, delta_lat = standardize_coordinates(
        sample_coordinates_df, x_offset=0.1, y_offset=0.1
    )

    # Check that offsets are applied
    assert "x" in result_df.columns, "x column should be created"
    assert "y" in result_df.columns, "y column should be created"

    # Check that values are numeric
    assert all(
        isinstance(x, (int, float, np.number)) for x in result_df["x"]
    ), "x values should be numeric"
    assert all(
        isinstance(y, (int, float, np.number)) for y in result_df["y"]
    ), "y values should be numeric"


def test_standardize_coordinates_with_reference(sample_coordinates_df, sample_reference_df):
    """Test coordinate standardization with reference DataFrame."""
    result_df, delta_lon, delta_lat = standardize_coordinates(
        sample_coordinates_df, reference_df=sample_reference_df
    )

    # Check that x and y columns are created
    assert "x" in result_df.columns, "x column should be created"
    assert "y" in result_df.columns, "y column should be created"

    # Check that delta values are returned
    assert delta_lon is not None, "delta_longitude should be returned"
    assert delta_lat is not None, "delta_latitude should be returned"


def test_standardize_coordinates_with_deltas(sample_coordinates_df):
    """Test coordinate standardization with provided deltas."""
    result_df, delta_lon, delta_lat = standardize_coordinates(
        sample_coordinates_df, delta_x=1.0, delta_y=1.0
    )

    # Check that provided deltas are returned
    assert delta_lon == 1.0, "delta_longitude should match input"
    assert delta_lat == 1.0, "delta_latitude should match input"


# ==================================================================================================
# Test lag_distance_matrix
# ------------------------
def test_lag_distance_matrix_basic(sample_coordinates_array):
    """Test basic lag distance matrix calculation."""
    distance_matrix, azimuth_matrix = lag_distance_matrix(
        sample_coordinates_array, self=True, azimuth_matrix=False
    )

    # Check that distance matrix is returned
    assert isinstance(distance_matrix, np.ndarray), "Distance matrix should be numpy array"

    # Check that azimuth matrix is empty when not requested
    assert azimuth_matrix.size == 0, "Azimuth matrix should be empty when not requested"

    # Check that distance matrix has some reasonable properties
    assert distance_matrix.ndim >= 2, "Distance matrix should be at least 2D"


def test_lag_distance_matrix_dataframe(sample_coordinates_df):
    """Test lag distance matrix with DataFrame input."""
    distance_matrix, azimuth_matrix = lag_distance_matrix(
        sample_coordinates_df,
        coordinate_names=("longitude", "latitude"),
        self=True,
        azimuth_matrix=False,
    )

    # Check that distance matrix is returned
    assert isinstance(distance_matrix, np.ndarray), "Distance matrix should be numpy array"
    assert distance_matrix.ndim >= 2, "Distance matrix should be at least 2D"


def test_lag_distance_matrix_two_sets(sample_coordinates_array):
    """Test lag distance matrix between two different coordinate sets."""
    coords_2 = np.array([[2.0, 2.0], [3.0, 3.0]])

    distance_matrix, azimuth_matrix = lag_distance_matrix(
        sample_coordinates_array, coordinates_2=coords_2, self=False, azimuth_matrix=False
    )

    # Check that distance matrix is returned
    assert isinstance(distance_matrix, np.ndarray), "Distance matrix should be numpy array"
    assert distance_matrix.ndim >= 2, "Distance matrix should be at least 2D"


# ==================================================================================================
# Test filter_lag_matrix
# ----------------------
def test_filter_lag_matrix_basic(sample_lag_matrix, sample_mask_matrix):
    """Test basic lag matrix filtering."""
    filtered_matrix = filter_lag_matrix(sample_lag_matrix, sample_mask_matrix)

    # Check that result is 1D array
    assert filtered_matrix.ndim == 1, "Filtered matrix should be 1D"

    # Check that only True mask values are retained
    expected_count = np.sum(sample_mask_matrix)
    assert len(filtered_matrix) == expected_count, f"Should have {expected_count} elements"


def test_filter_lag_matrix_with_azimuth(
    sample_lag_matrix, sample_mask_matrix, sample_azimuth_matrix
):
    """Test lag matrix filtering with azimuth threshold."""
    filtered_matrix = filter_lag_matrix(
        sample_lag_matrix,
        sample_mask_matrix,
        azimuth_matrix=sample_azimuth_matrix,
        azimuth_angle_threshold=90.0,
    )

    # Check that result is 1D array
    assert filtered_matrix.ndim == 1, "Filtered matrix should be 1D"

    # Should have fewer elements due to azimuth filtering
    assert len(filtered_matrix) <= np.sum(
        sample_mask_matrix
    ), "Azimuth filtering should reduce elements"


def test_filter_lag_matrix_1d_input(sample_mask_matrix):
    """Test lag matrix filtering with 1D input."""
    data_1d = np.array([1, 2, 3, 4])

    filtered_matrix = filter_lag_matrix(data_1d, sample_mask_matrix)

    # Check that result is 1D array
    assert filtered_matrix.ndim == 1, "Filtered matrix should be 1D"


def test_filter_lag_matrix_no_azimuth(sample_lag_matrix, sample_mask_matrix):
    """Test lag matrix filtering without azimuth matrix."""
    filtered_matrix = filter_lag_matrix(
        sample_lag_matrix,
        sample_mask_matrix,
        azimuth_matrix=np.array([]),  # Empty azimuth matrix
        azimuth_angle_threshold=None,
    )

    # Check that result is 1D array
    assert filtered_matrix.ndim == 1, "Filtered matrix should be 1D"

    # Should have same count as basic filtering
    expected_count = np.sum(sample_mask_matrix)
    assert len(filtered_matrix) == expected_count, f"Should have {expected_count} elements"
