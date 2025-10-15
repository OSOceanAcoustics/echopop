import numpy as np
import pytest

from echopop.geostatistics.variogram import empirical_variogram, quantize_lags, semivariance


# ==================================================================================================
# Test quantize_lags
# ------------------
def test_quantize_lags_basic(
    sample_estimates,
    sample_variogram_lag_matrix,
    sample_variogram_mask_matrix,
    sample_variogram_azimuth_matrix,
):
    """Test basic quantize_lags functionality."""
    n_lags = 6

    lag_counts, lag_estimates, lag_estimates_squared, lag_deviations = quantize_lags(
        sample_estimates,
        sample_variogram_lag_matrix,
        sample_variogram_mask_matrix,
        sample_variogram_azimuth_matrix,
        n_lags,
    )

    # Check that all outputs are arrays
    assert isinstance(lag_counts, np.ndarray), "lag_counts should be numpy array"
    assert isinstance(lag_estimates, np.ndarray), "lag_estimates should be numpy array"
    assert isinstance(
        lag_estimates_squared, np.ndarray
    ), "lag_estimates_squared should be numpy array"
    assert isinstance(lag_deviations, np.ndarray), "lag_deviations should be numpy array"

    # Check that lag_counts are integers
    assert lag_counts.dtype.kind in "iu", "lag_counts should be integers"

    # Check that other outputs are floats
    assert lag_estimates.dtype.kind == "f", "lag_estimates should be floats"
    assert lag_estimates_squared.dtype.kind == "f", "lag_estimates_squared should be floats"
    assert lag_deviations.dtype.kind == "f", "lag_deviations should be floats"

    # Check that all arrays have the same length
    assert len(lag_counts) == len(lag_estimates), "All outputs should have same length"
    assert len(lag_counts) == len(lag_estimates_squared), "All outputs should have same length"
    assert len(lag_counts) == len(lag_deviations), "All outputs should have same length"


def test_quantize_lags_with_azimuth_threshold(
    sample_estimates,
    sample_variogram_lag_matrix,
    sample_variogram_mask_matrix,
    sample_variogram_azimuth_matrix,
):
    """Test quantize_lags with azimuth angle threshold."""
    n_lags = 6
    azimuth_threshold = 45.0

    lag_counts, lag_estimates, lag_estimates_squared, lag_deviations = quantize_lags(
        sample_estimates,
        sample_variogram_lag_matrix,
        sample_variogram_mask_matrix,
        sample_variogram_azimuth_matrix,
        n_lags,
        azimuth_angle_threshold=azimuth_threshold,
    )

    # Check that all outputs are arrays (length may vary based on actual data)
    assert isinstance(lag_counts, np.ndarray), "lag_counts should be numpy array"
    assert isinstance(lag_estimates, np.ndarray), "lag_estimates should be numpy array"
    assert isinstance(
        lag_estimates_squared, np.ndarray
    ), "lag_estimates_squared should be numpy array"
    assert isinstance(lag_deviations, np.ndarray), "lag_deviations should be numpy array"

    # Check that all arrays have the same length
    assert len(lag_counts) == len(lag_estimates), "All outputs should have same length"
    assert len(lag_counts) == len(lag_estimates_squared), "All outputs should have same length"
    assert len(lag_counts) == len(lag_deviations), "All outputs should have same length"


def test_quantize_lags_invalid_estimates():
    """Test quantize_lags with invalid estimates array."""
    estimates_2d = np.array([[1, 2], [3, 4]])
    lag_matrix = np.array([[0, 1], [1, 0]])
    mask_matrix = np.array([[True, False], [False, True]])
    azimuth_matrix = np.array([[0, 45], [225, 0]])

    with pytest.raises(ValueError, match="Estimates array.*must be a 1D array"):
        quantize_lags(estimates_2d, lag_matrix, mask_matrix, azimuth_matrix, 3)


def test_quantize_lags_invalid_dimensions():
    """Test quantize_lags with invalid matrix dimensions."""
    estimates = np.array([1, 2, 3])
    lag_matrix_1d = np.array([0, 1, 2])  # Invalid 1D array
    mask_matrix = np.array([[True, False], [False, True]])
    azimuth_matrix = np.array([[0, 45], [225, 0]])

    with pytest.raises(ValueError, match="requires arrays to be 2D"):
        quantize_lags(estimates, lag_matrix_1d, mask_matrix, azimuth_matrix, 3)


# ==================================================================================================
# Test semivariance
# ----------------
def test_semivariance_basic(
    sample_estimates,
    sample_lag_estimates,
    sample_lag_estimates_squared,
    sample_lag_counts,
    sample_lag_deviations,
    sample_head_index,
):
    """Test basic semivariance calculation."""
    gamma_h, lag_covariance = semivariance(
        sample_estimates,
        sample_lag_estimates,
        sample_lag_estimates_squared,
        sample_lag_counts,
        sample_lag_deviations,
        sample_head_index,
    )

    # Check that outputs are arrays of correct length
    assert len(gamma_h) == len(sample_lag_counts), "gamma_h should match lag_counts length"
    assert isinstance(lag_covariance, (float, np.floating)), "lag_covariance should be a scalar"

    # Check that semivariance values are non-negative
    assert all(gamma >= 0 for gamma in gamma_h), "Semivariance values should be non-negative"


def test_semivariance_output_types(
    sample_estimates,
    sample_lag_estimates,
    sample_lag_estimates_squared,
    sample_lag_counts,
    sample_lag_deviations,
    sample_head_index,
):
    """Test that semivariance returns correct data types."""
    gamma_h, lag_covariance = semivariance(
        sample_estimates,
        sample_lag_estimates,
        sample_lag_estimates_squared,
        sample_lag_counts,
        sample_lag_deviations,
        sample_head_index,
    )

    # Check data types
    assert isinstance(gamma_h, np.ndarray), "gamma_h should be numpy array"
    assert gamma_h.dtype == np.float64, "gamma_h should be float64"
    assert isinstance(lag_covariance, (float, np.floating)), "lag_covariance should be float"


# ==================================================================================================
# Test empirical_variogram
# ------------------------
def test_empirical_variogram_basic(sample_transect_df):
    """Test basic empirical variogram calculation."""
    n_lags = 10
    lag_resolution = 0.02

    lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
        sample_transect_df,
        n_lags=n_lags,
        lag_resolution=lag_resolution,
        azimuth_filter=False,
        azimuth_angle_threshold=180.0,
    )

    # Check that outputs are arrays of correct length
    assert len(lags) == n_lags, f"lags should have {n_lags} elements"
    assert len(gamma_h) == n_lags, f"gamma_h should have {n_lags} elements"
    assert len(lag_counts) == n_lags, f"lag_counts should have {n_lags} elements"
    assert isinstance(lag_covariance, (float, np.floating)), "lag_covariance should be a scalar"

    # Check that first lag is 0 (due to force_lag_zero=True)
    assert lags[0] == 0.0, "First lag should be 0.0"
    assert gamma_h[0] == 0.0, "First gamma_h should be 0.0"


def test_empirical_variogram_with_azimuth(sample_transect_df):
    """Test empirical variogram with azimuth filtering."""
    n_lags = 10
    lag_resolution = 0.02

    lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
        sample_transect_df,
        n_lags=n_lags,
        lag_resolution=lag_resolution,
        azimuth_filter=True,
        azimuth_angle_threshold=45.0,
    )

    # Check that outputs are arrays of correct length
    assert len(lags) == n_lags, f"lags should have {n_lags} elements"
    assert len(gamma_h) == n_lags, f"gamma_h should have {n_lags} elements"
    assert len(lag_counts) == n_lags, f"lag_counts should have {n_lags} elements"
    assert isinstance(lag_covariance, (float, np.floating)), "lag_covariance should be a scalar"


def test_empirical_variogram_different_variable(sample_transect_df):
    """Test empirical variogram with different variable name."""
    # Add temperature column to the fixture data
    test_df = sample_transect_df.copy()
    test_df["temperature"] = [15.0, 17.2, 14.8, 18.1, 16.3, 19.7, 13.9, 15.4]

    n_lags = 10
    lag_resolution = 0.02

    lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
        test_df,
        n_lags=n_lags,
        lag_resolution=lag_resolution,
        azimuth_filter=False,
        azimuth_angle_threshold=180.0,
        variable="temperature",
    )

    # Check that outputs are arrays of correct length
    assert len(lags) == n_lags, f"lags should have {n_lags} elements"
    assert len(gamma_h) == n_lags, f"gamma_h should have {n_lags} elements"
    assert len(lag_counts) == n_lags, f"lag_counts should have {n_lags} elements"


def test_empirical_variogram_different_coordinates(sample_transect_df):
    """Test empirical variogram with different coordinate names."""
    # Rename coordinates in the fixture data
    test_df = sample_transect_df.copy()
    test_df = test_df.rename(columns={"x": "easting", "y": "northing"})

    n_lags = 10
    lag_resolution = 0.02

    lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
        test_df,
        n_lags=n_lags,
        lag_resolution=lag_resolution,
        azimuth_filter=False,
        azimuth_angle_threshold=180.0,
        coordinate_names=("easting", "northing"),  # Fixed parameter name
    )

    # Check that outputs are arrays of correct length
    assert len(lags) == n_lags, f"lags should have {n_lags} elements"
    assert len(gamma_h) == n_lags, f"gamma_h should have {n_lags} elements"
    assert len(lag_counts) == n_lags, f"lag_counts should have {n_lags} elements"


def test_empirical_variogram_no_force_lag_zero(sample_transect_df):
    """Test empirical variogram without forcing lag zero."""
    n_lags = 10
    lag_resolution = 0.02

    lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
        sample_transect_df,
        n_lags=n_lags,
        lag_resolution=lag_resolution,
        azimuth_filter=False,
        azimuth_angle_threshold=180.0,
        force_lag_zero=False,
    )

    # Check that outputs are arrays of correct length (should be n_lags - 1)
    assert len(lags) == n_lags - 1, f"lags should have {n_lags - 1} elements"
    assert len(gamma_h) == n_lags - 1, f"gamma_h should have {n_lags - 1} elements"
    assert len(lag_counts) == n_lags - 1, f"lag_counts should have {n_lags - 1} elements"

    # Check that first lag is not 0 (since force_lag_zero=False)
    assert lags[0] != 0.0, "First lag should not be 0.0 when force_lag_zero=False"
