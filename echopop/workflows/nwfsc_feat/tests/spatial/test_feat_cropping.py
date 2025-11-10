import numpy as np
import pandas as pd

from echopop.workflows.nwfsc_feat import (
    get_survey_western_extents,
    transect_ends_crop,
    western_boundary_search_strategy,
)
from echopop.workflows.nwfsc_feat.parameters import transect_mesh_region_2019


# ==================================================================================================
# Test western_boundary_search_strategy
# --------------------------------------
def test_western_boundary_search_strategy_basic(sample_western_extent):
    """Test basic functionality of western_boundary_search_strategy."""
    # Create mock data
    kriging_mesh = pd.DataFrame(
        {"x": [-1.0, 0.0, 1.0], "y": [0.0, 1.0, 2.0], "transect_num": [1, 2, 3]}
    )

    sparse_radii = np.array([0, 1, 2])
    valid_distances = np.array([2, 2, 2])
    local_points = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    distance_matrix_masked = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    nearby_indices = np.array([[0, 1], [1, 0], [2, 1]])
    k_min = 2
    k_max = 3
    search_radius = 1.5
    wr_indices = np.array([[0, 1, 2], [1, 0, 2], [2, 1, 0]], dtype=float)
    oos_indices = np.full((3, 2), np.nan)
    oos_weights = np.ones(3)

    result = western_boundary_search_strategy(
        kriging_mesh,
        sample_western_extent,
        coordinate_names=("x", "y"),
        sparse_radii=sparse_radii,
        valid_distances=valid_distances,
        local_points=local_points,
        distance_matrix_masked=distance_matrix_masked,
        nearby_indices=nearby_indices,
        k_min=k_min,
        k_max=k_max,
        search_radius=search_radius,
        wr_indices=wr_indices,
        oos_indices=oos_indices,
        oos_weights=oos_weights,
    )

    # Check that result is a tuple with 3 elements
    assert isinstance(result, tuple)
    assert len(result) == 3

    # Check that each element is a numpy array
    updated_wr_indices, updated_oos_indices, updated_oos_weights = result
    assert isinstance(updated_wr_indices, np.ndarray)
    assert isinstance(updated_oos_indices, np.ndarray)
    assert isinstance(updated_oos_weights, np.ndarray)

    # Check that weights are positive
    assert all(w > 0 for w in updated_oos_weights)  # All weights should be positive


# ==================================================================================================
# Test get_survey_western_extents
# --------------------------------
def test_get_survey_western_extents_basic(sample_mesh_coordinates):
    """Test basic functionality of get_survey_western_extents."""
    # Add latitude column to mesh coordinates
    mesh_with_lat = sample_mesh_coordinates.copy()
    mesh_with_lat["latitude"] = mesh_with_lat["y"] + 46.0  # Add realistic latitude

    result = get_survey_western_extents(
        mesh_with_lat, coordinate_names=("x", "y"), latitude_threshold=45.0
    )

    # Check that result is a DataFrame
    assert isinstance(result, pd.DataFrame)

    # Check that it has the expected columns
    expected_columns = {"transect_num", "x", "y"}
    assert set(result.columns) == expected_columns

    # Check that each transect has one entry
    transect_counts = result["transect_num"].value_counts()
    assert all(count == 1 for count in transect_counts.values)


def test_get_survey_western_extents_westernmost(sample_mesh_coordinates):
    """Test that get_survey_western_extents returns westernmost points."""
    # Add latitude column to mesh coordinates
    mesh_with_lat = sample_mesh_coordinates.copy()
    mesh_with_lat["latitude"] = mesh_with_lat["y"] + 46.0  # Add realistic latitude

    result = get_survey_western_extents(
        mesh_with_lat, coordinate_names=("x", "y"), latitude_threshold=45.0
    )

    # For each transect, the x-coordinate should be the minimum
    for transect_num in result["transect_num"].unique():
        original_transect = mesh_with_lat[mesh_with_lat["transect_num"] == transect_num]
        result_transect = result[result["transect_num"] == transect_num]

        min_x = original_transect["x"].min()
        result_x = result_transect["x"].iloc[0]

        assert result_x == min_x


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
