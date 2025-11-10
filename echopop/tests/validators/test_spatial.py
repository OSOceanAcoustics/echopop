import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError

from echopop.validators.spatial import (
    MeshDF,
    TransectsDF,
    ValidateHullCropArgs,
)


# ==================================================================================================
# Test MeshDF
# -----------
def test_mesh_df_valid_lon_lat():
    """Test MeshDF with valid longitude/latitude coordinates."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.5, -123.0],
            "latitude": [46.0, 46.5, 47.0],
            "area": [1.5, 2.0, 1.8],
        }
    )

    result = MeshDF.validate(df)
    pd.testing.assert_frame_equal(result, df)


def test_mesh_df_valid_x_y():
    """Test MeshDF with valid x/y coordinates."""

    df = pd.DataFrame({"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0], "fraction": [0.1, 0.15, 0.12]})

    result = MeshDF.validate(df)
    pd.testing.assert_frame_equal(result, df)


def test_mesh_df_valid_all_coordinates():
    """Test MeshDF with all coordinate types."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.5, -123.0],
            "latitude": [46.0, 46.5, 47.0],
            "x": [0.0, 1.0, 2.0],
            "y": [0.0, 1.0, 2.0],
            "area": [1.5, 2.0, 1.8],
            "fraction": [0.1, 0.15, 0.12],
        }
    )

    result = MeshDF.validate(df)
    pd.testing.assert_frame_equal(result, df)


def test_mesh_df_invalid_longitude():
    """Test MeshDF with invalid longitude values."""

    df = pd.DataFrame(
        {
            "longitude": [-200.0, -123.5, -123.0],  # Invalid: < -180
            "latitude": [46.0, 46.5, 47.0],
            "area": [1.5, 2.0, 1.8],
        }
    )

    with pytest.raises(ValueError):
        MeshDF.validate(df)


def test_mesh_df_invalid_latitude():
    """Test MeshDF with invalid latitude values."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.5, -123.0],
            "latitude": [100.0, 46.5, 47.0],  # Invalid: > 90
            "area": [1.5, 2.0, 1.8],
        }
    )

    with pytest.raises(ValueError):
        MeshDF.validate(df)


def test_mesh_df_missing_coordinate_pairs():
    """Test MeshDF with missing coordinate pairs."""

    # Missing latitude
    df = pd.DataFrame({"longitude": [-124.0, -123.5, -123.0], "area": [1.5, 2.0, 1.8]})

    with pytest.raises(ValueError, match="requires either paired"):
        MeshDF.validate(df)

    # Missing y
    df = pd.DataFrame({"x": [0.0, 1.0, 2.0], "area": [1.5, 2.0, 1.8]})

    with pytest.raises(ValueError, match="requires either paired"):
        MeshDF.validate(df)


def test_mesh_df_missing_area_fraction():
    """Test MeshDF with missing area and fraction columns."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.5, -123.0],
            "latitude": [46.0, 46.5, 47.0],
            # Missing both area and fraction
        }
    )

    with pytest.raises(ValueError, match="requires either an 'area' column"):
        MeshDF.validate(df)


def test_mesh_df_type_coercion():
    """Test MeshDF with type coercion."""

    df = pd.DataFrame(
        {
            "longitude": ["-124.0", "-123.5", "-123.0"],  # String values
            "latitude": ["46.0", "46.5", "47.0"],
            "area": ["1.5", "2.0", "1.8"],
        }
    )

    result = MeshDF.validate(df)
    assert result["longitude"].dtype == float
    assert result["latitude"].dtype == float
    assert result["area"].dtype == float


# ==================================================================================================
# Test TransectsDF
# ----------------
def test_transects_df_valid_lon_lat():
    """Test TransectsDF with valid longitude/latitude coordinates."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6, -123.4],
            "latitude": [46.0, 46.2, 46.4, 46.6],
            "biomass_density": [10.5, 15.2, 12.8, 18.1],
        }
    )

    result = TransectsDF.validate(df)
    pd.testing.assert_frame_equal(result, df)


def test_transects_df_valid_x_y():
    """Test TransectsDF with valid x/y coordinates."""

    df = pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": [0.0, 0.5, 1.0, 1.5],
            "biomass_density": [10.5, 15.2, 12.8, 18.1],
        }
    )

    result = TransectsDF.validate(df)
    pd.testing.assert_frame_equal(result, df)


def test_transects_df_invalid_longitude():
    """Test TransectsDF with invalid longitude values."""

    df = pd.DataFrame(
        {
            "longitude": [-200.0, -123.8, -123.6, -123.4],  # Invalid: < -180
            "latitude": [46.0, 46.2, 46.4, 46.6],
            "biomass_density": [10.5, 15.2, 12.8, 18.1],
        }
    )

    with pytest.raises(ValueError):
        TransectsDF.validate(df)


def test_transects_df_invalid_latitude():
    """Test TransectsDF with invalid latitude values."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6, -123.4],
            "latitude": [100.0, 46.2, 46.4, 46.6],  # Invalid: > 90
            "biomass_density": [10.5, 15.2, 12.8, 18.1],
        }
    )

    with pytest.raises(ValueError):
        TransectsDF.validate(df)


def test_transects_df_missing_coordinate_pairs():
    """Test TransectsDF with missing coordinate pairs."""

    # Missing latitude
    df = pd.DataFrame(
        {"longitude": [-124.0, -123.8, -123.6, -123.4], "biomass_density": [10.5, 15.2, 12.8, 18.1]}
    )

    with pytest.raises(ValueError, match="requires either paired"):
        TransectsDF.validate(df)

    # Missing y
    df = pd.DataFrame({"x": [0.0, 1.0, 2.0, 3.0], "biomass_density": [10.5, 15.2, 12.8, 18.1]})

    with pytest.raises(ValueError, match="requires either paired"):
        TransectsDF.validate(df)


def test_transects_df_type_coercion():
    """Test TransectsDF with type coercion."""

    df = pd.DataFrame(
        {
            "longitude": ["-124.0", "-123.8", "-123.6", "-123.4"],  # String values
            "latitude": ["46.0", "46.2", "46.4", "46.6"],
            "biomass_density": ["10.5", "15.2", "12.8", "18.1"],
        }
    )

    result = TransectsDF.validate(df)
    assert result["longitude"].dtype == float
    assert result["latitude"].dtype == float


# ==================================================================================================
# Test ValidateHullCropArgs
# -------------------------
def test_validate_hull_crop_args_valid():
    """Test ValidateHullCropArgs with valid arguments."""

    transects_df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6, -123.4],
            "latitude": [46.0, 46.2, 46.4, 46.6],
            "biomass_density": [10.5, 15.2, 12.8, 18.1],
        }
    )

    mesh_df = pd.DataFrame(
        {
            "longitude": [-124.5, -124.0, -123.5, -123.0],
            "latitude": [45.8, 46.0, 46.2, 46.4],
            "area": [1.5, 2.0, 1.8, 2.2],
        }
    )

    args = ValidateHullCropArgs(
        transects=transects_df,
        mesh=mesh_df,
        num_nearest_transects=3,
        mesh_buffer_distance=2.5,
        projection="epsg:4326",
        coordinate_names=("longitude", "latitude"),
    )

    assert args.num_nearest_transects == 3
    assert args.mesh_buffer_distance == 2.5
    assert args.projection == "epsg:4326"
    assert args.coordinate_names == ("longitude", "latitude")
    pd.testing.assert_frame_equal(args.transects, transects_df)
    pd.testing.assert_frame_equal(args.mesh, mesh_df)


def test_validate_hull_crop_args_invalid_num_nearest():
    """Test ValidateHullCropArgs with invalid num_nearest_transects."""

    transects_df = pd.DataFrame(
        {"longitude": [-124.0, -123.8], "latitude": [46.0, 46.2], "biomass_density": [10.5, 15.2]}
    )

    mesh_df = pd.DataFrame(
        {"longitude": [-124.5, -124.0], "latitude": [45.8, 46.0], "area": [1.5, 2.0]}
    )

    with pytest.raises(ValidationError):
        ValidateHullCropArgs(
            transects=transects_df,
            mesh=mesh_df,
            num_nearest_transects=0,  # Invalid: must be > 0
            mesh_buffer_distance=2.5,
            projection="epsg:4326",
            coordinate_names=("longitude", "latitude"),
        )


def test_validate_hull_crop_args_invalid_buffer_distance():
    """Test ValidateHullCropArgs with invalid mesh_buffer_distance."""

    transects_df = pd.DataFrame(
        {"longitude": [-124.0, -123.8], "latitude": [46.0, 46.2], "biomass_density": [10.5, 15.2]}
    )

    mesh_df = pd.DataFrame(
        {"longitude": [-124.5, -124.0], "latitude": [45.8, 46.0], "area": [1.5, 2.0]}
    )

    with pytest.raises(ValidationError):
        ValidateHullCropArgs(
            transects=transects_df,
            mesh=mesh_df,
            num_nearest_transects=3,
            mesh_buffer_distance=-1.0,  # Invalid: must be >= 0
            projection="epsg:4326",
            coordinate_names=("longitude", "latitude"),
        )


def test_validate_hull_crop_args_projection_validation():
    """Test ValidateHullCropArgs projection field validation."""

    transects_df = pd.DataFrame(
        {"longitude": [-124.0, -123.8], "latitude": [46.0, 46.2], "biomass_density": [10.5, 15.2]}
    )

    mesh_df = pd.DataFrame(
        {"longitude": [-124.5, -124.0], "latitude": [45.8, 46.0], "area": [1.5, 2.0]}
    )

    # Test different valid projection formats
    valid_projections = ["4326", "epsg:4326", "EPSG:4326", 4326]

    for proj in valid_projections:
        args = ValidateHullCropArgs(
            transects=transects_df,
            mesh=mesh_df,
            num_nearest_transects=3,
            mesh_buffer_distance=2.5,
            projection=proj,
            coordinate_names=("longitude", "latitude"),
        )
        assert args.projection == "epsg:4326"


def test_validate_hull_crop_args_invalid_projection():
    """Test ValidateHullCropArgs with invalid projection format."""

    transects_df = pd.DataFrame(
        {"longitude": [-124.0, -123.8], "latitude": [46.0, 46.2], "biomass_density": [10.5, 15.2]}
    )

    mesh_df = pd.DataFrame(
        {"longitude": [-124.5, -124.0], "latitude": [45.8, 46.0], "area": [1.5, 2.0]}
    )

    with pytest.raises(ValidationError):
        ValidateHullCropArgs(
            transects=transects_df,
            mesh=mesh_df,
            num_nearest_transects=3,
            mesh_buffer_distance=2.5,
            projection="invalid_projection",
            coordinate_names=("longitude", "latitude"),
        )


def test_validate_hull_crop_args_invalid_transects():
    """Test ValidateHullCropArgs with invalid transects DataFrame."""

    # Invalid transects (missing coordinates)
    transects_df = pd.DataFrame({"biomass_density": [10.5, 15.2]})

    mesh_df = pd.DataFrame(
        {"longitude": [-124.5, -124.0], "latitude": [45.8, 46.0], "area": [1.5, 2.0]}
    )

    with pytest.raises(ValueError):
        ValidateHullCropArgs(
            transects=transects_df,
            mesh=mesh_df,
            num_nearest_transects=3,
            mesh_buffer_distance=2.5,
            projection="epsg:4326",
            coordinate_names=("longitude", "latitude"),
        )


def test_validate_hull_crop_args_invalid_mesh():
    """Test ValidateHullCropArgs with invalid mesh DataFrame."""

    transects_df = pd.DataFrame(
        {"longitude": [-124.0, -123.8], "latitude": [46.0, 46.2], "biomass_density": [10.5, 15.2]}
    )

    # Invalid mesh (missing coordinates and area/fraction)
    mesh_df = pd.DataFrame({"other_column": [1, 2]})

    with pytest.raises(ValueError):
        ValidateHullCropArgs(
            transects=transects_df,
            mesh=mesh_df,
            num_nearest_transects=3,
            mesh_buffer_distance=2.5,
            projection="epsg:4326",
            coordinate_names=("longitude", "latitude"),
        )


# ==================================================================================================
# Integration tests
# -----------------
def test_spatial_validators_integration():
    """Test integration of spatial validators."""

    # Create valid test data
    transects_df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6, -123.4, -123.2],
            "latitude": [46.0, 46.2, 46.4, 46.6, 46.8],
            "x": [0.0, 0.2, 0.4, 0.6, 0.8],
            "y": [0.0, 0.2, 0.4, 0.6, 0.8],
            "biomass_density": [10.5, 15.2, 12.8, 18.1, 14.3],
            "transect_num": [1, 1, 2, 2, 3],
        }
    )

    mesh_df = pd.DataFrame(
        {
            "longitude": np.linspace(-124.5, -123.0, 10),
            "latitude": np.linspace(45.8, 47.0, 10),
            "x": np.linspace(-0.5, 1.0, 10),
            "y": np.linspace(-0.2, 1.2, 10),
            "area": np.ones(10) * 1.5,
            "fraction": np.ones(10) * 0.1,
        }
    )

    # Test individual DataFrame validations
    validated_transects = TransectsDF.validate(transects_df)
    validated_mesh = MeshDF.validate(mesh_df)

    # Test hull crop args validation
    hull_args = ValidateHullCropArgs(
        transects=validated_transects,
        mesh=validated_mesh,
        num_nearest_transects=3,
        mesh_buffer_distance=2.5,
        projection="epsg:4326",
        coordinate_names=("longitude", "latitude"),
    )

    assert hull_args.num_nearest_transects == 3
    assert len(hull_args.transects) == 5
    assert len(hull_args.mesh) == 10
