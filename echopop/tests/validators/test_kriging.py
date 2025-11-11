import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError

from echopop.validators.kriging import (
    KrigingParameters,
    ValidateKrigingClass,
    ValidateMeshCropArgs,
    VariogramKrigeModelParameters,
)


# ==================================================================================================
# Test KrigingParameters
# ----------------------
def test_kriging_parameters_valid():
    """Test KrigingParameters with valid parameters."""

    params = KrigingParameters(aspect_ratio=0.5, k_min=5, k_max=20, search_radius=10.0)

    assert params.aspect_ratio == 0.5
    assert params.k_min == 5
    assert params.k_max == 20
    assert params.search_radius == 10.0


def test_kriging_parameters_defaults():
    """Test KrigingParameters with default values."""

    params = KrigingParameters(search_radius=10.0)

    assert params.aspect_ratio == 1e-3  # Default value
    assert params.k_min == 5  # Default value
    assert params.k_max == 20  # Default value
    assert params.search_radius == 10.0


def test_kriging_parameters_invalid_aspect_ratio():
    """Test KrigingParameters with invalid aspect_ratio."""

    # aspect_ratio must be > 0
    with pytest.raises(ValidationError):
        KrigingParameters(aspect_ratio=0.0, search_radius=10.0)

    # aspect_ratio must be <= 1
    with pytest.raises(ValidationError):
        KrigingParameters(aspect_ratio=1.5, search_radius=10.0)


def test_kriging_parameters_invalid_k_min():
    """Test KrigingParameters with invalid k_min."""

    # k_min must be >= 2
    with pytest.raises(ValidationError):
        KrigingParameters(k_min=1, search_radius=10.0)


def test_kriging_parameters_invalid_search_radius():
    """Test KrigingParameters with invalid search_radius."""

    # search_radius must be > 0
    with pytest.raises(ValidationError):
        KrigingParameters(search_radius=0.0)

    # search_radius cannot be infinite
    with pytest.raises(ValidationError):
        KrigingParameters(search_radius=np.inf)


def test_kriging_parameters_k_interval_validation():
    """Test k_min <= k_max validation."""

    # Valid: k_min < k_max
    params = KrigingParameters(k_min=5, k_max=20, search_radius=10.0)
    assert params.k_min == 5
    assert params.k_max == 20

    # Valid: k_min = k_max
    params = KrigingParameters(k_min=10, k_max=10, search_radius=10.0)
    assert params.k_min == 10
    assert params.k_max == 10

    # Invalid: k_min > k_max
    with pytest.raises(ValidationError):
        KrigingParameters(k_min=20, k_max=5, search_radius=10.0)


# ==================================================================================================
# Test VariogramKrigeModelParameters
# ----------------------------------
def test_variogram_krige_model_parameters_valid():
    """Test VariogramKrigeModelParameters with valid parameters."""

    params = VariogramKrigeModelParameters(
        model="exponential", sill=1.0, nugget=0.1, correlation_range=0.5
    )

    assert params.model == "exponential"
    assert params.sill == 1.0
    assert params.nugget == 0.1
    assert params.correlation_range == 0.5


def test_variogram_krige_model_parameters_composite_model():
    """Test VariogramKrigeModelParameters with composite model."""

    params = VariogramKrigeModelParameters(
        model=["bessel", "exponential"],
        sill=1.0,
        nugget=0.1,
        correlation_range=0.5,
        hole_effect_range=0.3,
        decay_power=1.5,
    )

    assert params.model == ["bessel", "exponential"]
    assert params.sill == 1.0
    assert params.nugget == 0.1
    assert params.correlation_range == 0.5
    assert params.hole_effect_range == 0.3
    assert params.decay_power == 1.5


def test_variogram_krige_model_parameters_invalid_values():
    """Test VariogramKrigeModelParameters with invalid parameter values."""

    # Invalid sill (must be > 0)
    with pytest.raises(ValidationError):
        VariogramKrigeModelParameters(model="exponential", sill=0.0)

    # Invalid nugget (must be >= 0)
    with pytest.raises(ValidationError):
        VariogramKrigeModelParameters(model="exponential", nugget=-0.1)

    # Invalid correlation_range (must be > 0)
    with pytest.raises(ValidationError):
        VariogramKrigeModelParameters(model="exponential", correlation_range=0.0)

    # Invalid decay_power (must be > 0)
    with pytest.raises(ValidationError):
        VariogramKrigeModelParameters(model="exponential", decay_power=0.0)


# ==================================================================================================
# Test ValidateKrigingClass
# -------------------------
def test_validate_kriging_class_valid():
    """Test ValidateKrigingClass with valid arguments."""

    mesh_df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6],
            "latitude": [46.0, 46.2, 46.4],
            "area": [1.5, 2.0, 1.8],
        }
    )

    kriging_params = KrigingParameters(aspect_ratio=0.1, k_min=3, k_max=10, search_radius=5.0)

    variogram_params = VariogramKrigeModelParameters(
        model="exponential", sill=1.0, nugget=0.1, correlation_range=0.5
    )

    args = ValidateKrigingClass(
        mesh=mesh_df,
        kriging_params=kriging_params,
        variogram_params=variogram_params,
        coordinate_names=("longitude", "latitude"),
    )

    pd.testing.assert_frame_equal(args.mesh, mesh_df)
    assert args.kriging_params == kriging_params
    assert args.variogram_params == variogram_params
    assert args.coordinate_names == ("longitude", "latitude")


def test_validate_kriging_class_invalid_mesh():
    """Test ValidateKrigingClass with invalid mesh DataFrame."""

    # Invalid mesh (missing coordinates and area/fraction)
    mesh_df = pd.DataFrame({"other_column": [1, 2, 3]})

    kriging_params = KrigingParameters(search_radius=5.0)
    variogram_params = VariogramKrigeModelParameters(model="exponential")

    with pytest.raises(ValueError):
        ValidateKrigingClass(
            mesh=mesh_df,
            kriging_params=kriging_params,
            variogram_params=variogram_params,
            coordinate_names=("longitude", "latitude"),
        )


# ==================================================================================================
# Test ValidateMeshCropArgs
# -------------------------
def test_validate_mesh_crop_args_valid():
    """Test ValidateMeshCropArgs with valid arguments."""

    def dummy_crop_function(mesh, **kwargs):
        return mesh

    args = ValidateMeshCropArgs(
        crop_function=dummy_crop_function,
        coordinate_names=("longitude", "latitude"),
        kwargs={"param1": "value1", "param2": 42},
    )

    assert args.crop_function == dummy_crop_function
    assert args.coordinate_names == ("longitude", "latitude")
    assert args.kwargs == {"param1": "value1", "param2": 42}


def test_validate_mesh_crop_args_empty_kwargs():
    """Test ValidateMeshCropArgs with empty kwargs."""

    def dummy_crop_function(mesh, **kwargs):
        return mesh

    args = ValidateMeshCropArgs(
        crop_function=dummy_crop_function, coordinate_names=("x", "y"), kwargs={}
    )

    assert args.crop_function == dummy_crop_function
    assert args.coordinate_names == ("x", "y")
    assert args.kwargs == {}


# ==================================================================================================
# Integration tests
# -----------------
def test_kriging_validators_integration():
    """Test integration of kriging validators."""

    # Create test mesh
    mesh_df = pd.DataFrame(
        {
            "longitude": np.linspace(-124.5, -123.0, 20),
            "latitude": np.linspace(45.8, 47.0, 20),
            "x": np.linspace(-0.5, 1.0, 20),
            "y": np.linspace(-0.2, 1.2, 20),
            "area": np.ones(20) * 1.5,
        }
    )

    # Create kriging parameters
    kriging_params = KrigingParameters(aspect_ratio=0.001, k_min=3, k_max=8, search_radius=5.0)

    # Create variogram parameters
    variogram_params = VariogramKrigeModelParameters(
        model=["bessel", "exponential"],
        sill=0.91,
        nugget=0.0,
        correlation_range=0.007,
        hole_effect_range=0.001,  # Changed from 0.0 to 0.001 to satisfy validator constraints
        decay_power=1.5,
    )

    # Test kriging class validation
    kriging_class = ValidateKrigingClass(
        mesh=mesh_df,
        kriging_params=kriging_params,
        variogram_params=variogram_params,
        coordinate_names=("longitude", "latitude"),
    )

    assert len(kriging_class.mesh) == 20
    assert kriging_class.kriging_params.k_min == 3
    assert kriging_class.variogram_params.sill == 0.91

    # Test mesh crop args validation
    def sample_crop_function(mesh, transects, **kwargs):
        return mesh[:10]  # Return subset of mesh

    crop_args = ValidateMeshCropArgs(
        crop_function=sample_crop_function,
        coordinate_names=("longitude", "latitude"),
        kwargs={"num_nearest_transects": 3, "mesh_buffer_distance": 2.5},
    )

    assert crop_args.crop_function == sample_crop_function
    assert crop_args.kwargs["num_nearest_transects"] == 3


def test_kriging_validators_realistic_parameters():
    """Test kriging validators with realistic FEAT hake parameters."""

    # Realistic mesh for hake survey
    n_points = 1000
    mesh_df = pd.DataFrame(
        {
            "longitude": np.random.uniform(-125.0, -122.0, n_points),
            "latitude": np.random.uniform(45.0, 49.0, n_points),
            "x": np.random.uniform(-1.0, 1.0, n_points),
            "y": np.random.uniform(-1.0, 1.0, n_points),
            "area": np.ones(n_points) * 6.25,  # 2.5 nmi x 2.5 nmi
        }
    )

    # Realistic kriging parameters from FEAT workflow
    kriging_params = KrigingParameters(
        aspect_ratio=0.001, k_min=3, k_max=10, search_radius=0.021  # 3 * correlation_range
    )

    # Realistic variogram parameters
    variogram_params = VariogramKrigeModelParameters(
        model=["exponential", "bessel"],
        sill=0.91,
        nugget=0.0,
        correlation_range=0.007,
        hole_effect_range=0.0,
        decay_power=1.5,
        enhance_semivariance=False,
    )

    # Validate the complete kriging setup
    kriging_setup = ValidateKrigingClass(
        mesh=mesh_df,
        kriging_params=kriging_params,
        variogram_params=variogram_params,
        coordinate_names=("x", "y"),
    )

    assert len(kriging_setup.mesh) == n_points
    assert kriging_setup.kriging_params.search_radius == 0.021
    assert kriging_setup.variogram_params.model == ["exponential", "bessel"]
