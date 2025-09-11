import pandas as pd
import pytest
from lmfit import Parameters
from pydantic import ValidationError

from echopop.validators.variogram import (
    ValidateEmpiricalVariogramArgs,
    ValidateFitVariogramArgs,
    ValidateVariogramClass,
    VariogramModelParameters,
)


# ==================================================================================================
# Test VariogramModelParameters
# -----------------------------
def test_variogram_model_parameters_valid():
    """Test VariogramModelParameters with valid parameters."""

    params = VariogramModelParameters(
        sill=1.0,
        nugget=0.1,
        correlation_range=0.5,
        hole_effect_range=0.3,
        decay_power=1.5,
        smoothness_parameter=2.5,
        shape_parameter=5.0,
        power_exponent=1.2,
        enhance_semivariance=True,
    )

    assert params.sill == 1.0
    assert params.nugget == 0.1
    assert params.correlation_range == 0.5
    assert params.hole_effect_range == 0.3
    assert params.decay_power == 1.5
    assert params.smoothness_parameter == 2.5
    assert params.shape_parameter == 5.0
    assert params.power_exponent == 1.2
    assert params.enhance_semivariance is True


def test_variogram_model_parameters_optional():
    """Test VariogramModelParameters with optional parameters."""

    params = VariogramModelParameters(sill=1.0, nugget=0.1)

    assert params.sill == 1.0
    assert params.nugget == 0.1
    assert params.correlation_range is None


def test_variogram_model_parameters_negative_values():
    """Test VariogramModelParameters with invalid negative values."""

    with pytest.raises(ValidationError):
        VariogramModelParameters(sill=-1.0, nugget=0.1)

    with pytest.raises(ValidationError):
        VariogramModelParameters(sill=1.0, nugget=-0.1)

    with pytest.raises(ValidationError):
        VariogramModelParameters(sill=1.0, nugget=0.1, correlation_range=-0.5)


def test_variogram_model_parameters_decay_power_validation():
    """Test decay_power validation."""

    # Valid values
    params = VariogramModelParameters(decay_power=1.5)
    assert params.decay_power == 1.5

    params = VariogramModelParameters(decay_power=2.0)
    assert params.decay_power == 2.0

    # Invalid values
    with pytest.raises(ValidationError):
        VariogramModelParameters(decay_power=0.0)

    with pytest.raises(ValidationError):
        VariogramModelParameters(decay_power=2.5)


def test_variogram_model_parameters_power_exponent_validation():
    """Test power_exponent validation."""

    # Valid values
    params = VariogramModelParameters(power_exponent=1.5)
    assert params.power_exponent == 1.5

    params = VariogramModelParameters(power_exponent=1.99)
    assert params.power_exponent == 1.99

    # Invalid values
    with pytest.raises(ValidationError):
        VariogramModelParameters(power_exponent=0.0)

    with pytest.raises(ValidationError):
        VariogramModelParameters(power_exponent=2.0)

    with pytest.raises(ValidationError):
        VariogramModelParameters(power_exponent=2.5)


def test_variogram_model_parameters_smoothness_parameter_validation():
    """Test smoothness_parameter validation."""

    # Valid values
    params = VariogramModelParameters(smoothness_parameter=0.5)
    assert params.smoothness_parameter == 0.5

    params = VariogramModelParameters(smoothness_parameter=10.0)
    assert params.smoothness_parameter == 10.0

    # Invalid values
    with pytest.raises(ValidationError):
        VariogramModelParameters(smoothness_parameter=0.0)

    with pytest.raises(ValidationError):
        VariogramModelParameters(smoothness_parameter=10.5)


def test_variogram_model_parameters_shape_parameter_validation():
    """Test shape_parameter validation."""

    # Valid values
    params = VariogramModelParameters(shape_parameter=1.0)
    assert params.shape_parameter == 1.0

    params = VariogramModelParameters(shape_parameter=100.0)
    assert params.shape_parameter == 100.0

    # Invalid values
    with pytest.raises(ValidationError):
        VariogramModelParameters(shape_parameter=0.0)

    with pytest.raises(ValidationError):
        VariogramModelParameters(shape_parameter=100.5)


def test_variogram_model_parameters_sill_nugget_relationship():
    """Test sill-nugget relationship validation."""

    # Valid: sill > nugget
    params = VariogramModelParameters(sill=1.0, nugget=0.5)
    assert params.sill == 1.0
    assert params.nugget == 0.5

    # Invalid: sill <= nugget
    with pytest.raises(ValidationError):
        VariogramModelParameters(sill=0.5, nugget=1.0)

    with pytest.raises(ValidationError):
        VariogramModelParameters(sill=0.5, nugget=0.5)


# ==================================================================================================
# Test ValidateVariogramClass
# ---------------------------
def test_validate_variogram_class():
    """Test ValidateVariogramClass with valid parameters."""

    params = ValidateVariogramClass(
        coordinate_names=("longitude", "latitude"), lag_resolution=0.5, n_lags=20
    )

    assert params.coordinate_names == ("longitude", "latitude")
    assert params.lag_resolution == 0.5
    assert params.n_lags == 20


def test_validate_variogram_class_invalid_lag_resolution():
    """Test ValidateVariogramClass with invalid lag_resolution."""

    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("longitude", "latitude"), lag_resolution=0.0, n_lags=20
        )

    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("longitude", "latitude"), lag_resolution=-0.5, n_lags=20
        )


def test_validate_variogram_class_invalid_n_lags():
    """Test ValidateVariogramClass with invalid n_lags."""

    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("longitude", "latitude"), lag_resolution=0.5, n_lags=0
        )

    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("longitude", "latitude"), lag_resolution=0.5, n_lags=-5
        )


# ==================================================================================================
# Test ValidateEmpiricalVariogramArgs
# -----------------------------------
def test_validate_empirical_variogram_args_valid():
    """Test ValidateEmpiricalVariogramArgs with valid data."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6],
            "latitude": [46.0, 46.2, 46.4],
            "biomass_density": [10.5, 15.2, 12.8],
        }
    )

    params = ValidateEmpiricalVariogramArgs(
        azimuth_angle_threshold=90.0,
        azimuth_filter=True,
        coordinate_names=("longitude", "latitude"),
        data=df,
        force_lag_zero=True,
        variable="biomass_density",
    )

    assert params.azimuth_angle_threshold == 90.0
    assert params.azimuth_filter is True
    assert params.coordinate_names == ("longitude", "latitude")
    assert params.force_lag_zero is True
    assert params.variable == "biomass_density"
    pd.testing.assert_frame_equal(params.data, df)


def test_validate_empirical_variogram_args_invalid_azimuth():
    """Test ValidateEmpiricalVariogramArgs with invalid azimuth_angle_threshold."""

    df = pd.DataFrame(
        {"longitude": [-124.0, -123.8], "latitude": [46.0, 46.2], "biomass_density": [10.5, 15.2]}
    )

    with pytest.raises(ValidationError):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=-10.0,
            azimuth_filter=True,
            coordinate_names=("longitude", "latitude"),
            data=df,
            force_lag_zero=True,
            variable="biomass_density",
        )

    with pytest.raises(ValidationError):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=200.0,
            azimuth_filter=True,
            coordinate_names=("longitude", "latitude"),
            data=df,
            force_lag_zero=True,
            variable="biomass_density",
        )


def test_validate_empirical_variogram_args_missing_coordinates():
    """Test ValidateEmpiricalVariogramArgs with missing coordinate columns."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8],
            # Missing latitude column
            "biomass_density": [10.5, 15.2],
        }
    )

    with pytest.raises(ValidationError, match="Transect DataFrame requires either paired"):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=90.0,
            azimuth_filter=True,
            coordinate_names=("longitude", "latitude"),
            data=df,
            force_lag_zero=True,
            variable="biomass_density",
        )


def test_validate_empirical_variogram_args_missing_variable():
    """Test ValidateEmpiricalVariogramArgs with missing variable column."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8],
            "latitude": [46.0, 46.2],
            # Missing biomass_density column
        }
    )

    with pytest.raises(KeyError):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=90.0,
            azimuth_filter=True,
            coordinate_names=("longitude", "latitude"),
            data=df,
            force_lag_zero=True,
            variable="biomass_density",
        )


# ==================================================================================================
# Test ValidateFitVariogramArgs
# -----------------------------
def test_validate_fit_variogram_args_single_model():
    """Test ValidateFitVariogramArgs with single model."""

    params = Parameters()
    params.add("nugget", value=0.1, min=0.0)
    params.add("sill", value=1.0, min=0.0)
    params.add("correlation_range", value=0.5, min=0.0)

    args = ValidateFitVariogramArgs(
        model="exponential", model_parameters=params, optimizer_kwargs={"max_nfev": 500}
    )

    assert args.model == "exponential"
    assert isinstance(args.model_parameters, Parameters)
    assert args.optimizer_kwargs == {"max_nfev": 500}


def test_validate_fit_variogram_args_composite_model():
    """Test ValidateFitVariogramArgs with composite model."""

    params = Parameters()
    params.add("nugget", value=0.1, min=0.0)
    params.add("sill", value=1.0, min=0.0)
    params.add("correlation_range", value=0.5, min=0.0)
    params.add("hole_effect_range", value=0.3, min=0.0)

    args = ValidateFitVariogramArgs(
        model=["bessel", "exponential"], model_parameters=params, optimizer_kwargs={}
    )

    assert args.model == ["bessel", "exponential"]
    assert isinstance(args.model_parameters, Parameters)
    assert args.optimizer_kwargs == {}


def test_validate_fit_variogram_args_invalid_model():
    """Test ValidateFitVariogramArgs with invalid model."""

    params = Parameters()
    params.add("nugget", value=0.1, min=0.0)

    with pytest.raises(LookupError):
        ValidateFitVariogramArgs(
            model="invalid_model", model_parameters=params, optimizer_kwargs={}
        )


# ==================================================================================================
# Integration tests
# -----------------
def test_variogram_validators_integration():
    """Test integration of variogram validators."""

    # Create test data
    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6, -123.4],
            "latitude": [46.0, 46.2, 46.4, 46.6],
            "biomass_density": [10.5, 15.2, 12.8, 18.1],
        }
    )

    # Test empirical variogram validation
    empirical_args = ValidateEmpiricalVariogramArgs(
        azimuth_angle_threshold=180.0,
        azimuth_filter=False,
        coordinate_names=("longitude", "latitude"),
        data=df,
        force_lag_zero=False,
        variable="biomass_density",
    )

    # Test variogram class validation
    class_args = ValidateVariogramClass(
        coordinate_names=("longitude", "latitude"), lag_resolution=0.002, n_lags=30
    )

    # Test model parameters
    model_params = VariogramModelParameters(sill=1.0, nugget=0.1, correlation_range=0.5)

    assert empirical_args.variable == "biomass_density"
    assert class_args.n_lags == 30
    assert model_params.sill == 1.0
