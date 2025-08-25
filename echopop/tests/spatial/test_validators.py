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
    """Test valid variogram model parameters."""
    params = VariogramModelParameters(
        correlation_range=0.5,
        sill=1.0,
        nugget=0.1,
        smoothness_parameter=1.5,
        shape_parameter=2.0,
        power_exponent=1.5,
    )

    assert params.correlation_range == 0.5
    assert params.sill == 1.0
    assert params.nugget == 0.1
    assert params.smoothness_parameter == 1.5
    assert params.shape_parameter == 2.0
    assert params.power_exponent == 1.5


def test_variogram_model_parameters_sill_nugget_validation():
    """Test sill-nugget relationship validation."""
    # Valid case: sill > nugget
    params = VariogramModelParameters(sill=1.0, nugget=0.1)
    assert params.sill == 1.0
    assert params.nugget == 0.1

    # Invalid case: sill <= nugget
    with pytest.raises(ValidationError, match="sill.*must be greater than nugget"):
        VariogramModelParameters(sill=0.5, nugget=0.6)

    with pytest.raises(ValidationError, match="sill.*must be greater than nugget"):
        VariogramModelParameters(sill=0.5, nugget=0.5)


def test_variogram_model_parameters_decay_power_validation():
    """Test decay_power parameter validation."""
    # Valid values
    VariogramModelParameters(decay_power=0.5)
    VariogramModelParameters(decay_power=1.0)
    VariogramModelParameters(decay_power=2.0)

    # Invalid values
    with pytest.raises(ValidationError, match="decay_power must be in interval \\(0, 2\\]"):
        VariogramModelParameters(decay_power=0.0)

    with pytest.raises(ValidationError, match="decay_power must be in interval \\(0, 2\\]"):
        VariogramModelParameters(decay_power=2.5)


def test_variogram_model_parameters_power_exponent_validation():
    """Test power_exponent parameter validation."""
    # Valid values
    VariogramModelParameters(power_exponent=0.5)
    VariogramModelParameters(power_exponent=1.0)
    VariogramModelParameters(power_exponent=1.9)

    # Invalid values
    with pytest.raises(ValidationError, match="power_exponent must be in interval \\(0, 2\\)"):
        VariogramModelParameters(power_exponent=0.0)

    with pytest.raises(ValidationError, match="power_exponent must be in interval \\(0, 2\\)"):
        VariogramModelParameters(power_exponent=2.0)

    with pytest.raises(ValidationError, match="power_exponent must be in interval \\(0, 2\\)"):
        VariogramModelParameters(power_exponent=2.5)


def test_variogram_model_parameters_smoothness_parameter_validation():
    """Test smoothness_parameter validation."""
    # Valid values
    VariogramModelParameters(smoothness_parameter=0.5)
    VariogramModelParameters(smoothness_parameter=1.5)
    VariogramModelParameters(smoothness_parameter=10.0)

    # Invalid values
    with pytest.raises(
        ValidationError, match="smoothness_parameter must be in interval \\(0, 10\\]"
    ):
        VariogramModelParameters(smoothness_parameter=0.0)

    with pytest.raises(
        ValidationError, match="smoothness_parameter must be in interval \\(0, 10\\]"
    ):
        VariogramModelParameters(smoothness_parameter=15.0)


def test_variogram_model_parameters_shape_parameter_validation():
    """Test shape_parameter validation."""
    # Valid values
    VariogramModelParameters(shape_parameter=1.0)
    VariogramModelParameters(shape_parameter=50.0)
    VariogramModelParameters(shape_parameter=100.0)

    # Invalid values
    with pytest.raises(ValidationError, match="shape_parameter must be in interval \\(0, 100\\]"):
        VariogramModelParameters(shape_parameter=0.0)

    with pytest.raises(ValidationError, match="shape_parameter must be in interval \\(0, 100\\]"):
        VariogramModelParameters(shape_parameter=150.0)


def test_variogram_model_parameters_negative_values():
    """Test validation of negative values."""
    # correlation_range must be positive
    with pytest.raises(ValidationError):
        VariogramModelParameters(correlation_range=-0.5)

    # sill must be positive
    with pytest.raises(ValidationError):
        VariogramModelParameters(sill=-1.0)

    # nugget cannot be negative
    with pytest.raises(ValidationError):
        VariogramModelParameters(nugget=-0.1)

    # hole_effect_range must be positive
    with pytest.raises(ValidationError):
        VariogramModelParameters(hole_effect_range=-0.3)


def test_variogram_model_parameters_optional_fields():
    """Test that all fields are optional."""
    # Should work with no parameters
    params = VariogramModelParameters()
    assert params.correlation_range is None
    assert params.sill is None
    assert params.nugget is None


# ==================================================================================================
# Test ValidateVariogramClass
# ---------------------------
def test_validate_variogram_class_valid():
    """Test valid variogram class parameters."""
    params = ValidateVariogramClass(
        coordinate_names=("x", "y"),
        lag_resolution=0.1,
        n_lags=30,
    )

    assert params.coordinate_names == ("x", "y")
    assert params.lag_resolution == 0.1
    assert params.n_lags == 30


def test_validate_variogram_class_invalid_lag_resolution():
    """Test invalid lag_resolution values."""
    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("x", "y"),
            lag_resolution=0.0,
            n_lags=30,
        )

    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("x", "y"),
            lag_resolution=-0.1,
            n_lags=30,
        )


def test_validate_variogram_class_invalid_n_lags():
    """Test invalid n_lags values."""
    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("x", "y"),
            lag_resolution=0.1,
            n_lags=0,
        )

    with pytest.raises(ValidationError):
        ValidateVariogramClass(
            coordinate_names=("x", "y"),
            lag_resolution=0.1,
            n_lags=-5,
        )


# ==================================================================================================
# Test ValidateEmpiricalVariogramArgs
# -----------------------------------
def test_validate_empirical_variogram_args_valid(sample_transect_df):
    """Test valid empirical variogram arguments."""
    params = ValidateEmpiricalVariogramArgs(
        azimuth_angle_threshold=90.0,
        azimuth_filter=True,
        coordinate_names=("x", "y"),
        data=sample_transect_df,
        force_lag_zero=True,
        variable="biomass_density",
    )

    assert params.azimuth_angle_threshold == 90.0
    assert params.azimuth_filter is True
    assert params.coordinate_names == ("x", "y")
    assert params.variable == "biomass_density"


def test_validate_empirical_variogram_args_missing_columns():
    """Test validation with missing columns."""
    # Create DataFrame without biomass_density column
    data = pd.DataFrame(
        {
            "x": [1, 2, 3],
            "y": [1, 2, 3],
            "other_var": [10, 20, 30],
        }
    )

    with pytest.raises(KeyError, match="missing the defined column"):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=90.0,
            azimuth_filter=False,
            coordinate_names=("x", "y"),
            data=data,
            force_lag_zero=True,
            variable="biomass_density",
        )


def test_validate_empirical_variogram_args_missing_coordinates():
    """Test validation with missing coordinate columns."""
    # Create DataFrame without x column
    data = pd.DataFrame(
        {
            "y": [1, 2, 3],
            "biomass_density": [10, 20, 30],
        }
    )

    with pytest.raises(ValidationError, match="Transect DataFrame requires either paired"):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=90.0,
            azimuth_filter=False,
            coordinate_names=("x", "y"),
            data=data,
            force_lag_zero=True,
            variable="biomass_density",
        )


def test_validate_empirical_variogram_args_invalid_azimuth_threshold():
    """Test invalid azimuth angle threshold values."""
    data = pd.DataFrame(
        {
            "x": [1, 2, 3],
            "y": [1, 2, 3],
            "biomass_density": [10, 20, 30],
        }
    )

    with pytest.raises(ValidationError):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=-10.0,
            azimuth_filter=True,
            coordinate_names=("x", "y"),
            data=data,
            force_lag_zero=True,
            variable="biomass_density",
        )

    with pytest.raises(ValidationError):
        ValidateEmpiricalVariogramArgs(
            azimuth_angle_threshold=190.0,
            azimuth_filter=True,
            coordinate_names=("x", "y"),
            data=data,
            force_lag_zero=True,
            variable="biomass_density",
        )


# ==================================================================================================
# Test ValidateFitVariogramArgs
# -----------------------------
def test_validate_fit_variogram_args_single_model():
    """Test validation of fit variogram args with single model."""
    params = Parameters()
    params.add_many(
        ("nugget", 0.0, True, 0.0, None),
        ("sill", 1.0, True, 0.0, None),
        ("correlation_range", 0.5, True, 0.0, None),
    )

    args = ValidateFitVariogramArgs(
        model="exponential",
        model_parameters=params,
        optimizer_kwargs={},
    )

    assert args.model == "exponential"
    assert isinstance(args.model_parameters, Parameters)
    assert args.optimizer_kwargs == {}


def test_validate_fit_variogram_args_composite_model():
    """Test validation of fit variogram args with composite model."""
    params = Parameters()
    params.add_many(
        ("nugget", 0.0, True, 0.0, None),
        ("sill", 1.0, True, 0.0, None),
        ("correlation_range", 0.5, True, 0.0, None),
        ("hole_effect_range", 0.3, True, 0.0, None),
    )

    args = ValidateFitVariogramArgs(
        model=["bessel", "exponential"],
        model_parameters=params,
        optimizer_kwargs={"max_nfev": 500},
    )

    assert args.model == ["bessel", "exponential"]
    assert isinstance(args.model_parameters, Parameters)
    assert args.optimizer_kwargs == {"max_nfev": 500}


def test_validate_fit_variogram_args_invalid_model():
    """Test validation with invalid model name."""
    params = Parameters()
    params.add("nugget", 0.0)

    with pytest.raises(LookupError):
        ValidateFitVariogramArgs(
            model="invalid_model",
            model_parameters=params,
            optimizer_kwargs={},
        )
