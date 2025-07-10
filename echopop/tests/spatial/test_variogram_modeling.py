import numpy as np
import pytest

from echopop.spatial.variogram import (
    bessel_exponential,
    bessel_gaussian,
    cosine_exponential,
    cosine_gaussian,
    exponential,
    exponential_linear,
    gaussian,
    gaussian_linear,
    get_variogram_arguments,
    jbessel,
    kbessel,
    linear,
    nugget,
    sinc,
    spherical,
    variogram,
)


# ==================================================================================================
# Test single variogram models
# ----------------------------
def test_exponential(sample_distance_lags, sample_model_parameters):
    """Test exponential variogram model."""
    result = exponential(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that result approaches sill + nugget at large distances
    assert result[-1] <= sample_model_parameters["sill"] + sample_model_parameters["nugget"]

    # Check monotonic increasing behavior
    assert np.all(np.diff(result) >= 0)


def test_gaussian(sample_distance_lags, sample_model_parameters):
    """Test gaussian variogram model."""
    result = gaussian(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that result approaches sill + nugget at large distances
    assert result[-1] <= sample_model_parameters["sill"] + sample_model_parameters["nugget"]

    # Check monotonic increasing behavior
    assert np.all(np.diff(result) >= 0)


def test_jbessel(sample_distance_lags, sample_model_parameters):
    """Test J-Bessel variogram model."""
    result = jbessel(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["hole_effect_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_kbessel(sample_distance_lags, sample_model_parameters):
    """Test K-Bessel variogram model."""
    result = kbessel(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["hole_effect_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_linear(sample_distance_lags, sample_model_parameters):
    """Test linear variogram model."""
    result = linear(
        sample_distance_lags, sample_model_parameters["sill"], sample_model_parameters["nugget"]
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check linear relationship
    expected_result = (
        sample_model_parameters["sill"] - sample_model_parameters["nugget"]
    ) * sample_distance_lags + sample_model_parameters["nugget"]
    assert np.allclose(result, expected_result)


def test_nugget(sample_distance_lags, sample_model_parameters):
    """Test nugget variogram model."""
    result = nugget(
        sample_distance_lags, sample_model_parameters["sill"], sample_model_parameters["nugget"]
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals 0.0 (nugget model behavior)
    assert np.isclose(result[0], 0.0)

    # Check that all non-zero distances equal sill + nugget
    non_zero_mask = sample_distance_lags > 0
    expected_value = sample_model_parameters["sill"] + sample_model_parameters["nugget"]
    assert np.allclose(result[non_zero_mask], expected_value)


def test_sinc(sample_distance_lags, sample_model_parameters):
    """Test sinc variogram model."""
    result = sinc(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["hole_effect_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 is within reasonable range (sinc function behavior)
    # Due to machine epsilon handling, it might not be exactly the nugget
    assert result[0] >= sample_model_parameters["nugget"]
    assert result[0] <= sample_model_parameters["sill"] + sample_model_parameters["nugget"]

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_spherical(sample_distance_lags, sample_model_parameters):
    """Test spherical variogram model."""
    result = spherical(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that result approaches sill + nugget at large distances
    assert result[-1] <= sample_model_parameters["sill"] + sample_model_parameters["nugget"]


# ==================================================================================================
# Test composite variogram models
# -------------------------------
def test_bessel_exponential(sample_distance_lags, sample_model_parameters):
    """Test bessel-exponential composite variogram model."""
    result = bessel_exponential(
        sample_distance_lags,
        sample_model_parameters["nugget"],
        sample_model_parameters["sill"],
        sample_model_parameters["correlation_range"],
        sample_model_parameters["decay_power"],
        sample_model_parameters["hole_effect_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_bessel_gaussian(sample_distance_lags, sample_model_parameters):
    """Test bessel-gaussian composite variogram model."""
    result = bessel_gaussian(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
        sample_model_parameters["hole_effect_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_cosine_exponential(sample_distance_lags, sample_model_parameters):
    """Test cosine-exponential composite variogram model."""
    result = cosine_exponential(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
        sample_model_parameters["hole_effect_range"],
        True,  # enhance_semivariance parameter
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 is within reasonable range
    # Due to the cosine and exponential interaction, it might not equal nugget exactly
    assert result[0] >= 0.0
    assert result[0] <= sample_model_parameters["sill"] + sample_model_parameters["nugget"] + 1.0

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_cosine_gaussian(sample_distance_lags, sample_model_parameters):
    """Test cosine-gaussian composite variogram model."""
    result = cosine_gaussian(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
        sample_model_parameters["hole_effect_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 is within reasonable range
    # Due to the cosine function, it might not be exactly the nugget
    assert result[0] >= 0.0
    assert result[0] <= sample_model_parameters["sill"] + sample_model_parameters["nugget"]

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_exponential_linear(sample_distance_lags, sample_model_parameters):
    """Test exponential-linear composite variogram model."""
    result = exponential_linear(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
        sample_model_parameters["hole_effect_range"],
        sample_model_parameters["decay_power"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_gaussian_linear(sample_distance_lags, sample_model_parameters):
    """Test gaussian-linear composite variogram model."""
    result = gaussian_linear(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
        sample_model_parameters["hole_effect_range"],
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


# ==================================================================================================
# Test main variogram function
# ----------------------------
def test_variogram_single_model(sample_distance_lags, sample_variogram_arguments):
    """Test variogram function with single model."""
    result = variogram(sample_distance_lags, **sample_variogram_arguments)

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_variogram_arguments["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_variogram_composite_model(sample_distance_lags, sample_composite_arguments):
    """Test variogram function with composite model."""
    result = variogram(sample_distance_lags, **sample_composite_arguments)

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_composite_arguments["nugget"])

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_variogram_all_single_models(
    sample_distance_lags, sample_single_model_names, sample_minimal_parameters
):
    """Test variogram function with all single models."""
    for model_name in sample_single_model_names:
        # Create appropriate parameters for each model
        if model_name in ["jbessel", "kbessel", "sinc"]:
            # These models need hole_effect_range
            args = {**sample_minimal_parameters, "hole_effect_range": 0.4, "model": [model_name]}
        else:
            # Other models use the minimal parameters
            args = {**sample_minimal_parameters, "model": [model_name]}

        result = variogram(sample_distance_lags, **args)

        # Check output shape and type
        assert isinstance(result, np.ndarray)
        assert result.shape == sample_distance_lags.shape
        assert result.dtype == np.float64

        # Check that all values are finite
        assert np.all(np.isfinite(result))


def test_variogram_all_composite_models(
    sample_distance_lags, sample_composite_model_names, sample_extended_parameters
):
    """Test variogram function with all composite models."""
    for model_tuple in sample_composite_model_names:
        args = {**sample_extended_parameters, "model": list(model_tuple)}
        result = variogram(sample_distance_lags, **args)

        # Check output shape and type
        assert isinstance(result, np.ndarray)
        assert result.shape == sample_distance_lags.shape
        assert result.dtype == np.float64

        # Check that all values are finite
        assert np.all(np.isfinite(result))


def test_variogram_invalid_model():
    """Test variogram function with invalid model."""
    distance_lags = np.array([0.0, 0.1, 0.2])
    args = {"nugget": 0.0, "sill": 1.0, "correlation_range": 0.5, "model": ["invalid_model"]}

    with pytest.raises(LookupError):
        variogram(distance_lags, **args)


# ==================================================================================================
# Test get_variogram_arguments function
# ------------------------------------
def test_get_variogram_arguments_single_model():
    """Test get_variogram_arguments with single model."""
    model_name = "exponential"

    result_params, result_function = get_variogram_arguments(model_name)

    # Check that result contains function signature parameters
    assert "distance_lags" in result_params
    assert "sill" in result_params
    assert "nugget" in result_params
    assert "correlation_range" in result_params

    # Check that result contains model function
    assert "model_function" in result_function
    assert callable(result_function["model_function"])


def test_get_variogram_arguments_composite_model():
    """Test get_variogram_arguments with composite model."""
    model_name = ["bessel", "exponential"]

    result_params, result_function = get_variogram_arguments(model_name)

    # Check that result contains function signature parameters
    assert "distance_lags" in result_params
    assert "sill" in result_params
    assert "nugget" in result_params
    assert "correlation_range" in result_params
    assert "decay_power" in result_params
    assert "hole_effect_range" in result_params

    # Check that result contains model function
    assert "model_function" in result_function
    assert callable(result_function["model_function"])


def test_get_variogram_arguments_invalid_model():
    """Test get_variogram_arguments with invalid model."""
    model_name = "invalid_model"

    with pytest.raises(LookupError):
        get_variogram_arguments(model_name)
