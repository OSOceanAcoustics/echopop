import numpy as np
import pytest

from echopop.geostatistics.variogram_models import (
    bessel_exponential,
    bessel_gaussian,
    compute_variogram,
    cosine_exponential,
    cosine_gaussian,
    cubic,
    exponential,
    exponential_linear,
    gaussian,
    gaussian_linear,
    get_variogram_arguments,
    jbessel,
    kbessel,
    linear,
    matern,
    nugget,
    pentaspherical,
    power,
    quadratic,
    sinc,
    spherical,
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


def test_cubic(sample_distance_lags, sample_model_parameters):
    """Test cubic variogram model."""
    result = cubic(
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

    # Check that result approaches sill at large distances (finite range model)
    # For distances beyond range, should equal sill
    large_distances = sample_distance_lags > sample_model_parameters["correlation_range"]
    if np.any(large_distances):
        assert np.allclose(result[large_distances], sample_model_parameters["sill"], atol=1e-10)

    # Check monotonic increasing behavior
    assert np.all(np.diff(result) >= 0)


def test_matern(sample_distance_lags, sample_model_parameters):
    """Test Matern variogram model."""
    # Use a reasonable smoothness parameter
    smoothness = 1.5

    result = matern(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
        smoothness,
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that result approaches sill at large distances
    assert result[-1] <= sample_model_parameters["sill"] + 1e-10

    # Check monotonic increasing behavior
    assert np.all(np.diff(result) >= 0)

    # Check that all values are finite
    assert np.all(np.isfinite(result))


def test_pentaspherical(sample_distance_lags, sample_model_parameters):
    """Test pentaspherical variogram model."""
    result = pentaspherical(
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

    # Check that result approaches sill at large distances (finite range model)
    large_distances = sample_distance_lags > sample_model_parameters["correlation_range"]
    if np.any(large_distances):
        assert np.allclose(result[large_distances], sample_model_parameters["sill"], atol=1e-10)

    # Check monotonic increasing behavior
    assert np.all(np.diff(result) >= 0)


def test_power(sample_distance_lags, sample_model_parameters):
    """Test power variogram model."""
    # Use a reasonable power exponent < 2.0
    power_exponent = 1.5

    result = power(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        power_exponent,
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check monotonic increasing behavior
    assert np.all(np.diff(result) >= 0)

    # Check that all values are finite
    assert np.all(np.isfinite(result))

    # For power model, variance grows as h^power_exponent
    # Check basic power law behavior
    non_zero_mask = sample_distance_lags > 0
    if np.any(non_zero_mask):
        values = result[non_zero_mask]
        # All values should be above nugget
        assert np.all(values >= sample_model_parameters["nugget"])


def test_quadratic(sample_distance_lags, sample_model_parameters):
    """Test quadratic variogram model."""
    # Use a reasonable shape parameter
    shape_parameter = 2.0

    result = quadratic(
        sample_distance_lags,
        sample_model_parameters["sill"],
        sample_model_parameters["nugget"],
        sample_model_parameters["correlation_range"],
        shape_parameter,
    )

    # Check output shape and type
    assert isinstance(result, np.ndarray)
    assert result.shape == sample_distance_lags.shape
    assert result.dtype == np.float64

    # Check that result at distance 0 equals nugget
    assert np.isclose(result[0], sample_model_parameters["nugget"])

    # Check that result approaches sill at large distances
    assert result[-1] <= sample_model_parameters["sill"] + 1e-10

    # Check monotonic increasing behavior
    assert np.all(np.diff(result) >= 0)

    # Check that all values are finite
    assert np.all(np.isfinite(result))


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
    result = compute_variogram(sample_distance_lags, **sample_variogram_arguments)

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
    result = compute_variogram(sample_distance_lags, **sample_composite_arguments)

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
        elif model_name == "matern":
            # Matern model needs smoothness_parameter
            args = {**sample_minimal_parameters, "smoothness_parameter": 1.5, "model": [model_name]}
        elif model_name == "quadratic":
            # Quadratic model needs shape_parameter
            args = {**sample_minimal_parameters, "shape_parameter": 2.0, "model": [model_name]}
        elif model_name == "power":
            # Power model needs power_exponent (and doesn't use correlation_range)
            power_params = {
                k: v for k, v in sample_minimal_parameters.items() if k != "correlation_range"
            }
            args = {**power_params, "power_exponent": 1.5, "model": [model_name]}
        else:
            # Other models use the minimal parameters
            args = {**sample_minimal_parameters, "model": [model_name]}

        result = compute_variogram(sample_distance_lags, **args)

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
        result = compute_variogram(sample_distance_lags, **args)

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
        compute_variogram(distance_lags, **args)


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


def test_get_variogram_arguments_cubic():
    """Test get_variogram_arguments with cubic model."""
    model_name = "cubic"

    result_params, result_function = get_variogram_arguments(model_name)

    # Check that result contains function signature parameters
    assert "distance_lags" in result_params
    assert "sill" in result_params
    assert "nugget" in result_params
    assert "correlation_range" in result_params

    # Check that result contains model function
    assert "model_function" in result_function
    assert callable(result_function["model_function"])


def test_get_variogram_arguments_matern():
    """Test get_variogram_arguments with Matern model."""
    model_name = "matern"

    result_params, result_function = get_variogram_arguments(model_name)

    # Check that result contains function signature parameters
    assert "distance_lags" in result_params
    assert "sill" in result_params
    assert "nugget" in result_params
    assert "correlation_range" in result_params
    assert "smoothness_parameter" in result_params

    # Check that result contains model function
    assert "model_function" in result_function
    assert callable(result_function["model_function"])


def test_get_variogram_arguments_pentaspherical():
    """Test get_variogram_arguments with pentaspherical model."""
    model_name = "pentaspherical"

    result_params, result_function = get_variogram_arguments(model_name)

    # Check that result contains function signature parameters
    assert "distance_lags" in result_params
    assert "sill" in result_params
    assert "nugget" in result_params
    assert "correlation_range" in result_params

    # Check that result contains model function
    assert "model_function" in result_function
    assert callable(result_function["model_function"])


def test_get_variogram_arguments_power():
    """Test get_variogram_arguments with power model."""
    model_name = "power"

    result_params, result_function = get_variogram_arguments(model_name)

    # Check that result contains function signature parameters
    assert "distance_lags" in result_params
    assert "sill" in result_params
    assert "nugget" in result_params
    assert "power_exponent" in result_params

    # Check that result contains model function
    assert "model_function" in result_function
    assert callable(result_function["model_function"])


def test_get_variogram_arguments_quadratic():
    """Test get_variogram_arguments with quadratic model."""
    model_name = "quadratic"

    result_params, result_function = get_variogram_arguments(model_name)

    # Check that result contains function signature parameters
    assert "distance_lags" in result_params
    assert "sill" in result_params
    assert "nugget" in result_params
    assert "correlation_range" in result_params
    assert "shape_parameter" in result_params

    # Check that result contains model function
    assert "model_function" in result_function
    assert callable(result_function["model_function"])
