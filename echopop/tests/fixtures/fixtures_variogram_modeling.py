import numpy as np
import pytest


@pytest.fixture
def sample_distance_lags():
    """Create sample distance lags for variogram model testing."""
    return np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])


@pytest.fixture
def sample_model_parameters():
    """Create sample model parameters for variogram testing."""
    return {
        "nugget": 0.1,
        "sill": 0.8,
        "correlation_range": 0.3,
        "hole_effect_range": 0.4,
        "decay_power": 1.5,
        "smoothness_parameter": 1.5,
        "shape_parameter": 2.0,
        "power_exponent": 1.5,
    }


@pytest.fixture
def sample_variogram_arguments():
    """Create sample arguments for variogram function testing."""
    return {
        "nugget": 0.05,
        "sill": 0.9,
        "correlation_range": 0.25,
        "hole_effect_range": 0.35,
        "decay_power": 2.0,
        "model": ["exponential"],
    }


@pytest.fixture
def sample_composite_arguments():
    """Create sample arguments for composite variogram models."""
    return {
        "nugget": 0.1,
        "sill": 0.7,
        "correlation_range": 0.3,
        "hole_effect_range": 0.4,
        "decay_power": 1.8,
        "model": ["bessel", "exponential"],
    }


@pytest.fixture
def sample_zero_distance_lags():
    """Create distance lags including zero for testing nugget effect."""
    return np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])


@pytest.fixture
def sample_large_distance_lags():
    """Create large distance lags for testing model behavior at long distances."""
    return np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0])


@pytest.fixture
def sample_single_model_names():
    """Create list of single model names for testing."""
    return [
        "cubic",
        "exponential",
        "gaussian",
        "jbessel",
        "kbessel",
        "linear",
        "matern",
        "nugget",
        "pentaspherical",
        "power",
        "quadratic",
        "sinc",
        "spherical",
    ]


@pytest.fixture
def sample_composite_model_names():
    """Create list of composite model names for testing."""
    return [
        ("bessel", "exponential"),
        ("bessel", "gaussian"),
        ("cosine", "exponential"),
        ("cosine", "gaussian"),
        ("exponential", "linear"),
        ("gaussian", "linear"),
    ]


@pytest.fixture
def sample_minimal_parameters():
    """Create minimal parameters for simple models."""
    return {"nugget": 0.0, "sill": 1.0, "correlation_range": 0.5}


@pytest.fixture
def sample_extended_parameters():
    """Create extended parameters for complex models."""
    return {
        "nugget": 0.1,
        "sill": 0.9,
        "correlation_range": 0.3,
        "hole_effect_range": 0.4,
        "decay_power": 1.5,
        "ytox_ratio": 1.2,
        "ztox_ratio": 0.8,
        "dim": 2.0,
        "enhance_semivariance": True,
    }
