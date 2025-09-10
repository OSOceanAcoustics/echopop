import numpy as np
import pytest
from lmfit import Parameters

from echopop.nwfsc_feat.spatial import fit_variogram


# ==================================================================================================
# Test fit_variogram
# -----------------
def test_fit_variogram_basic(sample_simple_variogram_data, sample_simple_parameters):
    """Test basic variogram fitting with simple data."""
    lags, lag_counts, gamma = sample_simple_variogram_data

    # Simple optimization parameters
    opt_params = {"max_nfev": 100, "ftol": 1e-06}

    # Test with exponential model
    model = ["exponential"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, sample_simple_parameters, model, opt_params
    )

    # Check that we get fitted parameters back
    assert isinstance(fitted_params, dict)
    assert "nugget" in fitted_params
    assert "sill" in fitted_params
    assert "correlation_range" in fitted_params

    # Check that costs are numeric
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))

    # Check that optimization improved the fit (lower cost)
    assert final_cost <= initial_cost


def test_fit_variogram_composite_model(sample_simple_variogram_data, sample_composite_parameters):
    """Test variogram fitting with composite model."""
    lags, lag_counts, gamma = sample_simple_variogram_data

    # Simple optimization parameters
    opt_params = {"max_nfev": 100, "ftol": 1e-06}

    # Test with composite model
    model = ["exponential", "linear"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, sample_composite_parameters, model, opt_params
    )

    # Check that we get fitted parameters back
    assert isinstance(fitted_params, dict)
    assert "nugget" in fitted_params
    assert "sill" in fitted_params
    assert "correlation_range" in fitted_params

    # Check that costs are numeric
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))


def test_fit_variogram_realistic_data(
    sample_variogram_data, sample_variogram_parameters, sample_optimization_parameters
):
    """Test variogram fitting with realistic data."""
    lags, lag_counts, gamma = sample_variogram_data

    # Test with bessel exponential model
    model = ["exponential", "bessel"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, sample_variogram_parameters, model, sample_optimization_parameters
    )

    # Check that we get fitted parameters back
    assert isinstance(fitted_params, dict)
    assert "nugget" in fitted_params
    assert "sill" in fitted_params
    assert "correlation_range" in fitted_params
    assert "hole_effect_range" in fitted_params
    assert "decay_power" in fitted_params

    # Check that costs are numeric
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))

    # Check that optimization improved the fit (lower cost)
    assert final_cost <= initial_cost


def test_fit_variogram_single_models(sample_simple_variogram_data, sample_simple_parameters):
    """Test variogram fitting with different single models."""
    lags, lag_counts, gamma = sample_simple_variogram_data

    # Simple optimization parameters
    opt_params = {"max_nfev": 50, "ftol": 1e-06}

    # Test different single models
    single_models = ["exponential", "gaussian", "linear", "nugget"]

    for model_name in single_models:
        model = [model_name]

        fitted_params, initial_cost, final_cost = fit_variogram(
            lags, lag_counts, gamma, sample_simple_parameters, model, opt_params
        )

        # Check that we get fitted parameters back
        assert isinstance(fitted_params, dict)
        assert isinstance(initial_cost, (int, float))
        assert isinstance(final_cost, (int, float))


def test_fit_variogram_parameter_bounds(sample_simple_variogram_data):
    """Test that fitted parameters respect bounds."""
    lags, lag_counts, gamma = sample_simple_variogram_data

    # Create parameters with specific bounds
    params = Parameters()
    params.add_many(
        ("nugget", 0.1, True, 0.0, 0.5),  # Bounded between 0 and 0.5
        ("sill", 0.8, True, 0.5, 1.5),  # Bounded between 0.5 and 1.5
        ("correlation_range", 0.3, True, 0.1, 1.0),  # Bounded between 0.1 and 1.0
    )

    # Simple optimization parameters
    opt_params = {"max_nfev": 50, "ftol": 1e-06}

    model = ["exponential"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, params, model, opt_params
    )

    # Check that fitted parameters respect bounds
    assert 0.0 <= fitted_params["nugget"] <= 0.5
    assert 0.5 <= fitted_params["sill"] <= 1.5
    assert 0.1 <= fitted_params["correlation_range"] <= 1.0


def test_fit_variogram_invalid_inputs():
    """Test error handling for invalid inputs."""
    lags = np.array([0.0, 0.1, 0.2])
    lag_counts = np.array([100, 150, 200])
    gamma = np.array([0.0, 0.2, 0.4])

    params = Parameters()
    params.add("nugget", 0.0, True, 0.0, None)

    opt_params = {"max_nfev": 50}

    # Test with empty model
    with pytest.raises(LookupError):
        fit_variogram(lags, lag_counts, gamma, params, [], opt_params)

    # Test with mismatched array lengths
    short_gamma = np.array([0.0, 0.2])
    with pytest.raises((ValueError, IndexError)):
        fit_variogram(lags, lag_counts, short_gamma, params, ["exponential"], opt_params)


def test_fit_variogram_data_types(sample_simple_variogram_data, sample_simple_parameters):
    """Test that fit_variogram returns correct data types."""
    lags, lag_counts, gamma = sample_simple_variogram_data

    opt_params = {"max_nfev": 50, "ftol": 1e-06}
    model = ["exponential"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, sample_simple_parameters, model, opt_params
    )

    # Check return types
    assert isinstance(fitted_params, dict)
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))

    # Check that all fitted parameters are numeric
    for key, value in fitted_params.items():
        assert isinstance(value, (int, float, np.number))
        assert np.isfinite(value)
