import numpy as np
from lmfit import Parameters

from echopop.geostatistics import fit_variogram


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


def test_fit_variogram_new_single_models(sample_simple_variogram_data, sample_simple_parameters):
    """Test variogram fitting with new single models."""
    lags, lag_counts, gamma = sample_simple_variogram_data

    # Simple optimization parameters
    opt_params = {"max_nfev": 50, "ftol": 1e-06}

    # Test new single models
    new_models = ["cubic", "pentaspherical"]  # Remove quadratic - it needs shape_parameter

    for model_name in new_models:
        model = [model_name]

        fitted_params, initial_cost, final_cost = fit_variogram(
            lags, lag_counts, gamma, sample_simple_parameters, model, opt_params
        )

        # Check that we get fitted parameters back
        assert isinstance(fitted_params, dict)
        assert "nugget" in fitted_params
        assert "sill" in fitted_params
        assert "correlation_range" in fitted_params
        assert isinstance(initial_cost, (int, float))
        assert isinstance(final_cost, (int, float))

    # Test quadratic separately with required shape_parameter
    quadratic_params = Parameters()
    quadratic_params.add_many(
        ("nugget", 0.0, True, 0.0, None),
        ("sill", 1.0, True, 0.0, None),
        ("correlation_range", 0.5, True, 0.0, None),
        ("shape_parameter", 2.0, True, 0.0, None),  # Required for quadratic
    )

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, quadratic_params, ["quadratic"], opt_params
    )

    # Check that we get fitted parameters back
    assert isinstance(fitted_params, dict)
    assert "nugget" in fitted_params
    assert "sill" in fitted_params
    assert "correlation_range" in fitted_params
    assert "shape_parameter" in fitted_params
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))


def test_fit_variogram_matern_model():
    """Test variogram fitting with Matern model (requires smoothness parameter)."""
    # Create simple test data
    lags = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    lag_counts = np.array([100, 150, 200, 180, 160, 140])
    gamma = np.array([0.0, 0.2, 0.35, 0.45, 0.52, 0.58])

    # Create parameters with smoothness parameter for Matern
    params = Parameters()
    params.add_many(
        ("nugget", 0.1, True, 0.0, None),
        ("sill", 0.6, True, 0.0, None),
        ("correlation_range", 0.3, True, 0.0, None),
        ("smoothness_parameter", 1.5, True, 0.1, 10.0),
    )

    opt_params = {"max_nfev": 50, "ftol": 1e-06}
    model = ["matern"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, params, model, opt_params
    )

    # Check that we get fitted parameters back including smoothness
    assert isinstance(fitted_params, dict)
    assert "nugget" in fitted_params
    assert "sill" in fitted_params
    assert "correlation_range" in fitted_params
    assert "smoothness_parameter" in fitted_params
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))


def test_fit_variogram_power_model():
    """Test variogram fitting with Power model (requires power exponent)."""
    # Create simple test data
    lags = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    lag_counts = np.array([100, 150, 200, 180, 160, 140])
    gamma = np.array([0.0, 0.15, 0.25, 0.32, 0.38, 0.43])

    # Create parameters with power exponent for Power model
    params = Parameters()
    params.add_many(
        ("nugget", 0.05, True, 0.0, None),
        ("sill", 0.5, True, 0.0, None),
        ("power_exponent", 1.2, True, 0.1, 2.0),
    )

    opt_params = {"max_nfev": 50, "ftol": 1e-06}
    model = ["power"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, params, model, opt_params
    )

    # Check that we get fitted parameters back including power exponent
    assert isinstance(fitted_params, dict)
    assert "nugget" in fitted_params
    assert "sill" in fitted_params
    assert "power_exponent" in fitted_params
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))


def test_fit_variogram_quadratic_model():
    """Test variogram fitting with Quadratic model (requires shape parameter)."""
    # Create simple test data
    lags = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    lag_counts = np.array([100, 150, 200, 180, 160, 140])
    gamma = np.array([0.0, 0.18, 0.32, 0.42, 0.48, 0.52])

    # Create parameters with shape parameter for Quadratic model
    params = Parameters()
    params.add_many(
        ("nugget", 0.08, True, 0.0, None),
        ("sill", 0.55, True, 0.0, None),
        ("correlation_range", 0.4, True, 0.0, None),
        ("shape_parameter", 0.8, True, 0.1, 2.0),
    )

    opt_params = {"max_nfev": 50, "ftol": 1e-06}
    model = ["quadratic"]

    fitted_params, initial_cost, final_cost = fit_variogram(
        lags, lag_counts, gamma, params, model, opt_params
    )

    # Check that we get fitted parameters back including shape parameter
    assert isinstance(fitted_params, dict)
    assert "nugget" in fitted_params
    assert "sill" in fitted_params
    assert "correlation_range" in fitted_params
    assert "shape_parameter" in fitted_params
    assert isinstance(initial_cost, (int, float))
    assert isinstance(final_cost, (int, float))


def test_fit_variogram_parameter_validation_new_models():
    """Test parameter validation for new models with specific parameter requirements."""
    lags = np.array([0.0, 0.1, 0.2, 0.3])
    lag_counts = np.array([100, 150, 200, 180])
    gamma = np.array([0.0, 0.2, 0.35, 0.45])

    opt_params = {"max_nfev": 25, "ftol": 1e-06}

    # Test Matern model requires smoothness_parameter
    params_basic = Parameters()
    params_basic.add_many(
        ("nugget", 0.1, True, 0.0, None),
        ("sill", 0.6, True, 0.0, None),
        ("correlation_range", 0.3, True, 0.0, None),
    )

    # This should work with Matern if smoothness_parameter is added
    params_matern = Parameters()
    params_matern.add_many(
        ("nugget", 0.1, True, 0.0, None),
        ("sill", 0.6, True, 0.0, None),
        ("correlation_range", 0.3, True, 0.0, None),
        ("smoothness_parameter", 1.0, True, 0.1, 10.0),
    )

    # Test Power model requires power_exponent
    params_power = Parameters()
    params_power.add_many(
        ("nugget", 0.1, True, 0.0, None),
        ("sill", 0.6, True, 0.0, None),
        ("power_exponent", 1.5, True, 0.1, 2.0),
    )

    # Test Quadratic model requires shape_parameter
    params_quadratic = Parameters()
    params_quadratic.add_many(
        ("nugget", 0.1, True, 0.0, None),
        ("sill", 0.6, True, 0.0, None),
        ("correlation_range", 0.3, True, 0.0, None),
        ("shape_parameter", 1.0, True, 0.1, 2.0),
    )

    # These should work without issues
    for model_name, params in [
        ("matern", params_matern),
        ("power", params_power),
        ("quadratic", params_quadratic),
    ]:
        fitted_params, initial_cost, final_cost = fit_variogram(
            lags, lag_counts, gamma, params, [model_name], opt_params
        )
        assert isinstance(fitted_params, dict)
        assert isinstance(initial_cost, (int, float))
        assert isinstance(final_cost, (int, float))


def test_fit_variogram_optimization_convergence():
    """Test that variogram fitting converges properly with different optimization settings."""
    # Create realistic test data
    lags = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
    lag_counts = np.array([500, 450, 400, 380, 350, 320, 300, 280, 260])
    gamma = np.array([0.0, 0.08, 0.18, 0.28, 0.36, 0.42, 0.47, 0.51, 0.54])

    params = Parameters()
    params.add_many(
        ("nugget", 0.05, True, 0.0, 0.2),
        ("sill", 0.5, True, 0.3, 0.8),
        ("correlation_range", 0.2, True, 0.1, 0.5),
    )

    # Test different optimization settings
    opt_settings = [
        {"max_nfev": 200, "ftol": 1e-08},
        {"max_nfev": 100, "ftol": 1e-06},
        {"max_nfev": 50, "ftol": 1e-04},
    ]

    for opt_params in opt_settings:
        fitted_params, initial_cost, final_cost = fit_variogram(
            lags, lag_counts, gamma, params, ["exponential"], opt_params
        )

        # Check convergence
        assert isinstance(fitted_params, dict)
        assert final_cost <= initial_cost  # Should improve
        assert np.isfinite(final_cost)

        # Check parameter reasonableness
        assert 0.0 <= fitted_params["nugget"] <= 0.2
        assert 0.3 <= fitted_params["sill"] <= 0.8
        assert 0.1 <= fitted_params["correlation_range"] <= 0.5
