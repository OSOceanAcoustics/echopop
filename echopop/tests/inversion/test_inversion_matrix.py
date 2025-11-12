import numpy as np
import pandas as pd
import pytest
from lmfit import Parameters

import echopop.inversion.inversion_matrix as im
from echopop import inversion
from echopop.inversion import InversionMatrix, InvParameters, estimate_population


# Test base InvParameters
def test_InvParameters(sample_InversionMatrix_parameters):
    """Test InvParameters structure and methods"""

    # Helper function
    def minmaxer(parameter_vals):
        return {
            key: (
                {
                    **kwargs,
                    "value": (
                        (kwargs["value"] - kwargs["min"]) / (kwargs["max"] - kwargs["min"])
                        if kwargs["max"] != kwargs["min"]
                        else 0.0
                    ),
                    "min": 0.0,
                    "max": 1.0,
                }
                if isinstance(kwargs, dict)
                else kwargs
            )
            for key, kwargs in parameter_vals.items()
        }

    # Initialize
    params = InvParameters(sample_InversionMatrix_parameters)

    # Test typing
    assert isinstance(params, InvParameters)

    # Test initial attributes
    assert set(
        ["_unscaled_parameters", "_scaled_parameters", "parameters", "parameter_bounds", "_scaled"]
    ) <= set(params.__dir__())
    # ---- Check that initialized values are correct
    assert params._scaled is False
    assert params.parameters == params._unscaled_parameters
    parameter_limits = {
        key: {"min": value["min"], "max": value["max"]}
        for key, value in sample_InversionMatrix_parameters.items()
    }
    assert params.parameter_bounds == parameter_limits

    # Check properties
    assert params.is_scaled is False
    assert params.bounds == parameter_limits
    parameter_values = {
        key: param["value"] for key, param in sample_InversionMatrix_parameters.items()
    }
    assert params.values == parameter_values

    # Check methods
    assert set(
        ["scale", "unscale", "inverse_transform", "from_series", "update_bounds", "to_lmfit"]
    ) <= set(params.__dir__())

    # Test scale toggling
    # ---- Scale
    params.scale()
    assert params.is_scaled is True
    assert params.parameters == minmaxer(sample_InversionMatrix_parameters)
    scaled_parameter_values = {
        k: v["value"] for k, v in minmaxer(sample_InversionMatrix_parameters).items()
    }
    assert params.values == scaled_parameter_values
    # ---- Unscale
    params.unscale()
    assert params.is_scaled is False
    assert params.parameters == sample_InversionMatrix_parameters
    assert params.values == {k: v["value"] for k, v in sample_InversionMatrix_parameters.items()}

    # Test external unscaling
    assert params.inverse_transform(scaled_parameter_values) == params.values

    # Test type conversion from Series to Dictionary
    s = pd.Series({"length_mean": 1, "g": 1})
    Inv_s = InvParameters.from_series(s)
    assert isinstance(Inv_s, InvParameters)
    assert Inv_s.values == {"length_mean": 1, "g": 1}
    assert Inv_s.bounds == {
        "length_mean": {"min": -np.inf, "max": np.inf},
        "g": {"min": -np.inf, "max": np.inf},
    }

    # Test bound updating
    bounds = {"length_mean": {"min": 0.5, "max": 1.5}}
    Inv_s.update_bounds(bounds)
    assert Inv_s.values == {"length_mean": 1, "g": 1}
    assert Inv_s.bounds == {
        "length_mean": {"min": 0.5, "max": 1.5},
        "g": {"min": -np.inf, "max": np.inf},
    }

    # Test `lmfit.Parameters` generation
    assert isinstance(params.to_lmfit(), Parameters)


# Test functions for mininizer_print_cb
def test_callback_prints_at_key_iterations(capsys, inv_parameters):
    """Test callback prints output at specified iterations."""

    # Get the parameters
    params = inv_parameters.to_lmfit()

    # Iteration: 1
    im.mininizer_print_cb(
        params,
        1,
        np.array([0.5])[0],
        Sv_measured=None,
        center_frequencies=None,
        inv_params=inv_parameters,
        model_settings=None,
    )
    captured = capsys.readouterr()
    assert "Iter: 1" in captured.out
    assert "Q abs [pred. - meas.]:" in captured.out

    # Iteration: 5 [silent]
    im.mininizer_print_cb(
        params,
        5,
        np.array([0.5])[0],
        Sv_measured=None,
        center_frequencies=None,
        inv_params=inv_parameters,
        model_settings=None,
    )
    captured = capsys.readouterr()
    assert captured.out == ""

    # Iteration: 25
    im.mininizer_print_cb(
        params,
        50,
        np.array([0.3])[0],
        Sv_measured=None,
        center_frequencies=None,
        inv_params=inv_parameters,
        model_settings=None,
    )
    captured = capsys.readouterr()
    assert "Iter: 50" in captured.out


# Test functions for monte_carlo_initialize
def test_monte_carlo_basic_functionality(
    inv_parameters, sample_frequencies, model_InversionMatrix_settings
):
    """Test Monte Carlo initialization basic functionality."""

    # Create MC parameters with small realizations for testing
    inv_parameters.simulate_parameter_sets(mc_realizations=3, rng=np.random.default_rng(999))
    lmfit_params = inv_parameters.realizations_to_lmfit()
    Sv_measured = pd.Series([-95.0, -85.0, -75.0, -70.0, -68.0]).to_numpy()

    # Test returns one of the parameter sets
    result = im.monte_carlo_initialize(
        lmfit_params,
        sample_frequencies,
        Sv_measured,
        inv_parameters,
        model_InversionMatrix_settings,
    )

    # Should return one of the parameter sets
    assert result is not None
    np.testing.assert_equal(
        result.valuesdict(), {k: v["value"] for k, v in inv_parameters.realizations[0].items()}
    )


def test_monte_carlo_empty_parameters(sample_frequencies, model_InversionMatrix_settings):
    """Test Monte Carlo with minimal parameter set."""
    # Create minimal valid parameters instead of empty
    minimal_params = InvParameters(
        {"number_density": {"value": 100.0, "min": 1.0, "max": 1000.0, "vary": True}}
    )
    minimal_params.simulate_parameter_sets(mc_realizations=1)
    lmfit_params = minimal_params.realizations_to_lmfit()
    sv_measured = pd.Series([]).to_numpy()

    # Expect an error
    with pytest.raises(TypeError):
        assert im.monte_carlo_initialize(
            lmfit_params,
            sample_frequencies,
            sv_measured,
            minimal_params,
            model_InversionMatrix_settings,
        )


# Test functions for perturb_parameters
def test_perturb_parameters_respects_bounds(inv_parameters):
    """Test parameter perturbation respects bounds."""
    lmfit_params = inv_parameters.to_lmfit()
    # original_values = {name: param.value for name, param in lmfit_params.items()}

    perturbed = im.perturb_parameters(lmfit_params, scale=0.1)

    for name, param in perturbed.items():
        assert param.min <= param.value <= param.max


def test_perturb_parameters_scale_effect(inv_parameters):
    """Test larger scale produces larger perturbations."""
    lmfit_params = inv_parameters.to_lmfit()

    pert_small = im.perturb_parameters(lmfit_params.copy(), scale=0.01)
    pert_large = im.perturb_parameters(lmfit_params.copy(), scale=0.5)

    # Calculate average perturbation magnitude
    small_mag = np.mean([abs(p.value - lmfit_params[name].value) for name, p in pert_small.items()])
    large_mag = np.mean([abs(p.value - lmfit_params[name].value) for name, p in pert_large.items()])

    assert large_mag >= small_mag


def test_perturb_parameters_fixed_unchanged():
    """Test fixed parameters (vary=False) are not perturbed."""
    params = Parameters()
    params.add("vary_param", value=1.0, min=0.5, max=1.5, vary=True)
    params.add("fixed_param", value=2.0, min=1.0, max=3.0, vary=False)

    perturbed = im.perturb_parameters(params, scale=0.5)
    assert perturbed["fixed_param"].value == 2.0


# Test functions for estimate_population
def test_estimate_population(inv_transect_info, inv_interval_info, inv_cells_info):
    """Test population estimation returns expected structure."""

    # Test transect-index
    transect_result = estimate_population(
        inv_transect_info["inverted"],
        inv_transect_info["coords"],
        density_sw=1026.0,
        reference_frequency=120e3,
        aggregate_method="transect",
    )
    # ---- Check: shape
    assert transect_result.shape == (5, 17)
    # ---- Get keys:
    params = list(inv_transect_info["inverted"]["parameters"].loc[1].values.keys())
    # ---- Check: columns
    expected_cols = ["longitude", "latitude", "nasc", "number_density", "biomass_density"] + params
    assert all([col in transect_result.columns for col in expected_cols])
    # ---- Check: index
    assert list(transect_result.index.names) == ["transect_num"]

    # Test interval-index
    interval_result = estimate_population(
        inv_interval_info["inverted"],
        inv_interval_info["coords"],
        density_sw=1026.0,
        reference_frequency=120e3,
        aggregate_method="interval",
    )
    # ---- Check: shape
    assert interval_result.shape == (20, 15)
    # ---- Get keys:
    params = list(inv_interval_info["inverted"]["parameters"].loc[1, 1, 1].values.keys())
    # ---- Check: columns
    expected_cols = ["nasc", "number_density", "biomass_density"] + params
    assert all([col in interval_result.columns for col in expected_cols])
    # ---- Check: index
    assert all(
        [col in ["longitude", "latitude", "interval"] for col in list(interval_result.index.names)]
    )

    # Test cells-index
    cells_result = estimate_population(
        inv_cells_info["inverted"],
        inv_cells_info["coords"],
        density_sw=1026.0,
        reference_frequency=120e3,
        aggregate_method="cells",
    )
    # ---- Check: shape
    assert cells_result.shape == (10, 15)
    # ---- Get keys:
    params = list(inv_cells_info["inverted"]["parameters"].loc[1, 1, -1, -1].values.keys())
    # ---- Check: columns
    expected_cols = ["nasc", "number_density", "biomass_density"] + params
    assert all([col in cells_result.columns for col in expected_cols])
    # ---- Check: index
    assert all(
        [
            col in ["longitude", "latitude", "interval", "layer"]
            for col in list(cells_result.index.names)
        ]
    )


# Test functions for scattering operations
def test_wavenumber():
    """Test wavenumber calculation."""

    # Shared
    SOUND_SPEED = 1500.0

    # Scalar input
    freq = 120e3

    result = inversion.wavenumber(freq, SOUND_SPEED)
    expected = 2 * np.pi * freq / SOUND_SPEED
    np.testing.assert_almost_equal(result, expected)

    # Array input
    freqs = np.array([38e3, 120e3, 200e3])
    result = inversion.wavenumber(freqs, SOUND_SPEED)
    expected = 2 * np.pi * freqs / SOUND_SPEED
    np.testing.assert_array_almost_equal(result, expected)


def test_reflection_coefficient():
    """Test reflection coefficient calculation."""

    # Scalar
    g = 1.02  # Density contrast
    h = 1.01  # Sound speed contrast
    result = inversion.reflection_coefficient(g, h)
    expected = (1 - g * h * h) / (g * h * h) - (g - 1) / g
    np.testing.assert_almost_equal(result, expected)

    # Array
    g = np.array([1.01, 1.02, 1.03])
    h = np.array([1.005, 1.01, 1.015])
    result = inversion.reflection_coefficient(g, h)
    assert len(result) == 3
    assert all(isinstance(val, (float, np.floating)) for val in result)


def test_orientation_average():
    """Test the orientation-averaging operation"""

    # Shared
    THETA_MEAN = 0.0
    THETA_SD = 10.0
    N_THETA = 30
    FREQUENCY_INT = 2e3
    LENGTH_SD_NORM = 0.10
    RNG = np.random.default_rng(987)

    # Generate orientation distribution
    theta_values = np.linspace(THETA_MEAN - 3.1 * THETA_SD, THETA_MEAN + 3.1 * THETA_SD, N_THETA)

    # Helper function
    def _simulate_fbs(center_freqs: np.ndarray):

        # Generate frequency interval
        freqs = inversion.generate_frequency_interval(center_freqs, LENGTH_SD_NORM, FREQUENCY_INT)

        # Sample values
        freq_samples = []
        for f in range(len(freqs)):
            # ---- Get non-NaN
            val = np.sum(~np.isnan(freqs[f]))
            # ---- Generate samples
            samples = RNG.uniform(1e-10, 1e-5, (val, N_THETA, 2))
            # ---- Simulate complex numbers
            freq_samples.extend([[samples[:, :, 0] + 1j * samples[:, :, 1]]])

        # Return the samples
        return freq_samples

    # Single-frequency
    fbs_single_freq = _simulate_fbs(np.array([120e3]))
    # => Test orientation average: Gaussian
    single_freq_gauss = inversion.orientation_average(
        theta_values, fbs_single_freq, THETA_MEAN, THETA_SD, "gaussian"
    )
    # Assert: typing
    assert isinstance(single_freq_gauss, list)
    assert isinstance(single_freq_gauss[0], np.ndarray)
    assert isinstance(single_freq_gauss[0][0], float)
    # Assert: shape
    assert len(single_freq_gauss) == 1
    assert single_freq_gauss[0].shape == (39,)
    # => Test orientation average: Uniform
    single_freq_unif = inversion.orientation_average(
        theta_values, fbs_single_freq, THETA_MEAN, THETA_SD, "uniform"
    )
    # Assert: typing
    assert isinstance(single_freq_unif, list)
    assert isinstance(single_freq_unif[0], np.ndarray)
    assert isinstance(single_freq_unif[0][0], float)
    # Assert: shape
    assert len(single_freq_unif) == 1
    assert single_freq_unif[0].shape == (39,)
    # Assert: values
    assert not all(single_freq_gauss[0] == single_freq_unif[0])

    # Multi-frequency
    fbs_multi_freq = _simulate_fbs(np.array([38e3, 120e3, 200e3]))
    # => Test orientation average: Gaussian
    multi_freq_gauss = inversion.orientation_average(
        theta_values, fbs_multi_freq, THETA_MEAN, THETA_SD, "gaussian"
    )
    # Assert: typing
    assert isinstance(multi_freq_gauss, list)
    assert isinstance(multi_freq_gauss[0], np.ndarray)
    assert isinstance(multi_freq_gauss[0][0], float)
    # Assert: shape
    assert len(multi_freq_gauss) == 3
    assert all([len(multi_freq_gauss[f]) == [13, 39, 63][f] for f in range(len(multi_freq_gauss))])
    # => Test orientation average: Uniform
    multi_freq_unif = inversion.orientation_average(
        theta_values, fbs_multi_freq, THETA_MEAN, THETA_SD, "uniform"
    )
    # Assert: typing
    assert isinstance(multi_freq_unif, list)
    assert isinstance(multi_freq_unif[0], np.ndarray)
    assert isinstance(multi_freq_unif[0][0], float)
    # Assert: shape
    assert len(multi_freq_unif) == 3
    assert all([len(multi_freq_unif[f]) == [13, 39, 63][f] for f in range(len(multi_freq_unif))])
    assert all(
        [
            not np.allclose(multi_freq_gauss[f], multi_freq_unif[f], atol=1e-13)
            for f in range(len(multi_freq_unif))
        ]
    )


def test_length_average():
    """
    Test the length-averaging operation
    """

    # Shared
    LENGTH_MEAN = 0.10
    LENGTH_SD_NORM = 0.10
    N_LENGTH = 100
    LENGTH_RADIUS_RATIO = 10.0
    FREQUENCY_INT = 2e3
    LENGTH_SD_NORM = 0.10
    RNG = np.random.default_rng(987)
    SOUND_SPEED_SW = 1500.0

    # Generate length distribution
    length_values = np.linspace(
        LENGTH_MEAN - 3 * (LENGTH_SD_NORM * LENGTH_MEAN),
        LENGTH_MEAN + 3 * (LENGTH_SD_NORM * LENGTH_MEAN),
        N_LENGTH,
    )

    # Helper function
    def _simulate_fbs(center_freqs: np.ndarray):

        # Generate frequency interval
        freqs = inversion.generate_frequency_interval(center_freqs, LENGTH_SD_NORM, FREQUENCY_INT)

        # Compute the acoustic wavenumbers weighted by target size
        # ---- Center frequencies
        k_c = inversion.wavenumber(center_freqs, SOUND_SPEED_SW)
        # ---- Compute ka (center frequencies)
        ka_c = k_c * LENGTH_MEAN / LENGTH_RADIUS_RATIO

        # Frequency intervals
        # ---- Just wavenumber (`k`)
        k_f = inversion.wavenumber(freqs, SOUND_SPEED_SW)
        # ---- Now `ka`
        ka_f = k_f * LENGTH_MEAN / LENGTH_RADIUS_RATIO

        # Sample values
        freq_samples = []
        for f in range(len(center_freqs)):
            # ---- Get non-NaN
            val = np.sum(~np.isnan(freqs[f]))
            # ---- Generate samples
            samples = RNG.uniform(1e-10, 1e-5, val)
            # ---- Simulate complex numbers
            freq_samples.extend([samples])

        # Return the samples
        return freq_samples, ka_c, ka_f

    # Single frequency entry
    fbs_single_freq, ka_c, ka_f = _simulate_fbs(np.array([38e3]))
    # => Test orientation average: Gaussian
    single_freq_gauss = inversion.length_average(
        length_values,
        ka_f,
        ka_c,
        fbs_single_freq,
        LENGTH_MEAN,
        LENGTH_MEAN * LENGTH_SD_NORM,
        "gaussian",
    )
    # Assert: typing
    assert isinstance(single_freq_gauss, np.ndarray)
    assert isinstance(single_freq_gauss[0], float)
    # Assert: shape
    assert len(single_freq_gauss) == 1
    # => Test orientation average: Uniform
    single_freq_unif = inversion.length_average(
        length_values,
        ka_f,
        ka_c,
        fbs_single_freq,
        LENGTH_MEAN,
        LENGTH_MEAN * LENGTH_SD_NORM,
        "uniform",
    )
    # Assert: typing
    assert isinstance(single_freq_unif, np.ndarray)
    assert isinstance(single_freq_unif[0], float)
    # Assert: shape
    assert len(single_freq_unif) == 1
    # Assert: values
    assert single_freq_gauss != single_freq_unif

    # Multi-frequency
    fbs_multi_freq, ka_c, ka_f = _simulate_fbs(np.array([38e3, 120e3, 200e3]))
    # => Test orientation average: Gaussian
    multi_freq_gauss = inversion.length_average(
        length_values,
        ka_f,
        ka_c,
        fbs_multi_freq,
        LENGTH_MEAN,
        LENGTH_MEAN * LENGTH_SD_NORM,
        "gaussian",
    )
    # Assert: typing
    assert isinstance(multi_freq_gauss, np.ndarray)
    assert all([isinstance(multi_freq_gauss[f], float) for f in range(3)])
    # Assert: shape
    assert len(multi_freq_gauss) == 3
    # => Test orientation average: Uniform
    multi_freq_unif = inversion.length_average(
        length_values,
        ka_f,
        ka_c,
        fbs_multi_freq,
        LENGTH_MEAN,
        LENGTH_MEAN * LENGTH_SD_NORM,
        "uniform",
    )
    # Assert: typing
    assert isinstance(multi_freq_unif, np.ndarray)
    assert all([isinstance(multi_freq_unif[f], float) for f in range(3)])
    # Assert: shape
    assert len(multi_freq_unif) == 3
    assert all(
        [
            not np.allclose(multi_freq_gauss[f], multi_freq_unif[f], atol=1e-13)
            for f in range(len(multi_freq_unif))
        ]
    )


# Integration tests
def test_inv_parameters_workflow(inv_parameters):
    """Test complete parameter handling workflow."""
    original_values = inv_parameters.values.copy()

    # Test scaling
    inv_parameters.scale()
    scaled_values = inv_parameters.values

    # Values should be in [0,1] range after scaling
    for val in scaled_values.values():
        assert 0 <= val <= 1

    # Test unscaling
    inv_parameters.unscale()
    unscaled_values = inv_parameters.values

    # Should recover original values
    for key in original_values:
        np.testing.assert_almost_equal(unscaled_values[key], original_values[key])


def test_mc_parameters_generation(inv_parameters):
    """Test Monte Carlo parameter generation."""

    # Test shape
    inv_parameters.simulate_parameter_sets(mc_realizations=5, rng=np.random.default_rng(999))
    assert len(inv_parameters.realizations) == 5

    # All samples should have same parameter names
    param_names = set(inv_parameters.realizations[0].keys())
    for i in range(1, 5):
        assert set(inv_parameters.realizations[i].keys()) == param_names


def test_inversion_matrix_creation():
    """Test InversionMatrix creation and basic functionality."""
    # Create mock data
    frequencies = [38e3, 120e3]
    n_obs = 5

    # Create MultiIndex columns
    cols = pd.MultiIndex.from_product(
        [["sv_mean", "nasc", "thickness_mean"], frequencies], names=[None, "frequency"]
    )

    # Generate random data
    np.random.seed(131)
    data = np.random.randn(n_obs, len(cols))
    df = pd.DataFrame(data, columns=cols)

    # Simulation settings
    sim_settings = {
        "monte_carlo": True,
        "mc_realizations": 3,
        "scale_parameters": True,
        "minimum_frequency_count": 2,
    }

    # Create InversionMatrix
    inv_matrix = InversionMatrix(df, sim_settings, verbose=False)

    assert inv_matrix is not None
    assert inv_matrix.inversion_method == "scattering_model"
    assert hasattr(inv_matrix, "measurements")
    assert hasattr(inv_matrix, "simulation_settings")
