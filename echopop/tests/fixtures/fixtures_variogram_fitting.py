import numpy as np
import pytest
from lmfit import Parameters


@pytest.fixture
def sample_variogram_data():
    """Create sample variogram data for fitting tests."""
    lags = np.array(
        [
            0.0,
            0.002,
            0.004,
            0.006,
            0.008,
            0.01,
            0.012,
            0.014,
            0.016,
            0.018,
            0.02,
            0.022,
            0.024,
            0.026,
            0.028,
            0.03,
            0.032,
            0.034,
            0.036,
            0.038,
            0.04,
            0.042,
            0.044,
            0.046,
            0.048,
            0.05,
            0.052,
            0.054,
            0.056,
            0.058,
        ]
    )

    lag_counts = np.array(
        [
            9286,
            14493,
            29152,
            32618,
            25616,
            101514,
            98660,
            75886,
            68026,
            125668,
            124696,
            99036,
            89542,
            117678,
            122378,
            100952,
            94122,
            103910,
            113492,
            97696,
            91952,
            91318,
            99372,
            92484,
            87960,
            85726,
            90442,
            88866,
            87110,
            85014,
        ],
        dtype=np.int64,
    )

    gamma = np.array(
        [
            0.0,
            0.19503051,
            0.5907352,
            0.83452804,
            0.91885088,
            0.97740984,
            0.97675474,
            0.96503398,
            0.95147839,
            0.95787024,
            0.93243558,
            0.93820508,
            0.98178919,
            0.96986873,
            0.96917296,
            0.98028802,
            0.96361323,
            0.94339072,
            0.95201984,
            0.95366344,
            0.97300715,
            0.98992196,
            0.98765544,
            0.98094379,
            0.97437998,
            0.95234293,
            0.91041791,
            0.92755462,
            0.92370374,
            0.96271239,
        ]
    )

    return lags, lag_counts, gamma


@pytest.fixture
def sample_variogram_parameters():
    """Create sample variogram parameters for fitting tests."""
    dict_variogram_params = {
        "range": 0.006,
        "lag_resolution": 0.002,
        "vario": 1.0,
        "corr": 0.0,
        "nugget": 0.0,
        "sill": 0.91,
        "correlation_range": 0.007,
        "decay_power": 1.5,
        "hole_effect_range": 0.0,
        "model": 13.0,
        "ytox_ratio": 1.0,
        "ztox_ratio": 1.0,
        "dim": 1.0,
    }

    # Set up `lmfit` parameters
    variogram_parameters = Parameters()
    variogram_parameters.add_many(
        ("nugget", dict_variogram_params["nugget"], True, 0.0, None),
        ("sill", dict_variogram_params["sill"], True, 0.0, None),
        ("correlation_range", dict_variogram_params["correlation_range"], True, 0.0, None),
        ("hole_effect_range", dict_variogram_params["hole_effect_range"], True, 0.0, None),
        ("decay_power", dict_variogram_params["decay_power"], True, 0.0, None),
    )

    return variogram_parameters


@pytest.fixture
def sample_optimization_parameters():
    """Create sample optimization parameters for fitting tests."""
    return {
        "max_nfev": 500,
        "ftol": 1e-06,
        "gtol": 0.0001,
        "xtol": 1e-06,
        "diff_step": 1e-08,
        "tr_solver": "exact",
        "x_scale": "jac",
        "jac": "3-point",
    }


@pytest.fixture
def sample_simple_variogram_data():
    """Create simple synthetic variogram data for basic tests."""
    lags = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    lag_counts = np.array([100, 150, 200, 180, 160, 140])
    gamma = np.array([0.0, 0.2, 0.4, 0.6, 0.7, 0.8])

    return lags, lag_counts, gamma


@pytest.fixture
def sample_simple_parameters():
    """Create simple parameters for basic variogram fitting tests."""
    params = Parameters()
    params.add_many(
        ("nugget", 0.0, True, 0.0, None),
        ("sill", 1.0, True, 0.0, None),
        ("correlation_range", 0.5, True, 0.0, None),
    )
    return params


@pytest.fixture
def sample_composite_parameters():
    """Create parameters for composite variogram models."""
    params = Parameters()
    params.add_many(
        ("nugget", 0.0, True, 0.0, None),
        ("sill", 1.0, True, 0.0, None),
        ("correlation_range", 0.5, True, 0.0, None),
        ("hole_effect_range", 0.3, True, 0.0, None),
        ("decay_power", 1.5, True, 0.0, None),
    )
    return params
