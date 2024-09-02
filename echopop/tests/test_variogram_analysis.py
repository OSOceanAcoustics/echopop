import echopop.spatial.variogram as esv
from echopop.utils.validate import (
    VariogramOptimize
)
import pytest
import re

@pytest.mark.parametrize(
    "model, expected, exception",
    [
        (
            "exponential",
            ({"distance_lags", "sill", "nugget", "correlation_range"},
             esv.exponential), None
        ),
        (
             "gaussian",
             ({"distance_lags", "sill", "nugget", "correlation_range"},
             esv.gaussian), None
        ),
        (
             "jbessel",
             ({"distance_lags", "sill", "nugget", "hole_effect_range"},
             esv.jbessel), None
        ),
        (
             "kbessel",
             ({"distance_lags", "sill", "nugget", "hole_effect_range"}, 
             esv.kbessel), None
        ),
        (
             "linear",
             ({"distance_lags", "sill", "nugget"}, 
             esv.linear), None
        ),
        (
             "nugget",
             ({"distance_lags", "sill", "nugget"}, 
             esv.nugget), None
        ),
        (
             "sinc",
             ({"distance_lags", "sill", "nugget", "hole_effect_range"}, 
             esv.sinc), None
        ),
        (
             "spherical",
             ({"distance_lags", "sill", "nugget", "correlation_range"}, 
             esv.spherical), None
        ),
        (
             ("bessel", "exponential"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range", 
              "decay_power"}, esv.bessel_exponential), None
        ),
        (
             ("bessel", "gaussian"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"}, 
              esv.bessel_gaussian), None
        ),
        (
             ("cosine", "exponential"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range", 
              "enhance_semivariance"}, esv.cosine_exponential), None
        ),
        (
             ("cosine", "gaussian"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"}, 
              esv.cosine_gaussian), None
        ),
        (
             ("exponential", "linear"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range", 
              "decay_power"}, esv.exponential_linear), None
        ),
        (
             ("gaussian", "linear"),
             ({"distance_lags", "sill", "nugget", "correlation_range", "hole_effect_range"}, 
              esv.gaussian_linear), None
        ),
        (
             "invalid",
             (None, None), "could not be matched to an existing variogram method."
        ),
        (
             ("invalid", "tuple"),
             (None, None), "could not be matched to an existing variogram method."
        ),
        (
             ["invalid", "list"],
             (None, None), "could not be matched to an existing variogram method."
        ),
    ],
    ids=[
        "Exponential arguments",
        "Gaussian arguments",
        "J-Bessel arguments",
        "K-Bessel arguments",
        "Linear arguments",
        "Nugget arguments",
        "Sinc arguments",
        "Spherical arguments",
        "Bessel-Exponential arguments",
        "Bessel-Gaussian arguments",
        "Cosine-Exponential arguments",
        "Cosine-Gaussian arguments",
        "Exponential-Linear arguments",
        "Gaussian-Linear arguments",
        "Unknown single-family variogram model (str)",
        "Unknown composite variogram model (tuple)",
        "Invalid composite variogram model (list)"
    ],
)
def test_get_variogram_arguments(model, expected, exception):

    if exception is not None:
        with pytest.raises(LookupError, match=re.escape(exception)):
            assert esv.get_variogram_arguments(model)
    else:
        # Get arguments
        args, func = esv.get_variogram_arguments(model)
        # Assert argument equality        
        assert set(args.keys()) == set(expected[0])
        # Assert function equality
        assert func["model_function"] == expected[1]

@pytest.mark.parametrize(
    "description",
    ["Test `initialize_optimization_config` translation"],
    ids=["Test `initialize_optimization_config` translation"],
)
def test_initialize_optimization_config(description):

    # -------------------------
    # Mock data [ TEST 1 ]: ASSESS CORRECT `lmfit` PARAMETER/ARGUMENT CONVERSION
    # ---- Create input dictionary
    ECHOPOP_OPTIMIZATION_INPUTS = {
        "max_fun_evaluations": 500,
        "cost_fun_tolerance": 1e-8,
        "solution_tolerance": 1e-8,
        "gradient_tolerance": 1e-8,
        "finite_step_size": 1e-8,
        "trust_region_solver": "exact",
        "x_scale": "jacobian",
        "jacobian_approx": "forward",        
    }

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- Create output
    output = esv.initialize_optimization_config(ECHOPOP_OPTIMIZATION_INPUTS)
    # ---- Create translation dictionary
    TRANSLATION_DICT = {
        "max_fun_evaluations": {
            "name": "max_nfev",
            "value": 500
        },
        "cost_fun_tolerance": {
            "name": "ftol",
            "value": 1e-8
        },
        "solution_tolerance": {
            "name": "xtol",
            "value": 1e-8
        },
        "gradient_tolerance": {
            "name": "gtol",
            "value": 1e-8
        },
        "finite_step_size": {
            "name": "diff_step",
            "value": 1e-8
        },
        "trust_region_solver": {
            "name": "tr_solver",
            "value": "exact"
        },
        "x_scale": {
            "name": "x_scale",
            "value": "jac"
        },
        "jacobian_approx": {
            "name": "jac",
            "value": "2-point"
        },
    }
    # ---- Iterate through keys to assess equality across names [ EQUIVALENT KEYS ]
    assert all([TRANSLATION_DICT[key]["name"] in output.keys()
                for key in ECHOPOP_OPTIMIZATION_INPUTS.keys()])
    # ---- Iterate through keys to assess equality across values [ EQUIVALENT VALUES ]
    assert all([output[TRANSLATION_DICT[key]["name"]] == TRANSLATION_DICT[key]["value"] 
                for key in ECHOPOP_OPTIMIZATION_INPUTS.keys()])
