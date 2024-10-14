import re

import numpy as np
import pytest

from echopop.tests.conftest import assert_dictionary_values_equal
from echopop.utils.validate import posfloat, posint, realcircle, realposfloat

from ..utils.validate_dict import (
    InitialValues,
    VariogramBase,
    VariogramEmpirical,
    VariogramInitial,
    VariogramOptimize,
)


def test_posint():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert posint(int(0)) == 0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(float(0.0))

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert posint(int(5)) == 5

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(float(5.1))

    # -------------------------
    # [ TEST 5 ]: ASSERT: str("test") FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(str("test"))

    # -------------------------
    # [ TEST 6 ]: ASSERT: numpy.inf FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(np.inf)

    # -------------------------
    # [ TEST 7 ]: ASSERT: numpy.nan FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(np.nan)

    # -------------------------
    # [ TEST 8 ]: ASSERT: int(-1) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(int(-1))

    # -------------------------
    # [ TEST 9 ]: ASSERT: float(-1.0) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(float(-1.0))


def test_posfloat():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert posfloat(int(0)) == 0.0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0)
    assert posfloat(float(0.0)) == 0.0

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert posfloat(int(5)) == 5.0

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1))
    assert posfloat(float(5.1)) == 5.1

    # -------------------------
    # [ TEST 5 ]: ASSERT: str("test") FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(str("test"))

    # -------------------------
    # [ TEST 6 ]: ASSERT: np.inf
    assert posfloat(np.inf) == np.inf

    # -------------------------
    # [ TEST 7 ]: ASSERT: np.inf
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-np.inf)

    # -------------------------
    # [ TEST 8 ]: ASSERT: np.nan
    assert np.isnan(posfloat(np.nan))

    # -------------------------
    # [ TEST 9 ]: ASSERT: int(-2) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-2)

    # -------------------------
    # [ TEST 10 ]: ASSERT: int(-2.5) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-2.5)


def test_realposfloat():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert realposfloat(int(0)) == 0.0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0)
    assert realposfloat(float(0.0)) == 0.0

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert realposfloat(int(5)) == 5.0

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1))
    assert realposfloat(float(5.1)) == 5.1

    # -------------------------
    # [ TEST 5 ]: ASSERT: str("test") FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative real number."):
        assert realposfloat(str("test"))

    # -------------------------
    # [ TEST 6 ]: ASSERT: np.inf FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative real number."):
        assert realposfloat(np.inf)

    # -------------------------
    # [ TEST 7 ]: ASSERT: -np.inf FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative real number."):
        assert realposfloat(-np.inf)

    # -------------------------
    # [ TEST 8 ]: ASSERT: np.nan
    assert np.isnan(realposfloat(np.nan))

    # -------------------------
    # [ TEST 9 ]: ASSERT: int(-2) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert realposfloat(-2)

    # -------------------------
    # [ TEST 10 ]: ASSERT: float(-2.5) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-2.5)


def test_realcircle():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert realcircle(int(0)) == 0.0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0)
    assert realcircle(float(0.0)) == 0.0

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert realcircle(int(5)) == 5.0

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1))
    assert realcircle(float(5.1)) == 5.1

    # -------------------------
    # [ TEST 5 ]: ASSERT: int(360)
    assert realcircle(int(360)) == 360

    # -------------------------
    # [ TEST 6 ]: ASSERT: float(360.0))
    assert realcircle(float(360.0)) == 360.0

    # -------------------------
    # [ TEST 7 ]: ASSERT: int(361) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(int(361))

    # -------------------------
    # [ TEST 8 ]: ASSERT: float(361.1)) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(float(361.1))

    # -------------------------
    # [ TEST 9 ]: ASSERT: str("test") FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(str("test"))

    # -------------------------
    # [ TEST 10 ]: ASSERT: np.inf FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(np.inf)

    # -------------------------
    # [ TEST 11 ]: ASSERT: -np.inf FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(-np.inf)

    # -------------------------
    # [ TEST 12 ]: ASSERT: np.nan
    assert np.isnan(realcircle(np.nan))

    # -------------------------
    # [ TEST 13 ]: ASSERT: int(-2) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(int(-2))

    # -------------------------
    # [ TEST 14 ]: ASSERT: float(-2.5) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(float(-2.5))


@pytest.mark.parametrize(
    "description",
    ["Test `VariogramEmpirical` class structure"],
    ids=["Test `VariogramEmpirical` class structure"],
)
def test_VariogramEmpirical_class(description):

    # -------------------------
    # Attributes [ TEST 1 ]
    # ---- Get the primary attributes
    attrs = dir(VariogramEmpirical)

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- ATTRIBUTES
    assert set(["DEFAULT_VALUES", "EXPECTED_DTYPES"]).issubset(attrs)
    # ---- METHODS
    assert set(["create", "validate"]).issubset(attrs)

    # -------------------------
    # Annotations [ TEST 2 ]
    # ---- Get the key arguments
    key_args = VariogramInitial.__annotations__

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    assert set(["fit_parameters"]).issubset(key_args)


@pytest.mark.parametrize(
    "description",
    ["Test `VariogramBase` class structure"],
    ids=["Test `VariogramBase` class structure"],
)
def test_VariogramBase_class(description):

    # -------------------------
    # Attributes [ TEST 1 ]
    # ---- Get the primary attributes
    attrs = dir(VariogramBase)

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- ATTRIBUTES
    assert set(["DEFAULT_VALUES", "EXPECTED_DTYPES", "_ROOT_DEFAULT"]).issubset(attrs)
    # ---- METHODS
    assert set(["create", "validate", "update_defaults", "restore_defaults"]).issubset(attrs)

    # -------------------------
    # Annotations [ TEST 2 ]
    # ---- Get the key arguments
    key_args = VariogramBase.__annotations__

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    assert set(
        [
            "model",
            "n_lags",
            "lag_resolution",
            "max_range",
            "sill",
            "nugget",
            "hole_effect_range",
            "correlation_range",
            "enhance_semivariance",
            "decay_power",
        ]
    ).issubset(key_args)

    # -------------------------
    # Mock new defaults [ TEST 3 ]: UPDATE VALID DEFAULTS + NO INPUTS
    # ---- Mock defaults
    MOCK_NEW_DEFAULTS = {"model": "exponential"}
    # ---- Get original defaults
    ORIGINAL_DEFAULTS = VariogramBase.DEFAULT_VALUES.copy()

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- Update the values
    VariogramBase.update_defaults(MOCK_NEW_DEFAULTS)
    # ---- Assert inequality in the updated default keys [ EXPECT FAILURE ]
    with pytest.raises(AssertionError, match=re.escape("Values for key 'model' are not the same.")):
        assert_dictionary_values_equal(VariogramBase.DEFAULT_VALUES, ORIGINAL_DEFAULTS)
    # ---- Create typed dictionary
    params = VariogramBase.create()
    # ---- Assert equality of `model`
    assert params["model"] == "exponential"

    # -------------------------
    # Mock new defaults [ TEST 4 ]: RESTORE DEFAULTS
    # ---- Restore the original defaults
    VariogramBase.restore_defaults()

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    assert_dictionary_values_equal(VariogramBase.DEFAULT_VALUES, ORIGINAL_DEFAULTS)

    # -------------------------
    # Mock new defaults [ TEST 5 ]: UPDATE VALID DEFAULTS + VALID INPUTS
    # ---- Restore the original defaults
    # ---- Mock defaults
    MOCK_NEW_DEFAULTS = {"n_lags": 50, "nugget": 0.1, "sill": 1.0}
    # ---- Mock inputs
    MOCK_INPUTS = {"lag_resolution": 0.02}

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- Update the values
    VariogramBase.update_defaults(MOCK_NEW_DEFAULTS)
    # ---- Assert inequality in the updated default keys [ EXPECT FAILURE ]
    with pytest.raises(AssertionError, match=re.escape("Values for key 'n_lags' are not close.")):
        assert_dictionary_values_equal(VariogramBase.DEFAULT_VALUES, ORIGINAL_DEFAULTS)
    # ---- Create typed dictionary
    params = VariogramBase.create(**MOCK_INPUTS)
    # ---- Assert inequality of new default parameters
    assert {
        params[key] != ORIGINAL_DEFAULTS[key] and params[key] == MOCK_NEW_DEFAULTS[key]
        for key in MOCK_NEW_DEFAULTS.keys()
    }
    # ---- Assert equality of input parameters
    assert {
        params[key] != ORIGINAL_DEFAULTS[key] and params[key] == MOCK_INPUTS[key]
        for key in MOCK_INPUTS.keys()
    }

    # -------------------------
    # Mock new defaults [ TEST 6 ]: UPDATE DEFAULTS WITH INVALID VALUE
    # ---- Mock defaults
    MOCK_NEW_DEFAULTS = {"n_lags": -50}

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- Update the values
    with pytest.raises(
        TypeError,
        match=re.escape("Value for 'n_lags' (-50, type: 'int') must be a non-negative integer."),
    ):
        VariogramBase.update_defaults(MOCK_NEW_DEFAULTS)

    # -------------------------
    # Mock new defaults [ TEST 7 ]: UPDATE DEFAULTS WITH INVALID KEY
    # ---- Mock defaults [ EXPECTED: Key is omitted]
    MOCK_NEW_DEFAULTS = {"invalid_parameter": "should not get added"}

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- Update the values
    VariogramBase.update_defaults(MOCK_NEW_DEFAULTS)
    # ---- Search for 'invalid_parameter' in defaults attribute
    assert "invalid_parameter" not in VariogramBase.DEFAULT_VALUES

    # -------------------------
    VariogramBase.restore_defaults()


@pytest.mark.parametrize(
    "description",
    ["Test `VariogramOptimize` class structure"],
    ids=["Test `VariogramOptimize` class structure"],
)
def test_VariogramOptimize_class(description):

    # -------------------------
    # Attributes [ TEST 1 ]
    # ---- Get the primary attributes
    attrs = dir(VariogramOptimize)

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- ATTRIBUTES
    assert set(["DEFAULT_VALUES", "EXPECTED_DTYPES"]).issubset(attrs)
    # ---- METHODS
    assert set(["create", "validate"]).issubset(attrs)


@pytest.mark.parametrize(
    "description",
    ["Test `VariogramInitial` class structure"],
    ids=["Test `VariogramInitial` class structure"],
)
def test_VariogramInitial_class(description):

    # -------------------------
    # Attributes [ TEST 1 ]
    # ---- Get the primary attributes
    attrs = dir(VariogramInitial)

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- ATTRIBUTES
    assert set(["VALID_PARAMETERS"]).issubset(attrs)
    # ---- METHODS
    assert set(["create", "validate"]).issubset(attrs)


@pytest.mark.parametrize(
    "description",
    ["Test `InitialValues` class structure"],
    ids=["Test `InitialValues` class structure"],
)
def test_InitialValues_class(description):

    # -------------------------
    # Attributes [ TEST 1 ]
    # ---- Get the primary attributes
    attrs = dir(InitialValues)

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    # ---- ATTRIBUTES
    assert set(["DEFAULT_STRUCTURE", "EXPECTED_DTYPES"]).issubset(attrs)

    # -------------------------
    # Attributes [ TEST 2 ]
    # ---- Validate the structure
    defaults = InitialValues.DEFAULT_STRUCTURE

    # -------------------------
    # Evaluate [ DICTIONARY ] AND Assert
    assert_dictionary_values_equal(defaults, {"min": 0.0, "value": 0.0, "max": np.inf})


@pytest.mark.parametrize(
    "params, expected, exception",
    [
        (
            VariogramBase.DEFAULT_VALUES,
            {
                "model": ["bessel", "exponential"],
                "n_lags": 30,
                "lag_resolution": None,
                "max_range": None,
                "sill": None,
                "nugget": None,
                "hole_effect_range": None,
                "correlation_range": None,
                "enhance_semivariance": None,
                "decay_power": None,
            },
            None,
        ),
        ({}, VariogramBase.DEFAULT_VALUES, None),
        (
            {
                "model": ["bessel", "exponential"],
                "n_lags": 30,
                "lag_resolution": None,
                "max_range": None,
                "sill": None,
                "nugget": None,
                "hole_effect_range": None,
                "correlation_range": None,
                "enhance_semivariance": None,
                "decay_power": None,
            },
            VariogramBase.DEFAULT_VALUES,
            None,
        ),
        ({"sill": 1.00}, {**VariogramBase.DEFAULT_VALUES, **{"sill": 1.00}}, None),
        (
            {"sill": 1.00, "nugget": 0.1, "correlation_range": 0.02},
            {
                **VariogramBase.DEFAULT_VALUES,
                **{"sill": 1.00, "nugget": 0.1, "correlation_range": 0.02},
            },
            None,
        ),
        ({"model": 1.0}, None, "does not match expected type"),
        ({"decay_power": -1.2}, None, "must be a non-negative real number"),
        ({"model": ["bessel", "exponential"]}, VariogramBase.DEFAULT_VALUES, None),
    ],
    ids=[
        "Defaults",
        "Empty input",
        "Full input (match defaults)",
        "Switch 1 input",
        "Switch 3 inputs",
        "Invalid `model`",
        "Invalid `decay_power`",
        "Input that is same as default",
    ],
)
def test_VariogramBase_validation(params, expected, exception):

    if exception is not None:
        with pytest.raises(TypeError, match=exception):
            assert VariogramBase.create(**params)
    else:
        # Test creation with various parameters
        result = VariogramBase.create(**params)
        assert_dictionary_values_equal(result, expected)


@pytest.mark.parametrize(
    "params, expected, exception",
    [
        (
            VariogramEmpirical.DEFAULT_VALUES,
            {"force_lag_zero": True, "standardize_coordinates": True, "azimuth_range": 360.0},
            None,
        ),
        ({}, VariogramEmpirical.DEFAULT_VALUES, None),
        (
            {"force_lag_zero": True, "standardize_coordinates": True, "azimuth_range": 360.0},
            VariogramEmpirical.DEFAULT_VALUES,
            None,
        ),
        (
            {"force_lag_zero": False},
            {**VariogramEmpirical.DEFAULT_VALUES, **{"force_lag_zero": False}},
            None,
        ),
        (
            {"force_lag_zero": False, "azimuth_range": 180.0},
            {
                **VariogramEmpirical.DEFAULT_VALUES,
                **{"force_lag_zero": False, "azimuth_range": 180.0},
            },
            None,
        ),
        ({"force_lag_zero": 1.0}, None, "does not match expected type"),
        ({"azimuth_range": -300}, None, "must be a non-negative real angle"),
        ({"standardize_coordinates": True}, VariogramEmpirical.DEFAULT_VALUES, None),
    ],
    ids=[
        "Defaults",
        "Empty input",
        "Full input (match defaults)",
        "Switch 1 input",
        "Switch 2 inputs",
        "Invalid `force_lag_zero`",
        "Invalid `azimuth_range`",
        "Input that is same as default",
    ],
)
def test_VariogramEmpirical_validation(params, expected, exception):

    if exception is not None:
        with pytest.raises(TypeError, match=exception):
            assert VariogramEmpirical.create(**params)
    else:
        # Test creation with various parameters
        result = VariogramEmpirical.create(**params)
        assert_dictionary_values_equal(result, expected)


@pytest.mark.parametrize(
    "params, expected, exception",
    [
        (
            VariogramOptimize.DEFAULT_VALUES,
            {
                "max_fun_evaluations": 500,
                "cost_fun_tolerance": 1e-6,
                "solution_tolerance": 1e-6,
                "gradient_tolerance": 1e-4,
                "finite_step_size": 1e-8,
                "trust_region_solver": "exact",
                "x_scale": "jacobian",
                "jacobian_approx": "central",
            },
            None,
        ),
        ({}, VariogramOptimize.DEFAULT_VALUES, None),
        (
            {
                "max_fun_evaluations": 500,
                "cost_fun_tolerance": 1e-6,
                "solution_tolerance": 1e-6,
                "gradient_tolerance": 1e-4,
                "finite_step_size": 1e-8,
                "trust_region_solver": "exact",
                "x_scale": "jacobian",
                "jacobian_approx": "central",
            },
            VariogramOptimize.DEFAULT_VALUES,
            None,
        ),
        (
            {"max_fun_evaluations": 400},
            {**VariogramOptimize.DEFAULT_VALUES, **{"max_fun_evaluations": 400}},
            None,
        ),
        (
            {"gradient_tolerance": 1e-6, "trust_region_solver": "base", "x_scale": np.array([1.0])},
            {
                **VariogramOptimize.DEFAULT_VALUES,
                **{
                    "gradient_tolerance": 1e-6,
                    "trust_region_solver": "base",
                    "x_scale": np.array([1.0]),
                },
            },
            None,
        ),
        ({"trust_region_solver": 1.0}, None, "does not match expected type"),
        ({"finite_step_size": -300}, None, "must be a non-negative real number"),
        ({"jacobian_approx": "central"}, VariogramOptimize.DEFAULT_VALUES, None),
    ],
    ids=[
        "Defaults",
        "Empty input",
        "Default input (type: list)",
        "Switch 1 input",
        "Switch 3 inputs",
        "Invalid `trust_region_solver`",
        "Invalid `finite_step_size`",
        "Input that is same as default",
    ],
)
def test_VariogramOptimize_validation(params, expected, exception):

    if exception:
        with pytest.raises(TypeError, match=exception):
            assert VariogramOptimize.create(**params)
    else:
        # Test creation with various parameters
        result = VariogramOptimize.create(**params)
        assert_dictionary_values_equal(result, expected)


@pytest.mark.parametrize(
    "params, expected, exception",
    [
        ({}, {}, None),
        (
            VariogramInitial.VALID_PARAMETERS,
            {
                "correlation_range": {"min": 0.0, "value": 0.0, "max": np.inf},
                "decay_power": {"min": 0.0, "value": 0.0, "max": np.inf},
                "hole_effect_range": {"min": 0.0, "value": 0.0, "max": np.inf},
                "nugget": {"min": 0.0, "value": 0.0, "max": np.inf},
                "sill": {"min": 0.0, "value": 0.0, "max": np.inf},
            },
            None,
        ),
        (
            ["correlation_range", "sill"],
            {
                "correlation_range": {"min": 0.0, "value": 0.0, "max": np.inf},
                "sill": {"min": 0.0, "value": 0.0, "max": np.inf},
            },
            None,
        ),
        ({"sill": {}}, {"sill": {"min": 0.0, "value": 0.0, "max": np.inf}}, None),
        ({"sill": {"value": 0.5}}, {"sill": {"min": 0.0, "value": 0.5, "max": np.inf}}, None),
        (
            {"sill": {"min": 0.1, "value": 0.5}},
            {"sill": {"min": 0.1, "value": 0.5, "max": np.inf}},
            None,
        ),
        (
            {"sill": {"min": 0.1, "value": 0.5, "max": 1.0}},
            {"sill": {"min": 0.1, "value": 0.5, "max": 1.0}},
            None,
        ),
        (
            {
                "sill": {"min": 0.1, "value": 0.5, "max": 1.0},
                "nugget": {},
                "decay_power": {"value": 1.5},
                "correlation_range": {"min": 0.1, "value": 0.2, "max": np.inf},
            },
            {
                "sill": {"min": 0.1, "value": 0.5, "max": 1.0},
                "nugget": {"min": 0.0, "value": 0.0, "max": np.inf},
                "decay_power": {"min": 0.0, "value": 1.5, "max": np.inf},
                "correlation_range": {"min": 0.1, "value": 0.2, "max": np.inf},
            },
            None,
        ),
        (["invalid_parameter"], None, "Unexpected parameter(s): 'invalid_parameter'"),
        ({"invalid_parameter": {}}, None, "Unexpected parameter(s): 'invalid_parameter'"),
        (["enhance_semivariance"], None, "Unexpected parameter(s): 'enhance_semivariance'"),
        (
            {"enhance_semivariance": True},
            None,
            "Value for 'enhance_semivariance' must be a dictionary.",
        ),
        (
            {"enhance_semivariance": {"value": True}},
            None,
            "Unexpected parameter(s): 'enhance_semivariance'",
        ),
        ({"sill": {"min": 0.2}}, None, "Invalid initial values for: 'sill'."),
        ({"nugget": {"min": 0.2, "value": 0.0}}, None, "Invalid initial values for: 'nugget'."),
        (
            {"decay_power": {"value": 0.2, "max": 0.1}},
            None,
            "Invalid initial values for: 'decay_power'.",
        ),
        ({"decay_power": {"value": -0.1, "max": 0.1}}, None, "must be a non-negative real number"),
        ({"decay_power": {"value": 0.0, "max": -0.1}}, None, "must be a non-negative float"),
        ({"decay_power": {"min": -0.2, "value": 0.0}}, None, "must be a non-negative float"),
    ],
    ids=[
        "Empty input",
        "Default output for valid parameters",
        "User input (type: list)",
        "User input (type: dictionary: empty)",
        "User input (type: dictionary: value)",
        "User input (type: dictionary: min+value)",
        "User input (type: dictionary: min+value+max)",
        "User input (type: dictionary: mixed parameters)",
        "Invalid input (type: list)",
        "Invalid input (type: dictionary)",
        "Unexpected input (type: list)",
        "Unexpected input and wrong type (type: dictionary)",
        "Unexpected input (type: dictionary)",
        "Invalid case: 'min' > 'value' (no 'value' input)",
        "Invalid case: 'min' > 'value' (provided 'value' input)",
        "Invalid case: 'value' > 'max' (provided 'value' input)",
        "Invalid 'value' datatype",
        "Invalid 'max' datatype",
        "Invalid 'min' datatype",
    ],
)
def test_VariogramInitial_validation(params, expected, exception):

    if exception is not None:
        with pytest.raises((TypeError, ValueError), match=re.escape(exception)):
            assert VariogramInitial.create(params)
    else:
        # Test creation with various parameters
        result = VariogramInitial.create(params)
        assert_dictionary_values_equal(result, expected)
