import re
from collections import OrderedDict
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
import pytest
from pydantic import RootModel, ValidationError

from echopop.utils.validate_dict import (
    InitialValues,
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
    VariogramInitial,
    VariogramModel,
    VariogramOptimize,
)

from ..utils.validate import posfloat, posint, realcircle, realposfloat


@pytest.fixture
def InitialValues_fields() -> Dict[str, Any]:

    return {
        "min": {
            "default": None,
            "annotation": Optional[realposfloat],
        },
        "value": {
            "default": 0.0,
            "annotation": realposfloat,
            "allow_inf_nan": False,
        },
        "max": {
            "default": None,
            "annotation": Optional[posfloat],
        },
        "vary": {"default": False, "annotation": bool},
    }


@pytest.mark.parametrize(
    "description",
    ["Assess `InitialValues` pydantic model structure"],
    ids=["Assess `InitialValues` pydantic model structure"],
)
def test_InitialValues_model_structure(description, InitialValues_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(InitialValues_fields).issubset(InitialValues.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            InitialValues.model_fields[param]._attributes_set == InitialValues_fields[param]
            for param in InitialValues.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert InitialValues.__base__ == VariogramModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(InitialValues)
    # ---- Validate correct setting
    assert InitialValues.model_config == dict(arbitrary_types_allowed=True)

    # -------------------------
    # ASSERT: 'create' method
    # ---- Check existence
    assert "create" in dir(InitialValues)
    # ---- Check default output
    assert InitialValues.create() == {"value": 0.0, "vary": False}

    # -------------------------
    # ASSERT: private methods
    # ---- Internal class methods for filling default values
    assert set(["_DEFAULT_STRUCTURE_EMPTY", "_DEFAULT_STRUCTURE_OPTIMIZE"]).issubset(
        InitialValues.__dict__
    )
    # ---- Output for '._DEFAULT_STRUCTURE_EMPTY'
    assert InitialValues._DEFAULT_STRUCTURE_EMPTY() == {"vary": False}
    # ---- Output for '._DEFAULT_STRUCTURE_OPTIMIZE'
    assert InitialValues._DEFAULT_STRUCTURE_OPTIMIZE() == {
        "min": 0.0,
        "value": 0.0,
        "max": np.inf,
        "vary": True,
    }

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = InitialValues.__pydantic_decorators__.field_validators
    # ---- 'model_validators' 'Decorator' object
    model_decor = InitialValues.__pydantic_decorators__.model_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert set(["validate_realposfloat", "validate_posfloat"]).issubset(InitialValues.__dict__)

    # -------------------------
    # ASSERT: 'valid_realposfloat'
    # ---- Check output for 'validate_realposfloat' [VALID]
    assert InitialValues.validate_realposfloat(2) == 2.0
    # ---- Check output for 'validate_realposfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert InitialValues.validate_realposfloat(-2.0)
    # ---- Check output for 'validate_realposfloat' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real number.")):
        assert InitialValues.validate_realposfloat(np.inf)
    # ---- Check the applicable fields
    assert field_decor["validate_realposfloat"].info.fields == ("min", "value")
    # ---- Check the applicable validation mode
    assert field_decor["validate_realposfloat"].info.mode == "before"

    # -------------------------
    # ASSERT: 'valid_posfloat'
    # ---- Checkout output for 'validate_posfloat' [VALID]
    assert InitialValues.validate_posfloat(1) == 1.0
    # ---- Check output for 'validate_posfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert InitialValues.validate_posfloat(-1.0)
    # ---- Check the applicable fields
    assert field_decor["validate_posfloat"].info.fields == ("max",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_posfloat"].info.mode == "before"

    # -------------------------
    # ASSERT: model validator decorators
    # ---- Check whether the field exists
    assert "validate_value_sort" in InitialValues.__dict__
    # MOCK: parameterized model
    MOCK_SELF = InitialValues(**{"max": 2.0, "min": 1.0, "value": 1.5, "vary": True})
    # ASSERT: mode validator correctly sorts the outputs
    # ---- [NOTE]: This is not required as an 'OrderedDict' -- it is more for console-related
    # ---- messages/reports
    # ---- Check that key-value pairs are identical
    assert OrderedDict(InitialValues.validate_value_sort(MOCK_SELF).model_dump()) == OrderedDict(
        {"min": 1.0, "value": 1.5, "max": 2.0, "vary": True}
    )
    # ---- Check the applicable validation mode
    assert model_decor["validate_value_sort"].info.mode == "after"


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            {},
            {"value": 0.0, "vary": False},
            None,
        ),
        (
            {"vary": False},
            {"value": 0.0, "vary": False},
            None,
        ),
        (
            {"vary": True},
            {"min": 0.0, "value": 0.0, "max": np.inf, "vary": True},
            None,
        ),
        (
            {"min": 0.0, "value": 0.5, "max": 1.0, "vary": False},
            {"value": 0.5, "vary": False},
            None,
        ),
        (
            {"min": 0.0, "max": 1.0, "vary": True},
            {"min": 0.0, "value": 0.0, "max": 1.0, "vary": True},
            None,
        ),
        (
            {"max": np.inf, "vary": True},
            {"min": 0.0, "value": 0.0, "max": np.inf, "vary": True},
            None,
        ),
        (
            {"min": 1, "value": 2, "max": 3, "vary": True},
            {"min": 1.0, "value": 2.0, "max": 3.0, "vary": True},
            None,
        ),
        (
            {"min": -1.0, "value": -2.0, "max": -3.0, "vary": "dummy"},
            None,
            [
                "min\n  Value error, Value must be a non-negative float",
                "value\n  Value error, Value must be a non-negative float.",
                "max\n  Value error, Value must be a non-negative float.",
                "vary\n  Input should be a valid boolean",
            ],
        ),
        (
            {"min": 2.0, "value": 1.0, "max": 3.0, "vary": True},
            None,
            "Value error, Optimization minimum, starting, and maximum values  must satisfy the "
            "logic:",
        ),
        (
            {"min": 1.0, "value": 2.0, "max": 1.0, "vary": True},
            None,
            "Value error, Optimization minimum, starting, and maximum values  must satisfy the "
            "logic:",
        ),
        (
            {"min": 3.0, "max": 1.0, "vary": True},
            None,
            "Value error, Optimization minimum, starting, and maximum values  must satisfy the "
            "logic:",
        ),
        (
            {"min": 3.0, "value": 4.0, "vary": True, "extra": True},
            {"min": 3.0, "value": 4.0, "max": np.inf, "vary": True},
            None,
        ),
    ],
    ids=[
        "Empty dictionary [fill with default values]",
        "Default output ['vary'=False]",
        "Default output ['vary'=True]",
        "Override 'min' and 'max' inputs when ['vary'=False]",
        "Fill default 'value' when 'vary' = True",
        "'max' set to 'np.inf' with no other inputs ['vary'=True]",
        "Valid coercion of 'int' to 'float'",
        "Invalid datatyping for each input",
        "Value for 'min' exceeds 'value'",
        "Value for 'value' exceeds 'max'",
        "Value for 'min' exceeds 'max'",
        "Automatic pruning of function arguments",
    ],
)
def test_InitialValues_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert InitialValues.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert InitialValues.create(**input)
    else:
        result = InitialValues.create(**input)
        assert result == expected


@pytest.fixture
def VariogramBase_fields() -> Dict[str, Any]:

    return {
        "model": {
            "annotation": Union[str, List[str]],
            "union_mode": "left_to_right",
        },
        "n_lags": {
            "annotation": posint,
            "ge": 1,
            "allow_inf_nan": False,
        },
        "lag_resolution": {
            "default": None,
            "annotation": Optional[realposfloat],
            "gt": 0.0,
            "allow_inf_nan": False,
        },
        "sill": {
            "default": None,
            "annotation": Optional[realposfloat],
            "ge": 0.0,
            "allow_inf_nan": False,
        },
        "nugget": {
            "default": None,
            "annotation": Optional[realposfloat],
            "ge": 0.0,
            "allow_inf_nan": False,
        },
        "hole_effect_range": {
            "default": None,
            "annotation": Optional[realposfloat],
            "ge": 0.0,
            "allow_inf_nan": False,
        },
        "correlation_range": {
            "default": None,
            "annotation": Optional[realposfloat],
            "ge": 0.0,
            "allow_inf_nan": False,
        },
        "enhance_semivariance": {
            "default": None,
            "annotation": Optional[bool],
        },
        "decay_power": {
            "default": None,
            "annotation": Optional[realposfloat],
            "ge": 0.0,
            "allow_inf_nan": False,
        },
    }


@pytest.mark.parametrize(
    "description",
    ["Assess `VariogramBase` pydantic model structure"],
    ids=["Assess `VariogramBase` pydantic model structure"],
)
def test_VariogramBase_model_structure(description, VariogramBase_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(VariogramBase_fields).issubset(VariogramBase.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            VariogramBase.model_fields[param]._attributes_set == VariogramBase_fields[param]
            for param in VariogramBase.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert VariogramBase.__base__ == VariogramModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(VariogramBase)
    # ---- Validate correct setting
    assert VariogramBase.model_config == dict(arbitrary_types_allowed=True)

    # -------------------------
    # ASSERT: 'create' method
    # ---- Check existence
    assert "create" in dir(VariogramBase)
    # ---- Check default output
    with pytest.raises(ValidationError, match=re.escape("2 validation errors for VariogramBase")):
        assert VariogramBase.create()

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = VariogramBase.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert set(["validate_posint", "validate_realposfloat"]).issubset(VariogramBase.__dict__)

    # -------------------------
    # ASSERT: 'valid_realposfloat'
    # ---- Check output for 'validate_realposfloat' [VALID]
    assert VariogramBase.validate_realposfloat(2) == 2.0
    # ---- Check output for 'validate_realposfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert VariogramBase.validate_realposfloat(-2.0)
    # ---- Check output for 'validate_realposfloat' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real number.")):
        assert VariogramBase.validate_realposfloat(np.inf)
    # ---- Check the applicable fields
    assert field_decor["validate_realposfloat"].info.fields == (
        "lag_resolution",
        "sill",
        "nugget",
        "hole_effect_range",
        "correlation_range",
        "decay_power",
    )
    # ---- Check the applicable validation mode
    assert field_decor["validate_realposfloat"].info.mode == "before"

    # -------------------------
    # ASSERT: 'valid_posint'
    # ---- Checkout output for 'validate_posfloat' [VALID]
    assert VariogramBase.validate_posint(1) == 1
    # ---- Checkout output for 'validate_posfloat' [INVALID: float]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert VariogramBase.validate_posint(1.0)
    # ---- Check output for 'validate_posint' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert VariogramBase.validate_posint(-1)
    # ---- Check output for 'validate_posint' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert VariogramBase.validate_posint(np.nan)
    # ---- Check the applicable fields
    assert field_decor["validate_posint"].info.fields == ("n_lags",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_posint"].info.mode == "before"


@pytest.fixture
def VariogramEmpirical_fields() -> Dict[str, Any]:

    return {
        "azimuth_range": {
            "default": 360.0,
            "annotation": realcircle,
            "ge": 0.0,
            "le": 360.0,
            "allow_inf_nan": False,
        },
        "force_lag_zero": {
            "default": True,
            "annotation": bool,
        },
        "standardize_coordinates": {
            "default": True,
            "annotation": bool,
        },
    }


@pytest.mark.parametrize(
    "description",
    ["Assess `VariogramEmpirical` pydantic model structure"],
    ids=["Assess `VariogramEmpirical` pydantic model structure"],
)
def test_VariogramEmpirical_model_structure(description, VariogramEmpirical_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(VariogramEmpirical_fields).issubset(VariogramEmpirical.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            VariogramEmpirical.model_fields[param]._attributes_set
            == VariogramEmpirical_fields[param]
            for param in VariogramEmpirical.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert VariogramEmpirical.__base__ == VariogramModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(VariogramEmpirical)
    # ---- Validate correct setting
    assert VariogramEmpirical.model_config == dict(arbitrary_types_allowed=True)

    # -------------------------
    # ASSERT: 'create' method
    # ---- Check existence
    assert "create" in dir(VariogramEmpirical)
    # ---- Check default output
    assert VariogramEmpirical.create() == {
        "azimuth_range": 360.0,
        "force_lag_zero": True,
        "standardize_coordinates": True,
    }

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = VariogramEmpirical.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert "validate_realcircle" in VariogramEmpirical.__dict__
    # -------------------------
    # ASSERT: 'validate_realcircle'
    # ---- Check output for 'validate_realcircle' [VALID]
    assert VariogramEmpirical.validate_realcircle(2) == 2.0
    # ---- Check output for 'validate_realcircle' [INVALID: < 0.0]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real angle")):
        assert VariogramEmpirical.validate_realcircle(-2.0)
    # ---- Check output for 'validate_realcircle' [INVALID: > 360.0]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real angle")):
        assert VariogramEmpirical.validate_realcircle(361.0)
    # ---- Check output for 'validate_realposfloat' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real angle")):
        assert VariogramEmpirical.validate_realcircle(np.inf)
    # ---- Check the applicable fields
    assert field_decor["validate_realcircle"].info.fields == ("azimuth_range",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_realcircle"].info.mode == "before"


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            {"azimuth_range": 360.0, "force_lag_zero": True, "standardize_coordinates": True},
            {"azimuth_range": 360.0, "force_lag_zero": True, "standardize_coordinates": True},
            None,
        ),
        (
            {},
            {"azimuth_range": 360.0, "force_lag_zero": True, "standardize_coordinates": True},
            None,
        ),
        (
            {"azimuth_range": 180.0, "force_lag_zero": False, "standardize_coordinates": False},
            {"azimuth_range": 180.0, "force_lag_zero": False, "standardize_coordinates": False},
            None,
        ),
        (
            {"azimuth_range": 180.0},
            {"azimuth_range": 180.0, "force_lag_zero": True, "standardize_coordinates": True},
            None,
        ),
        (
            {"force_lag_zero": False},
            {"azimuth_range": 360.0, "force_lag_zero": False, "standardize_coordinates": True},
            None,
        ),
        (
            {"standardize_coordinates": False},
            {"azimuth_range": 360.0, "force_lag_zero": True, "standardize_coordinates": False},
            None,
        ),
        (
            {"azimuth_range": int(360)},
            {"azimuth_range": 360.0, "force_lag_zero": True, "standardize_coordinates": True},
            None,
        ),
        (
            {
                "azimuth_range": 361.0,
                "force_lag_zero": "Invalid",
                "standardize_coordinates": "Invalid",
            },
            None,
            [
                "Value must be a non-negative real angle",
                "Input should be a valid boolean",
                "Input should be a valid boolean",
            ],
        ),
        (
            {"extra": True},
            {"azimuth_range": 360.0, "force_lag_zero": True, "standardize_coordinates": True},
            None,
        ),
    ],
    ids=[
        "Default values",
        "Empty dictionary [fill with default values]",
        "Change parameter values for each",
        "Define only 'azimuth_range'",
        "Define only 'force_lag_zero'",
        "Define only 'standardize_coordinates'",
        "Acceptable coercion of 'azimuth_range' ['int' to 'float']",
        "Error handling when input datatypes are all incorrect",
        "Automatic pruning of erroneous arguments",
    ],
)
def test_VariogramEmpirical_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert VariogramEmpirical.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert VariogramEmpirical.create(**input)
    else:
        result = VariogramEmpirical.create(**input)
        assert result == expected


@pytest.fixture
def VariogramInitial_fields() -> Dict[str, Any]:

    return {
        "root": {
            "annotation": Dict[str, InitialValues],
            "frozen": None,
        },
    }


@pytest.mark.parametrize(
    "description",
    ["Assess `VariogramInitial` pydantic model structure"],
    ids=["Assess `VariogramInitial` pydantic model structure"],
)
def test_VariogramInitial_model_structure(description, VariogramInitial_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(VariogramInitial_fields).issubset(VariogramInitial.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            VariogramInitial.model_fields[param]._attributes_set == VariogramInitial_fields[param]
            for param in VariogramInitial.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert VariogramInitial.__base__ == RootModel[InitialValues]

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(VariogramInitial)
    # ---- Validate correct setting
    assert VariogramInitial.model_config == {}

    # -------------------------
    # ASSERT: private methods
    # ---- Internal class methods for filling default values
    assert "_VALID_PARAMETERS" in VariogramInitial.__dict__
    # ---- Output for '._VALID_PARAMETERS'
    assert VariogramInitial._VALID_PARAMETERS() == [
        "correlation_range",
        "decay_power",
        "hole_effect_range",
        "nugget",
        "sill",
    ]
    # ---- Check existence for 'create'
    assert "create" in dir(VariogramInitial)

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'model_validators' 'Decorator' object
    model_decor = VariogramInitial.__pydantic_decorators__.model_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert "validate_model_params" in VariogramInitial.__dict__
    # MOCK: parameterized values input from the pydantic model
    MOCK_V = {
        "correlation_range": {},
        "decay_power": {},
        "hole_effect_range": {},
        "nugget": {},
        "sill": {},
    }
    # ASSERT: mode validator correctly sorts the outputs
    assert VariogramInitial.validate_model_params(MOCK_V) == MOCK_V
    # ---- Check the applicable validation mode
    assert model_decor["validate_model_params"].info.mode == "before"


# input = {"dummy1": {}, "dummy2": {}}
# exception = "Unexpected optimization parameters: ['dummy1', 'dummy2']."
# from echopop.utils.validate_dict import VariogramInitial
# VariogramInitial.create(**input)


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            {"root": {}},
            {},
            None,
        ),
        (
            {"correlation_range": {"vary": False}},
            {"correlation_range": {"value": 0.0, "vary": False}},
            None,
        ),
        (
            {"correlation_range": {"vary": True}},
            {"correlation_range": {"min": 0.0, "value": 0.0, "max": np.inf, "vary": True}},
            None,
        ),
        (
            {"nugget": {"vary": False}, "sill": {"vary": False}},
            {"nugget": {"value": 0.0, "vary": False}, "sill": {"value": 0.0, "vary": False}},
            None,
        ),
        (
            {"nugget": {"vary": True}, "sill": {"vary": True}},
            {
                "nugget": {"min": 0.0, "value": 0.0, "max": np.inf, "vary": True},
                "sill": {"min": 0.0, "value": 0.0, "max": np.inf, "vary": True},
            },
            None,
        ),
        (
            {"nugget": {"vary": False}, "sill": {"vary": True}},
            {
                "nugget": {"value": 0.0, "vary": False},
                "sill": {"min": 0.0, "value": 0.0, "max": np.inf, "vary": True},
            },
            None,
        ),
        (
            {
                "correlation_range": {},
                "decay_power": {},
                "hole_effect_range": {},
                "nugget": {},
                "sill": {},
            },
            {
                "correlation_range": {"value": 0.0, "vary": False},
                "decay_power": {"value": 0.0, "vary": False},
                "hole_effect_range": {"value": 0.0, "vary": False},
                "nugget": {"value": 0.0, "vary": False},
                "sill": {"value": 0.0, "vary": False},
            },
            None,
        ),
        (
            {"dummy": {}},
            None,
            "Value error, Unexpected optimization parameters: ['dummy']",
        ),
        (
            {"dummy1": {}, "dummy2": {}},
            None,
            "Unexpected optimization parameters",
        ),
        (
            {"dummy1": {}, "dummy2": {}, "nugget": {}},
            None,
            "Unexpected optimization parameters",
        ),
    ],
    ids=[
        "Empty input [root]",
        "Single input ['vary'=False]",
        "Single input['vary'=True]",
        "Multiple inputs ['vary'=False]",
        "Multiple inputs ['vary'=True]",
        "Multiple inputs ['vary'=True/False]",
        "All valid parameters [Empty entries]",
        "Single invalid input",
        "Multiple invalid inputs",
        "Multiple invalid inputs with valid parameters",
    ],
)
def test_VariogramInitial_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert VariogramInitial.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert VariogramInitial.create(**input)
    else:
        result = VariogramInitial.create(**input)
        assert result == expected


@pytest.mark.parametrize(
    "description",
    ["Assess `VariogramModel` pydantic model structure"],
    ids=["Assess `VariogramModel` pydantic model structure"],
)
def test_VariogramModel_model_structure(description):

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(VariogramModel)
    # ---- Validate correct setting
    assert VariogramModel.model_config == dict(arbitrary_types_allowed=True)

    # -------------------------
    # ASSERT: 'create' method
    # ---- Check existence
    assert "create" in dir(VariogramModel)
    # ---- Validate correct output from empty DataFrame
    assert VariogramModel.create(**dict()) == {}
    # ---- Validate that any arbitrary output yields an empty DataFrame
    assert VariogramModel.create(**dict(dum1="dum2")) == {}


@pytest.fixture
def VariogramOptimize_fields() -> Dict[str, Any]:

    return {
        "max_fun_evaluations": {
            "default": 500,
            "annotation": posint,
            "gt": 0,
            "allow_inf_nan": False,
        },
        "cost_fun_tolerance": {
            "default": 1e-6,
            "annotation": realposfloat,
            "gt": 0.0,
            "allow_inf_nan": False,
        },
        "gradient_tolerance": {
            "default": 1e-4,
            "annotation": realposfloat,
            "gt": 0.0,
            "allow_inf_nan": False,
        },
        "solution_tolerance": {
            "default": 1e-6,
            "annotation": realposfloat,
            "gt": 0.0,
            "allow_inf_nan": False,
        },
        "finite_step_size": {
            "default": 1e-8,
            "annotation": realposfloat,
            "gt": 0.0,
            "allow_inf_nan": False,
        },
        "trust_region_solver": {
            "default": "exact",
            "annotation": Literal["base", "exact"],
        },
        "x_scale": {
            "default": "jacobian",
            "annotation": Union[Literal["jacobian"], np.ndarray[realposfloat]],
        },
        "jacobian_approx": {
            "default": "central",
            "annotation": Literal["forward", "central"],
        },
    }


@pytest.mark.parametrize(
    "description",
    ["Assess `VariogramOptimize` pydantic model structure"],
    ids=["Assess `VariogramOptimize` pydantic model structure"],
)
def test_VariogramOptimize_model_structure(description, VariogramOptimize_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(VariogramOptimize_fields).issubset(VariogramOptimize.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            VariogramOptimize.model_fields[param]._attributes_set == VariogramOptimize_fields[param]
            for param in VariogramOptimize.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert VariogramOptimize.__base__ == VariogramModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(VariogramOptimize)
    # ---- Validate correct setting
    assert VariogramOptimize.model_config == dict(arbitrary_types_allowed=True)

    # -------------------------
    # ASSERT: 'create' method
    # ---- Check existence
    assert "create" in dir(VariogramOptimize)
    # ---- Check default output
    assert VariogramOptimize.create() == {
        "max_fun_evaluations": 500,
        "cost_fun_tolerance": 1e-06,
        "gradient_tolerance": 1e-04,
        "solution_tolerance": 1e-06,
        "finite_step_size": 1e-08,
        "trust_region_solver": "exact",
        "x_scale": "jacobian",
        "jacobian_approx": "central",
    }

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = VariogramOptimize.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert set(["validate_posint", "validate_realposfloat", "validate_xscale"]).issubset(
        VariogramOptimize.__dict__
    )

    # -------------------------
    # ASSERT: 'valid_realposfloat'
    # ---- Check output for 'validate_realposfloat' [VALID]
    assert VariogramOptimize.validate_realposfloat(2) == 2.0
    # ---- Check output for 'validate_realposfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert VariogramOptimize.validate_realposfloat(-2.0)
    # ---- Check output for 'validate_realposfloat' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real number.")):
        assert VariogramOptimize.validate_realposfloat(np.inf)
    # ---- Check the applicable fields
    assert field_decor["validate_realposfloat"].info.fields == (
        "cost_fun_tolerance",
        "gradient_tolerance",
        "finite_step_size",
        "solution_tolerance",
    )
    # ---- Check the applicable validation mode
    assert field_decor["validate_realposfloat"].info.mode == "before"

    # -------------------------
    # ASSERT: 'valid_posint'
    # ---- Checkout output for 'validate_posfloat' [VALID]
    assert VariogramOptimize.validate_posint(1) == 1
    # ---- Checkout output for 'validate_posfloat' [INVALID: float]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert VariogramOptimize.validate_posint(1.0)
    # ---- Check output for 'validate_posint' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert VariogramOptimize.validate_posint(-1)
    # ---- Check output for 'validate_posint' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert VariogramOptimize.validate_posint(np.nan)
    # ---- Check the applicable fields
    assert field_decor["validate_posint"].info.fields == ("max_fun_evaluations",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_posint"].info.mode == "before"

    # -------------------------
    # ASSERT: 'valid_xscale'
    # ---- Checkout output for 'validate_xscale' [VALID: Literal['jacobian']]
    assert VariogramOptimize.validate_xscale("jacobian")
    # ---- Checkout output for 'validate_xscale' [INVALID: Literal['other']]
    with pytest.raises(
        ValueError,
        match=re.escape("Input should be either the Literal " "'jacobian' or a NumPy array"),
    ):
        assert VariogramOptimize.validate_xscale("other")
    # ---- Checkout output for 'validate_xscale' [VALID: np.ndarray[1.0]]
    assert VariogramOptimize.validate_xscale(np.array([1.0])) == np.array([1.0])
    # ---- Checkout output for 'validate_xscale' [VALID: np.ndarray[1]]
    assert VariogramOptimize.validate_xscale(np.array([1])) == np.array([1.0])
    # ---- Checkout output for 'validate_xscale' [VALID: np.ndarray[1.0, 2.0]]
    assert all(VariogramOptimize.validate_xscale(np.array([1.0, 2.0])) == np.array([1.0, 2.0]))
    # ---- Checkout output for 'validate_xscale' [VALID: np.ndarray[1, 2]]
    assert all(VariogramOptimize.validate_xscale(np.array([1, 2])) == np.array([1.0, 2.0]))
    # ---- Checkout output for 'validate_xscale' [INVALID: np.ndarray[-1.0]]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert VariogramOptimize.validate_xscale(np.array([-1.0]))
    # ---- Checkout output for 'validate_xscale' [INVALID: 1.0]
    with pytest.raises(
        ValueError,
        match=re.escape("Input should be either the Literal " "'jacobian' or a NumPy array"),
    ):
        assert VariogramOptimize.validate_xscale(1.0)
    # ---- Check output for 'validate_posint' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real number.")):
        assert VariogramOptimize.validate_xscale(np.array([np.inf]))
    # ---- Check the applicable fields
    assert field_decor["validate_xscale"].info.fields == ("x_scale",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_xscale"].info.mode == "before"


def VariogramOptimize_DEFAULT():
    return VariogramOptimize.create(**{})


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            VariogramOptimize_DEFAULT(),
            VariogramOptimize_DEFAULT(),
            None,
        ),
        (
            {},
            VariogramOptimize_DEFAULT(),
            None,
        ),
        (
            {"max_fun_evaluations": 100},
            {**VariogramOptimize_DEFAULT(), **{"max_fun_evaluations": 100}},
            None,
        ),
        (
            {"max_fun_evaluations": 100, "solution_tolerance": 1e-1},
            {
                **VariogramOptimize_DEFAULT(),
                **{"max_fun_evaluations": 100, "solution_tolerance": 1e-1},
            },
            None,
        ),
        (
            {
                "max_fun_evaluations": 100,
                "cost_fun_tolerance": 1e-1,
                "solution_tolerance": 1e-1,
                "gradient_tolerance": 1e-1,
                "finite_step_size": 1e-1,
                "trust_region_solver": "base",
                "x_scale": np.array([1.0]),
                "jacobian_approx": "forward",
            },
            {
                "max_fun_evaluations": 100,
                "cost_fun_tolerance": 1e-1,
                "solution_tolerance": 1e-1,
                "gradient_tolerance": 1e-1,
                "finite_step_size": 1e-1,
                "trust_region_solver": "base",
                "x_scale": np.array([1.0]),
                "jacobian_approx": "forward",
            },
            None,
        ),
        (
            {
                "max_fun_evaluations": -1,
                "cost_fun_tolerance": -1,
                "solution_tolerance": -1,
                "gradient_tolerance": -1,
                "finite_step_size": -1,
                "trust_region_solver": "invalid",
                "x_scale": "invalid",
                "jacobian_approx": "invalid",
            },
            None,
            [
                "max_fun_evaluations\n  Value error, Value must be a non-negative integer",
                "cost_fun_tolerance\n  Value error, Value must be a non-negative float",
                "solution_tolerance\n  Value error, Value must be a non-negative float",
                "gradient_tolerance\n  Value error, Value must be a non-negative float",
                "finite_step_size\n  Value error, Value must be a non-negative float",
                "trust_region_solver\n  Input should be 'base' or 'exact'",
                "x_scale\n  Value error, Input should be either the Literal "
                "'jacobian' or a NumPy array",
                "jacobian_approx\n  Input should be 'forward' or 'central'",
            ],
        ),
    ],
    ids=[
        "Default input",
        "Empty input [return defaults]",
        "Change single input",
        "Change multiple inputs",
        "Change all inputs",
        "All invalid inputs",
    ],
)
def test_VariogramOptimize_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert VariogramOptimize.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert VariogramOptimize.create(**input)
    else:
        result = VariogramOptimize.create(**input)
        assert result == expected


@pytest.mark.parametrize(
    "description",
    ["Assess `KrigingParameterInputs` pydantic model structure"],
    ids=["Assess `KrigingParameterInputs` pydantic model structure"],
)
def test_KrigingParameterInputs_model_structure(description):

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(KrigingParameterInputs)
    # ---- Validate correct setting
    assert KrigingParameterInputs.model_config == dict(
        arbitrary_types_allowed=True, title="kriging model parameters ('kriging_parameters')"
    )

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = KrigingParameterInputs.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert set(["validate_realposfloat"]).issubset(KrigingParameterInputs.__dict__)

    # -------------------------
    # ASSERT: 'valid_realposfloat'
    # ---- Check output for 'validate_realposfloat' [VALID]
    assert KrigingParameterInputs.validate_realposfloat(2) == 2.0
    # ---- Check output for 'validate_realposfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert KrigingParameterInputs.validate_realposfloat(-2.0)
    # ---- Check output for 'validate_realposfloat' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real number.")):
        assert KrigingParameterInputs.validate_realposfloat(np.inf)
    # ---- Check output for 'validate_realposfloat' [VALID: None]
    assert KrigingParameterInputs.validate_realposfloat(None) is None
    # ---- Check the applicable fields
    assert field_decor["validate_realposfloat"].info.fields == (
        "anisotropy",
        "correlation_range",
        "search_radius",
    )
    # ---- Check the applicable validation mode
    assert field_decor["validate_realposfloat"].info.mode == "before"


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (dict(), None, "Both 'correlation_range' and 'search_radius' arguments are missing"),
        (
            dict(correlation_range=1.0),
            dict(anisotropy=1e-3, kmin=3, kmax=10, correlation_range=1.0, search_radius=3.0),
            None,
        ),
        (
            dict(correlation_range=1.0, search_radius=5.0),
            dict(anisotropy=1e-3, kmin=3, kmax=10, correlation_range=1.0, search_radius=5.0),
            None,
        ),
        (
            dict(anisotropy=1, correlation_range=2, search_radius=3),
            dict(anisotropy=1.0, kmin=3, kmax=10, correlation_range=2.0, search_radius=3.0),
            None,
        ),
        (
            dict(kmin=3.0, kmax=10.0, correlation_range=1.0),
            None,
            ["Value must be a non-negative integer", "Value must be a non-negative integer"],
        ),
        (
            dict(kmin=10, kmax=3, correlation_range=1.0),
            None,
            "Defined 'kmax' (3) must be greater than or equal to 'kmin' (10)",
        ),
        (
            dict(anisotropy=0.00, kmin=1, kmax=2, correlation_range=-1.0, search_radius=-1.0),
            None,
            [
                "Input should be greater than or equal to 3",
                "Input should be greater than or equal to 3",
                "Value must be a non-negative float",
                "Value must be a non-negative float",
            ],
        ),
        (
            dict(
                anisotropy=np.nan,
                kmin=np.nan,
                kmax=np.nan,
                correlation_range=np.nan,
                search_radius=np.nan,
            ),
            None,
            [
                "Input should be a finite number",
                "Value must be a non-negative integer",
                "Value must be a non-negative integer",
                "Input should be greater than 0",
                "Input should be greater than 0",
            ],
        ),
        (
            dict(
                anisotropy=np.inf,
                kmin=np.inf,
                kmax=np.inf,
                correlation_range=np.inf,
                search_radius=np.inf,
            ),
            None,
            [
                "Value must be a non-negative real number",
                "Value must be a non-negative integer",
                "Value must be a non-negative integer",
                "Value must be a non-negative real number",
                "Value must be a non-negative real number",
            ],
        ),
    ],
    ids=[
        "Empty inputs [invalid]",
        "Produce valid 'search_radius' based on valid 'correlation_range' input",
        "Define both 'search_radius' and 'correlation_range'",
        "Coerce integer inputs for 'anisotropy', 'correlation_range', and 'search_radius'",
        "Invalid float datatyping for 'kmin' and 'kmax'",
        "Enforcing 'kmax' > 'kmin'",
        "Invalid values for numerics [lower limits]",
        "NaN inputs",
        "Inf inputs",
    ],
)
def test_KrigingParameters_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert KrigingParameterInputs.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert KrigingParameterInputs.create(**input)
    else:
        result = KrigingParameterInputs.create(**input)
        assert result == expected


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            dict(),
            dict(
                best_fit_variogram=False,
                coordinate_transform=True,
                extrapolate=False,
                variable="biomass",
                verbose=True,
            ),
            None,
        ),
        (
            dict(best_fit_variogram=3, coordinate_transform=3, extrapolate=3, verbose=3),
            None,
            [
                "Input should be a valid boolean",
                "Input should be a valid boolean",
                "Input should be a valid boolean",
                "Input should be a valid boolean",
            ],
        ),
        (dict(variable="krakens"), None, "Input should be 'biomass'"),
    ],
    ids=[
        "Default values [no inputs, empty dictionary]",
        "Invalid boolean inputs [integers, not bool/str]",
        "Invalid Literal input for 'variable'",
    ],
)
def test_KrigingAnalysis_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert KrigingAnalysis.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert KrigingAnalysis.create(**input)
    else:
        result = KrigingAnalysis.create(**input)
        assert result == expected


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            dict(),
            dict(
                crop_method="transect_ends",
                num_nearest_transects=4,
                mesh_buffer_distance=1.25,
                latitude_resolution=1.25,
                bearing_tolerance=15.0,
            ),
            None,
        ),
        (dict(crop_method="invalid"), None, "Input should be 'transect_ends' or 'convex_hull'"),
        (
            dict(mesh_buffer_distance=1.0, latitude_resolution=1.0, bearing_tolerance=15),
            dict(
                crop_method="transect_ends",
                num_nearest_transects=4,
                mesh_buffer_distance=1.0,
                latitude_resolution=1.0,
                bearing_tolerance=15.0,
            ),
            None,
        ),
        (
            dict(num_nearest_transects=1.0),
            None,
            "Value error, Value must be a non-negative integer.",
        ),
        (dict(bearing_tolerance="a"), None, "Value must be a non-negative real angle"),
        (
            dict(mesh_buffer_distance=-1.0, latitude_resolution=-1.0, bearing_tolerance=-1.0),
            None,
            [
                "Value must be a non-negative float",
                "Value must be a non-negative float",
                "Value must be a non-negative real angle",
            ],
        ),
        (
            dict(
                num_nearest_transects=np.nan,
                mesh_buffer_distance=np.nan,
                latitude_resolution=np.nan,
                bearing_tolerance=np.nan,
            ),
            None,
            [
                "Value must be a non-negative integer.",
                "Input should be greater than 0",
                "Input should be greater than 0",
                "Input should be greater than 0",
            ],
        ),
        (
            dict(
                num_nearest_transects=np.inf,
                mesh_buffer_distance=np.inf,
                latitude_resolution=np.inf,
                bearing_tolerance=np.inf,
            ),
            None,
            [
                "Value must be a non-negative integer.",
                "Value must be a non-negative real number",
                "Value must be a non-negative real number",
                "Value must be a non-negative real angle",
            ],
        ),
        (dict(bearing_tolerance=181.0), None, "Input should be less than or equal to 180"),
    ],
    ids=[
        "Default values [no inputs, empty dictionary]",
        "Invalid Literal for 'crop_method'",
        "Valid int-to-float coercion",
        "Invalid floats where value should be int",
        "Invalid 'bearing_tolerance' [str]",
        "Invalid values below limits",
        "All NaN values for numeric inputs",
        "All Inf values for numeric inputs",
        "Invalid 'bearing_tolerance' input [upper limit]",
    ],
)
def test_MeshCrop_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert MeshCrop.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert MeshCrop.create(**input)
    else:
        result = MeshCrop.create(**input)
        assert result == expected
