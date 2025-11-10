from typing import Optional

import pandas as pd
import pandera.pandas as pa
import pytest
from pydantic import ValidationError, field_validator

from echopop.core.validators import BaseDataFrame, BaseDictionary


# ==================================================================================================
# Test BaseDictionary
# -------------------
def test_base_dictionary_judge_success():
    """Test BaseDictionary judge method with valid data."""

    class TestDict(BaseDictionary):
        value: int
        name: str

    result = TestDict.judge(value=42, name="test")
    assert result.value == 42
    assert result.name == "test"


def test_base_dictionary_judge_validation_error():
    """Test BaseDictionary judge method with invalid data."""

    class TestDict(BaseDictionary):
        value: int
        name: str

    with pytest.raises(ValidationError):
        TestDict.judge(value="not_an_int", name="test")


def test_base_dictionary_create_success():
    """Test BaseDictionary create factory method."""

    class TestDict(BaseDictionary):
        value: int
        name: str = "default"

    result = TestDict.create(value=42)
    expected = {"value": 42, "name": "default"}
    assert result == expected


def test_base_dictionary_create_exclude_none():
    """Test BaseDictionary create method excludes None values."""

    class TestDict(BaseDictionary):
        value: int
        optional_value: Optional[int] = None

    result = TestDict.create(value=42)
    expected = {"value": 42}
    assert result == expected


# ==================================================================================================
# Test BaseDataFrame
# ------------------
def test_base_dataframe_pre_validate():
    """Test BaseDataFrame pre_validate method."""

    class TestDF(BaseDataFrame):
        value: int

    df = pd.DataFrame({"value": [1, 2, 3]})
    result = TestDF.pre_validate(df)

    # Should return the same DataFrame
    pd.testing.assert_frame_equal(result, df)


def test_base_dataframe_validate_success():
    """Test BaseDataFrame validate method with valid DataFrame."""

    class TestDF(BaseDataFrame):
        value: int = pa.Field(ge=0)
        name: str

    df = pd.DataFrame({"value": [1, 2, 3], "name": ["a", "b", "c"]})

    result = TestDF.validate(df)
    pd.testing.assert_frame_equal(result, df)


def test_base_dataframe_validate_coercion():
    """Test BaseDataFrame validation with type coercion."""

    class TestDF(BaseDataFrame):
        value: float
        name: str

    df = pd.DataFrame(
        {
            "value": ["1.0", "2.0", "3.0"],  # String values that can be coerced
            "name": ["a", "b", "c"],
        }
    )

    result = TestDF.validate(df)
    assert result["value"].dtype == float


def test_base_dataframe_validate_failure():
    """Test BaseDataFrame validation failure."""

    class TestDF(BaseDataFrame):
        value: int = pa.Field(ge=0)
        name: str

    df = pd.DataFrame({"value": [-1, 2, 3], "name": ["a", "b", "c"]})  # Invalid: negative value

    with pytest.raises(ValueError):
        TestDF.validate(df)


# ==================================================================================================
# Integration tests
# -----------------
def test_base_integration_custom_model():
    """Test integration of BaseDictionary with custom validation."""

    class CustomModel(BaseDictionary):
        positive_number: float
        optional_string: Optional[str] = None

        @field_validator("positive_number")
        def validate_positive(cls, v):
            if v <= 0:
                raise ValueError("Must be positive")
            return v

    # Test valid data
    result = CustomModel.create(positive_number=5.5, optional_string="test")
    expected = {"positive_number": 5.5, "optional_string": "test"}
    assert result == expected

    # Test excluding None
    result_no_optional = CustomModel.create(positive_number=3.14)
    expected_no_optional = {"positive_number": 3.14}
    assert result_no_optional == expected_no_optional
