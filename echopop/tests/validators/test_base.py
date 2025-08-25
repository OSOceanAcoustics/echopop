import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError
import pandera as pa

from echopop.validators.base import BaseDictionary, BaseDataFrame, EchopopValidationError


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
        optional_value: int = None
    
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
    
    df = pd.DataFrame({
        "value": [1, 2, 3],
        "name": ["a", "b", "c"]
    })
    
    result = TestDF.validate(df)
    pd.testing.assert_frame_equal(result, df)


def test_base_dataframe_validate_coercion():
    """Test BaseDataFrame validation with type coercion."""
    
    class TestDF(BaseDataFrame):
        value: float
        name: str
    
    df = pd.DataFrame({
        "value": ["1.0", "2.0", "3.0"],  # String values that can be coerced
        "name": ["a", "b", "c"]
    })
    
    result = TestDF.validate(df)
    assert result["value"].dtype == float


def test_base_dataframe_validate_failure():
    """Test BaseDataFrame validation failure."""
    
    class TestDF(BaseDataFrame):
        value: int = pa.Field(ge=0)
        name: str
    
    df = pd.DataFrame({
        "value": [-1, 2, 3],  # Invalid: negative value
        "name": ["a", "b", "c"]
    })
    
    with pytest.raises(pa.errors.SchemaError):
        TestDF.validate(df)


# ==================================================================================================
# Test EchopopValidationError
# ---------------------------
def test_echopop_validation_error_with_exception():
    """Test EchopopValidationError with an exception."""
    
    original_error = ValueError("Original error message")
    error = EchopopValidationError(original_error)
    
    assert error.exception == original_error
    assert str(error) == "Original error message"


def test_echopop_validation_error_without_exception():
    """Test EchopopValidationError without an exception."""
    
    error = EchopopValidationError()
    
    assert error.exception is None
    assert str(error) == ""


def test_echopop_validation_error_inheritance():
    """Test that EchopopValidationError inherits from Exception."""
    
    error = EchopopValidationError()
    assert isinstance(error, Exception)


# ==================================================================================================
# Integration tests
# -----------------
def test_base_integration_custom_model():
    """Test integration of BaseDictionary with custom validation."""
    
    class CustomModel(BaseDictionary):
        positive_number: float
        optional_string: str = None
        
        @pytest.fixture  # This would be a field_validator in real code
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
