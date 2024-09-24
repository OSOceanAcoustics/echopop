import re

# import numpy as np
import pandas as pd
import pytest
from pandera.errors import SchemaError

from echopop.tests.conftest import assert_dataframe_equal
from echopop.utils.validate_df import BaseDataFrame, IsobathData, KrigedMesh


@pytest.mark.parametrize(
    "description",
    ["Assess `BaseDataFrame` pandera model structure"],
    ids=["Assess `BaseDataFrame` pandera model structure"],
)
def test_BaseDataFrame_model_structure(description):

    # -------------------------
    # ASSERT: 'Config' attribute
    # ---- Check existence
    assert "Config" in dir(BaseDataFrame)
    # ---- Check that the attribute comprises the correct entries
    assert set(["metadata", "strict"]).issubset(BaseDataFrame.Config.__dict__)
    # ---- Verify that 'metadata' is a dictionary
    assert isinstance(BaseDataFrame.Config.__dict__["metadata"], dict)
    # ---- Verify that 'strict' is set to 'False'
    assert not BaseDataFrame.Config.__dict__["strict"]

    # -------------------------
    # ASSERT: '_DTYPE_TESTS' and '_DTYPE_COERCION' attributes
    # ---- Check existence
    assert set(["_DTYPE_TESTS", "_DTYPE_COERCION"]).issubset(dir(BaseDataFrame))
    # ---- Verify that both attributes are dictionaries
    assert all(
        [
            isinstance(getattr(BaseDataFrame, attr), dict)
            for attr in ["_DTYPE_TESTS", "_DTYPE_COERCION"]
        ]
    )
    # ---- Check that the dictionary keys comprise the 'int', 'float', and 'str' datatypes
    assert all(
        [
            set(getattr(BaseDataFrame, attr)).issubset(set([int, float, str]))
            for attr in ["_DTYPE_TESTS", "_DTYPE_COERCION"]
        ]
    )

    # -------------------------
    # ASSERT: 'BaseDataFrame' methods
    # ---- Check that all necessary methods exist
    assert set(["get_column_types", "judge", "validate_df"]).issubset(dir(BaseDataFrame))


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            pd.DataFrame(dict(latitude=[0.0], longitude=[0.0])),
            pd.DataFrame(dict(latitude=[0.0], longitude=[0.0])),
            None,
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0])),
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0])),
            None,
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0], longitude=[-181.0, 0.0, 1.0])),
            None,
            "greater_than_or_equal_to(-180.0)",
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 181.0])),
            None,
            "less_than_or_equal_to(180.0)",
        ),
        (
            pd.DataFrame(dict(latitude=[-91.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0])),
            None,
            "greater_than_or_equal_to(-90.0)",
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 91.0], longitude=[1.0, 0.0, 1.0])),
            None,
            "less_than_or_equal_to(90.0)",
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0])),
            None,
            "'.*longitude.*' did not match any columns in the dataframe",
        ),
        (
            pd.DataFrame(dict(longitude=[1.0, 0.0, 1.0])),
            None,
            "'.*latitude.*' did not match any columns in the dataframe",
        ),
        (
            pd.DataFrame(dict()),
            None,
            [
                "'.*longitude.*' did not match any columns in the dataframe",
                "'.*latitude.*' did not match any columns in the dataframe",
            ],
        ),
        (
            pd.DataFrame(dict(latitude=[-1, 0, 1], longitude=[-1, 0, 1])),
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0])),
            None,
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0], longitude=["a", "b", "c"])),
            None,
            "Longitude column must be a Series of 'float64' values",
        ),
        (
            pd.DataFrame(dict(latitude=["a", "b", "c"], longitude=[-1.0, 0.0, 1.0])),
            None,
            "Latitude column must be a Series of 'float64' values",
        ),
        (
            pd.DataFrame(dict(latitude=["a", "b", "c"], longitude=["a", "b", "c"])),
            None,
            [
                "Latitude column must be a Series of 'float64' values",
                "Longitude column must be a Series of 'float64' values",
            ],
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, "c"], longitude=["a", 0.0, 1.0])),
            None,
            [
                "Latitude column must be a Series of 'float64' values",
                "Longitude column must be a Series of 'float64' values",
            ],
        ),
        (
            pd.DataFrame(dict(central_latitude=[-1.0, 0.0, 1.0], longitude_funny=[-1.0, 0.0, 1.0])),
            pd.DataFrame(dict(central_latitude=[-1.0, 0.0, 1.0], longitude_funny=[-1.0, 0.0, 1.0])),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid longitude [lower limit]",
        "Invalid longitude [upper limit]",
        "Invalid latitude [lower limit]",
        "Invalid latitude [upper limit]",
        "Missing column [longitude]",
        "Missing column [latitude]",
        "Missing all columns",
        "Incorrect datatyping but coercible",
        "Incorrect 'Longitude' datatyping and not coercible",
        "Incorrect 'Latitude' datatyping and not coercible",
        "Incorrect 'Latitude' and 'Longitude' datatyping and not coercible",
        "Partially incorrect 'Latitude' and 'Longitude' datatyping and not coercible",
        "Coerced column names based on regex",
    ],
)
def test_IsobathData_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert IsobathData.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert IsobathData.validate_df(input)
    else:
        # Test creation with various parameters
        result = IsobathData.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            pd.DataFrame(dict(latitude=[0.0], longitude=[0.0], fraction=[1.0])),
            pd.DataFrame(dict(latitude=[0.0], longitude=[0.0], fraction=[1.0])),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0], fraction=[0.0, 0.5, 1.0]
                )
            ),
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0], fraction=[0.0, 0.5, 1.0]
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0], fraction=[-0.1, 0.5, 1.0]
                )
            ),
            None,
            "greater_than_or_equal_to(0.0)",
        ),
        (
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0], fraction=[0.0, 0.5, 1.1]
                )
            ),
            None,
            "less_than_or_equal_to(1.0)",
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0])),
            None,
            "'.*fraction.*' did not match any columns in the dataframe",
        ),
        (
            pd.DataFrame(dict()),
            None,
            [
                "'.*longitude.*' did not match any columns in the dataframe",
                "'.*latitude.*' did not match any columns in the dataframe",
                "'.*fraction.*' did not match any columns in the dataframe",
            ],
        ),
        (
            pd.DataFrame(dict(latitude=[-1, 0, 1], longitude=[-1, 0, 1], fraction=[1, 1, 0])),
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0], fraction=[1.0, 1.0, 0.0]
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0], longitude=[-1.0, 0.0, 1.0], fraction=["a", "b", "c"]
                )
            ),
            None,
            "Fraction column must be a Series of 'float64' values",
        ),
        (
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0],
                    longitude=[-1.0, 0.0, 1.0],
                    fraction=["a", 0.25, 0.50],
                )
            ),
            None,
            "Fraction column must be a Series of 'float64' values",
        ),
        (
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0],
                    longitude=[-1.0, 0.0, 1.0],
                    the_entire_or_maybe_part_of_fraction=[1.0, 1.0, 0.0],
                )
            ),
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0],
                    longitude=[-1.0, 0.0, 1.0],
                    the_entire_or_maybe_part_of_fraction=[1.0, 1.0, 0.0],
                )
            ),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid fraction [lower limit]",
        "Invalid fraction [upper limit]",
        "Missing column [fraction]",
        "Missing all columns",
        "Incorrect datatyping but coercible",
        "Incorrect 'Fraction' datatyping and not coercible",
        "Partially incorrect 'Fraction' datatyping and not coercible",
        "Coerced column names based on regex",
    ],
)
def test_KrigedMesh_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert KrigedMesh.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert KrigedMesh.validate_df(input)
    else:
        # Test creation with various parameters
        result = KrigedMesh.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)
