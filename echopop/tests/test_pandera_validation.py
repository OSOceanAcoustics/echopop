import re

import numpy as np
import pandas as pd
import pytest
from pandera.errors import SchemaError
from pandera.typing import Series

from ..utils.validate_df import (
    AcousticData,
    BaseDataFrame,
    CatchBiodata,
    GeoStrata,
    IsobathData,
    KrigedMesh,
    KSStrata,
    LengthBiodata,
    SpecimenBiodata,
)
from .conftest import assert_dataframe_equal

pytestmark = pytest.mark.skip(reason="Temporarily disable this module")


####################################################################################################
# MOCK MODEL FIXTURE
# --------------------------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def mock_model():

    class MOCKModel(BaseDataFrame):
        dum1: Series[int]
        dum2: Series[float]

    return MOCKModel


@pytest.fixture
def MOCK_DF() -> pd.DataFrame:

    return pd.DataFrame(dict(dum1=[1, 2, 3], dum2=[1.0, 2.0, 3.0]))


####################################################################################################
# BASE DATAFRAME MODEL STRUCTURE
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "description",
    ["Assess `BaseDataFrame` pandera model structure"],
    ids=["Assess `BaseDataFrame` pandera model structure"],
)
def test_BaseDataFrame_model_structure(description, MOCK_DF, mock_model):

    # -------------------------
    # ASSERT: 'Config' attribute
    # ---- Check existence
    assert "Config" in dir(BaseDataFrame)
    # ---- Check that the attribute comprises the correct entries
    assert set(["metadata", "strict"]).issubset(BaseDataFrame.Config.__dict__)
    # ---- Verify that 'metadata' is a dictionary
    assert isinstance(BaseDataFrame.Config.__dict__["metadata"], dict)
    # ---- Verify that 'strict' is set to 'filter'
    assert BaseDataFrame.Config.__dict__["strict"] == "filter"

    # -------------------------
    # ASSERT: Base Model Inheritance
    assert mock_model.__base__ == BaseDataFrame

    # -------------------------
    # ASSERT: 'BaseDataFrame' methods
    # ---- Check that all necessary methods exist
    assert set(["pre_validate", "validate_df"]).issubset(dir(BaseDataFrame))
    # ---- Validate 'pre_validate'
    assert mock_model.pre_validate(MOCK_DF).equals(
        pd.DataFrame(dict(dum1=[1, 2, 3], dum2=[1.0, 2.0, 3.0]))
    )
    # ---- Validate 'validate_df'
    assert mock_model.validate_df(MOCK_DF).equals(
        pd.DataFrame(dict(dum1=[1, 2, 3], dum2=[1.0, 2.0, 3.0]))
    )


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
            "'longitude':  No valid column name found matching the accepted regex pattern "
            + "[.*longitude.*]",
        ),
        (
            pd.DataFrame(dict(longitude=[1.0, 0.0, 1.0])),
            None,
            "'latitude':  No valid column name found matching the accepted regex pattern "
            + "[.*latitude.*]",
        ),
        (
            pd.DataFrame(dict()),
            None,
            [
                "- 'latitude':  No valid column name found matching the accepted regex pattern "
                + "[.*latitude.*]",
                "- 'longitude':  No valid column name found matching the accepted regex pattern "
                + "[.*longitude.*]",
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
            [
                "- '.*longitude.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'longitude': Expected series 'longitude' to have type float64, got object.",
                "- 'longitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'longitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
        ),
        (
            pd.DataFrame(dict(latitude=["a", "b", "c"], longitude=[-1.0, 0.0, 1.0])),
            None,
            [
                "- '.*latitude.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'latitude': Expected series 'latitude' to have type float64, got object.",
                "- 'latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
        ),
        (
            pd.DataFrame(dict(latitude=["a", "b", "c"], longitude=["a", "b", "c"])),
            None,
            [
                "- '.*longitude.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'longitude': Expected series 'longitude' to have type float64, got object.",
                "- 'longitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'longitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- '.*latitude.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'latitude': Expected series 'latitude' to have type float64, got object.",
                "- 'latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
        ),
        (
            pd.DataFrame(dict(latitude=[-1.0, 0.0, "c"], longitude=["a", 0.0, 1.0])),
            None,
            [
                "- '.*longitude.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'longitude': Expected series 'longitude' to have type float64, got object.",
                "- 'longitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'longitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- '.*latitude.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'latitude': Expected series 'latitude' to have type float64, got object.",
                "- 'latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
        ),
        (
            pd.DataFrame(dict(central_latitude=[-1.0, 0.0, 1.0], longitude_funny=[-1.0, 0.0, 1.0])),
            pd.DataFrame(dict(central_latitude=[-1.0, 0.0, 1.0], longitude_funny=[-1.0, 0.0, 1.0])),
            None,
        ),
        (
            pd.DataFrame(
                dict(central_latitude=[np.nan, 0.0, 1.0], longitude_funny=[-1.0, np.nan, 1.0])
            ),
            pd.DataFrame(dict(central_latitude=[1.0], longitude_funny=[1.0])),
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
        "Drop invalid NaN rows",
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
            "- 'fraction':  No valid column name found matching the accepted regex pattern "
            + "[.*fraction.*]",
        ),
        (
            pd.DataFrame(dict()),
            None,
            [
                "- 'fraction':  No valid column name found matching the accepted regex "
                + "pattern [.*fraction.*]",
                "- 'latitude':  No valid column name found matching the accepted regex "
                + "pattern [.*latitude.*]",
                "'longitude':  No valid column name found matching the accepted regex "
                + "pattern [.*longitude.*]",
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
            [
                "- '.*fraction.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'fraction': Expected series 'fraction' to have type float64, got object.",
                "- 'fraction': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'fraction': TypeError(\"'<=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
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
            [
                "- '.*fraction.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'fraction': Expected series 'fraction' to have type float64, got object.",
                "- 'fraction': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'fraction': TypeError(\"'<=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
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
        (
            pd.DataFrame(
                dict(
                    latitude=[-1.0, 0.0, 1.0],
                    longitude=[-1.0, 0.0, 1.0],
                    fraction=[np.nan, 0.5, 1.0],
                )
            ),
            pd.DataFrame(dict(latitude=[0.0, 1.0], longitude=[0.0, 1.0], fraction=[0.5, 1.0])),
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
        "Drop invalid NaN rows",
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


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            pd.DataFrame(dict(haul=[1], northlimit_latitude=[1.0], stratum=[1])),
            pd.DataFrame(dict(haul=[1], northlimit_latitude=[1.0], stratum=[1])),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-91.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            None,
            "greater_than_or_equal_to(-90.0)",
        ),
        (
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 91.0], stratum=[1, 2, 3])
            ),
            None,
            "less_than_or_equal_to(90.0)",
        ),
        (
            pd.DataFrame(dict(northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
            "'haul':  No valid column name found matching the accepted regex pattern "
            + "[(haul|haul_num|haul_start|haul_end)]",
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], stratum=[1, 2, 3])),
            None,
            "The following columns were missing from the DataFrame: ['northlimit_latitude']",
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0])),
            None,
            "'stratum':  No valid column name found matching the accepted regex pattern "
            + "[stratum(_(?:num|inpfc))?$]",
        ),
        (
            pd.DataFrame(dict()),
            None,
            [
                "The following columns were missing from the DataFrame: ['northlimit_latitude']",
                "'haul':  No valid column name found matching the accepted regex pattern "
                + "[(haul|haul_num|haul_start|haul_end)]",
                "'stratum':  No valid column name found matching the accepted regex pattern "
                + "[stratum(_(?:num|inpfc))?$]",
            ],
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], northlimit_latitude=[-1, 0, 1], stratum=[1, 2, 3])),
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul=[1.0, 2.0, 3.0], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul=[1.0, 2.0, 3.0],
                    northlimit_latitude=[-1.0, 0.0, 1.0],
                    stratum=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul=[1.0, 2.0, 3.0],
                    northlimit_latitude=[-1.0, 0.0, 1.0],
                    stratum=["1", "2", "3"],
                )
            ),
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul=[1.0, 2.0, 3.0],
                    northlimit_latitude=[-1.0, 0.0, 1.0],
                    stratum=[1, "2a", "3b"],
                )
            ),
            pd.DataFrame(
                dict(haul=[1, 2, 3], northlimit_latitude=[-1.0, 0.0, 1.0], stratum=[1, "2a", "3b"])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul=[1.0, 2.0, 3.0], northlimit_latitude=["a", "b", "c"], stratum=[1, 2, 3])
            ),
            None,
            (
                "The following DataFrame validation errors were flagged for None: \n"
                "    - 'northlimit_latitude': Values within column are not coercible to the "
                + "expected datatype 'float64'.\n"
                "    - 'northlimit_latitude': Expected series 'northlimit_latitude' to have "
                + "type float64, got object.\n"
                "    - 'northlimit_latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").\n"
                "    - 'northlimit_latitude': TypeError(\"'<=' not supported between instances "
                + "of 'str' and 'float'\")."
            ),
        ),
        (
            pd.DataFrame(
                dict(haul=[1.0, 2.0, 3.0], northlimit_latitude=["a", 0.0, 1.0], stratum=[1, 2, 3])
            ),
            None,
            (
                "The following DataFrame validation errors were flagged for None: \n"
                "    - 'northlimit_latitude': Values within column are not coercible to the "
                + "expected datatype 'float64'.\n"
                "    - 'northlimit_latitude': Expected series 'northlimit_latitude' to have "
                + "type float64, got object.\n"
                "    - 'northlimit_latitude': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").\n"
                "    - 'northlimit_latitude': TypeError(\"'<=' not supported between instances "
                + "of 'str' and 'float'\")."
            ),
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1.0, 2.0, 3.0],
                    northlimit_latitude=[-1.0, 0.0, 1.0],
                    stratum_num=[1, 2, 3],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1.0, 2.0, 3.0],
                    northlimit_latitude=[-1.0, 0.0, 1.0],
                    stratum_num=[1, 2, 3],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[np.nan, 2, 3, 4],
                    northlimit_latitude=[-1.0, 0.0, 1.0, 2.0],
                    stratum_num=[1, 2, np.nan, 4],
                )
            ),
            pd.DataFrame(dict(haul_num=[2, 4], northlimit_latitude=[0.0, 2.0], stratum_num=[2, 4])),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid 'northlimit_latitude' [lower limit]",
        "Invalid 'northlimit_latitude' [upper limit]",
        "Missing column [haul]",
        "Missing column [northlimit_latitude]",
        "Missing column [stratum]",
        "Missing all columns",
        "Incorrect 'northlimit_latitude' datatyping but coercible",
        "Haul column datatyping [float]",
        "Stratum column datatyping [float]",
        "Stratum column datatyping [str]",
        "Stratum column datatyping [mixed]",
        "Incorrect 'northlimit_latitude' datatyping and not coercible",
        "Partially incorrect 'northlimit_latitude' datatyping and not coercible",
        "Coerced column names based on [haul, stratum]",
        "Drop invalid NaN rows",
    ],
)
def test_GeoStrata_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert GeoStrata.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert GeoStrata.validate_df(input)
    else:
        # Test creation with various parameters
        result = GeoStrata.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            pd.DataFrame(dict(haul=[1], fraction=[1.0], stratum=[1])),
            pd.DataFrame(dict(fraction=[1.0], haul=[1], stratum=[1])),
            None,
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[-1.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
            "greater_than_or_equal_to(0.0)",
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 2.0], stratum=[1, 2, 3])),
            None,
            "less_than_or_equal_to(1.0)",
        ),
        (
            pd.DataFrame(dict(fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
            " - 'haul':  No valid column name found matching the accepted regex pattern "
            + "[(haul|haul_num|haul_start|haul_end)]",
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], stratum=[1, 2, 3])),
            None,
            "- 'fraction':  No valid column name found matching the accepted regex pattern "
            + "[.*fraction.*]",
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0])),
            None,
            "- 'stratum':  No valid column name found matching the accepted regex pattern "
            + "[stratum(_(?:num|inpfc))?$]",
        ),
        (
            pd.DataFrame(dict()),
            None,
            [
                "- 'fraction':  No valid column name found matching the accepted regex pattern "
                + "[.*fraction.*]",
                "- 'haul':  No valid column name found matching the accepted regex pattern "
                + "[(haul|haul_num|haul_start|haul_end)]",
                "- 'stratum':  No valid column name found matching the accepted regex pattern "
                + "[stratum(_(?:num|inpfc))?$]",
            ],
        ),
        (
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0, 0, 1], stratum=[1, 2, 3])),
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
        ),
        (
            pd.DataFrame(dict(haul=[1.0, 2.0, 3.0], fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul=[1.0, 2.0, 3.0], fraction=[0.0, 0.0, 1.0], stratum=[1.0, 2.0, 3.0])
            ),
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul=[1.0, 2.0, 3.0], fraction=[0.0, 0.0, 1.0], stratum=["1", "2", "3"])
            ),
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0], stratum=[1, 2, 3])),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul=[1.0, 2.0, 3.0], fraction=[0.0, 0.0, 1.0], stratum=[1, "2a", "3b"])
            ),
            pd.DataFrame(dict(haul=[1, 2, 3], fraction=[0.0, 0.0, 1.0], stratum=[1, "2a", "3b"])),
            None,
        ),
        (
            pd.DataFrame(dict(haul=[1.0, 2.0, 3.0], fraction=["a", "b", "c"], stratum=[1, 2, 3])),
            None,
            [
                "- '.*fraction.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'fraction': Expected series 'fraction' to have type float64, got object.",
                "- 'fraction': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'fraction': TypeError(\"'<=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
        ),
        (
            pd.DataFrame(dict(haul=[1.0, 2.0, 3.0], fraction=["a", 0.0, 1.0], stratum=[1, 2, 3])),
            None,
            [
                "- '.*fraction.*': Values within column are not coercible to the expected "
                + "datatype 'float64'.",
                "- 'fraction': Expected series 'fraction' to have type float64, got object.",
                "- 'fraction': TypeError(\"'>=' not supported between instances "
                + "of 'str' and 'float'\").",
                "- 'fraction': TypeError(\"'<=' not supported between instances "
                + "of 'str' and 'float'\").",
            ],
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1.0, 2.0, 3.0], fraction_catch=[0.0, 0.0, 1.0], stratum_num=[1, 2, 3]
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1.0, 2.0, 3.0], fraction_catch=[0.0, 0.0, 1.0], stratum_num=[1, 2, 3]
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[np.nan, 2, 3, 4],
                    catch_fraction=[0.0, np.nan, 0.5, 1.0],
                    stratum_num=[1, 2, np.nan, 4],
                )
            ),
            pd.DataFrame(dict(haul_num=[4], catch_fraction=[1.0], stratum_num=[4])),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid 'fraction' [lower limit]",
        "Invalid 'fraction' [upper limit]",
        "Missing column [haul]",
        "Missing column [fraction]",
        "Missing column [stratum]",
        "Missing all columns",
        "Incorrect 'fraction' datatyping but coercible",
        "Haul column datatyping [float]",
        "Stratum column datatyping [float]",
        "Stratum column datatyping [str]",
        "Stratum column datatyping [mixed]",
        "Incorrect 'fraction' datatyping and not coercible",
        "Partially incorrect 'fraction' datatyping and not coercible",
        "Coerced column names based on [haul, stratum, fraction]",
        "Drop invalid NaN rows",
    ],
)
def test_KSStrata_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert KSStrata.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert KSStrata.validate_df(input)
    else:
        # Test creation with various parameters
        result = KSStrata.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            pd.DataFrame(dict(haul_num=[1], haul_weight=[1.0], species_id=[1])),
            pd.DataFrame(dict(haul_num=[1], haul_weight=[1.0], species_id=[1])),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=[1, 2, 3])
            ),
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=[1, 2, 3])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[-1.0, 2.0, 3.0], species_id=[1, 2, 3])
            ),
            None,
            "greater_than_or_equal_to(0.0)",
        ),
        (
            pd.DataFrame(dict(haul_weight=[1.0, 2.0, 3.0], species_id=[1, 2, 3])),
            None,
            "The following columns were missing from the DataFrame: ['haul_num']",
        ),
        (
            pd.DataFrame(dict(haul_num=[1, 2, 3], species_id=[1, 2, 3])),
            None,
            "The following columns were missing from the DataFrame: ['haul_weight']",
        ),
        (
            pd.DataFrame(dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0])),
            None,
            "The following columns were missing from the DataFrame: ['species_id']",
        ),
        (
            pd.DataFrame(dict()),
            None,
            "The following columns were missing from the DataFrame: "
            + "['haul_num', 'haul_weight', 'species_id']",
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1.0, 2.0, 3.0], haul_weight=[1.0, 2.0, 3.0], species_id=[1, 2, 3])
            ),
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=[1, 2, 3])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=[1, 2, 3])
            ),
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=[1.0, 2.0, 3.0])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=["1a", "2b", "3c"])
            ),
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=["1a", "2b", "3c"])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=[1, "2b", "3c"])
            ),
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, 2.0, 3.0], species_id=[1, "2b", "3c"])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=["a", "b", "c"], species_id=[1, 2, 3])
            ),
            None,
            ["Expected series 'haul_weight' to have type float64, got object", "TypeError"],
        ),
        (
            pd.DataFrame(
                dict(haul_num=[1, 2, 3], haul_weight=[1.0, "b", "c"], species_id=[1, 2, 3])
            ),
            None,
            ["Expected series 'haul_weight' to have type float64, got object", "TypeError"],
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[np.nan, 2, 3, 4],
                    haul_weight=[1.0, 2.0, 3.0, 4.0],
                    species_id=[1, 2, np.nan, 4],
                )
            ),
            pd.DataFrame(dict(haul_num=[2, 4], haul_weight=[2.0, 4.0], species_id=[2, 4])),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid 'haul_weight' [lower limit]",
        "Missing column [haul]",
        "Missing column [haul_weight]",
        "Missing column [species_id]",
        "Missing all columns",
        "Haul column datatyping [float]",
        "Species_id column datatyping [float]",
        "Species_id column datatyping [str]",
        "Stratum column datatyping [mixed]",
        "Incorrect 'haul_weight' datatyping and not coercible",
        "Partially incorrect 'haul_weight' datatyping and not coercible",
        "Drop invalid NaN rows",
    ],
)
def test_CatchBiodata_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert CatchBiodata.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert CatchBiodata.validate_df(input)
    else:
        # Test creation with various parameters
        result = CatchBiodata.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (  # 1
            pd.DataFrame(
                dict(haul_num=[1], length=[1.0], length_count=[1], sex=[1], species_id=[1])
            ),
            pd.DataFrame(
                dict(haul_num=[1], length=[1.0], length_count=[1], sex=[1], species_id=[1])
            ),
            None,
        ),
        (  # 2
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            None,
        ),
        (  # 3
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[-1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            None,
            "greater_than(0.0)",
        ),
        (  # 4
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[-1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            None,
            "greater_than_or_equal_to(0)",
        ),
        (  # 5
            pd.DataFrame(dict(length=[1.0], length_count=[1], sex=[1], species_id=[1])),
            None,
            "The following columns were missing from the DataFrame: ['haul_num']",
        ),
        (  # 6
            pd.DataFrame(dict(haul_num=[1], length_count=[1], sex=[1], species_id=[1])),
            None,
            "The following columns were missing from the DataFrame: ['length']",
        ),
        (  # 1
            pd.DataFrame(dict(haul_num=[1], length=[1.0], sex=[1], species_id=[1])),
            None,
            "The following columns were missing from the DataFrame: ['length_count']",
        ),
        (
            pd.DataFrame(dict(haul_num=[1], length=[1.0], length_count=[1], species_id=[1])),
            None,
            "The following columns were missing from the DataFrame: ['sex']",
        ),
        (
            pd.DataFrame(dict(haul_num=[1], length=[1.0], length_count=[1], sex=[1])),
            None,
            "The following columns were missing from the DataFrame: ['species_id']",
        ),
        (
            pd.DataFrame(dict()),
            None,
            "The following columns were missing from the DataFrame: "
            + "['haul_num', 'length', 'length_count', 'sex', 'species_id']",
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1.0, 2.0, 3.0],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=["1a", "2b", "3c"],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=["1a", "2b", "3c"],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, "2b", "3c"],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=[1, 2, 3],
                    species_id=[1, "2b", "3c"],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=["male", "female", "unsexed"],
                    species_id=[1, 2, 3],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=["male", "female", "unsexed"],
                    species_id=[1, 2, 3],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=["m", "f", "u"],
                    species_id=[1, 2, 3],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=["m", "f", "u"],
                    species_id=[1, 2, 3],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=["sex1", "sex2", "sex3"],
                    species_id=[1, 2, 3],
                )
            ),
            None,
            "datatype should either be all 'int', or all 'str' contained within",
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    length_count=[1, 2, 3],
                    sex=["male", "female", "sex3"],
                    species_id=[1, 2, 3],
                )
            ),
            None,
            "datatype should either be all 'int', or all 'str' contained within",
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[np.nan, 2, 3, np.nan, 5],
                    length=[1.0, np.nan, 3.0, 4.0, 5.0],
                    length_count=[1, 2, 3, 4, 5],
                    sex=[1, 2, 3, 4, 5],
                    species_id=[1, np.nan, 3, 4, 5],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[3, 5],
                    length=[3.0, 5.0],
                    length_count=[3, 5],
                    sex=[3, 5],
                    species_id=[3, 5],
                )
            ),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid 'length' [lower limit]",
        "Invalid 'length_count' [lower limit]",
        "Missing column [haul_num]",
        "Missing column [length]",
        "Missing column [length_count]",
        "Missing column [sex]",
        "Missing column [species_id]",
        "Missing all columns",
        "Haul column datatyping [float]",
        "Length_count column datatyping [float]",
        "Species_id column datatyping [float]",
        "Species_id column datatyping [str]",
        "Species_id column datatyping [mixed]",
        "Sex column Literal input ['male', 'female', 'unsexed']",
        "Sex column Literal input ['m', 'f', 'u']",
        "Invalid 'sex' column Literal input ['sex1', 'sex2', 'sex3']",
        "Erroneous 'sex' Literal",
        "Drop invalid NaN rows",
    ],
)
def test_LengthBiodata_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert LengthBiodata.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert LengthBiodata.validate_df(input)
    else:
        # Test creation with various parameters
        result = LengthBiodata.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            pd.DataFrame(
                dict(age=[0], haul_num=[1], length=[1.0], sex=[1], species_id=[1], weight=[1.0])
            ),
            pd.DataFrame(
                dict(age=[0], haul_num=[1], length=[1.0], sex=[1], species_id=[1], weight=[1.0])
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[-1, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
            "greater_than_or_equal_to(0)",
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[-1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
            "greater_than(0.0)",
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[-1.0, 2.0, 3.0],
                )
            ),
            None,
            "greater_than(0.0)",
        ),
        (
            pd.DataFrame(dict(haul_num=[1], length=[1.0], sex=[1], species_id=[1], weight=[1.0])),
            None,
            "- The following columns were missing from the DataFrame: ['age'].",
        ),
        (
            pd.DataFrame(dict(age=[0], length=[1.0], sex=[1], species_id=[1], weight=[1.0])),
            None,
            "- The following columns were missing from the DataFrame: ['haul_num'].",
        ),
        (
            pd.DataFrame(dict(age=[0], haul_num=[1], sex=[1], species_id=[1], weight=[1.0])),
            None,
            "- The following columns were missing from the DataFrame: ['length'].",
        ),
        (
            pd.DataFrame(dict(age=[0], haul_num=[1], length=[1.0], species_id=[1], weight=[1.0])),
            None,
            "- The following columns were missing from the DataFrame: ['sex'].",
        ),
        (
            pd.DataFrame(dict(age=[0], haul_num=[1], length=[1.0], sex=[1], weight=[1.0])),
            None,
            "- The following columns were missing from the DataFrame: ['species_id'].",
        ),
        (
            pd.DataFrame(dict(age=[0], haul_num=[1], length=[1.0], sex=[1], species_id=[1])),
            None,
            "- The following columns were missing from the DataFrame: ['weight'].",
        ),
        (
            pd.DataFrame(dict()),
            None,
            "- The following columns were missing from the DataFrame: "
            + "['age', 'haul_num', 'length', 'sex', 'species_id', 'weight'].",
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0.0, 1.0, 2.0],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1.0, 2.0, 3.0],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1.0, 2.0, 3.0],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=["1a", "2b", "3c"],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=["1a", "2b", "3c"],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, "2b", "3c"],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=[1, 2, 3],
                    species_id=[1, "2b", "3c"],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=["male", "female", "unsexed"],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=["male", "female", "unsexed"],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=["m", "f", "u"],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=["m", "f", "u"],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=["sex1", "sex2", "sex3"],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
            "- 'sex': Column datatype should either be all 'int', or all 'str' contained within "
            + "['male', 'm', 'female', 'f', 'unsexed', 'u'].",
        ),
        (
            pd.DataFrame(
                dict(
                    age=[0, 1, 2],
                    haul_num=[1, 2, 3],
                    length=[1.0, 2.0, 3.0],
                    sex=["male", "female", "sex3"],
                    species_id=[1, 2, 3],
                    weight=[1.0, 2.0, 3.0],
                )
            ),
            None,
            "- 'sex': Column datatype should either be all 'int', or all 'str' contained within "
            + "['male', 'm', 'female', 'f', 'unsexed', 'u'].",
        ),
        (
            pd.DataFrame(
                dict(
                    age=[np.nan, 2, 3, 4, 5, 6],
                    haul_num=[1, np.nan, 3, 4, 5, 6],
                    length=[1.0, 2.0, np.nan, np.nan, 5.0, 6.0],
                    sex=[1, 2, 3, 4, 5, 6],
                    species_id=[1, 2, 3, 4, np.nan, 6],
                    weight=[1.0, 2.0, 3.0, 4.0, 5.0, np.nan],
                )
            ),
            pd.DataFrame(
                dict(
                    age=[np.nan, 6],
                    haul_num=[1, 6],
                    length=[1.0, 6.0],
                    sex=[1, 6],
                    species_id=[1.0, 6.0],
                    weight=[1.0, np.nan],
                )
            ),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid 'age' [lower limit]",
        "Invalid 'length' [lower limit]",
        "Invalid 'weight' [lower limit]",
        "Missing column [age]",
        "Missing column [haul_num]",
        "Missing column [length]",
        "Missing column [sex]",
        "Missing column [species_id]",
        "Missing column [weight]",
        "Missing all columns",
        "Age column datatyping [float]",
        "Haul column datatyping [float]",
        "Species_id column datatyping [float]",
        "Species_id column datatyping [str]",
        "Species_id column datatyping [mixed]",
        "Sex column Literal input ['male', 'female', 'unsexed']",
        "Sex column Literal input ['m', 'f', 'u']",
        "Invalid 'sex' column Literal input ['sex1', 'sex2', 'sex3']",
        "Erroneous 'sex' Literal",
        "Drop invalid NaN rows",
    ],
)
def test_SpecimenBiodata_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert SpecimenBiodata.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert SpecimenBiodata.validate_df(input)
    else:
        # Test creation with various parameters
        result = SpecimenBiodata.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            pd.DataFrame(
                dict(
                    haul=[1],
                    latitude=[1.0],
                    longitude=[1.0],
                    nasc=[1.0],
                    transect_num=[1],
                    transect_spacing=[1.0],
                    vessel_log_start=[0.0],
                    vessel_log_end=[1.0],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1],
                    latitude=[1.0],
                    longitude=[1.0],
                    nasc=[1.0],
                    transect_num=[1],
                    transect_spacing=[1.0],
                    distance_s=[0.0],
                    distance_e=[1.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[-91.0, 2.0, 91.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
            ["less_than_or_equal_to(90.0)", "greater_than_or_equal_to(-90.0)"],
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[-181.0, 2.0, 181.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
            ["less_than_or_equal_to(180.0)", "greater_than_or_equal_to(-180.0)"],
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[-1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
            "greater_than_or_equal_to(0.0)",
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[-3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
            "greater_than_or_equal_to(0.0)",
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[-1.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
            "greater_than_or_equal_to(0.0)",
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[-1.0, 2.0, 3.0],
                )
            ),
            None,
            "greater_than_or_equal_to(0.0)",
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1.0, 2.0, 3.0],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1.0, 2.0, 3.0],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_name=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_name=[1, 2, 3],
                    latitude=[1.0, 2.0, 3.0],
                    longitude=[1.0, 2.0, 3.0],
                    nasc=[1.0, 2.0, 3.0],
                    transect_num=[1, 2, 3],
                    transect_spacing=[3.0, 3.0, 3.0],
                    distance_s=[0.0, 1.0, 2.0],
                    distance_e=[1.0, 2.0, 3.0],
                )
            ),
            None,
        ),
        (
            pd.DataFrame(
                dict(
                    haul_num=[
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        7,
                        np.nan,
                        np.nan,
                        10,
                    ],
                    latitude=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    longitude=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    nasc=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    transect_num=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    transect_spacing=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    distance_s=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    distance_e=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                )
            ),
            pd.DataFrame(
                dict(
                    haul_num=[7, 10],
                    latitude=[7.0, 10.0],
                    longitude=[7.0, 10.0],
                    nasc=[7.0, 10.0],
                    transect_num=[7, 10],
                    transect_spacing=[7.0, 10.0],
                    distance_s=[7.0, 10.0],
                    distance_e=[7.0, 10.0],
                )
            ),
            None,
        ),
    ],
    ids=[
        "Simple DataFrame input [single row]",
        "Simple DataFrame input [multiple rows]",
        "Invalid 'latitude' [lower and upper limit]",
        "Invalid 'longitude' [lower and upper limit]",
        "Invalid 'nasc' [lower limit]",
        "Invalid 'transect_spacing' [lower limit]",
        "Invalid 'vessel_log_start' [lower limit]",
        "Invalid 'vessel_log_end' [lower limit]",
        "Haul column datatyping [float]",
        "Transect column datatyping [float]",
        "Coerced column names based on [haul, transect]",
        "Drop invalid NaN rows",
    ],
)
def test_AcousticData_model(input, expected, exception):

    if exception:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(SchemaError, match=re.escape(e)):
                    assert AcousticData.validate_df(input)
        else:
            with pytest.raises(SchemaError, match=re.escape(exception)):
                assert AcousticData.validate_df(input)
    else:
        # Test creation with various parameters
        result = AcousticData.validate_df(input)
        assert_dataframe_equal(result, result.dtypes, expected)
