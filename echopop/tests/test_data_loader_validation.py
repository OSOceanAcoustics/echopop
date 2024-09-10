import pytest
from pydantic import ValidationError

from ..utils.validate import *


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"filename": "blurgh/blargh", "sheetname": "sheet1"}, None),
        ({"sheetname": "sheet1"}, ValidationError),
        ({"filename": "blurgh/blargh"}, ValidationError),
        ({}, ValidationError),
        ({"filename": None, "sheetname": "sheet1"}, ValidationError),
        ({"filename": "blurgh/blargh", "sheetname": None}, ValidationError),
        ({"filename": 1, "sheetname": "sheet1"}, ValidationError),
        ({"filename": "blurgh/blargh", "sheetname": 1}, ValidationError),
        ({"filename": "blurgh/blarg", "sheetname": "sheet1", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `FileSettings`",
        "Missing 'filename'",
        "Missing 'sheetname'",
        "Empty dictionary",
        "filename key missing value",
        "Sheetname key missing value",
        "filename value not a string",
        "Sheetname value not a string",
        "Excess keys (valid)",
    ],
)
def test_FileSettings(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert FileSettings(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert FileSettings(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(FileSettings(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(FileSettings(**input).model_dump(exclude_none=True))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"pattern": "pretty", "label": "ponies"}, None),
        ({"pattern": "pretty"}, ValidationError),
        ({"label": "ponies"}, ValidationError),
        ({}, ValidationError),
        ({"pattern": None, "label": "ponies"}, ValidationError),
        ({"pattern": "pretty", "label": None}, ValidationError),
        ({"pattern": 1, "label": "ponies"}, ValidationError),
        ({"pattern": "pretty", "label": 1}, ValidationError),
        ({"pattern": "pretty", "label": "ponies", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `PatternParts`",
        "Missing 'label'",
        "Missing 'pattern'",
        "Empty dictionary",
        "Pattern key missing value",
        "Label key missing value",
        "filename value not a string",
        "Sheetname value not a string",
        "Excess keys (valid)",
    ],
)
def test_PatternParts(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert PatternParts(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert PatternParts(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(PatternParts(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(PatternParts(**input).model_dump(exclude_none=True))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"filename": "blargh", "sheetname": "sheet1"}, None),
        ({"sheetname": "sheet1"}, ValidationError),
        ({"filename": "blargh"}, ValidationError),
        ({}, ValidationError),
        ({"filename": None, "sheetname": "sheet1"}, ValidationError),
        ({"filename": "blargh", "sheetname": None}, ValidationError),
        ({"filename": 1, "sheetname": "sheet1"}, ValidationError),
        ({"filename": "blargh", "sheetname": 1}, ValidationError),
        ({"filename": "blargh", "sheetname": "sheet1", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `XLSXFiles`",
        "Missing 'filename'",
        "Missing 'sheetname'",
        "Empty dictionary",
        "Filename key missing value",
        "Sheetname key missing value",
        "Filename value not a string",
        "Sheetname value not a string",
        "Excess keys (valid)",
    ],
)
def test_XLSXFiles(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert XLSXFiles(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert XLSXFiles(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(XLSXFiles(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(XLSXFiles(**input).model_dump(exclude_none=True))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": "km"},
            None,
        ),
        ({"TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": "km"}, ValidationError),
        ({"number_code": 12345, "TS_L_intercept": 2.0, "length_units": "km"}, ValidationError),
        ({"number_code": 12345, "TS_L_slope": 1.0, "length_units": "km"}, ValidationError),
        ({"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0}, ValidationError),
        ({}, ValidationError),
        (
            {"number_code": None, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": "km"},
            ValidationError,
        ),
        (
            {"number_code": 12345, "TS_L_slope": None, "TS_L_intercept": 2.0, "length_units": "km"},
            ValidationError,
        ),
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": None, "length_units": "km"},
            ValidationError,
        ),
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": None},
            ValidationError,
        ),
        (
            {
                "number_code": "12345",
                "TS_L_slope": "1.0",
                "TS_L_intercept": "2.0",
                "length_units": "km",
            },
            None,
        ),
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": 1},
            ValidationError,
        ),
        (
            {"number_code": 12345.0, "TS_L_slope": 1, "TS_L_intercept": 2, "length_units": "km"},
            None,
        ),
        (
            {
                "number_code": 12345.0,
                "TS_L_slope": 1,
                "TS_L_intercept": 2,
                "length_units": "km",
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `TSLRegressionParameters`",
        "Missing 'number_code'",
        "Missing 'TS_L_slope'",
        "Missing 'TS_L_intercept'",
        "Missing 'length_units'",
        "Empty dictionary",
        "Number_code key missing value",
        "TS_L_slope key missing value",
        "TS_L_intercept key missing value",
        "Length_units key missing value",
        "All values as strings (valid)",
        "Length_units as numeric (invalid)",
        "Number_code as float, TS_L_slope/TS_L_intercept as integers (valid)",
        "Excess keys (valid)",
    ],
)
def test_TSLRegressionParameters(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert TSLRegressionParameters(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert TSLRegressionParameters(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(TSLRegressionParameters(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (
                    set(input) - set(TSLRegressionParameters(**input).model_dump(exclude_none=True))
                ) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"init": "epsg:1000"}, None),
        ({"init": "EPSG:1000"}, None),
        ({"init": "epsg1000"}, None),
        ({"init": "EPSG1000"}, None),
        ({"init": "1000"}, None),
        ({"init": 1000}, None),
        ({"init": 1000.5}, ValidationError),
        ({"init": "ABCD:1000"}, ValidationError),
        ({}, ValidationError),
        ({"init": None}, AttributeError),
        ({"init": "epsg:1000", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `Geospatial`",
        "Uppercase EPSG (valid)",
        "Missing ':' between 'epsg' and number code (valid)",
        "Missing ':' between 'EPSG' and number code (valid)",
        "Missing 'epsg' from number code (valid)",
        "Number code integer input (valid)",
        "Number code float input (invalid)",
        "Unknown coordinate system type (non-EPSG)",
        "Empty dictionary",
        "Init key missing value",
        "Excess keys (valid)",
    ],
)
def test_Geospatial(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert Geospatial(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert Geospatial(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(Geospatial(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(Geospatial(**input).model_dump(exclude_none=True))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            None,
        ),
        ({"num_replicates": 10, "mesh_transects_per_latitude": 1}, ValidationError),
        ({"strata_transect_proportion": 1.0, "mesh_transects_per_latitude": 1}, ValidationError),
        ({"strata_transect_proportion": 1.0, "num_replicates": 10}, ValidationError),
        ({}, ValidationError),
        (
            {
                "strata_transect_proportion": None,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": None,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": 10,
                "mesh_transects_per_latitude": None,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            None,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10.0,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": "1",
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": "10",
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10,
                "mesh_transects_per_latitude": "1",
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `StratifiedSurveyMeanParameters`",
        "Missing 'strata_transect_proportion'",
        "Missing 'num_replicates'",
        "Missing 'mesh_transects_per_latitude'",
        "Empty dictionary",
        "Strata_transect_proportion key missing value",
        "Num_replicates key missing value",
        "Mesh_transects_per_latitude key missing value",
        "Strata_transect_proportion as int (valid)",
        "Num_replicates as float (invalid)",
        "Mesh_transects_per_latitude as float (invalid)",
        "Strata_transect_proportion as str (invalid)",
        "Num_replicates as str (invalid)",
        "Mesh_transects_per_latitude as str (invalid)",
        "Excess keys (valid)",
    ],
)
def test_StratifiedSurveyMeanParameters(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert StratifiedSurveyMeanParameters(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert StratifiedSurveyMeanParameters(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(
                StratifiedSurveyMeanParameters(**input).model_dump(exclude_none=True)
            ) == set(input)
        except AssertionError:
            try:
                assert (
                    set(input)
                    - set(StratifiedSurveyMeanParameters(**input).model_dump(exclude_none=True))
                ) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")
