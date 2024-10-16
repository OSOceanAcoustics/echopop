from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from ..utils.validate_dict import (
    CONFIG_DATA_MODEL,
    CONFIG_INIT_MODEL,
    BiologicalFiles,
    FileSettings,
    Geospatial,
    HaulTransectMap,
    KrigingFiles,
    KrigingParameters,
    NASCExports,
    PatternParts,
    StratificationFiles,
    StratifiedSurveyMeanParameters,
    TransectRegionMap,
    TSLRegressionParameters,
    XLSXFiles,
)


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"directory": "blurgh/blargh", "sheetname": "sheet1"}, None),
        ({"sheetname": "sheet1"}, ValidationError),
        ({"directory": "blurgh/blargh"}, ValidationError),
        ({}, ValidationError),
        ({"directory": None, "sheetname": "sheet1"}, ValidationError),
        ({"directory": "blurgh/blargh", "sheetname": None}, ValidationError),
        ({"directory": 1, "sheetname": "sheet1"}, ValidationError),
        ({"directory": "blurgh/blargh", "sheetname": 1}, ValidationError),
        ({"directory": "blurgh/blarg", "sheetname": "sheet1", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `FileSettings`",
        "Missing 'directory'",
        "Missing 'sheetname'",
        "Empty dictionary",
        "Directory key missing value",
        "Sheetname key missing value",
        "Directory value not a string",
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
        ({"filename": "blargh", "sheetname": ["sheet1", "sheet2"]}, None),
        ({"filename": ["blargh", "blurgh"], "sheetname": "sheet1"}, ValidationError),
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
        "List of sheets (valid)",
        "List of filenames (invalid)",
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


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "strata": {"filename": "blurgh", "sheetname": "sheet1"},
                "geo_strata": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            None,
        ),
        ({"strata": {"filename": "blurgh", "sheetname": "sheet1"}}, ValidationError),
        ({"geo_strata": {"filename": "blurgh", "sheetname": "sheet1"}}, ValidationError),
        ({}, ValidationError),
        (
            {
                "strata": {"filename": "blurgh", "sheetname": "sheet1"},
                "geo_strata": {"filename": "blurgh", "sheetname": "sheet1"},
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `StratificationFiles`",
        "Missing 'geo_strata'",
        "Missing 'strata'",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_StratificationFiles(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert StratificationFiles(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert StratificationFiles(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(StratificationFiles(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (
                    set(input) - set(StratificationFiles(**input).model_dump(exclude_none=True))
                ) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            None,
        ),
        (
            {
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            ValidationError,
        ),
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            ValidationError,
        ),
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            ValidationError,
        ),
        ({}, ValidationError),
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `KrigingFiles`",
        "Missing 'vario_krig_para'",
        "Missing 'isobath_200m'",
        "Missing 'mesh'",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_KrigingFiles(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert KrigingFiles(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert KrigingFiles(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(KrigingFiles(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(KrigingFiles(**input).model_dump(exclude_none=True))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
            },
            None,
        ),
        (
            {"longitude_reference": 0.0, "longitude_offset": -1.0, "latitude_offset": 1.0},
            ValidationError,
        ),
        ({"A0": 1.0, "longitude_offset": -1.0, "latitude_offset": 1.0}, ValidationError),
        ({"A0": 1.0, "longitude_reference": 0.0, "latitude_offset": 1.0}, ValidationError),
        ({"A0": 1.0, "longitude_reference": 0.0, "longitude_offset": -1.0}, ValidationError),
        (
            {
                "A0": None,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": None,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": None,
                "latitude_offset": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": None,
            },
            ValidationError,
        ),
        ({}, ValidationError),
        ({"A0": 1, "longitude_reference": 0, "longitude_offset": -1, "latitude_offset": 1}, None),
        (
            {"A0": -1.0, "longitude_reference": 0, "longitude_offset": -1, "latitude_offset": 1},
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `KrigingParameters`",
        "Missing 'A0'",
        "Missing 'longitude_reference'",
        "Missing 'longitude_offset'",
        "Missing 'latitude_offset'",
        "A0 key value missing",
        "Longitude_reference key value missing",
        "Longitude_offset key value missing",
        "Latitude_offset key value missing",
        "Empty dictionary",
        "All key values as integers (valid)",
        "Negative A0 (invalid)",
        "Excess keys (valid)",
    ],
)
def test_KrigingParameters(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert KrigingParameters(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert KrigingParameters(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(KrigingParameters(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (
                    set(input) - set(KrigingParameters(**input).model_dump(exclude_none=True))
                ) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            None,
        ),
        (
            {
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            ValidationError,
        ),
        (
            {"save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx", "country_code": ["A", "B"]},
            ValidationError,
        ),
        (
            {
                "save_file_template": None,
                "country_code": [["A", "B"]],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            TypeError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": None,
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            TypeError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {},
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            None,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}_{ERRONEOUS}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            (TypeError, ValueError),
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            ValueError,
        ),
        ({}, ValidationError),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `HaulTransectMap`",
        "Missing key 'save_file_template'",
        "Missing key 'country_code'",
        "Missing key 'file_settings'",
        "Save_file_template key value missing",
        "Country_code key value missing",
        "File_settings key value missing",
        "Single filename ID in template (valid)",
        "Erroneous filename ID in template (invalid)",
        "Mismatched 'file_setting' keys and items in 'country_code' list",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_HaulTransectMap(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert HaulTransectMap(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert HaulTransectMap(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(HaulTransectMap(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (
                    set(input) - set(HaulTransectMap(**input).model_dump(exclude_none=True))
                ) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            None,
        ),
        (
            {
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            None,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}_{ERRONEOUS}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        ({}, ValidationError),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `NASCExports`",
        "Missing key 'save_file_template'",
        "Missing key 'nasc_export_directory'",
        "Missing key 'export_file_directory'",
        "Missing key 'save_file_template'",
        "Missing key 'save_file_sheetname'",
        "Missing key 'regions'",
        "Missing key 'file_columns'",
        "Single filename ID in template (valid)",
        "Erroneous filename ID in template (invalid)",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_NASCExports(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert NASCExports(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert NASCExports(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(NASCExports(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(NASCExports(**input).model_dump(exclude_none=True))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",  #
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            None,
        ),
        (
            {
                "parts": {  #
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            ValidationError,
        ),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",  #
            },
            ValidationError,
        ),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",  #
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            None,
        ),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",  #
                "parts": {
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            ValidationError,
        ),
        (
            {
                "pattern": "{REGION_CLASS}",  #
                "parts": {"REGION_CLASS": [{"pattern": "a", "label": "A"}]},
            },
            None,
        ),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}{ERRONEOUS}",
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            ValidationError,
        ),
        (
            {
                "pattern": "{REGION_CLASS}",
                "parts": {"REGION_CLASS": {"pattern": "a", "label": "A"}},
            },
            ValidationError,
        ),
        ({}, ValidationError),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `TransectRegionMap`",
        "Missing key 'pattern'",
        "Missing key 'parts'",
        "Single filename ID in template (valid)",
        "Mismatch between 'parts' keys and those in pattern",
        "Single pattern ID in template (valid)",
        "Erroneous pattern ID in template (invalid)",
        "Parts components as a list and not a dictionary (invalid)",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_TransectRegionMap(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert TransectRegionMap(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert TransectRegionMap(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(TransectRegionMap(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (
                    set(input) - set(TransectRegionMap(**input).model_dump(exclude_none=True))
                ) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            ValidationError,
        ),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            ValidationError,
        ),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            ValidationError,
        ),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
            },
            ValidationError,
        ),
        (
            {
                "length": {"A": {"filename": "blurgh", "sheetname": "sheet1"}},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "length": {
                    "A": {"filename": "blurgh", "sheetname": "sheet1"},
                    "B": {"filename": "blaurrgh", "sheetname": "sheet5"},
                },
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "length": {"A": {"filename": "blurgh", "sheetname": "sheet1"}},
                "specimen": {"A": {"filename": "blargh", "sheetname": "sheet2"}},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "length": {
                    "A": {"filename": "blurgh", "sheetname": "sheet1"},
                    "B": {"filename": "blaurrgh", "sheetname": "sheet5"},
                },
                "specimen": {
                    "A": {"filename": "blargh", "sheetname": "sheet2"},
                    "B": {"filename": "bleurgh", "sheetname": "sheet6"},
                },
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        ({}, ValidationError),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `KrigingFiles` (unnested dictionary)",
        "Missing key 'length'",
        "Missing 'specimen'",
        "Missing 'catch'",
        "Missing 'haul_to_transect'",
        "Single nested directory with one key (valid)",
        "Single nested directory with two keys (valid)",
        "Multiple nested directories with one key (valid)",
        "Multiple nested directories with two keys (valid)",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_BiologicalFiles(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert BiologicalFiles(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert BiologicalFiles(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(BiologicalFiles(**input).model_dump(exclude_none=True)) == set(input)
        except AssertionError:
            try:
                assert (
                    set(input) - set(BiologicalFiles(**input).model_dump(exclude_none=True))
                ) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "description",
    ["Test `CONFIG_INIT_MODEL` and `CONFIG_DATA_MODEL` validation"],
    ids=["Test `CONFIG_INIT_MODEL` and `CONFIG_DATA_MODEL` validation"],
)
def test_config_pydantic(test_path, description):

    # -------------------------
    # Read in the initialization and file configuration
    # ---- Initialization
    init_config = yaml.safe_load(Path(test_path["CONFIG"] / "config_init.yml").read_text())
    # ---- Files
    files_config = yaml.safe_load(Path(test_path["CONFIG"] / "config_survey.yml").read_text())

    # -------------------------
    # [ TEST 1 ]: FIRST CHECK OF PRE-VALIDATED DATA - `CONFIG_INIT_MODEL`
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 2 ]: ADD 'transect_region_mapping' KEY
    # ---- Add key values
    init_config.update(
        {
            "transect_region_mapping": {
                "save_file_template": "blurgh_{COUNTRY}_{YEAR}_{GROUP}.xlsx",
                "save_file_directory": "blurgh/blargh",
                "pattern": "{REGION_CLASS}",
                "save_file_sheetname": "sheet1",
                "parts": {"REGION_CLASS": [{"pattern": "a", "label": "A"}]},
            }
        }
    )
    # ---- ASSERT
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 3 ]: ADD 'nasc_exports' KEY
    # ---- Add key values
    init_config.update(
        {
            "nasc_exports": {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            }
        }
    )
    # ---- ASSERT
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 4 ]: ADD 'haul_to_transect_mapping' KEY
    # ---- Add key values
    init_config.update(
        {
            "haul_to_transect_mapping": {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            }
        }
    )
    # ---- ASSERT
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 5 ]: FIRST CHECK OF PRE-VALIDATED DATA - `CONFIG_DATA_MODEL`
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 6 ]: ADD 'ship_id' [int] KEY
    # ---- Add key values
    files_config.update({"ship_id": 1})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 7 ]: ADD 'ship_id' [float] KEY
    # ---- Add key values
    files_config.update({"ship_id": 1.1})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 8 ]: ADD 'ship_id' [str] KEY
    # ---- Add key values
    files_config.update({"ship_id": "1"})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 9 ]: ADD 'export_regions' KEY
    # ---- Add key values
    files_config.update({"export_regions": {"A": {"filename": "blee", "sheetname": "blah"}}})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)
