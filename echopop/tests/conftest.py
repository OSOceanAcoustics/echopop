from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import pytest
from _pytest.assertion.util import assertrepr_compare

from echopop import Survey

# Set up path to the `test_data` folder
HERE = Path(__file__).parent.absolute()
TEST_DATA_ROOT = HERE.parent / "test_data"


# Fixtures
# ---- Test root/config/input file paths
@pytest.fixture(scope="session")
def test_path():

    return {
        "ROOT": TEST_DATA_ROOT,
        "CONFIG": TEST_DATA_ROOT / "config_files",
        "INPUT": TEST_DATA_ROOT / "input_files",
    }


# ---- Mock `Survey` class object
@pytest.fixture(scope="session")
def mock_survey(test_path) -> Survey:

    return Survey(
        init_config_path=Path(test_path["CONFIG"] / "config_init.yml"),
        survey_year_config_path=Path(test_path["CONFIG"] / "config_survey.yml"),
    )


# Hook functions
def pytest_assertrepr_compare(config, op, left, right):
    """
    Hook function that always shows the full `diff` on assertion
    failures by increasing the verbosity (`config.option.verbose`)
    """

    # Adjust configuration `diff` verbosity
    config.option.verbose = 2

    return assertrepr_compare(config, op, left, right)


# Utility functions
# ---- DICTIONARY
# ++++ Shape and structure
def dictionary_shape(dictionary: dict):
    """
    A utility test function that extracts the shape of a nested dictionary
    """

    if isinstance(dictionary, dict):
        return {i: dictionary_shape(dictionary[i]) for i in dictionary}
    else:
        return None


# ---- DATAFRAME
# ++++ Shape
def dataframe_shape(input: Union[pd.DataFrame, dict]):

    # DataFrame
    if isinstance(input, pd.DataFrame):

        return input.shape

    # Dictionary (bundled dataframes)
    elif isinstance(input, dict):
        dataframe_shapes = {}

        for key, value in input.items():
            if isinstance(value, pd.DataFrame):
                dataframe_shapes[key] = value.shape
            elif isinstance(value, dict):
                dataframe_shapes[key] = dataframe_shape(value)

        return dataframe_shapes


# Assertion functions
# ---- DICTIONARY
# ---- Shape and dimensions
def assert_dictionary_structure_equal(dictionary1: dict, dictionary2: dict):
    """
    Tests equality between the shapes of two nested dictionaries
    """

    result = dictionary_shape(dictionary1) == dictionary_shape(dictionary2)

    if result:
        assert result
    else:
        if set(dictionary_shape(dictionary1)) <= set(dictionary_shape(dictionary2)):
            tracked_true = []

            for j in dictionary2.keys():
                test = set(dictionary1[j].keys()) <= (dictionary2[j].keys())
                tracked_true.append(test)

            if np.all(tracked_true):
                assert True
            else:
                assert result
        else:
            assert result


# ---- dtypes
def assert_dictionary_dtypes_equal(dictionary, reference_dictionary):

    for key in reference_dictionary:
        if isinstance(reference_dictionary[key], dict):
            assert isinstance(
                dictionary[key], dict
            ), f"Key '{key}' has different types in the dictionaries."
            assert_dictionary_dtypes_equal(dictionary[key], reference_dictionary[key])
        elif isinstance(dictionary[key], type):
            assert np.issubdtype(
                type(dictionary[key]), reference_dictionary[key]
            ), f"Datatype for key '{key}' is not a subdtype of the reference datatype."
        elif isinstance(reference_dictionary[key], np.ndarray):
            assert isinstance(
                dictionary[key], np.ndarray
            ), f"Datatype for key '{key}' is not the same as in reference dictionary."
            assert np.issubdtype(
                dictionary[key].dtype, reference_dictionary[key].dtype
            ), f"Dtype for key '{key}' is not a subdtype of the reference dtype."


# ---- Values
def assert_dictionary_values_equal(dictionary, reference_dictionary):
    for key in dictionary:
        if isinstance(dictionary[key], dict):
            assert isinstance(
                reference_dictionary[key], dict
            ), f"Key '{key}' has different types in the dictionaries."
            assert_dictionary_values_equal(dictionary[key], reference_dictionary[key])
        elif isinstance(dictionary[key], np.ndarray):
            assert np.allclose(
                dictionary[key], reference_dictionary[key]
            ), f"Arrays for key '{key}' are not close."
        else:
            assert np.isclose(
                dictionary[key], reference_dictionary[key]
            ), f"Values for key '{key}' are not close."


# ---- DATAFRAME
# ---- Shape and dimensions
def assert_dataframe_shape_equal(input: Union[pd.DataFrame, dict], reference: Union[tuple, dict]):

    # DataFrame
    if (isinstance(input, pd.DataFrame)) & (isinstance(reference, tuple)):
        assert input.shape == reference

    # Dictionary
    elif (isinstance(input, dict)) & (isinstance(reference, dict)):
        assert dataframe_shape(input) == dataframe_shape(reference)


# ---- dtypes
# ~~~~ !!!! ATTN: this is a nested function within `assert_dataframe_dtypes_equal`!
def _assert_dataframe_dtypes_equal(dataframe: pd.DataFrame, reference_dictionary: dict):

    # Separate evaluation for categorical-type
    # ---- Parse expected categorical variables
    categorical_columns = [
        k for k, v in reference_dictionary.items() if isinstance(v, pd.CategoricalDtype)
    ]

    # ---- Assert that all categorical columns in the reference dictionary match the categorical
    # ----- columns in the tested dataframe
    assert np.all(dataframe.select_dtypes(include=["category"]).columns.isin(categorical_columns))

    # ---- Remove categorical columns from the dataframe
    dataframe = dataframe.copy().drop(categorical_columns, axis=1)

    # Loop through columns to assert that dtypes from the tested dataframe
    # match those expected in a reference dictionary
    for column, dtype in dataframe.dtypes.items():
        assert np.issubdtype(
            dtype, reference_dictionary.get(column, object)
        ), f"Data type mismatch for column '{column}'"


# ~~~~ dtypes --> compatible with direct DataFrame or bundled DataFrames within a dictionary
def assert_dataframe_dtypes_equal(input: Union[pd.DataFrame, dict], reference: dict):

    # DataFrame
    if isinstance(input, pd.DataFrame):
        _assert_dataframe_dtypes_equal(input, reference)

    # Dictionary
    elif isinstance(input, dict):
        for category, data in reference.items():

            # ---- Single Dictionary layer
            if isinstance(input[category], pd.DataFrame):
                _assert_dataframe_dtypes_equal(input[category], reference[category])

            # ---- Nested Dictionary layers
            else:
                for df_name, _ in data.items():
                    _assert_dataframe_dtypes_equal(
                        input[category][df_name], reference[category][df_name]
                    )


# ---- Values
# ~~~~ !!!! ATTN: this is a nested function within `assert_dataframe_equal`!
def _aassert_dataframe_values_equal(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame):

    # Evaluate equality between numerical values
    assert np.allclose(
        dataframe1.select_dtypes(include=["number"]),
        dataframe2.select_dtypes(include=["number"]),
        equal_nan=True,
    )

    # Evaluate equality between non-numerical values
    # ---- Mask out "NaN"
    dataframe1_nan_mask = dataframe1.isna().any(axis=1)
    dataframe2_nan_mask = dataframe2.isna().any(axis=1)
    # ---- Evaluate equality
    dataframe1_nan_mask == dataframe2_nan_mask
    # ---- Evaluate equality among "real" values
    dataframe1_masked = dataframe1[~dataframe1_nan_mask]
    dataframe2_masked = dataframe2[~dataframe2_nan_mask]
    assert np.all(
        dataframe1_masked.select_dtypes(exclude=["number"])
        == dataframe2_masked.select_dtypes(exclude=["number"])
    )


# ~~~~ Values --> compatible with direct DataFrame or bundled DataFrames within a dictionary
def assert_dataframe_values_equal(
    input: Union[pd.DataFrame, dict], reference: Union[pd.DataFrame, dict]
):

    # Direct DataFrame
    if isinstance(input, pd.DataFrame) & (isinstance(reference, pd.DataFrame)):
        _aassert_dataframe_values_equal(input, reference)

    # Iterate through nested DataFrames within each dictionary
    else:
        for key, expected_df in reference.items():

            if isinstance(input[key], pd.DataFrame):
                _aassert_dataframe_values_equal(input[key], expected_df)

            else:
                for sub_key, _ in reference[key].items():
                    _aassert_dataframe_values_equal(input[key][sub_key], expected_df[sub_key])


# ++++ DICTIONARY + DATAFRAME BUNDLING
# ---> Dictionary
def assert_dictionary_equal(
    input: dict,
    reference_dtypes: dict,
    reference_values: dict,
):

    # Shape
    assert_dictionary_structure_equal(input, reference_values)
    # dtypes
    assert_dictionary_dtypes_equal(input, reference_dtypes)
    # Values
    assert_dictionary_values_equal(input, reference_values)


# ---> DataFrame
def assert_dataframe_equal(
    input: Union[pd.DataFrame, dict],
    reference_dtypes: dict,
    reference_values: Union[pd.DataFrame, dict],
):
    # Shape
    assert_dataframe_shape_equal(input, reference_values)
    # dtypes
    assert_dataframe_dtypes_equal(input, reference_dtypes)
    # Values
    assert_dataframe_values_equal(input, reference_values)
