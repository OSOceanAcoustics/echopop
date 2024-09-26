import json
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import pytest
from _pytest.assertion.util import assertrepr_compare

# Set up path to the `test_data` folder
from echopop import Survey

# Set up path to the `test_data` folder
HERE = Path(__file__).parent.absolute()
TEST_DATA_ROOT = HERE.parent / "test_data"


# FIXTURES
# ---- Test root/config/input file paths
@pytest.fixture(scope="session")
def test_path():

    return {
        "ROOT": TEST_DATA_ROOT,
        "CONFIG": TEST_DATA_ROOT / "config_files",
        "INPUT": TEST_DATA_ROOT / "input_files",
        "EXPECTED": TEST_DATA_ROOT / "expected_outputs",
    }


# ---- Mock `Survey` class object
@pytest.fixture(scope="session")
def mock_survey(test_path) -> Survey:

    return Survey(
        init_config_path=Path(test_path["CONFIG"] / "config_init.yml"),
        survey_year_config_path=Path(test_path["CONFIG"] / "config_survey.yml"),
    )


# HOOK FUNCTIONS
def pytest_assertrepr_compare(config, op, left, right):
    """
    Hook function that always shows the full `diff` on assertion
    failures by increasing the verbosity (`config.option.verbose`)
    """

    # Adjust configuration `diff` verbosity
    # Adjust configuration `diff` verbosity
    config.option.verbose = 2

    return assertrepr_compare(config, op, left, right)


# UTILITY FUNCTIONS
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
        elif isinstance(dictionary[key], (int, float)):
            assert np.isclose(
                dictionary[key], reference_dictionary[key]
            ), f"Values for key '{key}' are not close."
        elif dictionary[key] == reference_dictionary[key]:
            continue
        else:
            raise AssertionError(
                f"Values for key '{key}' are not the same. Got: {dictionary[key]}, expected: "
                f"{reference_dictionary[key]}."
            )


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

                for df_name, _ in data.items():
                    _assert_dataframe_dtypes_equal(
                        input[category][df_name], reference[category][df_name]
                    )


# ---- Values
# ~~~~ !!!! ATTN: this is a nested function within `assert_dataframe_equal`!
def _assert_dataframe_values_equal(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame):

    # Evaluate equality between numerical values
    assert np.allclose(
        dataframe1.select_dtypes(include=["number"]),
        dataframe2.select_dtypes(include=["number"]),
        equal_nan=True,
    )

    # Evaluate equality between non-numerical values
    # ---- Mask out "NaN"
    dataframe1_nan_mask = dataframe1.isna().any(axis=1).reset_index(drop=True)
    dataframe2_nan_mask = dataframe2.isna().any(axis=1).reset_index(drop=True)
    # ---- Evaluate equality
    assert np.all(dataframe1_nan_mask == dataframe2_nan_mask)
    # ---- Evaluate equality among "real" values
    dataframe1_masked = dataframe1.loc[~dataframe1_nan_mask].reset_index(drop=True)
    dataframe2_masked = dataframe2.loc[~dataframe2_nan_mask].reset_index(drop=True)
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
        _assert_dataframe_values_equal(input, reference)

    # Iterate through nested DataFrames within each dictionary
    else:
        for key, expected_df in reference.items():

            if isinstance(input[key], pd.DataFrame):
                _assert_dataframe_values_equal(input[key], expected_df)

            else:
                for sub_key, _ in reference[key].items():
                    _assert_dataframe_values_equal(input[key][sub_key], expected_df[sub_key])


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


# Utility functions (JSON)
def update_json_file(filename: str, test_name: str, new_entries: dict):
    """
    Updates an existing JSON file with new entries, overwriting existing entries with the same key.

    Parameters
    ----------
    filename: str
        The name of the JSON file.
    test_name: str
        The name of the associated `pytest`.
    new_entries: dict
        A dictionary containing the new entries to add or update.

    Returns
    ----------
    None

    Notes
    ----------
    This is primarily used for updating associated .JSON files that store values for expected
    test results.
    """

    # Determine whether the associated test-key already exists
    try:
        with open(filename, "r") as f:
            existing_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        existing_data = {}

    # Helper function for converting Python into JSON data types
    def convert_to_json(data):
        if isinstance(data, dict):
            return {k: convert_to_json(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [convert_to_json(item) for item in data]
        elif isinstance(data, np.ndarray):
            # Handle NaN values by converting them to None for JSON serialization
            return data.tolist()
        elif isinstance(data, pd.DataFrame):
            return data.to_dict(orient="records")
        else:
            return data

    # Convert the data
    converted_entries = convert_to_json(new_entries)

    # Helper function for gathering data types
    def gather_metadata(data, is_numpy, is_pandas, prefix=""):
        if isinstance(data, dict):
            for k, v in data.items():
                gather_metadata(v, is_numpy, is_pandas, f"{prefix}{k}.")
        elif isinstance(data, list):
            for item in data:
                gather_metadata(item, is_numpy, is_pandas, prefix)
        elif isinstance(data, np.ndarray):
            is_numpy.append(prefix.rstrip("."))
        elif isinstance(data, pd.DataFrame):
            is_pandas.append(prefix.rstrip("."))

    # Gather metadata to determine whether datatypes are numpy.ndarray or pandas.DataFrame
    # ---- Initialize lists
    is_numpy = []
    is_pandas = []
    # ---- Populate the lists
    gather_metadata(new_entries, is_numpy, is_pandas)

    # Among new entries, convert any arrays to a JSON-friendly format
    # ---- Add a metadata key `is_numpy`: a list of keys
    converted_entries["is_numpy"] = list(set([key.split(".")[-1] for key in is_numpy]))
    # ---- Add a metadata key `is_pandas`: a list of keys
    converted_entries["is_pandas"] = list(set([key.split(".")[-1] for key in is_pandas]))

    # Update the test-key
    existing_data[test_name] = converted_entries

    # Write/update the JSON file
    with open(filename, "w") as f:
        json.dump(existing_data, f, indent=4)


def load_json_data(filename: str, default_value=None, test_name: Optional[str] = None):
    """
    Loads JSON data from a file. If the file doesn't exist or an error occurs, returns the default
    value.

    Parameters
    ----------
    filename: str
        The name of the JSON file.
    default_value: dict
        The default value to return if loading fails.
    test_name: Optional[str]
        The name of the associated `pytest`.

    Returns
    ----------
    The loaded JSON data or the default value.
    """

    # Extract the expected results
    try:
        with open(filename, "r") as f:
            data = json.load(f)
            if test_name is not None:
                # ---- Get expected data
                expected_data = data.get(test_name, default_value)
            else:
                expected_data = data
            # ---- Get the key names that are numpy arrays
            is_numpy = expected_data.get("is_numpy", [])
            # ---- Get the key names that are pandas DataFrames
            is_pandas = expected_data.get("is_pandas", [])

            # Helper function for converting JSON to Python data types
            def convert_from_json(data, is_numpy, is_pandas):
                if isinstance(data, dict):
                    for k, v in data.items():
                        # Check if the current key is in `is_numpy` or `is_pandas`
                        if k in is_numpy:
                            data[k] = np.array(v)
                        elif k in is_pandas:
                            data[k] = pd.DataFrame(v)
                        else:
                            # Recursively process nested dictionaries or lists
                            data[k] = convert_from_json(v, is_numpy, is_pandas)
                elif isinstance(data, list):
                    for i in range(len(data)):
                        data[i] = convert_from_json(data[i], is_numpy, is_pandas)
                return data

            # Return the converted data
            return convert_from_json(expected_data, is_numpy, is_pandas)

    except (FileNotFoundError, json.JSONDecodeError):
        return default_value
