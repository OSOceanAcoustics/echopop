from functools import reduce
from typing import Any, List

import pandas as pd


def pull_nested_dict(nested_dictionary, variable_tree_list):
    """
    Utility function that pulls target values nested within ``nested_dictionary``
    using variable/data layer names provided via ``variable_tree_list`` as a
    list. The primary purpose of this function is to specifically pull dataframes
    to evaluate whether or not they exist, should be concatenated with a proposal
    dataframe, or overwritten by a proposal dataframe.
    """
    return reduce(dict.get, variable_tree_list, nested_dictionary)


def push_nested_dict(nested_dictionary, variable_tree_list, data):
    """
    Utility function that pushes a proposal dataframe, ``data``, to a specific
    nested directory within the ``Survey`` class object via variable/data layer names
    provided within ``variable_tree_list`` as a list. This iterately generates a new
    dictionary tree/branch that appends itself to already existing branches, or from
    scratch.
    """

    for key in variable_tree_list[:-1]:
        nested_dictionary = nested_dictionary.setdefault(key, {})

    nested_dictionary[variable_tree_list[-1]] = data


def map_imported_datasets(dictionary: dict) -> List[str]:
    """
    Utility function for mapping datasets that have been successfully imported.
    """

    # Helper function for detecting whether a dictionary key is filled/empty
    def is_present(value: Any) -> bool:
        if isinstance(value, pd.DataFrame):
            return not value.empty
        elif isinstance(value, dict):
            # ---- Check if any value in the dictionary is non-empty
            return any(is_present(v) for v in value.values())
        elif isinstance(value, list):
            # ---- Non-empty if there are elements
            return bool(value)
        else:
            # ---- Non-empty if not false
            return bool(value)

    # Apply the check to the top-level dictionary
    keys_with_data = filter(lambda item: is_present(item[1]), dictionary.items())

    # Extract the keys and return the subsequent list
    return list(map(lambda item: item[0], keys_with_data))
