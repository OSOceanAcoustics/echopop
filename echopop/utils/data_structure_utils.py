from functools import reduce


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
