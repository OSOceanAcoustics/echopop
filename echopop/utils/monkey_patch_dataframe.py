from typing import Callable

import pandas as pd


def patch_method_to_DataFrame(cls: Callable = pd.DataFrame) -> Callable[[Callable], Callable]:
    """
    A monkey patch method for adding attributes to the pandas.DataFrame object

    Parameters
    ----------
    cls : Callable
        The pandas.DataFrame object

    Notes
    -----
    This is primarily a support function that modifies the pandas.DataFrame object to
    extend usage for various utility functions used throughout the calculations embedded
    within the echopop.survey object. This only modifies the pandas.DataFrame class
    definitions within the instance and therefore does not affect the actual installation.
    This allows for modifying the class during runtime when only necessary and enables a
    '@patch_method_to_DataFrame' decorator to be associated with defined methods that are
    added to the pandas.DataFrame attributes.
    """

    def decorator(function) -> Callable:

        setattr(cls, function.__name__, function)

        return function

    return decorator
