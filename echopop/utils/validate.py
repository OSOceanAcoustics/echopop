"""
Validation functions.
"""

# TODO: Compile all package validators here since they may not belong elsewhere
import numpy as np


# CLASS-SPECIFIC CORE API
class posint(int):
    """Positive-only integer (includes 0)"""

    __failstate__ = "must be a non-negative integer"
    __origin__ = "posint"

    def __new__(cls, value):
        if not isinstance(value, int) or value < 0:
            raise ValueError("Value must be a non-negative integer.")
        return super().__new__(cls, value)


class posfloat(float):
    """Positive-only float (includes 0.0)"""

    __failstate__ = "must be a non-negative float"
    __origin__ = "posfloat"

    def __new__(cls, value):
        if not isinstance(value, (float, int)) or value < 0:
            raise ValueError("Value must be a non-negative float.")
        return super().__new__(cls, value)


class realposfloat(posfloat):
    """Real number positive-only float (includes 0.0)"""

    __failstate__ = "must be a non-negative real number"
    __origin__ = "realposfloat"

    def __new__(cls, value):
        if not isinstance(value, (float, int)) or np.isinf(value):  # Check if value is infinity
            raise ValueError(f"Value {cls.__failstate__}.")
        return super().__new__(cls, value)


class realcircle(realposfloat):
    """Real number in a unit circle"""

    __failstate__ = "must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees"
    __origin__ = "realcircle"

    def __new__(cls, value):
        if not isinstance(value, (float, int)) or (value < 0.0 or value > 360.0):
            raise ValueError(f"Value {cls.__failstate__}.")
        return super().__new__(cls, value)
