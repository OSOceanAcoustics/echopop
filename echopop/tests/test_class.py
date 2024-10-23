import re

import numpy as np
import pytest

from echopop.utils.validate import posfloat, posint, realcircle, realposfloat


def test_posint():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert posint(int(0)) == 0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(float(0.0))

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert posint(int(5)) == 5

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(float(5.1))

    # -------------------------
    # [ TEST 5 ]: ASSERT: str("test") FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(str("test"))

    # -------------------------
    # [ TEST 6 ]: ASSERT: numpy.inf FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(np.inf)

    # -------------------------
    # [ TEST 7 ]: ASSERT: numpy.nan FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(np.nan)

    # -------------------------
    # [ TEST 8 ]: ASSERT: int(-1) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(int(-1))

    # -------------------------
    # [ TEST 9 ]: ASSERT: float(-1.0) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative integer"):
        assert posint(float(-1.0))


def test_posfloat():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert posfloat(int(0)) == 0.0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0)
    assert posfloat(float(0.0)) == 0.0

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert posfloat(int(5)) == 5.0

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1))
    assert posfloat(float(5.1)) == 5.1

    # -------------------------
    # [ TEST 5 ]: ASSERT: str("test") FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(str("test"))

    # -------------------------
    # [ TEST 6 ]: ASSERT: np.inf
    assert posfloat(np.inf) == np.inf

    # -------------------------
    # [ TEST 7 ]: ASSERT: np.inf
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-np.inf)

    # -------------------------
    # [ TEST 8 ]: ASSERT: np.nan
    assert np.isnan(posfloat(np.nan))

    # -------------------------
    # [ TEST 9 ]: ASSERT: int(-2) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-2)

    # -------------------------
    # [ TEST 10 ]: ASSERT: int(-2.5) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-2.5)


def test_realposfloat():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert realposfloat(int(0)) == 0.0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0)
    assert realposfloat(float(0.0)) == 0.0

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert realposfloat(int(5)) == 5.0

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1))
    assert realposfloat(float(5.1)) == 5.1

    # -------------------------
    # [ TEST 5 ]: ASSERT: str("test") FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative real number."):
        assert realposfloat(str("test"))

    # -------------------------
    # [ TEST 6 ]: ASSERT: np.inf FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative real number."):
        assert realposfloat(np.inf)

    # -------------------------
    # [ TEST 7 ]: ASSERT: -np.inf FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative real number."):
        assert realposfloat(-np.inf)

    # -------------------------
    # [ TEST 8 ]: ASSERT: np.nan
    assert np.isnan(realposfloat(np.nan))

    # -------------------------
    # [ TEST 9 ]: ASSERT: int(-2) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert realposfloat(-2)

    # -------------------------
    # [ TEST 10 ]: ASSERT: float(-2.5) FAIL
    with pytest.raises(ValueError, match="Value must be a non-negative float."):
        assert posfloat(-2.5)


def test_realcircle():

    # -------------------------
    # [ TEST 1 ]: ASSERT: int(0)
    assert realcircle(int(0)) == 0.0

    # -------------------------
    # [ TEST 2 ]: ASSERT: float(0.0)
    assert realcircle(float(0.0)) == 0.0

    # -------------------------
    # [ TEST 3 ]: ASSERT: int(5)
    assert realcircle(int(5)) == 5.0

    # -------------------------
    # [ TEST 4 ]: ASSERT: float(5.1))
    assert realcircle(float(5.1)) == 5.1

    # -------------------------
    # [ TEST 5 ]: ASSERT: int(360)
    assert realcircle(int(360)) == 360

    # -------------------------
    # [ TEST 6 ]: ASSERT: float(360.0))
    assert realcircle(float(360.0)) == 360.0

    # -------------------------
    # [ TEST 7 ]: ASSERT: int(361) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(int(361))

    # -------------------------
    # [ TEST 8 ]: ASSERT: float(361.1)) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(float(361.1))

    # -------------------------
    # [ TEST 9 ]: ASSERT: str("test") FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(str("test"))

    # -------------------------
    # [ TEST 10 ]: ASSERT: np.inf FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(np.inf)

    # -------------------------
    # [ TEST 11 ]: ASSERT: -np.inf FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(-np.inf)

    # -------------------------
    # [ TEST 12 ]: ASSERT: np.nan
    assert np.isnan(realcircle(np.nan))

    # -------------------------
    # [ TEST 13 ]: ASSERT: int(-2) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(int(-2))

    # -------------------------
    # [ TEST 14 ]: ASSERT: float(-2.5) FAIL
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value must be a non-negative real angle (as a 'float') between 0.0 and 360.0 degrees."
        ),
    ):
        assert realcircle(float(-2.5))
