import numpy as np

from echopop.spatial.variogram import (
    VARIOGRAM_MODELS,
    bessel_exponential,
    bessel_gaussian,
    cosine_exponential,
    cosine_gaussian,
    exponential,
    exponential_linear,
    gaussian,
    gaussian_linear,
    get_variogram_arguments,
    jbessel,
    kbessel,
    linear,
    nugget,
    sinc,
    spherical,
    variogram,
)


def test_exponential():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.25

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        exponential(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
        np.array([0.00000000, 0.55067104, 0.79810348, 0.90928205, 0.9592378, 0.98168436]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        exponential(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
        np.array([0.20000000, 0.64053683, 0.83848279, 0.92742564, 0.96739024, 0.98534749]),
    )


def test_gaussian():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        gaussian(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
        np.array([0.00000000, 0.63212056, 0.98168436, 0.99987659, 0.99999989, 1.00000000]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        gaussian(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
        np.array([0.20000000, 0.70569645, 0.98534749, 0.99990127, 0.99999991, 1.00000000]),
    )


def test_jbessel():

    # -------------------------
    # Mock values [ TEST 1 ]: NO HOLE EFFECT + NO BASE NUGGET EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.array_equal(
        jbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO HOLE EFFECT + BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.array_equal(
        jbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2]),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: HOLE EFFECT + NO BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.50

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        jbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.0000000, 0.00249844, 0.00997503, 0.02237375, 0.03960177, 0.06153019]),
    )

    # -------------------------
    # Mock values [ TEST 4 ]: HOLE EFFECT + BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        jbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.1000000, 0.10224859, 0.10897752, 0.12013638, 0.1356416, 0.15537717]),
    )


def test_kbessel():

    # -------------------------
    # Mock values [ TEST 1 ]: NO HOLE EFFECT + NO BASE NUGGET EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.array_equal(
        kbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO HOLE EFFECT + BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.array_equal(
        kbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.0, 0.2, 0.2, 0.2, 0.2, 0.2]),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: HOLE EFFECT + NO BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.50

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        kbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.0000000, 0.12625823, 0.31057469, 0.47848913, 0.61498574, 0.72026824]),
    )

    # -------------------------
    # Mock values [ TEST 4 ]: HOLE EFFECT + BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        kbessel(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.1000000, 0.21363241, 0.37951722, 0.53064022, 0.65348717, 0.74824141]),
    )


def test_linear():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        linear(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET), np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    )

    # -------------------------
    # Mock values [ TEST 2 ]: BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        linear(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET), np.array([0.20, 0.36, 0.52, 0.68, 0.84, 1.00])
    )


def test_nugget():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        nugget(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET), np.array([0.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    )

    # -------------------------
    # Mock values [ TEST 2 ]: BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        nugget(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET), np.array([0.0, 1.2, 1.2, 1.2, 1.2, 1.2])
    )


def test_sinc():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT RANGE
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.array_equal(
        sinc(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: BASE NUGGET EFFECT + NO HOLE EFFECT RANGE
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.array_equal(
        sinc(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: BASE NUGGET EFFECT + HOLE EFFECT RANGE
    # ---- Nugget
    MOCK_HOLE_EFFECT_RANGE = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        sinc(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_HOLE_EFFECT_RANGE),
        np.array([0.84000000, 0.84004266, 0.84017061, 0.84038372, 0.84068179, 0.84106454]),
    )


def test_spherical():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        spherical(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
        np.array([0.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.20

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        spherical(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
        np.array([0.2, 1.2, 1.2, 1.2, 1.2, 1.2]),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: BASE NUGGET EFFECT + INCREASE CORRELATION RANGE
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.50

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        spherical(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
        np.array([0.2000, 0.6544, 0.9552, 1.2000, 1.2000, 1.2000]),
    )


def test_bessel_exponential():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT + NO DECAY POWER
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.2
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00
    # ---- Decay power
    MOCK_DECAY_POWER = 0.0

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        bessel_exponential(
            MOCK_LAGS,
            MOCK_NUGGET,
            MOCK_SILL,
            MOCK_CORRELATION_RANGE,
            MOCK_DECAY_POWER,
            MOCK_HOLE_EFFECT_RANGE,
        ),
        np.array([0.63212056, 0.63212056, 0.63212056, 0.63212056, 0.63212056, 0.63212056]),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT
    # ---- Decay power
    MOCK_DECAY_POWER = 1.0

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        bessel_exponential(
            MOCK_LAGS,
            MOCK_NUGGET,
            MOCK_SILL,
            MOCK_CORRELATION_RANGE,
            MOCK_DECAY_POWER,
            MOCK_HOLE_EFFECT_RANGE,
        ),
        np.array([0.00000000, 0.63212056, 0.86466472, 0.95021293, 0.98168436, 0.99326205]),
    )
    # Evaluate [ ARRAY ] AND Assert convergence with base model lacking hole effect [ EXPONENTIAL ]
    assert np.allclose(
        bessel_exponential(
            MOCK_LAGS,
            MOCK_NUGGET,
            MOCK_SILL,
            MOCK_CORRELATION_RANGE,
            MOCK_DECAY_POWER,
            MOCK_HOLE_EFFECT_RANGE,
        ),
        exponential(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT + CHANGE DECAY POWER ORDER
    # ---- Decay power
    MOCK_DECAY_POWER = 2.0

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        bessel_exponential(
            MOCK_LAGS,
            MOCK_NUGGET,
            MOCK_SILL,
            MOCK_CORRELATION_RANGE,
            MOCK_DECAY_POWER,
            MOCK_HOLE_EFFECT_RANGE,
        ),
        np.array([0.00000000, 0.63212056, 0.98168436, 0.99987659, 0.99999989, 1.00000000]),
    )
    # Evaluate [ ARRAY ] AND Assert convergence with base model lacking hole effect [ GAUSSIAN ]
    assert np.allclose(
        bessel_exponential(
            MOCK_LAGS,
            MOCK_NUGGET,
            MOCK_SILL,
            MOCK_CORRELATION_RANGE,
            MOCK_DECAY_POWER,
            MOCK_HOLE_EFFECT_RANGE,
        ),
        gaussian(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
    )

    # -------------------------
    # Mock values [ TEST 4 ]: BASE NUGGET EFFECT
    # ---- Nugget
    MOCK_NUGGET = 0.1

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        bessel_exponential(
            MOCK_LAGS,
            MOCK_NUGGET,
            MOCK_SILL,
            MOCK_CORRELATION_RANGE,
            MOCK_DECAY_POWER,
            MOCK_HOLE_EFFECT_RANGE,
        ),
        np.array([0.10000000, 0.6689085, 0.98351593, 0.99988893, 0.9999999, 1.00000000]),
    )


def test_bessel_gaussian():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.2
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        bessel_gaussian(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array([0.00000000, 0.63212056, 0.98168436, 0.99987659, 0.99999989, 1.00000000]),
    )
    # Evaluate [ ARRAY ] AND Assert convergence with base model lacking hole effect [ GAUSSIAN ]
    assert np.allclose(
        bessel_gaussian(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        gaussian(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO BASE NUGGET EFFECT + HOLE EFFECT
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        bessel_gaussian(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array([0.00000000, 0.63205735, 0.98129173, 0.9989769, 0.99840053, 0.99750156]),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: BASE NUGGET EFFECT
    # ---- Hole effect range
    MOCK_NUGGET = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        bessel_gaussian(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array([0.10000000, 0.66885161, 0.98316255, 0.99907921, 0.99856047, 0.99775141]),
    )


def test_cosine_exponential():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT + NO ENHANCEMENT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.2
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00
    # ---- Semivariance enhancement
    MOCK_ENHANCEMENT = False

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        cosine_exponential(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_ENHANCEMENT,
        ),
        np.array([0.00000000, 0.63212056, 0.86466472, 0.95021293, 0.98168436, 0.99326205]),
    )
    # Evaluate [ ARRAY ] AND Assert convergence with base model lacking hole effect [ EXPONENTIAL ]
    assert np.allclose(
        cosine_exponential(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_ENHANCEMENT,
        ),
        exponential(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO BASE NUGGET EFFECT + HOLE EFFECT
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        cosine_exponential(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_ENHANCEMENT,
        ),
        np.array([0.00000000, 0.63219413, 0.86477297, 0.95030252, 0.98174294, 0.99329571]),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: BASE NUGGET EFFECT
    # ---- Hole effect range
    MOCK_NUGGET = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        cosine_exponential(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_ENHANCEMENT,
        ),
        np.array([0.10000000, 0.66897472, 0.87829567, 0.95527227, 0.98356865, 0.99396614]),
    )

    # -------------------------
    # Mock values [ TEST 4 ]: ENHANCE SEMIVARIANCE
    # ---- Semivariance enhancement
    MOCK_ENHANCEMENT = True

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        cosine_exponential(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_ENHANCEMENT,
        ),
        np.array([1.90000000, 1.33102528, 1.12170433, 1.04472773, 1.01643135, 1.00603386]),
    )


def test_cosine_gaussian():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT + NO ENHANCEMENT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.2
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.0

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        cosine_gaussian(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array(
            [
                1.00000000e00,
                3.67879441e-01,
                1.83156389e-02,
                1.23409804e-04,
                1.12535175e-07,
                1.38879439e-11,
            ]
        ),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO BASE NUGGET EFFECT + HOLE EFFECT
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        cosine_gaussian(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array(
            [
                1.00000000e00,
                3.67805868e-01,
                1.83009883e-02,
                1.23187733e-04,
                1.12175254e-07,
                1.38185620e-11,
            ]
        ),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: BASE NUGGET EFFECT
    # ---- Hole effect range
    MOCK_NUGGET = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        cosine_gaussian(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array([1.0000000, 0.43102528, 0.11647089, 0.10011087, 0.1000001, 0.1000000]),
    )


def test_exponential_linear():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT + ZEROTH DECAY POWER ORDER
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.2
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00
    # ---- Decay power
    MOCK_DECAY_POWER = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        exponential_linear(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_DECAY_POWER,
        ),
        np.array([0.63212056, 0.63212056, 0.63212056, 0.63212056, 0.63212056, 0.63212056]),
    )
    # Evaluate [ ARRAY ] AND Assert convergence with base model lacking hole effect
    # [ BESSEL-EXPONENTIAL ]
    assert np.allclose(
        exponential_linear(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_DECAY_POWER,
        ),
        bessel_exponential(MOCK_LAGS, MOCK_NUGGET, MOCK_SILL, MOCK_CORRELATION_RANGE, 0.00, 0.00),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO BASE NUGGET EFFECT + HOLE EFFECT
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        exponential_linear(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_DECAY_POWER,
        ),
        np.array([0.5689085, 0.5689085, 0.5689085, 0.5689085, 0.5689085, 0.5689085]),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: CHANGE DECAY POWER ORDER
    # ---- Decay power
    MOCK_DECAY_POWER = 1.50

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        exponential_linear(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_DECAY_POWER,
        ),
        np.array([0.00000000, 0.6264667, 0.9170913, 0.94824374, 0.92813437, 0.89998745]),
    )

    # -------------------------
    # Mock values [ TEST 4 ]: BASE NUGGET
    # ---- Nugget
    MOCK_NUGGET = 0.1

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        exponential_linear(
            MOCK_LAGS,
            MOCK_SILL,
            MOCK_NUGGET,
            MOCK_CORRELATION_RANGE,
            MOCK_HOLE_EFFECT_RANGE,
            MOCK_DECAY_POWER,
        ),
        np.array([0.10000000, 0.66382003, 0.92538217, 0.95341937, 0.93532093, 0.9099887]),
    )


def test_gaussian_linear():

    # -------------------------
    # Mock values [ TEST 1 ]: NO BASE NUGGET EFFECT + NO HOLE EFFECT
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- Sill
    MOCK_SILL = 1.00
    # ---- Nugget
    MOCK_NUGGET = 0.00
    # ---- Correlation range
    MOCK_CORRELATION_RANGE = 0.2
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.00

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        gaussian_linear(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array([0.00000000, 0.63212056, 0.98168436, 0.99987659, 0.99999989, 1.00000000]),
    )
    # Evaluate [ ARRAY ] AND Assert convergence with base model lacking hole effect [ GAUSSIAN ]
    assert np.allclose(
        gaussian_linear(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        gaussian(MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE),
    )

    # -------------------------
    # Mock values [ TEST 2 ]: NO BASE NUGGET EFFECT + HOLE EFFECT
    # ---- Hole effect range
    MOCK_HOLE_EFFECT_RANGE = 0.10

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        gaussian_linear(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array([0.00000000, 0.62959208, 0.96597741, 0.96388103, 0.93599989, 0.90000000]),
    )

    # -------------------------
    # Mock values [ TEST 3 ]: BASE NUGGET
    # ---- Nugget
    MOCK_NUGGET = 0.1

    # -------------------------
    # Evaluate [ ARRAY ] AND Assert
    assert np.allclose(
        gaussian_linear(
            MOCK_LAGS, MOCK_SILL, MOCK_NUGGET, MOCK_CORRELATION_RANGE, MOCK_HOLE_EFFECT_RANGE
        ),
        np.array([0.10000000, 0.66663287, 0.96937967, 0.96749293, 0.94239991, 0.91000000]),
    )


def test_variogram():

    # -------------------------
    # Mock values
    # ---- Distance lags
    MOCK_LAGS = np.linspace(0.00, 1.00, 6)
    # ---- PARAMETER dictionary
    PARAMETERS = {
        # "distance_lags": np.linspace(0.00, 1.00, 6),
        "sill": 1.00,
        "nugget": 0.10,
        "correlation_range": 0.20,
        "hole_effect_range": 0.10,
        "decay_power": 1.5,
        "enhance_semivariance": False,
    }
    # ----- MODELS list
    MODELS = list(VARIOGRAM_MODELS["single"].keys()) + list(VARIOGRAM_MODELS["composite"].keys())
    # -------- Update formatting of tuples into lists to match expected argument input
    MODELS_ARG = [
        list(model_name) if isinstance(model_name, tuple) else model_name for model_name in MODELS
    ]

    # -------------------------
    # Evaluate [ MODELS ] AND Assert
    # ---- Loop through each model to ensure that the wrapper function `variogram` is mapping
    # ---- the correct variogram model functions
    for model in MODELS_ARG:
        # ---- Get the associated arguments
        args, fun = get_variogram_arguments(model)
        # ---- Filter out the parameters
        SUB_PARAMETERS = {key: value for key, value in PARAMETERS.items() if key in list(args)}
        # ---- Compute accessing function directly
        direct_function = fun["model_function"](MOCK_LAGS, **SUB_PARAMETERS)
        # ---- Update the dictionary with the model name
        SUB_PARAMETERS.update({"model": model})
        # ---- Access via the `variogram` wrapper
        indirect_function = variogram(MOCK_LAGS, model=model, variogram_parameters=SUB_PARAMETERS)
        # ---- ASSERT
        assert np.allclose(direct_function, indirect_function)
