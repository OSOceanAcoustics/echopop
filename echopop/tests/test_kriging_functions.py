import numpy as np

from echopop.spatial.krige import kriging_lambda, kriging_matrix


def test_kriging_lambda():

    # Mock lagged semivariance values
    test_lagged_semivariogram = np.array([0.0, 0.50, 0.75, 0.875, 0.9375, 1.0])

    # Mock anisotropy value
    test_anisotropy = 0.001

    # Mock kriging matrix
    test_kriging_matrix = np.array(
        [
            [0.00, 0.11, 0.20, 0.24, 0.25, 1.00],
            [0.11, 0.00, 0.11, 0.20, 0.24, 1.00],
            [0.20, 0.11, 0.00, 0.11, 0.20, 1.00],
            [0.24, 0.20, 0.11, 0.00, 0.11, 1.00],
            [0.25, 0.24, 0.20, 0.11, 0.00, 1.00],
            [1.00, 1.00, 1.00, 1.00, 1.00, 0.00],
        ]
    )

    # Evaluate `kriging_lambda`
    eval_kriging_lambda = kriging_lambda(
        test_anisotropy, test_lagged_semivariogram, test_kriging_matrix
    )

    # -----------------
    # Test for equality
    # -----------------
    # Expected outcome
    expected_array = np.array(
        [3.2621651, -1.01812305, -0.14593900, -0.16842517, -0.9296778, 0.41402284]
    )
    # Test
    assert np.allclose(eval_kriging_lambda, expected_array)


def test_kriging_matrix():

    # Mock x-coordinates
    test_x_coordinates = np.linspace(0.0, 10.0, 5)

    # Mock y-coordinates
    test_y_coordinates = np.linspace(0.0, 10.0, 5)

    # Mock variogram parameters
    test_variogram_parameters = {
        "range": 10.0,
        "correlation_range": 2.0,
        "hole_effect_range": 0.0,
        "nugget": 0.0,
        "sill": 5.0,
        "decay_power": 1.0,
        "model": ["bessel", "exponential"],
    }

    # Evaluate `kriging_matrix`
    eval_kriging_matrix = kriging_matrix(
        test_x_coordinates, test_y_coordinates, test_variogram_parameters
    )

    # -----------------
    # Test for equality
    # -----------------
    # Expected outcome
    expected_matrix = np.array(
        [
            [0.00000000, 4.14643112, 4.85428403, 4.97512428, 4.99575337, 1.00000000],
            [4.14643112, 0.00000000, 4.14643112, 4.85428403, 4.97512428, 1.00000000],
            [4.85428403, 4.14643112, 0.00000000, 4.14643112, 4.85428403, 1.00000000],
            [4.97512428, 4.85428403, 4.14643112, 0.00000000, 4.14643112, 1.00000000],
            [4.99575337, 4.97512428, 4.85428403, 4.14643112, 0.00000000, 1.00000000],
            [1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 0.00000000],
        ]
    )

    # Test
    assert np.allclose(eval_kriging_matrix, expected_matrix)
