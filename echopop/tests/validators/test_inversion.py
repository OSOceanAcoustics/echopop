import numpy as np
import pytest
from pydantic import ValidationError

from echopop.validators.inversion import (
    TSLRegressionParameters,
    ValidateLengthTS,
)


# ==================================================================================================
# Test TSLRegressionParameters
# ----------------------------
def test_tsl_regression_parameters_valid():
    """Test TSLRegressionParameters with valid parameters."""

    params = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    assert params.slope == -20.5
    assert params.intercept == 67.3


def test_tsl_regression_parameters_zero_values():
    """Test TSLRegressionParameters with zero values."""

    params = TSLRegressionParameters(slope=0.0, intercept=0.0)

    assert params.slope == 0.0
    assert params.intercept == 0.0


def test_tsl_regression_parameters_positive_slope():
    """Test TSLRegressionParameters with positive slope."""

    params = TSLRegressionParameters(slope=15.2, intercept=-45.8)

    assert params.slope == 15.2
    assert params.intercept == -45.8


def test_tsl_regression_parameters_invalid_slope():
    """Test TSLRegressionParameters with invalid slope values."""

    # Test infinite slope
    with pytest.raises(ValidationError):
        TSLRegressionParameters(slope=np.inf, intercept=67.3)

    # Test NaN slope
    with pytest.raises(ValidationError):
        TSLRegressionParameters(slope=np.nan, intercept=67.3)


def test_tsl_regression_parameters_invalid_intercept():
    """Test TSLRegressionParameters with invalid intercept values."""

    # Test infinite intercept
    with pytest.raises(ValidationError):
        TSLRegressionParameters(slope=-20.5, intercept=np.inf)

    # Test NaN intercept
    with pytest.raises(ValidationError):
        TSLRegressionParameters(slope=-20.5, intercept=np.nan)


def test_tsl_regression_parameters_realistic_hake_values():
    """Test TSLRegressionParameters with realistic hake TS-length parameters."""

    # Typical hake TS-length regression parameters from literature
    params = TSLRegressionParameters(slope=-20.8, intercept=67.5)

    assert params.slope == -20.8
    assert params.intercept == 67.5


# ==================================================================================================
# Test ValidateLengthTS
# ---------------------
def test_validate_length_ts_valid():
    """Test ValidateLengthTS with valid parameters."""

    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    params = ValidateLengthTS(
        ts_length_regression=ts_regression,
        stratify_by=["stratum_ks"],
        expected_strata=np.array([1, 2, 3, 4, 5]),
        impute_missing_strata=True,
        haul_replicates=True,
    )

    assert params.ts_length_regression.slope == -20.5
    assert params.ts_length_regression.intercept == 67.3
    assert params.stratify_by == ["stratum_ks"]
    np.testing.assert_array_equal(params.expected_strata, np.array([1, 2, 3, 4, 5]))
    assert params.impute_missing_strata is True
    assert params.haul_replicates is True


def test_validate_length_ts_defaults():
    """Test ValidateLengthTS with default values."""

    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    params = ValidateLengthTS(ts_length_regression=ts_regression, stratify_by=["stratum_ks"])

    assert params.ts_length_regression.slope == -20.5
    assert params.expected_strata is None
    assert params.impute_missing_strata is True  # Default value
    assert params.haul_replicates is True  # Default value


def test_validate_length_ts_single_stratify_string():
    """Test ValidateLengthTS with single string for stratify_by."""

    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    params = ValidateLengthTS(
        ts_length_regression=ts_regression, stratify_by="stratum_ks"  # Single string
    )

    assert params.stratify_by == ["stratum_ks"]  # Should be converted to list


def test_validate_length_ts_multiple_stratify_by():
    """Test ValidateLengthTS with multiple stratification columns."""

    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    params = ValidateLengthTS(
        ts_length_regression=ts_regression, stratify_by=["stratum_ks", "region", "year"]
    )

    assert params.stratify_by == ["stratum_ks", "region", "year"]


def test_validate_length_ts_expected_strata_list():
    """Test ValidateLengthTS with expected_strata as list."""

    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    params = ValidateLengthTS(
        ts_length_regression=ts_regression,
        stratify_by=["stratum_ks"],
        expected_strata=[1, 2, 3, 4, 5],  # List instead of array
    )

    np.testing.assert_array_equal(params.expected_strata, np.array([1, 2, 3, 4, 5]))


def test_validate_length_ts_expected_strata_none():
    """Test ValidateLengthTS with expected_strata as None."""

    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    params = ValidateLengthTS(
        ts_length_regression=ts_regression, stratify_by=["stratum_ks"], expected_strata=None
    )

    assert params.expected_strata is None


def test_validate_length_ts_boolean_flags():
    """Test ValidateLengthTS with various boolean flag combinations."""

    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    # Test all combinations of boolean flags
    boolean_combinations = [(True, True), (True, False), (False, True), (False, False)]

    for impute, haul_rep in boolean_combinations:
        params = ValidateLengthTS(
            ts_length_regression=ts_regression,
            stratify_by=["stratum_ks"],
            impute_missing_strata=impute,
            haul_replicates=haul_rep,
        )

        assert params.impute_missing_strata == impute
        assert params.haul_replicates == haul_rep


def test_validate_length_ts_realistic_hake_scenario():
    """Test ValidateLengthTS with realistic hake survey parameters."""

    # Realistic hake TS-length regression
    ts_regression = TSLRegressionParameters(slope=-20.8, intercept=67.5)

    # Realistic hake survey stratification
    expected_strata = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])  # Typical hake strata

    params = ValidateLengthTS(
        ts_length_regression=ts_regression,
        stratify_by=["stratum_ks"],
        expected_strata=expected_strata,
        impute_missing_strata=True,
        haul_replicates=True,  # Recommended to avoid pseudoreplication
    )

    assert params.ts_length_regression.slope == -20.8
    assert params.ts_length_regression.intercept == 67.5
    assert params.stratify_by == ["stratum_ks"]
    np.testing.assert_array_equal(params.expected_strata, expected_strata)
    assert params.impute_missing_strata is True
    assert params.haul_replicates is True


def test_validate_length_ts_multiple_stratification():
    """Test ValidateLengthTS with complex stratification scenario."""

    ts_regression = TSLRegressionParameters(slope=-19.2, intercept=65.1)

    params = ValidateLengthTS(
        ts_length_regression=ts_regression,
        stratify_by=["stratum_ks", "transect_num", "year"],
        expected_strata=np.arange(1, 21),  # 20 strata
        impute_missing_strata=False,
        haul_replicates=True,
    )

    assert len(params.stratify_by) == 3
    assert "stratum_ks" in params.stratify_by
    assert "transect_num" in params.stratify_by
    assert "year" in params.stratify_by
    assert len(params.expected_strata) == 20
    assert params.impute_missing_strata is False


# ==================================================================================================
# Integration tests
# -----------------
def test_inversion_validators_integration():
    """Test integration of inversion validators."""

    # Create TS-length regression parameters
    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)

    # Test that we can use it in the main validator
    length_ts_config = ValidateLengthTS(
        ts_length_regression=ts_regression,
        stratify_by=["stratum_ks", "region"],
        expected_strata=np.array([1, 2, 3, 4, 5, 6]),
        impute_missing_strata=True,
        haul_replicates=True,
    )

    # Verify integration works
    assert isinstance(length_ts_config.ts_length_regression, TSLRegressionParameters)
    assert length_ts_config.ts_length_regression.slope == -20.5
    assert length_ts_config.ts_length_regression.intercept == 67.3
    assert len(length_ts_config.stratify_by) == 2
    assert len(length_ts_config.expected_strata) == 6


def test_inversion_validators_create_factory_method():
    """Test inversion validators with create factory method."""

    # Test TSLRegressionParameters create method
    ts_params_dict = TSLRegressionParameters.create(slope=-20.5, intercept=67.3)
    expected_ts = {"slope": -20.5, "intercept": 67.3}
    assert ts_params_dict == expected_ts

    # Test ValidateLengthTS create method
    ts_regression = TSLRegressionParameters(slope=-20.5, intercept=67.3)
    length_ts_dict = ValidateLengthTS.create(
        ts_length_regression=ts_regression,
        stratify_by=["stratum_ks"],
        impute_missing_strata=False,
        haul_replicates=True,
    )

    assert length_ts_dict["stratify_by"] == ["stratum_ks"]
    assert length_ts_dict["impute_missing_strata"] is False
    assert length_ts_dict["haul_replicates"] is True
    # expected_strata should be excluded (None values excluded)
    assert "expected_strata" not in length_ts_dict


def test_inversion_validators_edge_cases():
    """Test inversion validators with edge cases."""

    # Test with very small slope (near zero but not zero)
    ts_regression = TSLRegressionParameters(slope=-0.001, intercept=100.0)
    params = ValidateLengthTS(ts_length_regression=ts_regression, stratify_by=["stratum_ks"])
    assert params.ts_length_regression.slope == -0.001

    # Test with very large intercept
    ts_regression = TSLRegressionParameters(slope=-20.0, intercept=1000.0)
    params = ValidateLengthTS(ts_length_regression=ts_regression, stratify_by=["stratum_ks"])
    assert params.ts_length_regression.intercept == 1000.0

    # Test with single stratum
    params = ValidateLengthTS(
        ts_length_regression=ts_regression,
        stratify_by=["stratum_ks"],
        expected_strata=np.array([1]),
    )
    assert len(params.expected_strata) == 1
    assert params.expected_strata[0] == 1
