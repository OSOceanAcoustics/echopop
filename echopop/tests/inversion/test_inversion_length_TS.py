import pandas as pd
import pytest
from pydantic import ValidationError

from echopop import inversion


def test_inversion_length_ts_init(model_parameters):
    """Test InversionLengthTS initialization."""
    inverter = inversion.InversionLengthTS(model_parameters)

    assert inverter.inversion_method == "length_TS_regression"
    assert inverter.sigma_bs_haul is None
    assert "ts_length_regression" in inverter.model_params


def test_inversion_length_ts_init_missing_params():
    """Test initialization with missing required parameters."""
    bad_params = {
        "stratify_by": "stratum_ks",
        # Missing ts_length_regression
    }

    # Expected to raise a Validation Error
    with pytest.raises(ValidationError):
        assert inversion.InversionLengthTS(bad_params)


def test_get_stratified_sigma_bs_with_data(model_parameters, specimen_df):
    """Test get_stratified_sigma_bs with provided data."""
    inverter = inversion.InversionLengthTS(model_parameters)

    inverter.get_stratified_sigma_bs(specimen_df)

    # Check types and formatting
    # ---- sigma_bs_haul
    assert isinstance(inverter.sigma_bs_haul, pd.DataFrame)
    assert "sigma_bs" in inverter.sigma_bs_haul.columns
    assert len(inverter.sigma_bs_haul) > 0
    # ---- sigma_bs_strata
    assert isinstance(inverter.sigma_bs_strata, pd.DataFrame)
    assert "sigma_bs" in inverter.sigma_bs_strata.columns
    assert len(inverter.sigma_bs_strata) > 0


def test_invert_basic(model_parameters, nasc_df, specimen_df):
    """Test basic inversion functionality."""
    inverter = inversion.InversionLengthTS(model_parameters)

    result = inverter.invert(df_nasc=nasc_df, df_length=specimen_df)

    assert "number_density" in result.columns
    assert len(result) == len(nasc_df)
    assert not result["number_density"].isna().any()
    assert (result["number_density"] >= 0).all()


def test_invert_with_imputation(model_parameters, specimen_df):
    """Test inversion with missing strata imputation."""
    # Create NASC data with additional stratum not in specimen data
    nasc_df = pd.DataFrame(
        {
            "stratum_ks": [1, 2, 3, 4],  # 4 not in specimen_df
            "nasc": [1000, 1500, 1200, 800],
            "nasc_proportion": [1.0, 1.0, 1.0, 1.0],
        }
    )

    inverter = inversion.InversionLengthTS(model_parameters)
    result = inverter.invert(nasc_df, specimen_df)

    assert "number_density" in result.columns
    assert len(result) == len(nasc_df)
    # Should not have NaN values due to imputation
    assert not result["number_density"].isna().any()


def test_invert_no_imputation(model_parameters_no_impute, nasc_df, specimen_df):
    """Test inversion without imputation."""
    inverter = inversion.InversionLengthTS(model_parameters_no_impute)

    result = inverter.invert(nasc_df, specimen_df)

    assert "number_density" in result.columns
    assert len(result) == len(nasc_df)


def test_invert_preserves_original_data(model_parameters, nasc_df, specimen_df):
    """Test that inversion preserves original NASC data."""
    original_nasc = nasc_df.copy()
    inverter = inversion.InversionLengthTS(model_parameters)

    result = inverter.invert(nasc_df, specimen_df)

    # Original columns should be preserved (except for index changes)
    for col in ["nasc", "nasc_proportion"]:
        if col in original_nasc.columns:
            assert col in result.columns


def test_invert_large_dataset(model_parameters, large_specimen_df):
    """Test inversion performance with larger dataset."""
    # Create corresponding NASC data
    nasc_df = pd.DataFrame(
        {
            "stratum_ks": [1, 2, 3, 4, 5],
            "nasc": [2000, 2500, 1800, 2200, 1900],
            "nasc_proportion": [1.0, 1.0, 1.0, 1.0, 1.0],
        }
    )

    inverter = inversion.InversionLengthTS(model_parameters)

    # Should complete without errors
    result = inverter.invert(nasc_df, large_specimen_df)

    assert "number_density" in result.columns
    assert len(result) == len(nasc_df)


# ==============================================================================
# INTEGRATION TESTS
# ==============================================================================


def test_full_workflow(model_parameters, specimen_df, length_df, nasc_df):
    """Test complete workflow from data to inversion."""
    inverter = inversion.InversionLengthTS(model_parameters)

    # Perform inversion
    result = inverter.invert(nasc_df, [specimen_df, length_df])

    # Verify complete workflow
    assert "number_density" in result.columns
    assert len(result) == len(nasc_df)
    assert (result["number_density"] > 0).all()

    # Check that results are reasonable orders of magnitude
    assert (result["number_density"] < 1e10).all()  # Not absurdly large
    assert (result["number_density"] > 1e-10).all()  # Not absurdly small


def test_workflow_different_strata_coverage():
    """Test workflow when biological and acoustic data have different strata."""
    # Biological data has strata 1, 2, 3
    bio_df = pd.DataFrame(
        {
            "stratum_ks": [1, 1, 2, 2, 3, 3],
            "haul_num": [101, 102, 201, 202, 301, 302],
            "length": [20, 22, 24, 21, 23, 25],
        }
    )

    # NASC data has strata 1, 2, 3, 4, 5
    nasc_df = pd.DataFrame(
        {
            "stratum_ks": [1, 2, 3, 4, 5],
            "nasc": [1000, 1200, 800, 1500, 900],
            "nasc_proportion": [1.0, 1.0, 1.0, 1.0, 1.0],
        }
    )

    params = {
        "ts_length_regression": {"slope": 20.0, "intercept": -68.0},
        "stratify_by": "stratum_ks",
        "expected_strata": [1, 2, 3, 4, 5],
        "impute_missing_strata": True,
        "haul_replicates": True,
    }

    inverter = inversion.InversionLengthTS(params)
    result = inverter.invert(nasc_df, bio_df)

    # Should handle missing strata gracefully
    assert len(result) == len(nasc_df)
    assert "number_density" in result.columns


# ==============================================================================
# ERROR HANDLING TESTS
# ==============================================================================


def test_missing_required_columns():
    """Test error handling for missing required columns."""
    bad_nasc = pd.DataFrame({"wrong_col": [1, 2, 3]})
    specimen_df = pd.DataFrame({"length": [20, 22, 24]})

    params = {
        "ts_length_regression": {"slope": 20.0, "intercept": -68.0},
        "stratify_by": "stratum_ks",
    }

    inverter = inversion.InversionLengthTS(params)

    # Should handle missing columns gracefully or raise appropriate error
    try:
        result = inverter.invert(bad_nasc, specimen_df)
        # If it doesn't raise an error, it should at least not crash
        assert isinstance(result, pd.DataFrame)
    except (KeyError, ValueError):  # Expected behavior for missing required columns
        pass
