import pytest

from echopop import inversion


def test_inversion_base_abstract():
    """Test that InversionBase cannot be instantiated directly."""
    with pytest.raises(TypeError):
        inversion.InversionBase({})


def test_inversion_base_stratify_by_string(model_parameters):
    """Test that single string stratify_by gets converted to list."""
    params = model_parameters.copy()
    params["stratify_by"] = "stratum_ks"  # Single string

    # Create a concrete subclass for testing
    class TestInversion(inversion.InversionBase):
        def invert(self, df_nasc):
            return df_nasc

    inverter = TestInversion(params)
    assert isinstance(inverter.model_params["stratify_by"], list)
    assert inverter.model_params["stratify_by"] == ["stratum_ks"]


def test_inversion_base_stratify_by_list(model_parameters_multiple_strata):
    """Test that list stratify_by remains as list."""

    class TestInversion(inversion.InversionBase):
        def invert(self, df_nasc):
            return df_nasc

    inverter = TestInversion(model_parameters_multiple_strata)
    assert isinstance(inverter.model_params["stratify_by"], list)
    assert len(inverter.model_params["stratify_by"]) == 2
