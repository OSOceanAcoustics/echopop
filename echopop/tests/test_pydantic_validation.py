import re

import numpy as np
import pytest
from pydantic import ValidationError

from echopop.utils.validate_dict import KrigingAnalysis, KrigingParameterInputs, MeshCrop


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (dict(), None, "Both 'correlation_range' and 'search_radius' arguments are missing"),
        (
            dict(correlation_range=1.0),
            dict(anisotropy=0.0, kmin=3, kmax=10, correlation_range=1.0, search_radius=3.0),
            None,
        ),
        (
            dict(correlation_range=1.0, search_radius=5.0),
            dict(anisotropy=0.0, kmin=3, kmax=10, correlation_range=1.0, search_radius=5.0),
            None,
        ),
        (
            dict(anisotropy=1, correlation_range=2, search_radius=3),
            dict(anisotropy=1.0, kmin=3, kmax=10, correlation_range=2.0, search_radius=3.0),
            None,
        ),
        (
            dict(kmin=3.0, kmax=10.0, correlation_range=1.0),
            None,
            ["Value must be a non-negative integer", "Value must be a non-negative integer"],
        ),
        (
            dict(kmin=10, kmax=3, correlation_range=1.0),
            None,
            "Defined 'kmax' (3) must be greater than or equal to 'kmin' (10)",
        ),
        (
            dict(anisotropy=0.00, kmin=1, kmax=2, correlation_range=-1.0, search_radius=-1.0),
            None,
            [
                "Input should be greater than or equal to 3",
                "Input should be greater than or equal to 3",
                "Value must be a non-negative float",
                "Value must be a non-negative float",
            ],
        ),
        (
            dict(
                anisotropy=np.nan,
                kmin=np.nan,
                kmax=np.nan,
                correlation_range=np.nan,
                search_radius=np.nan,
            ),
            None,
            [
                "Input should be a finite number",
                "Value must be a non-negative integer",
                "Value must be a non-negative integer",
                "Input should be greater than 0",
                "Input should be greater than 0",
            ],
        ),
        (
            dict(
                anisotropy=np.inf,
                kmin=np.inf,
                kmax=np.inf,
                correlation_range=np.inf,
                search_radius=np.inf,
            ),
            None,
            [
                "Value must be a non-negative real number",
                "Value must be a non-negative integer",
                "Value must be a non-negative integer",
                "Value must be a non-negative real number",
                "Value must be a non-negative real number",
            ],
        ),
    ],
    ids=[
        "Empty inputs [invalid]",
        "Produce valid 'search_radius' based on valid 'correlation_range' input",
        "Define both 'search_radius' and 'correlation_range'",
        "Coerce integer inputs for 'anisotropy', 'correlation_range', and 'search_radius'",
        "Invalid float datatyping for 'kmin' and 'kmax'",
        "Enforcing 'kmax' > 'kmin'",
        "Invalid values for numerics [lower limits]",
        "NaN inputs",
        "Inf inputs",
    ],
)
def test_KrigingParameters_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert KrigingParameterInputs.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert KrigingParameterInputs.create(**input)
    else:
        result = KrigingParameterInputs.create(**input)
        assert result == expected


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            dict(),
            dict(
                best_fit_variogram=False,
                coordinate_transform=True,
                extrapolate=False,
                variable="biomass",
                verbose=True,
            ),
            None,
        ),
        (
            dict(best_fit_variogram=3, coordinate_transform=3, extrapolate=3, verbose=3),
            None,
            [
                "Input should be a valid boolean",
                "Input should be a valid boolean",
                "Input should be a valid boolean",
                "Input should be a valid boolean",
            ],
        ),
        (dict(variable="krakens"), None, "Input should be 'biomass'"),
    ],
    ids=[
        "Default values [no inputs, empty dictionary]",
        "Invalid boolean inputs [integers, not bool/str]",
        "Invalid Literal input for 'variable'",
    ],
)
def test_KrigingAnalysis_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert KrigingAnalysis.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert KrigingAnalysis.create(**input)
    else:
        result = KrigingAnalysis.create(**input)
        assert result == expected


@pytest.mark.parametrize(
    "input, expected, exception",
    [
        (
            dict(),
            dict(
                crop_method="transect_ends",
                num_nearest_transects=4,
                mesh_buffer_distance=1.25,
                latitude_resolution=1.25,
                bearing_tolerance=15.0,
            ),
            None,
        ),
        (dict(crop_method="invalid"), None, "Input should be 'transect_ends' or 'convex_hull'"),
        (
            dict(mesh_buffer_distance=1.0, latitude_resolution=1.0, bearing_tolerance=15),
            dict(
                crop_method="transect_ends",
                num_nearest_transects=4,
                mesh_buffer_distance=1.0,
                latitude_resolution=1.0,
                bearing_tolerance=15.0,
            ),
            None,
        ),
        (
            dict(num_nearest_transects=1.0),
            None,
            "Value error, Value must be a non-negative integer.",
        ),
        (dict(bearing_tolerance="a"), None, "Value must be a non-negative real angle"),
        (
            dict(mesh_buffer_distance=-1.0, latitude_resolution=-1.0, bearing_tolerance=-1.0),
            None,
            [
                "Value must be a non-negative float",
                "Value must be a non-negative float",
                "Value must be a non-negative real angle",
            ],
        ),
        (
            dict(
                num_nearest_transects=np.nan,
                mesh_buffer_distance=np.nan,
                latitude_resolution=np.nan,
                bearing_tolerance=np.nan,
            ),
            None,
            [
                "Value must be a non-negative integer.",
                "Input should be greater than 0",
                "Input should be greater than 0",
                "Input should be greater than 0",
            ],
        ),
        (
            dict(
                num_nearest_transects=np.inf,
                mesh_buffer_distance=np.inf,
                latitude_resolution=np.inf,
                bearing_tolerance=np.inf,
            ),
            None,
            [
                "Value must be a non-negative integer.",
                "Value must be a non-negative real number",
                "Value must be a non-negative real number",
                "Value must be a non-negative real angle",
            ],
        ),
        (dict(bearing_tolerance=181.0), None, "Input should be less than or equal to 180"),
    ],
    ids=[
        "Default values [no inputs, empty dictionary]",
        "Invalid Literal for 'crop_method'",
        "Valid int-to-float coercion",
        "Invalid floats where value should be int",
        "Invalid 'bearing_tolerance' [str]",
        "Invalid values below limits",
        "All NaN values for numeric inputs",
        "All Inf values for numeric inputs",
        "Invalid 'bearing_tolerance' input [upper limit]",
    ],
)
def test_MeshCrop_model(input, expected, exception):

    # -------------------------
    if exception is not None:
        if isinstance(exception, list):
            for e in exception:
                with pytest.raises(ValidationError, match=re.escape(e)):
                    assert MeshCrop.create(**input)
        else:
            with pytest.raises(ValidationError, match=re.escape(exception)):
                assert MeshCrop.create(**input)
    else:
        result = MeshCrop.create(**input)
        assert result == expected
