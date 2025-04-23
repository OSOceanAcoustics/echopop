import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pytest
import yaml
from pydantic import ValidationError

from ..utils.validate import posfloat, posint, realposfloat
from ..utils.validate_dict import (
    CONFIG_DATA_MODEL,
    CONFIG_INIT_MODEL,
    BiologicalFile,
    BiologicalFiles,
    FileSettings,
    Geospatial,
    HaulTransectMap,
    InputModel,
    KrigingFiles,
    KrigingParameters,
    NASCExports,
    PatternParts,
    SpeciesDefinition,
    StratificationFiles,
    StratifiedSurveyMeanParameters,
    TransectRegionMap,
    TSLRegressionParameters,
    XLSXFile,
)


####################################################################################################
# BASE MODEL STRUCTURE: TEST
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "description",
    ["Assess `InputModel` pydantic model structure"],
    ids=["Assess `InputModel` pydantic model structure"],
)
def test_InputModel_model_structure(description):

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(InputModel)
    # ---- Validate correct setting
    assert InputModel.model_config == dict()

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(InputModel)
    # ---- Validate correct output from empty input
    assert InputModel.judge(**dict()) == InputModel()
    # ---- Validate that any arbitrary output yields an empty output
    assert InputModel.judge(**dict(dum1="dum2")) == InputModel()
    # ---- Validate non-dictionary input
    with pytest.raises(
        TypeError,
        match=re.escape("judge() takes 1 positional argument but 2 were given"),
    ):
        assert InputModel.judge(list())
    # ---- Test empty
    assert InputModel.judge() == InputModel()

    # -------------------------
    # ASSERT: 'create' method
    # ---- Check existence
    assert "create" in dir(InputModel)
    # ---- Validate correct output from empty input
    assert InputModel.create(**dict()) == {}
    # ---- Validate that any arbitrary output yields an empty output
    assert InputModel.create(**dict(dum1="dum2")) == {}


####################################################################################################
# HIGHER MODEL STRUCTURE: FIXTURES
# --------------------------------------------------------------------------------------------------
@pytest.fixture
def BiologicalFiles_fields() -> Dict[str, Any]:

    return {
        "length": {
            "annotation": Union[Dict[str, XLSXFile], XLSXFile],
            "frozen": None,
        },
        "specimen": {
            "annotation": Union[Dict[str, XLSXFile], XLSXFile],
            "frozen": None,
        },
        "catch": {
            "annotation": Union[Dict[str, XLSXFile], XLSXFile],
            "frozen": None,
        },
        "haul_to_transect": {
            "annotation": Union[Dict[str, XLSXFile], XLSXFile, None],
            "default": None,
        },
    }


@pytest.fixture
def FileSettings_fields() -> Dict[str, Any]:

    return {
        "directory": {
            "annotation": str,
            "frozen": None,
        },
        "sheetname": {
            "annotation": str,
            "frozen": None,
        },
    }


@pytest.fixture
def Geospatial_fields() -> Dict[str, Any]:

    return {
        "init": {
            "annotation": str,
            "frozen": None,
        },
    }


@pytest.fixture
def HaulTransectMap_fields() -> Dict[str, Any]:

    return {
        "save_file_template": {
            "annotation": str,
            "frozen": None,
        },
        "country_code": {
            "annotation": List[str],
            "frozen": None,
        },
        "file_settings": {
            "annotation": Dict[str, FileSettings],
            "frozen": None,
        },
    }


@pytest.fixture
def KrigingFiles_fields() -> Dict[str, Any]:

    return {
        "vario_krig_para": {
            "annotation": XLSXFile,
            "frozen": None,
        },
        "isobath_200m": {
            "annotation": XLSXFile,
            "frozen": None,
        },
        "mesh": {
            "annotation": XLSXFile,
            "frozen": None,
        },
    }


@pytest.fixture
def KrigingParameters_fields() -> Dict[str, Any]:

    return {
        "A0": {
            "annotation": posfloat,
            "allow_inf_nan": False,
            "ge": 0.0,
        },
        "longitude_reference": {
            "annotation": float,
            "ge": -180.0,
            "le": 180.0,
            "allow_inf_nan": False,
        },
        "longitude_offset": {
            "annotation": float,
            "allow_inf_nan": False,
        },
        "latitude_offset": {
            "annotation": float,
            "allow_inf_nan": False,
        },
    }


@pytest.fixture
def NASCExports_fields() -> Dict[str, Any]:

    return {
        "export_file_directory": {
            "annotation": str,
            "frozen": None,
        },
        "nasc_export_directory": {
            "annotation": str,
            "frozen": None,
        },
        "save_file_template": {
            "annotation": str,
            "frozen": None,
        },
        "save_file_sheetname": {
            "annotation": str,
            "frozen": None,
        },
        "regions": {
            "annotation": Dict[str, List[str]],
            "frozen": None,
        },
        "max_transect_spacing": {
            "annotation": realposfloat,
            "ge": 0.0,
            "allow_inf_nan": False,
        },
        "file_columns": {
            "annotation": List[str],
            "frozen": None,
        },
    }


@pytest.fixture
def PatternParts_fields() -> Dict[str, Any]:

    return {
        "pattern": {
            "annotation": str,
            "frozen": None,
        },
        "label": {
            "annotation": str,
            "frozen": None,
        },
    }


@pytest.fixture
def StratificationFiles_fields() -> Dict[str, Any]:

    return {
        "geo_strata": {
            "annotation": XLSXFile,
            "frozen": None,
        },
        "strata": {
            "annotation": XLSXFile,
            "frozen": None,
        },
    }


@pytest.fixture
def StratifiedSurveyMeanParameters_fields() -> Dict[str, Any]:

    return {
        "strata_transect_proportion": {
            "annotation": posfloat,
            "frozen": None,
        },
        "num_replicates": {
            "annotation": posint,
            "frozen": None,
        },
        "mesh_transects_per_latitude": {
            "annotation": posint,
            "frozen": None,
        },
    }


@pytest.fixture
def TSLRegressionParameters_fields() -> Dict[str, Any]:

    return {
        "number_code": {
            "annotation": int,
            "frozen": None,
        },
        "TS_L_slope": {
            "annotation": float,
            "allow_inf_nan": False,
        },
        "TS_L_intercept": {
            "annotation": float,
            "allow_inf_nan": False,
        },
        "length_units": {
            "annotation": str,
            "frozen": None,
        },
    }


@pytest.fixture
def TransectRegionMap_fields() -> Dict[str, Any]:

    return {
        "pattern": {
            "annotation": str,
            "frozen": None,
        },
        "parts": {
            "annotation": Dict[str, List[PatternParts]],
            "frozen": None,
        },
    }


@pytest.fixture
def XLSXFile_fields() -> Dict[str, Any]:

    return {
        "filename": {
            "annotation": str,
            "frozen": None,
        },
        "sheetname": {
            "annotation": Union[str, List[str], Dict[str, str]],
            "frozen": None,
        },
    }


####################################################################################################
# HIGHER MODEL STRUCTURE: TESTS
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "description",
    ["Assess `BiologicalFiles` pydantic model structure"],
    ids=["Assess `BiologicalFiles` pydantic model structure"],
)
def test_BiologicalFiles_model_structure(description, BiologicalFiles_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(BiologicalFiles_fields).issubset(BiologicalFiles.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            BiologicalFiles.model_fields[param]._attributes_set == BiologicalFiles_fields[param]
            for param in BiologicalFiles.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert BiologicalFiles.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(BiologicalFiles)
    # ---- Validate correct setting
    assert BiologicalFiles.model_config == dict(title="biological file inputs")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(BiologicalFiles)
    # ---- Check default output
    with pytest.raises(ValidationError, match="3 validation errors for biological file inputs"):
        assert BiologicalFiles.judge(**dict())


@pytest.mark.parametrize(
    "description",
    ["Assess `FileSettings` pydantic model structure"],
    ids=["Assess `FileSettings` pydantic model structure"],
)
def test_FileSettings_model_structure(description, FileSettings_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(FileSettings_fields).issubset(FileSettings.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            FileSettings.model_fields[param]._attributes_set == FileSettings_fields[param]
            for param in FileSettings.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert FileSettings.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(FileSettings)
    # ---- Validate correct setting
    assert FileSettings.model_config == dict(title="parameter file settings")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(FileSettings)
    # ---- Check default output
    with pytest.raises(ValidationError, match="2 validation errors for parameter file settings"):
        assert FileSettings.judge(**dict())


@pytest.mark.parametrize(
    "description",
    ["Assess `Geospatial` pydantic model structure"],
    ids=["Assess `Geospatial` pydantic model structure"],
)
def test_Geospatial_model_structure(description, Geospatial_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(Geospatial_fields).issubset(Geospatial.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            Geospatial.model_fields[param]._attributes_set == Geospatial_fields[param]
            for param in Geospatial.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert Geospatial.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(Geospatial)
    # ---- Validate correct setting
    assert Geospatial.model_config == dict(title="EPSG code")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(Geospatial)
    # ---- Check default output
    with pytest.raises(ValidationError, match="1 validation error for EPSG code"):
        assert Geospatial.judge(**dict())

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = Geospatial.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert "validate_init" in Geospatial.__dict__

    # -------------------------
    # ASSERT: 'validate_init'
    # ---- Check output for 'validate_init' [VALID]
    assert Geospatial.validate_init("epsg:4326") == "epsg:4326"
    # ---- Check output for 'validate_init' [INVALID: float]
    with pytest.raises(ValueError, match=re.escape("Echopop cannot parse the defined EPSG code")):
        assert Geospatial.validate_init(4326.5)
    # ---- Check output for 'validate_init' [VALID: int]
    assert Geospatial.validate_init(4326) == "epsg:4326"
    # ---- Check output for 'validate_float' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Echopop cannot parse the defined EPSG code")):
        assert Geospatial.validate_init(np.nan)
    # ---- Check output for 'validate_float' [INVALID: inf]
    with pytest.raises(ValueError, match=re.escape("Echopop cannot parse the defined EPSG code")):
        assert Geospatial.validate_init(np.inf)
    # ---- Check output for 'validate_init' [VALID: capitalized string]
    assert Geospatial.validate_init("EPSG:4326") == "epsg:4326"
    # ---- Check output for 'validate_init' [VALID: imputed ':']
    assert Geospatial.validate_init("epsg4326") == "epsg:4326"
    # ---- Check the applicable fields
    assert field_decor["validate_init"].info.fields == ("init",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_init"].info.mode == "before"


@pytest.mark.parametrize(
    "description",
    ["Assess `HaulTransectMap` pydantic model structure"],
    ids=["Assess `HaulTransectMap` pydantic model structure"],
)
def test_HaulTransectMap_model_structure(description, HaulTransectMap_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(HaulTransectMap_fields).issubset(HaulTransectMap.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            HaulTransectMap.model_fields[param]._attributes_set == HaulTransectMap_fields[param]
            for param in HaulTransectMap.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert HaulTransectMap.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(HaulTransectMap)
    # ---- Validate correct setting
    assert HaulTransectMap.model_config == dict(
        arbitrary_types_allowed=True, title="haul-transect key mapping"
    )

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(HaulTransectMap)
    # ---- Check default output
    with pytest.raises(ValidationError, match="3 validation errors for haul-transect key mapping"):
        assert HaulTransectMap.judge(**dict())

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = HaulTransectMap.__pydantic_decorators__.field_validators
    # ---- 'model_validators' 'Decorator' object
    model_decor = HaulTransectMap.__pydantic_decorators__.model_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert "validate_save_file_template" in HaulTransectMap.__dict__

    # -------------------------
    # ASSERT: 'validate_save_file_template'
    # ---- Check output for 'validate_save_file_template' [VALID]
    assert HaulTransectMap.validate_save_file_template("") == ""
    # ---- Check output for 'validate_save_file_template' [VALID]
    assert HaulTransectMap.validate_save_file_template("{YEAR}{COUNTRY}") == "{YEAR}{COUNTRY}"
    # ---- Check output for 'validate_save_file_template' [VALID]
    assert HaulTransectMap.validate_save_file_template("{YEAR}") == "{YEAR}"
    # ---- Check output for 'validate_save_file_template' [INVALID: extra identifier]
    with pytest.raises(ValueError, match=re.escape("contains invalid identifiers")):
        assert HaulTransectMap.validate_save_file_template("{YEAR}{EXTRA}")
    # ---- Check the applicable fields
    assert field_decor["validate_save_file_template"].info.fields == ("save_file_template",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_save_file_template"].info.mode == "after"

    # -------------------------
    # ASSERT: model validator decorators
    # ---- Check whether the field exists
    assert "validate_country_files" in HaulTransectMap.__dict__
    # MOCK: parameterized model
    MOCK_VALUES = {
        "save_file_template": "mapping_{COUNTRY}_{YEAR}",
        "country_code": ["CANDY", "LAND"],
        "file_settings": {
            "CANDY": {
                "directory": "outer",
                "sheetname": "space",
            },
            "LAND": {
                "directory": "outer",
                "sheetname": "space",
            },
        },
    }
    # ASSERT: model validator
    assert HaulTransectMap.validate_country_files(MOCK_VALUES) == MOCK_VALUES
    # ---- Check the applicable validation mode
    assert model_decor["validate_country_files"].info.mode == "before"


@pytest.mark.parametrize(
    "description",
    ["Assess `KrigingFiles` pydantic model structure"],
    ids=["Assess `KrigingFiles` pydantic model structure"],
)
def test_KrigingFiles_model_structure(description, KrigingFiles_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(KrigingFiles_fields).issubset(KrigingFiles.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            KrigingFiles.model_fields[param]._attributes_set == KrigingFiles_fields[param]
            for param in KrigingFiles.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert KrigingFiles.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(KrigingFiles)
    # ---- Validate correct setting
    assert KrigingFiles.model_config == dict(title="kriging file inputs")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(KrigingFiles)
    # ---- Check default output
    with pytest.raises(ValidationError, match="3 validation errors for kriging file inputs"):
        assert KrigingFiles.judge(**dict())


@pytest.mark.parametrize(
    "description",
    ["Assess `KrigingParameters` pydantic model structure"],
    ids=["Assess `KrigingParameters` pydantic model structure"],
)
def test_KrigingParameters_model_structure(description, KrigingParameters_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(KrigingParameters_fields).issubset(KrigingParameters.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            KrigingParameters.model_fields[param]._attributes_set == KrigingParameters_fields[param]
            for param in KrigingParameters.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert KrigingParameters.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(KrigingParameters)
    # ---- Validate correct setting
    assert KrigingParameters.model_config == dict(
        arbitrary_types_allowed=True, title="kriging parameters"
    )

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(KrigingParameters)
    # ---- Check default output
    with pytest.raises(ValidationError, match="4 validation errors for kriging parameters"):
        assert KrigingParameters.judge(**dict())

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = KrigingParameters.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert "validate_posfloat" in KrigingParameters.__dict__

    # -------------------------
    # ASSERT: 'valid_posfloat'
    # ---- Check output for 'validate_posfloat' [VALID]
    assert KrigingParameters.validate_posfloat(2) == 2.0
    # ---- Check output for 'validate_posfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert KrigingParameters.validate_posfloat(-2.0)
    # ---- Check output for 'validate_float' [INVALID: NaN]
    assert KrigingParameters.validate_posfloat(np.inf) == np.inf
    # ---- Check the applicable fields
    assert field_decor["validate_posfloat"].info.fields == ("A0",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_posfloat"].info.mode == "before"


@pytest.mark.parametrize(
    "description",
    ["Assess `NASCExports` pydantic model structure"],
    ids=["Assess `NASCExports` pydantic model structure"],
)
def test_NASCExports_model_structure(description, NASCExports_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(NASCExports_fields).issubset(NASCExports.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            NASCExports.model_fields[param]._attributes_set == NASCExports_fields[param]
            for param in NASCExports.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert NASCExports.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(NASCExports)
    # ---- Validate correct setting
    assert NASCExports.model_config == dict(
        arbitrary_types_allowed=True, title="Echoview export processing parameters"
    )

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(NASCExports)
    # ---- Check default output
    with pytest.raises(
        ValidationError, match="7 validation errors for Echoview export processing parameters"
    ):
        assert NASCExports.judge(**dict())

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = NASCExports.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert set(["validate_realposfloat", "validate_save_file_template"]).issubset(
        NASCExports.__dict__
    )

    # -------------------------
    # ASSERT: 'valid_realposfloat'
    # ---- Check output for 'validate_realposfloat' [VALID]
    assert NASCExports.validate_realposfloat(2) == 2.0
    # ---- Check output for 'validate_realposfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert NASCExports.validate_realposfloat(-2.0)
    # ---- Check output for 'validate_realposfloat' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative real number.")):
        assert NASCExports.validate_realposfloat(np.inf)
    # ---- Check the applicable fields
    assert field_decor["validate_realposfloat"].info.fields == ("max_transect_spacing",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_realposfloat"].info.mode == "before"

    # -------------------------
    # ASSERT: 'validate_save_file_template'
    # ---- Check output for 'validate_save_file_template' [VALID]
    assert NASCExports.validate_save_file_template("") == ""
    # ---- Check output for 'validate_save_file_template' [VALID]
    assert NASCExports.validate_save_file_template("{YEAR}{GROUP}") == "{YEAR}{GROUP}"
    # ---- Check output for 'validate_save_file_template' [VALID]
    assert NASCExports.validate_save_file_template("{YEAR}") == "{YEAR}"
    # ---- Check output for 'validate_save_file_template' [INVALID: extra identifier]
    with pytest.raises(ValueError, match=re.escape("contains invalid identifiers")):
        assert NASCExports.validate_save_file_template("{YEAR}{EXTRA}")
    # ---- Check the applicable fields
    assert field_decor["validate_save_file_template"].info.fields == ("save_file_template",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_save_file_template"].info.mode == "after"


@pytest.mark.parametrize(
    "description",
    ["Assess `PatternParts` pydantic model structure"],
    ids=["Assess `PatternParts` pydantic model structure"],
)
def test_PatternParts_model_structure(description, PatternParts_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(PatternParts_fields).issubset(PatternParts.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            PatternParts.model_fields[param]._attributes_set == PatternParts_fields[param]
            for param in PatternParts.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert PatternParts.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(PatternParts)
    # ---- Validate correct setting
    assert PatternParts.model_config == dict(title="region name pattern")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(PatternParts)
    # ---- Check default output
    with pytest.raises(ValidationError, match="2 validation errors for region name pattern"):
        assert PatternParts.judge(**dict())


@pytest.mark.parametrize(
    "description",
    ["Assess `StratificationFiles` pydantic model structure"],
    ids=["Assess `StratificationFiles` pydantic model structure"],
)
def test_StratificationFiles_model_structure(description, StratificationFiles_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(StratificationFiles_fields).issubset(StratificationFiles.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            StratificationFiles.model_fields[param]._attributes_set
            == StratificationFiles_fields[param]
            for param in StratificationFiles.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert StratificationFiles.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(StratificationFiles)
    # ---- Validate correct setting
    assert StratificationFiles.model_config == dict(title="stratification file inputs")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(StratificationFiles)
    # ---- Check default output
    with pytest.raises(ValidationError, match="2 validation errors for stratification file inputs"):
        assert StratificationFiles.judge(**dict())


@pytest.mark.parametrize(
    "description",
    ["Assess `StratifiedSurveyMeanParameters` pydantic model structure"],
    ids=["Assess `StratifiedSurveyMeanParameters` pydantic model structure"],
)
def test_StratifiedSurveyMeanParameters_model_structure(
    description, StratifiedSurveyMeanParameters_fields
):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(StratifiedSurveyMeanParameters_fields).issubset(
        StratifiedSurveyMeanParameters.model_fields
    )
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            StratifiedSurveyMeanParameters.model_fields[param]._attributes_set
            == StratifiedSurveyMeanParameters_fields[param]
            for param in StratifiedSurveyMeanParameters.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert StratifiedSurveyMeanParameters.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(StratifiedSurveyMeanParameters)
    # ---- Validate correct setting
    assert StratifiedSurveyMeanParameters.model_config == dict(
        arbitrary_types_allowed=True, title="stratified survey parameters"
    )

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(StratifiedSurveyMeanParameters)
    # ---- Check default output
    with pytest.raises(
        ValidationError, match="3 validation errors for stratified survey parameters"
    ):
        assert StratifiedSurveyMeanParameters.judge(**dict())

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = StratifiedSurveyMeanParameters.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert set(["validate_posint", "validate_posfloat"]).issubset(
        StratifiedSurveyMeanParameters.__dict__
    )

    # -------------------------
    # ASSERT: 'valid_posfloat'
    # ---- Check output for 'validate_posfloat' [VALID]
    assert StratifiedSurveyMeanParameters.validate_posfloat(2) == 2.0
    # ---- Check output for 'validate_posfloat' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative float.")):
        assert StratifiedSurveyMeanParameters.validate_posfloat(-2.0)
    # ---- Check output for 'validate_float' [INVALID: NaN]
    assert StratifiedSurveyMeanParameters.validate_posfloat(np.inf) == np.inf

    # ---- Check the applicable fields
    assert field_decor["validate_posfloat"].info.fields == ("strata_transect_proportion",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_posfloat"].info.mode == "before"

    # -------------------------
    # ASSERT: 'valid_posint'
    # ---- Checkout output for 'validate_posfloat' [VALID]
    assert StratifiedSurveyMeanParameters.validate_posint(1) == 1
    # ---- Checkout output for 'validate_posfloat' [INVALID: float]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert StratifiedSurveyMeanParameters.validate_posint(1.0)
    # ---- Check output for 'validate_posint' [INVALID: negative]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert StratifiedSurveyMeanParameters.validate_posint(-1)
    # ---- Check output for 'validate_posint' [INVALID: NaN]
    with pytest.raises(ValueError, match=re.escape("Value must be a non-negative integer.")):
        assert StratifiedSurveyMeanParameters.validate_posint(np.nan)
    # ---- Check the applicable fields
    assert field_decor["validate_posint"].info.fields == (
        "num_replicates",
        "mesh_transects_per_latitude",
    )
    # ---- Check the applicable validation mode
    assert field_decor["validate_posint"].info.mode == "before"


@pytest.mark.parametrize(
    "description",
    ["Assess `TransectRegionMap` pydantic model structure"],
    ids=["Assess `TransectRegionMap` pydantic model structure"],
)
def test_TransectRegionMap_model_structure(description, TransectRegionMap_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(TransectRegionMap_fields).issubset(TransectRegionMap.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            TransectRegionMap.model_fields[param]._attributes_set == TransectRegionMap_fields[param]
            for param in TransectRegionMap.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert TransectRegionMap.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(TransectRegionMap)
    # ---- Validate correct setting
    assert TransectRegionMap.model_config == dict(
        arbitrary_types_allowed=True, title="transect-region mapping parameters"
    )

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(TransectRegionMap)
    # ---- Check default output
    with pytest.raises(
        ValidationError, match="2 validation errors for transect-region mapping parameters"
    ):
        assert TransectRegionMap.judge(**dict())

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = TransectRegionMap.__pydantic_decorators__.field_validators
    # ---- 'model_validators' 'Decorator' object
    model_decor = TransectRegionMap.__pydantic_decorators__.model_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert set(["validate_pattern"]).issubset(TransectRegionMap.__dict__)

    # -------------------------
    # ASSERT: 'validate_pattern'
    # ---- Check output for 'validate_pattern' [VALID]
    assert TransectRegionMap.validate_pattern("") == ""
    # ---- Check output for 'validate_pattern' [VALID]
    assert TransectRegionMap.validate_pattern("{HAUL_NUM}{COUNTRY}") == "{HAUL_NUM}{COUNTRY}"
    # ---- Check output for 'vvalidate_pattern' [VALID]
    assert TransectRegionMap.validate_pattern("{HAUL_NUM}") == "{HAUL_NUM}"
    # ---- Check output for 'validate_pattern' [INVALID: extra identifier]
    with pytest.raises(ValueError, match=re.escape("contains invalid identifiers")):
        assert TransectRegionMap.validate_pattern("{YEAR}{EXTRA}")
    # ---- Check the applicable fields
    assert field_decor["validate_pattern"].info.fields == ("pattern",)
    # ---- Check the applicable validation mode
    assert field_decor["validate_pattern"].info.mode == "after"

    # -------------------------
    # ASSERT: model validator decorators
    # ---- Check whether the field exists
    assert "validate_pattern_parts" in TransectRegionMap.__dict__
    # MOCK: parameterized model
    MOCK_VALUES = {
        "save_file_template": "mapping_{COUNTRY}_{YEAR}",
        "save_file_directory": "outer",
        "save_file_sheetname": "space",
        "pattern": "{YEAR}{COUNTRY}",
        "parts": {
            "YEAR": [
                {
                    "pattern": "[0-9]+",
                    "label": None,
                },
            ],
            "COUNTRY": [
                {
                    "pattern": "^[cC]*[ANDY]",
                    "label": "Candy",
                },
                {
                    "pattern": "^[lL]*[AND]",
                    "label": "Land",
                },
            ],
        },
    }
    # ASSERT: model validator
    assert TransectRegionMap.validate_pattern_parts(MOCK_VALUES) == MOCK_VALUES
    # ---- Check the applicable validation mode
    assert model_decor["validate_pattern_parts"].info.mode == "before"


@pytest.mark.parametrize(
    "description",
    ["Assess `TSLRegressionParameters` pydantic model structure"],
    ids=["Assess `TSLRegressionParameters` pydantic model structure"],
)
def test_TSLRegressionParameters_model_structure(description, TSLRegressionParameters_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(TSLRegressionParameters_fields).issubset(TSLRegressionParameters.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            TSLRegressionParameters.model_fields[param]._attributes_set
            == TSLRegressionParameters_fields[param]
            for param in TSLRegressionParameters.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert TSLRegressionParameters.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(TSLRegressionParameters)
    # ---- Validate correct setting
    assert TSLRegressionParameters.model_config == dict(title="TS-length regression parameters")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(TSLRegressionParameters)
    # ---- Check default output
    with pytest.raises(
        ValidationError, match="4 validation errors for TS-length regression parameters"
    ):
        assert TSLRegressionParameters.judge(**dict())


@pytest.mark.parametrize(
    "description",
    ["Assess `XLSXFile` pydantic model structure"],
    ids=["Assess `XLSXFile` pydantic model structure"],
)
def test_XLSXFile_model_structure(description, XLSXFile_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(XLSXFile_fields).issubset(XLSXFile.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            XLSXFile.model_fields[param]._attributes_set == XLSXFile_fields[param]
            for param in XLSXFile.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert XLSXFile.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(XLSXFile)
    # ---- Validate correct setting
    assert XLSXFile.model_config == dict(title="*.xlsx file tree")

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(XLSXFile)
    # ---- Check default output
    assert XLSXFile.judge(**dict(filename="dummy1", sheetname="dummy2")).model_dump() == {
        "filename": "dummy1",
        "sheetname": "dummy2",
    }


####################################################################################################
# HIGHER MODEL PARAMETERIZATION: TESTS
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            ValidationError,
        ),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            ValidationError,
        ),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            ValidationError,
        ),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
            },
            None,
        ),
        (
            {
                "length": {"A": {"filename": "blurgh", "sheetname": "sheet1"}},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "length": {
                    "A": {"filename": "blurgh", "sheetname": "sheet1"},
                    "B": {"filename": "blaurrgh", "sheetname": "sheet5"},
                },
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "length": {"A": {"filename": "blurgh", "sheetname": "sheet1"}},
                "specimen": {"A": {"filename": "blargh", "sheetname": "sheet2"}},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        (
            {
                "length": {
                    "A": {"filename": "blurgh", "sheetname": "sheet1"},
                    "B": {"filename": "blaurrgh", "sheetname": "sheet5"},
                },
                "specimen": {
                    "A": {"filename": "blargh", "sheetname": "sheet2"},
                    "B": {"filename": "bleurgh", "sheetname": "sheet6"},
                },
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
            },
            None,
        ),
        ({}, ValidationError),
        (
            {
                "length": {"filename": "blurgh", "sheetname": "sheet1"},
                "specimen": {"filename": "blargh", "sheetname": "sheet2"},
                "catch": {"filename": "blorgh", "sheetname": "sheet3"},
                "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `BiologicalFiles` (unnested dictionary)",
        "Missing key 'length'",
        "Missing 'specimen'",
        "Missing 'catch'",
        "Missing 'haul_to_transect'",
        "Single nested directory with one key (valid)",
        "Single nested directory with two keys (valid)",
        "Multiple nested directories with one key (valid)",
        "Multiple nested directories with two keys (valid)",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_BiologicalFiles(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert BiologicalFiles.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert BiologicalFiles.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(BiologicalFiles.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(BiologicalFiles.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"directory": "blurgh/blargh", "sheetname": "sheet1"}, None),
        ({"sheetname": "sheet1"}, ValidationError),
        ({"directory": "blurgh/blargh"}, ValidationError),
        ({}, ValidationError),
        ({"directory": None, "sheetname": "sheet1"}, ValidationError),
        ({"directory": "blurgh/blargh", "sheetname": None}, ValidationError),
        ({"directory": 1, "sheetname": "sheet1"}, ValidationError),
        ({"directory": "blurgh/blargh", "sheetname": 1}, ValidationError),
        ({"directory": "blurgh/blarg", "sheetname": "sheet1", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `FileSettings`",
        "Missing 'directory'",
        "Missing 'sheetname'",
        "Empty dictionary",
        "Directory key missing value",
        "Sheetname key missing value",
        "Directory value not a string",
        "Sheetname value not a string",
        "Excess keys (valid)",
    ],
)
def test_FileSettings(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert FileSettings.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert FileSettings.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(FileSettings.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(FileSettings(**input).model_dump(exclude_none=True))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"init": "epsg:1000"}, None),
        ({"init": "EPSG:1000"}, None),
        ({"init": "epsg1000"}, None),
        ({"init": "EPSG1000"}, None),
        ({"init": "1000"}, None),
        ({"init": 1000}, None),
        ({"init": 1000.5}, ValidationError),
        ({"init": "ABCD:1000"}, ValidationError),
        ({}, ValidationError),
        ({"init": None}, AttributeError),
        ({"init": "epsg:1000", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `Geospatial`",
        "Uppercase EPSG (valid)",
        "Missing ':' between 'epsg' and number code (valid)",
        "Missing ':' between 'EPSG' and number code (valid)",
        "Missing 'epsg' from number code (valid)",
        "Number code integer input (valid)",
        "Number code float input (invalid)",
        "Unknown coordinate system type (non-EPSG)",
        "Empty dictionary",
        "Init key missing value",
        "Excess keys (valid)",
    ],
)
def test_Geospatial(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert Geospatial.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert Geospatial.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(Geospatial.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(Geospatial.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            None,
        ),
        (
            {
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            ValidationError,
        ),
        (
            {"save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx", "country_code": ["A", "B"]},
            ValidationError,
        ),
        (
            {
                "save_file_template": None,
                "country_code": [["A", "B"]],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            TypeError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": None,
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            TypeError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {},
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            None,
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}_{ERRONEOUS}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            (TypeError, ValueError),
        ),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            },
            ValueError,
        ),
        ({}, ValidationError),
        (
            {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `HaulTransectMap`",
        "Missing key 'save_file_template'",
        "Missing key 'country_code'",
        "Missing key 'file_settings'",
        "Save_file_template key value missing",
        "Country_code key value missing",
        "File_settings key value missing",
        "Single filename ID in template (valid)",
        "Erroneous filename ID in template (invalid)",
        "Mismatched 'file_setting' keys and items in 'country_code' list",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_HaulTransectMap(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert HaulTransectMap.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert HaulTransectMap.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(HaulTransectMap.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(HaulTransectMap.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            None,
        ),
        (
            {
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            ValidationError,
        ),
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            ValidationError,
        ),
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            ValidationError,
        ),
        ({}, ValidationError),
        (
            {
                "vario_krig_para": {"filename": "blurgh", "sheetname": "sheet1"},
                "isobath_200m": {"filename": "blurgh", "sheetname": "sheet1"},
                "mesh": {"filename": "blurgh", "sheetname": "sheet1"},
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `KrigingFiles`",
        "Missing 'vario_krig_para'",
        "Missing 'isobath_200m'",
        "Missing 'mesh'",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_KrigingFiles(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert KrigingFiles.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert KrigingFiles.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(KrigingFiles.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(KrigingFiles.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
            },
            None,
        ),
        (
            {"longitude_reference": 0.0, "longitude_offset": -1.0, "latitude_offset": 1.0},
            ValidationError,
        ),
        ({"A0": 1.0, "longitude_offset": -1.0, "latitude_offset": 1.0}, ValidationError),
        ({"A0": 1.0, "longitude_reference": 0.0, "latitude_offset": 1.0}, ValidationError),
        ({"A0": 1.0, "longitude_reference": 0.0, "longitude_offset": -1.0}, ValidationError),
        (
            {
                "A0": None,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": None,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": None,
                "latitude_offset": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": None,
            },
            ValidationError,
        ),
        ({}, ValidationError),
        ({"A0": 1, "longitude_reference": 0, "longitude_offset": -1, "latitude_offset": 1}, None),
        (
            {"A0": -1.0, "longitude_reference": 0, "longitude_offset": -1, "latitude_offset": 1},
            ValidationError,
        ),
        (
            {
                "A0": 1.0,
                "longitude_reference": 0.0,
                "longitude_offset": -1.0,
                "latitude_offset": 1.0,
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `KrigingParameters`",
        "Missing 'A0'",
        "Missing 'longitude_reference'",
        "Missing 'longitude_offset'",
        "Missing 'latitude_offset'",
        "A0 key value missing",
        "Longitude_reference key value missing",
        "Longitude_offset key value missing",
        "Latitude_offset key value missing",
        "Empty dictionary",
        "All key values as integers (valid)",
        "Negative A0 (invalid)",
        "Excess keys (valid)",
    ],
)
def test_KrigingParameters(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert KrigingParameters.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert KrigingParameters.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(KrigingParameters.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(KrigingParameters.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            None,
        ),
        (
            {
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            None,
        ),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}_{ERRONEOUS}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            },
            ValidationError,
        ),
        ({}, ValidationError),
        (
            {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `NASCExports`",
        "Missing key 'save_file_template'",
        "Missing key 'nasc_export_directory'",
        "Missing key 'export_file_directory'",
        "Missing key 'save_file_template'",
        "Missing key 'save_file_sheetname'",
        "Missing key 'regions'",
        "Missing key 'file_columns'",
        "Single filename ID in template (valid)",
        "Erroneous filename ID in template (invalid)",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_NASCExports(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert NASCExports.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert NASCExports.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(NASCExports.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(NASCExports.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"pattern": "pretty", "label": "ponies"}, None),
        ({"pattern": "pretty"}, ValidationError),
        ({"label": "ponies"}, ValidationError),
        ({}, ValidationError),
        ({"pattern": None, "label": "ponies"}, ValidationError),
        ({"pattern": "pretty", "label": None}, ValidationError),
        ({"pattern": 1, "label": "ponies"}, ValidationError),
        ({"pattern": "pretty", "label": 1}, ValidationError),
        ({"pattern": "pretty", "label": "ponies", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `PatternParts`",
        "Missing 'label'",
        "Missing 'pattern'",
        "Empty dictionary",
        "Pattern key missing value",
        "Label key missing value",
        "filename value not a string",
        "Sheetname value not a string",
        "Excess keys (valid)",
    ],
)
def test_PatternParts(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert PatternParts.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert PatternParts.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(PatternParts.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(PatternParts.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "strata": {"filename": "blurgh", "sheetname": "sheet1"},
                "geo_strata": {"filename": "blurgh", "sheetname": "sheet1"},
            },
            None,
        ),
        ({"strata": {"filename": "blurgh", "sheetname": "sheet1"}}, ValidationError),
        ({"geo_strata": {"filename": "blurgh", "sheetname": "sheet1"}}, ValidationError),
        ({}, ValidationError),
        (
            {
                "strata": {"filename": "blurgh", "sheetname": "sheet1"},
                "geo_strata": {"filename": "blurgh", "sheetname": "sheet1"},
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `StratificationFiles`",
        "Missing 'geo_strata'",
        "Missing 'strata'",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_StratificationFiles(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert StratificationFiles.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert StratificationFiles.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(StratificationFiles.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(StratificationFiles.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            None,
        ),
        ({"num_replicates": 10, "mesh_transects_per_latitude": 1}, ValidationError),
        ({"strata_transect_proportion": 1.0, "mesh_transects_per_latitude": 1}, ValidationError),
        ({"strata_transect_proportion": 1.0, "num_replicates": 10}, ValidationError),
        ({}, ValidationError),
        (
            {
                "strata_transect_proportion": None,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": None,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": 10,
                "mesh_transects_per_latitude": None,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            None,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10.0,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1.0,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": "1",
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": "10",
                "mesh_transects_per_latitude": 1,
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1,
                "num_replicates": 10,
                "mesh_transects_per_latitude": "1",
            },
            ValidationError,
        ),
        (
            {
                "strata_transect_proportion": 1.0,
                "num_replicates": 10,
                "mesh_transects_per_latitude": 1,
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `StratifiedSurveyMeanParameters`",
        "Missing 'strata_transect_proportion'",
        "Missing 'num_replicates'",
        "Missing 'mesh_transects_per_latitude'",
        "Empty dictionary",
        "Strata_transect_proportion key missing value",
        "Num_replicates key missing value",
        "Mesh_transects_per_latitude key missing value",
        "Strata_transect_proportion as int (valid)",
        "Num_replicates as float (invalid)",
        "Mesh_transects_per_latitude as float (invalid)",
        "Strata_transect_proportion as str (invalid)",
        "Num_replicates as str (invalid)",
        "Mesh_transects_per_latitude as str (invalid)",
        "Excess keys (valid)",
    ],
)
def test_StratifiedSurveyMeanParameters(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert StratifiedSurveyMeanParameters.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert StratifiedSurveyMeanParameters.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(StratifiedSurveyMeanParameters.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(StratifiedSurveyMeanParameters.create(**input))) == {
                    "excess"
                }
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            None,
        ),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",
            },
            ValidationError,
        ),
        (
            {
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            ValidationError,
        ),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",
                "parts": {
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            ValidationError,
        ),
        (
            {
                "pattern": "{REGION_CLASS}",
                "parts": {"REGION_CLASS": [{"pattern": "a", "label": "A"}]},
            },
            None,
        ),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}{ERRONEOUS}",
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
            },
            ValidationError,
        ),
        ({}, ValidationError),
        (
            {
                "pattern": "{REGION_CLASS}{COUNTRY}{HAUL_NUM}",
                "parts": {
                    "REGION_CLASS": [{"pattern": "a", "label": "A"}],
                    "HAUL_NUM": [{"pattern": "b", "label": "B"}],
                    "COUNTRY": [{"pattern": "c", "label": "C"}],
                },
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `TransectRegionMap`",
        "Missing key 'parts'",
        "Missing key 'pattern'",
        "Mismatch between 'parts' keys and those in pattern",
        "Single pattern ID in template (valid)",
        "Parts components as a list and not a dictionary (invalid)",
        "Empty dictionary",
        "Excess keys (valid)",
    ],
)
def test_TransectRegionMap(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert TransectRegionMap.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert TransectRegionMap.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(TransectRegionMap.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(TransectRegionMap.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": "km"},
            None,
        ),
        ({"TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": "km"}, ValidationError),
        ({"number_code": 12345, "TS_L_intercept": 2.0, "length_units": "km"}, ValidationError),
        ({"number_code": 12345, "TS_L_slope": 1.0, "length_units": "km"}, ValidationError),
        ({"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0}, ValidationError),
        ({}, ValidationError),
        (
            {"number_code": None, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": "km"},
            ValidationError,
        ),
        (
            {"number_code": 12345, "TS_L_slope": None, "TS_L_intercept": 2.0, "length_units": "km"},
            ValidationError,
        ),
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": None, "length_units": "km"},
            ValidationError,
        ),
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": None},
            ValidationError,
        ),
        (
            {
                "number_code": "12345",
                "TS_L_slope": "1.0",
                "TS_L_intercept": "2.0",
                "length_units": "km",
            },
            None,
        ),
        (
            {"number_code": 12345, "TS_L_slope": 1.0, "TS_L_intercept": 2.0, "length_units": 1},
            ValidationError,
        ),
        (
            {"number_code": 12345.0, "TS_L_slope": 1, "TS_L_intercept": 2, "length_units": "km"},
            None,
        ),
        (
            {
                "number_code": 12345.0,
                "TS_L_slope": 1,
                "TS_L_intercept": 2,
                "length_units": "km",
                "excess": "erroneous",
            },
            None,
        ),
    ],
    ids=[
        "Valid `TSLRegressionParameters`",
        "Missing 'number_code'",
        "Missing 'TS_L_slope'",
        "Missing 'TS_L_intercept'",
        "Missing 'length_units'",
        "Empty dictionary",
        "Number_code key missing value",
        "TS_L_slope key missing value",
        "TS_L_intercept key missing value",
        "Length_units key missing value",
        "All values as strings (valid)",
        "Length_units as numeric (invalid)",
        "Number_code as float, TS_L_slope/TS_L_intercept as integers (valid)",
        "Excess keys (valid)",
    ],
)
def test_TSLRegressionParameters(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert TSLRegressionParameters.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert TSLRegressionParameters.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(TSLRegressionParameters.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(TSLRegressionParameters.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


@pytest.mark.parametrize(
    "input, exception",
    [
        ({"filename": "blargh", "sheetname": "sheet1"}, None),
        ({"sheetname": "sheet1"}, ValidationError),
        ({"filename": "blargh"}, ValidationError),
        ({}, ValidationError),
        ({"filename": None, "sheetname": "sheet1"}, ValidationError),
        ({"filename": "blargh", "sheetname": None}, ValidationError),
        ({"filename": 1, "sheetname": "sheet1"}, ValidationError),
        ({"filename": "blargh", "sheetname": 1}, ValidationError),
        ({"filename": "blargh", "sheetname": ["sheet1", "sheet2"]}, None),
        ({"filename": ["blargh", "blurgh"], "sheetname": "sheet1"}, ValidationError),
        ({"filename": "blargh", "sheetname": "sheet1", "excess": "erroneous"}, None),
    ],
    ids=[
        "Valid `XLSXFiles`",
        "Missing 'filename'",
        "Missing 'sheetname'",
        "Empty dictionary",
        "Filename key missing value",
        "Sheetname key missing value",
        "Filename value not a string",
        "Sheetname value not a string",
        "List of sheets (valid)",
        "List of filenames (invalid)",
        "Excess keys (valid)",
    ],
)
def test_XLSXFile(input, exception):

    # Test for exceptions
    if exception is not None:
        with pytest.raises(exception):
            assert XLSXFile.create(**input)
    # Assert valid entries
    else:
        # ---- Assert validity
        assert XLSXFile.create(**input)
        # ---- Comparison of result vs expected
        try:
            assert set(XLSXFile.create(**input)) == set(input)
        except AssertionError:
            try:
                assert (set(input) - set(XLSXFile.create(**input))) == {"excess"}
            except AssertionError as e:
                pytest.fail(f"Unexpected AssertionError: {e}")


####################################################################################################
# SUPER MODEL STRUCTURE: FIXTURES
# --------------------------------------------------------------------------------------------------
@pytest.fixture
def CONFIG_DATA_MODEL_fields() -> Dict[str, Any]:

    return {
        "survey_year": {
            "annotation": int,
            "frozen": None,
        },
        "biological": {
            "annotation": Union[BiologicalFile, BiologicalFiles],
            "frozen": None,
        },
        "stratification": {
            "annotation": StratificationFiles,
            "frozen": None,
        },
        "NASC": {
            "annotation": Dict[str, XLSXFile],
            "frozen": None,
        },
        "species": {
            "annotation": SpeciesDefinition,
            "frozen": None,
        },
        "kriging": {
            "annotation": KrigingFiles,
            "frozen": None,
        },
        "data_root_dir": {
            "annotation": Union[str, None],
            "default": None,
            "frozen": None,
        },
        "report_path": {
            "annotation": Union[str, None],
            "default": None,
            "frozen": None,
        },
        "CAN_haul_offset": {
            "annotation": Union[int, None],
            "default": None,
            "frozen": None,
        },
        "ship_id": {
            "annotation": Union[int, str, float, Dict[Any, Any], None],
            "default": None,
            "frozen": None,
        },
        "export_regions": {
            "annotation": Union[Dict[str, XLSXFile], None],
            "default": None,
            "frozen": None,
        },
    }


@pytest.fixture
def CONFIG_INIT_MODEL_fields() -> Dict[str, Any]:

    return {
        "stratified_survey_mean_parameters": {
            "annotation": StratifiedSurveyMeanParameters,
            "frozen": None,
        },
        "kriging_parameters": {
            "annotation": KrigingParameters,
            "frozen": None,
        },
        "bio_hake_age_bin": {
            "annotation": List[Union[posint, realposfloat]],
            "frozen": None,
        },
        "bio_hake_len_bin": {
            "annotation": List[Union[posint, realposfloat]],
            "frozen": None,
        },
        "TS_length_regression_parameters": {
            "annotation": Dict[str, TSLRegressionParameters],
            "frozen": None,
        },
        "geospatial": {
            "annotation": Geospatial,
            "frozen": None,
        },
        "nasc_exports": {
            "annotation": Optional[NASCExports],
            "default": None,
        },
        "haul_to_transect_mapping": {
            "annotation": Optional[HaulTransectMap],
            "default": None,
        },
        "transect_region_mapping": {
            "annotation": Optional[TransectRegionMap],
            "default": None,
        },
    }


####################################################################################################
# SUPER MODEL STRUCTURE: TESTS
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "description",
    ["Assess `CONFIG_DATA_MODEL` pydantic model structure"],
    ids=["Assess `CONFIG_DATA_MODEL` pydantic model structure"],
)
def test_CONFIG_DATA_MODEL_model_structure(description, CONFIG_DATA_MODEL_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(CONFIG_DATA_MODEL_fields).issubset(CONFIG_DATA_MODEL.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            CONFIG_DATA_MODEL.model_fields[param]._attributes_set == CONFIG_DATA_MODEL_fields[param]
            for param in CONFIG_DATA_MODEL.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert CONFIG_DATA_MODEL.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(CONFIG_DATA_MODEL)
    # ---- Validate correct setting
    assert CONFIG_DATA_MODEL.model_config == dict()

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(CONFIG_DATA_MODEL)
    # ---- Check default output
    # ASSERT: '__name__' attribute modification
    with pytest.raises(
        ValueError,
        match=re.escape("configured data files defined in dir/folder/file"),
    ):
        assert CONFIG_DATA_MODEL.judge(filename="dir/folder/file", **{})


@pytest.mark.parametrize(
    "description",
    ["Assess `CONFIG_INIT_MODEL` pydantic model structure"],
    ids=["Assess `CONFIG_INIT_MODEL` pydantic model structure"],
)
def test_CONFIG_INIT_MODEL_model_structure(description, CONFIG_INIT_MODEL_fields):

    # --------------------------
    # ASSERT: field annotations
    # ---- Check existence
    assert set(CONFIG_INIT_MODEL_fields).issubset(CONFIG_INIT_MODEL.model_fields)
    # ---- Check typing, defaults, and other expected attributes
    assert all(
        [
            CONFIG_INIT_MODEL.model_fields[param]._attributes_set == CONFIG_INIT_MODEL_fields[param]
            for param in CONFIG_INIT_MODEL.model_fields
        ]
    )

    # --------------------------
    # ASSERT: '__base__' class inheritance
    assert CONFIG_INIT_MODEL.__base__ == InputModel

    # --------------------------
    # ASSERT: 'model_config' parameterization
    # ---- Check existence
    assert "model_config" in dir(CONFIG_INIT_MODEL)
    # ---- Validate correct setting
    assert CONFIG_INIT_MODEL.model_config == dict(arbitrary_types_allowed=True)

    # --------------------------
    # MOCK: Dictionary for inputs
    MOCK_VALUES = {
        "stratified_survey_mean_parameters": {
            "strata_transect_proportion": 1.00,
            "num_replicates": 1,
            "mesh_transects_per_latitude": 1,
        },
        "kriging_parameters": {
            "A0": 1.0,
            "longitude_reference": 0.0,
            "longitude_offset": 0.0,
            "latitude_offset": 0.0,
        },
        "bio_hake_age_bin": [0, 1, 2],
        "bio_hake_len_bin": [0, 1, 2],
        "TS_length_regression_parameters": {
            "elephant": {
                "number_code": 1,
                "TS_L_slope": 1.0,
                "TS_L_intercept": 0.0,
                "length_units": "km",
            },
        },
        "geospatial": {
            "init": 4326,
        },
    }
    # ASSERT: Valid judgment
    assert CONFIG_INIT_MODEL.judge(filename="dir/folder/file", **MOCK_VALUES)

    # -------------------------
    # ASSERT: 'judge' method
    # ---- Check existence
    assert "judge" in dir(CONFIG_INIT_MODEL)
    # ---- Check default output
    # ASSERT: '__name__' attribute modification
    with pytest.raises(
        ValueError,
        match=re.escape("configured initialization parameters defined in dir/folder/file"),
    ):
        assert CONFIG_INIT_MODEL.judge(filename="dir/folder/file", **{})

    # --------------------------
    # EXTRACT: field validator decorator arguments
    # ---- 'field_validators' 'Decorator' object
    field_decor = CONFIG_INIT_MODEL.__pydantic_decorators__.field_validators

    # -------------------------
    # ASSERT: field validator decorators
    # ---- Check whether field decorators exist
    assert "validate_interval" in CONFIG_INIT_MODEL.__dict__

    # -------------------------
    # ASSERT: 'validate_intervalt'
    # ---- Check output for 'validate_interval' [VALID]
    assert CONFIG_INIT_MODEL.validate_interval([1.0, 2.0, 5]) == [1.0, 2.0, 5]
    # ---- Check output for 'validate_interval' [INVALID: negative]
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Invalid value detected within list. " "Every value must be a non-negative float."
        ),
    ):
        assert CONFIG_INIT_MODEL.validate_interval([1.0, -2.0, -5])
    # ---- Check output for 'validate_interval' [INVALID: NaN]
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Invalid value detected within list. " "Every value must be a non-negative real number."
        ),
    ):
        assert CONFIG_INIT_MODEL.validate_interval([1.0, 2.0, np.inf])
    # ---- Check output for 'validate_interval' [INVALID: < 3 inputs]
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Interval list must have a length of 3: ['starting_value', 'ending_value', 'number']."
        ),
    ):
        assert CONFIG_INIT_MODEL.validate_interval([1.0, 2.0])
    # ---- Check the applicable fields
    assert field_decor["validate_interval"].info.fields == (
        "bio_hake_age_bin",
        "bio_hake_len_bin",
    )
    # ---- Check the applicable validation mode
    assert field_decor["validate_interval"].info.mode == "before"


####################################################################################################
# CONFIGURATION: TEST
# --------------------------------------------------------------------------------------------------

test_path = {
    "CONFIG": Path("C:/Users/Brandyn/Documents/GitHub/echopop/echopop/test_data/config_files/")
}


@pytest.mark.parametrize(
    "description",
    ["Test `CONFIG_INIT_MODEL` and `CONFIG_DATA_MODEL` validation"],
    ids=["Test `CONFIG_INIT_MODEL` and `CONFIG_DATA_MODEL` validation"],
)
def test_config_pydantic(test_path, description):

    # -------------------------
    # Read in the initialization and file configuration
    # ---- Initialization
    init_config = yaml.safe_load(Path(test_path["CONFIG"] / "config_init.yml").read_text())
    # ---- Files
    files_config = yaml.safe_load(Path(test_path["CONFIG"] / "config_survey.yml").read_text())

    # -------------------------
    # [ TEST 1 ]: FIRST CHECK OF PRE-VALIDATED DATA - `CONFIG_INIT_MODEL`
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 2 ]: ADD 'transect_region_mapping' KEY
    # ---- Add key values
    init_config.update(
        {
            "transect_region_mapping": {
                "save_file_template": "blurgh_{COUNTRY}_{YEAR}_{GROUP}.xlsx",
                "save_file_directory": "blurgh/blargh",
                "pattern": "{REGION_CLASS}",
                "save_file_sheetname": "sheet1",
                "parts": {"REGION_CLASS": [{"pattern": "a", "label": "A"}]},
            }
        }
    )
    # ---- ASSERT
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 3 ]: ADD 'nasc_exports' KEY
    # ---- Add key values
    init_config.update(
        {
            "nasc_exports": {
                "save_file_template": "blurgh_{REGION}_{YEAR}_{GROUP}.xlsx",
                "nasc_export_directory": "blurgh/blargh",
                "export_file_directory": "inblurgh/inblargh",
                "save_file_sheetname": "sheet1",
                "regions": {"A": ["a", "b"], "B": ["c", "d"]},
                "max_transect_spacing": 1.0,
                "file_columns": ["dum1", "dum2", "dum3"],
            }
        }
    )
    # ---- ASSERT
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 4 ]: ADD 'haul_to_transect_mapping' KEY
    # ---- Add key values
    init_config.update(
        {
            "haul_to_transect_mapping": {
                "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
                "country_code": ["A", "B"],
                "file_settings": {
                    "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
                    "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
                },
            }
        }
    )
    # ---- ASSERT
    assert CONFIG_INIT_MODEL(test_path["CONFIG"] / "config_init.yml", **init_config)

    # -------------------------
    # [ TEST 5 ]: FIRST CHECK OF PRE-VALIDATED DATA - `CONFIG_DATA_MODEL`
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 6 ]: ADD 'ship_id' [int] KEY
    # ---- Add key values
    files_config.update({"ship_id": 1})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 7 ]: ADD 'ship_id' [float] KEY
    # ---- Add key values
    files_config.update({"ship_id": 1.1})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 8 ]: ADD 'ship_id' [str] KEY
    # ---- Add key values
    files_config.update({"ship_id": "1"})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)

    # -------------------------
    # [ TEST 9 ]: ADD 'export_regions' KEY
    # ---- Add key values
    files_config.update({"export_regions": {"A": {"filename": "blee", "sheetname": "blah"}}})
    # ---- ASSERT
    assert CONFIG_DATA_MODEL(test_path["CONFIG"] / "config_survey.yml", **files_config)
