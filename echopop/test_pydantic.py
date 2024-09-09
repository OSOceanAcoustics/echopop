import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Literal, Optional, Type, Union

import yaml
from pydantic import (
    BaseModel,
    Field,
    ValidationError,
    confloat,
    conint,
    conlist,
    field_validator,
    model_validator,
)

#
init_config_path = (
    "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
)
survey_year_config_path = (
    "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
)

# Load YAML files as dictionaries
init_config_params = yaml.safe_load(Path(init_config_path).read_text())
survey_year_config_params = yaml.safe_load(Path(survey_year_config_path).read_text())


class TSLengthRegressionParameters(BaseModel):
    number_code: int
    TS_L_slope: float
    TS_L_intercept: float
    length_units: str


class StratifiedSurveyMeanParameters(BaseModel):
    strata_transect_proportion: float
    num_replicates: int
    mesh_transects_per_latitude: int


class KrigingParameters(BaseModel):
    A0: float
    longitude_reference: float
    longitude_offset: float
    latitude_offset: float


class CONFIG_INIT_MODEL(BaseModel):
    stratified_survey_mean_parameters: StratifiedSurveyMeanParameters
    bio_hake_len_bin: List[confloat()]
    bio_hake_age_bin: List[confloat()]
    TS_length_regression_parameters: Dict[str, TSLengthRegressionParameters]
    geospatial: Dict[str, str]
    kriging_parameters: KrigingParameters


class BiologicalFile(BaseModel):
    filename: str
    sheetname: str


class Biological(BaseModel):
    length: Dict[str, BiologicalFile]
    specimen: Dict[str, BiologicalFile]
    catch: Dict[str, BiologicalFile]
    haul_to_transect: Dict[str, BiologicalFile]


class Stratification(BaseModel):
    strata: BiologicalFile
    geo_strata: Dict[str, List[str]]


class NASC(BaseModel):
    no_age1: BiologicalFile
    all_ages: BiologicalFile


class Kriging(BaseModel):
    mesh: BiologicalFile
    isobath_200m: BiologicalFile
    vario_krig_para: BiologicalFile


class CONFIG_DATA_MODEL(BaseModel):
    survey_year: int
    species: Dict[str, int]
    CAN_haul_offset: int
    data_root_dir: str
    biological: Biological
    stratification: Stratification
    NASC: NASC
    kriging: Kriging


from echopop.utils.validate import posfloat, posint, realcircle, realposfloat


class StratifiedSurveyMeanParameters(BaseModel, arbitrary_types_allowed=True):
    strata_transect_proportion: posfloat
    num_replicates: posint
    mesh_transects_per_latitude: posint

    @field_validator("num_replicates", "mesh_transects_per_latitude", mode="before")
    def validate_posint(cls, v):
        return posint(v)

    @field_validator("strata_transect_proportion", mode="before")
    def validate_posfloat(cls, v):
        return posfloat(v)


class KrigingParameters(BaseModel, arbitrary_types_allowed=True):
    A0: posfloat
    longitude_reference: float
    longitude_offset: float
    latitude_offset: float

    @field_validator("A0", mode="before")
    def validate_posfloat(cls, v):
        return posfloat(v)


class FileSettings(BaseModel):
    directory: str
    sheetname: str


class HaulTransectMap(BaseModel, arbitrary_types_allowed=True):
    save_file_template: str
    country_code: List[str]
    file_settings: Dict[str, FileSettings]

    @model_validator(mode="before")
    def validate_country_files(cls, values):
        # ---- Get the country code list
        country_codes = values.get("country_code", [])
        # ---- Get file settings keys
        file_settings_keys = list(values.get("file_settings", {}).keys())
        # ---- Keys within `file_settings` must match those defined in `country_code`
        if not set(file_settings_keys) == set(country_codes):
            # ---- Raise error
            raise ValueError(
                f"File settings keys {file_settings_keys} must match those defined in "
                f"'country_code' ({country_codes})."
            )
        # ---- Return values
        return values

    @field_validator("save_file_template", mode="after")
    def validate_save_file_template(cls, v):
        # ---- Find all strings contained within curly braces
        template_ids = re.findall(r"{(.*?)}", v)
        # ---- Evaluate valid id's
        if not set(template_ids).issubset(set(["YEAR", "COUNTRY"])):
            # ---- Get the unknown IDs
            unknown_ids = set(template_ids) - set(["YEAR", "COUNTRY"])
            # ---- Raise Error
            raise ValueError(
                f"Haul-to-transect mapping save file template ({v}) contains invalid identifiers "
                f"({list(unknown_ids)}). Valid identifiers within the filename template (bounded "
                f"by curly braces) include: ['YEAR', 'COUNTRY']."
            )


class PatternParts(BaseModel):
    pattern: str
    label: str


class TransectRegionMap(BaseModel, arbitrary_types_allowed=True):
    save_file_template: str
    save_file_directory: str
    save_file_sheetname: str
    pattern: str
    parts: List[PatternParts]

    @model_validator(mode="before")
    def validate_country_files(cls, values):
        # ---- Get the country code list
        country_codes = values.get("country_code", [])
        # ---- Get file settings keys
        file_settings_keys = list(values.get("file_settings", {}).keys())
        # ---- Keys within `file_settings` must match those defined in `country_code`
        if not set(file_settings_keys) == set(country_codes):
            # ---- Raise error
            raise ValueError(
                f"File settings keys {file_settings_keys} must match those defined in "
                f"'country_code' ({country_codes})."
            )
        # ---- Return values
        return values

    @field_validator("save_file_template", mode="after")
    def validate_save_file_template(cls, v):
        # ---- Find all strings contained within curly braces
        template_ids = re.findall(r"{(.*?)}", v)
        # ---- Evaluate valid id's
        if not set(template_ids).issubset(set(["YEAR", "COUNTRY", "GROUP"])):
            # ---- Get the unknown IDs
            unknown_ids = set(template_ids) - set(["YEAR", "COUNTRY", "GROUP"])
            # ---- Raise Error
            raise ValueError(
                f"Haul-to-transect mapping save file template ({v}) contains invalid identifiers "
                f"({list(unknown_ids)}). Valid identifiers within the filename template (bounded "
                f"by curly braces) include: ['YEAR', 'COUNTRY', 'GROUP']."
            )

    @field_validator("pattern", mode="after")
    def validate_pattern(cls, v):
        # ---- Find all strings contained within curly braces
        template_ids = re.findall(r"{(.*?)}", v)
        # ---- Evaluate valid id's
        if not set(template_ids).issubset(set(["REGION_CLASS", "HAUL_NUM", "COUNTRY"])):
            # ---- Get the unknown IDs
            unknown_ids = set(template_ids) - set(["REGION_CLASS", "HAUL_NUM", "COUNTRY"])
            # ---- Raise Error
            raise ValueError(
                f"Haul-to-transect mapping save file template ({v}) contains invalid identifiers "
                f"({list(unknown_ids)}). Valid identifiers within the filename template (bounded "
                f"by curly braces) include: ['REGION_CLASS', 'HAUL_NUM', 'GROUP']."
            )


FileSettings(**{"directory": "/Biological/US", "sheetname": "Sheet1"})
HaulTransectMap(
    init_config_path,
    **{
        "save_file_template": "A",
        "country_code": ["US", "CAN"],
        "file_settings": {
            "US": {"directory": "/Biological/US", "sheetname": "Sheet1"},
            "CAN": {"directory": "/Biological/CAN", "sheetname": "Sheet1"},
            "RUS": {"directory": "/Biological/CAN", "sheetname": "Sheet1"},
        },
    },
)


template = "haul_to_transect_mapping_{YEAR}_{COUNTRY}"
match = re.findall(r"{(.*?)}", template)


dir(cls)
values = {
    "save_file_template": "A",
    "country_code": ["US", "CAN"],
    "file_settings": {
        "US": {"direcdatory": "/Biological/US", "sheetname": "Sheet1"},
        "CAN": {"directory": "/Biological/CAN", "sheetname": "Sheet1"},
    },
}
cls.__config__.error_handler


class CONFIG_INIT_MODEL(BaseModel, arbitrary_types_allowed=True):
    stratified_survey_mean_parameters: StratifiedSurveyMeanParameters
    kriging_parameters: KrigingParameters
    bio_hake_age_bin: List[Union[posint, realposfloat]]
    bio_hake_len_bin: List[Union[posint, realposfloat]]
    haul_to_transect_mapping: Optional[HaulTransectMap] = None
    transect_region_mapping: Optional[TransectRegionMap] = None

    def __init__(self, filename, **kwargs):
        try:
            super().__init__(**kwargs)
        except ValidationError as e:
            # Customize error message
            new_message = str(e).replace(
                self.__class__.__name__, f"configuration parameters defined in {filename}"
            )
            raise ValueError(new_message) from e

    @field_validator("bio_hake_age_bin", "bio_hake_len_bin", mode="before")
    def validate_interval(cls, v):
        # ---- Check Union typing
        try:
            all(
                isinstance(value, (int, float))
                and (posint(value) if isinstance(value, int) else realposfloat(value))
                for value in v
            )
        except ValueError as e:
            raise ValueError(f"Invalid value detected within list. Every {str(e).lower()}")
        # ---- Check length
        if not len(v) == 3:
            raise ValueError(
                "Interval list must have a length of 3: "
                "['starting_value', 'ending_value', 'number']."
            )
        # ---- Check for any that may be 'realposfloat'
        any_posfloat = any(isinstance(value, (realposfloat, float)) for value in v)
        # ---- If true, then convert
        if any_posfloat:
            return [posint(value) if i == 2 else realposfloat(value) for i, value in enumerate(v)]
        else:
            return [posint(value) for value in v]


init_config_params["transect_region_mapping"]
CONFIG_INIT_MODEL(init_config_path, **init_config_params)
cls = CONFIG_INIT_MODEL
v = [1, 22, 1]
init_config_params.keys()

model_instance.model_validate()

try:
    CONFIG_INIT_MODEL(init_config_path, **init_config_params).model_validate()
except ValueError as e:
    print(e)
CONFIG_INIT_MODEL(
    init_config_path, **{**init_config_params, **{"dummy": {"A0": 6.25, "latitude_offset": -1}}}
)


class ExampleModel(BaseModel):
    value: Union[int, float]


# Example usage
try:
    ExampleModel(value=10)  # Valid, int
    ExampleModel(value=10.5)  # Valid, float
    ExampleModel(value="string")  # This will raise a validation error
except ValueError as e:
    print(e)

# Example usage
try:
    PairedKeys(key1=1, key2=2)  # This will work
    PairedKeys(key1=1)  # This will raise a validation error
except ValueError as e:
    print(e)


class Custom(BaseModel, arbitrary_types_allowed=True):
    strata_transect_proportion: posfloat = 1.00
    num_replicates: posint = 1
    mesh_transects_per_latitude: int = 5
    geospatial: str = "epsg:4326"
    dummy: Optional[str] = None
    dummy2: Optional[Literal["a", "b"]] = None

    @field_validator("num_replicates", "mesh_transects_per_latitude", mode="before")
    def validate_posint(cls, v):
        return posint(v)

    @field_validator("strata_transect_proportion", mode="before")
    def validate_posfloat(cls, v):
        return posfloat(v)

    def __init__(self, **kwargs):
        try:
            super().__init__(**kwargs)
        except ValidationError as e:
            # Customize error message
            new_message = str(e).replace(self.__class__.__name__, "Model parameter configuration")
            raise ValueError(new_message) from e

    @classmethod
    def create(cls, **kwargs):
        try:
            instance = cls(**kwargs)
            return instance.model_dump(exclude_none=True)
        except ValidationError as e:
            errors = str(e).replace(cls.__name__, "Model parameter configuration")
            raise ValueError(errors)


Custom.create(**{"strata_transect_proportion": 0.75, "num_replicates": 10, "dummy2": "b"})
Custom(**{"strata_transect_proportion": 0.75, "num_replicates": 10.1, "dummy2": "b"})


Custom.create(dummy="ed")
try:
    Custom.create({"num_replicates": "10", "dummy2": "c"})
except ValidationError as e:
    print(e)

KrigingParameters()
