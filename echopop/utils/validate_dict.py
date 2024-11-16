import re
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
from pydantic import BaseModel, Field, RootModel, ValidationError, field_validator, model_validator

from .validate import posfloat, posint, realcircle, realposfloat

####################################################################################################
# PYDANTIC VALIDATORS
# --------------------------------------------------------------------------------------------------


class InputModel(BaseModel):
    """
    Base Pydantic model for scrutinizing file inputs
    """

    # Validator method
    @classmethod
    def judge(cls, **kwargs):
        """
        Validator method
        """
        try:
            return cls(**kwargs)
        except ValidationError as e:
            e.__traceback__ = None
            raise e

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method

        Notes
        ----------
        This is for `pytest` testing.
        """

        return cls.judge(**kwargs).model_dump(exclude_none=True)


class XLSXFile(InputModel, title="*.xlsx file tree"):
    """
    .xlsx file tree structure

    Parameters
    ----------
    filename: str
        Filename (as a string).
    sheetname: Union[str, List[str]]
        Sheet name (or list of sheet names) of a *.xlsx file (as a string) that will be loaded.
    """

    filename: str
    sheetname: Union[str, List[str]]


class FileSettings(InputModel, title="parameter file settings"):
    """
    Parameter file settings

    Parameters
    ----------
    directory: str
        File directory path (as a string).
    sheetname: str
        Sheet name of a *.xlsx file (as a string) that will be loaded.
    """

    directory: str
    sheetname: str


class StratifiedSurveyMeanParameters(
    InputModel, title="stratified survey parameters", arbitrary_types_allowed=True
):
    """
    Stratified sampling parameters

    Parameters
    ----------
    strata_transect_proportion: posfloat
        Proportion of transects sampled from each stratum.
    num_replicates: posint
        The number of replicates for the stratified analysis.
    mesh_transects_per_latitude: posint
        The number of synthetic/virtual transect lines generated per latitude interval.
    """

    strata_transect_proportion: posfloat
    num_replicates: posint
    mesh_transects_per_latitude: posint

    @field_validator("num_replicates", "mesh_transects_per_latitude", mode="before")
    def validate_posint(cls, v):
        return posint(v)

    @field_validator("strata_transect_proportion", mode="before")
    def validate_posfloat(cls, v):
        return posfloat(v)


class KrigingParameters(InputModel, arbitrary_types_allowed=True, title="kriging parameters"):
    """
    Kriging model parameters

    Parameters
    ----------
    A0: posfloat
        Reference area (nmi^2) of a kriged mesh cell.
    longitude_reference: float
        Longitude reference for kriging mesh adjustment.
    longitude_offset: float
        Longitudinal offset for kriging mesh adjustment.
    latitude_offset: float
        Latitudinal offset for kriging mesh adjustment.
    """

    A0: posfloat = Field(ge=0.0, allow_inf_nan=False)
    longitude_reference: float = Field(ge=-180.0, le=180.0, allow_inf_nan=False)
    longitude_offset: float = Field(allow_inf_nan=False)
    latitude_offset: float = Field(allow_inf_nan=False)

    @field_validator("A0", mode="before")
    def validate_posfloat(cls, v):
        return posfloat(v)


class HaulTransectMap(InputModel, arbitrary_types_allowed=True, title="haul-transect key mapping"):
    """
    Haul-to-transect key mapping generation parameters

    Parameters
    ----------
    save_file_template: str
        Save file template name.
    country_code: List[str]
        List of country names, abbreviations, or codes.
    file_settings: Dict[str, FileSettings]
        A dictionary comprising directory names and sheetnames of associated files.
    """

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
        # ---- Return values
        return v


class PatternParts(InputModel, title="region name pattern"):
    """
    String pattern parts

    Parameters
    ----------
    pattern: str
        Pattern used for parsing region names.
    label: str
        Resulting label for each export group/region.
    """

    pattern: str
    label: str


class TransectRegionMap(
    InputModel, arbitrary_types_allowed=True, title="transect-region mapping parameters"
):
    """
    Transect-to-region mapping parameters

    Parameters
    ----------
    save_file_template: str
        Save file template name.
    save_file_directory: str
        File directory name.
    save_file_sheetname: str
        File sheetname.
    pattern: str
        Region map code/pattern.
    parts: Dict[str, List[PatternParts]]
        Dictionary of metadata-pattern paired codes.
    """

    pattern: str
    parts: Dict[str, List[PatternParts]]

    @model_validator(mode="before")
    def validate_pattern_parts(cls, values):
        # ---- Get the country code list
        pattern_codes = values.get("pattern", "")
        # ---- Extract the codes
        codes = re.findall(r"{(.*?)}", pattern_codes)
        # ---- Get file settings keys
        parts_keys = list(values.get("parts", {}).keys())
        # ---- Keys within `file_settings` must match those defined in `country_code`
        if not set(parts_keys) == set(codes):
            # ---- Raise error
            raise ValueError(
                f"Defined pattern dictionary keys {parts_keys} must match those defined in "
                f"'pattern' ({codes})."
            )
        # ---- Return values
        return values

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
                f"Transect-to-region mapping save file template ({v}) contains invalid identifiers "
                f"({list(unknown_ids)}). Valid identifiers within the filename template (bounded "
                f"by curly braces) include: ['REGION_CLASS', 'HAUL_NUM', 'COUNTRY']."
            )
            # ---- Return value
        return v


class TSLRegressionParameters(InputModel, title="TS-length regression parameters"):
    """
    Target strength - length regression parameters

    Parameters
    ----------
    number_code: int
        Numeric species code.
    TS_L_slope: float
        TS-length regression slope.
    TS_L_intercept: float
        TS-length regression intercept.
    length_units: str
        Length units for the TS-length regression.
    """

    number_code: int
    TS_L_slope: float = Field(allow_inf_nan=False)
    TS_L_intercept: float = Field(allow_inf_nan=False)
    length_units: str


class Geospatial(InputModel, title="EPSG code"):
    """
    Geospatial parameters

    Parameters
    ----------
    init: str
        EPSG projection code.
    """

    init: str

    @field_validator("init", mode="before")
    def validate_init(cls, v):
        # ---- Convert to a string if read in as an integer
        if isinstance(v, (int, float)):
            v = str(v)
        # ---- Convert to lowercase
        v = v.lower()
        # ---- Mold the entry into the expected format that includes a preceding 'epsg:'
        if not v.startswith("epsg"):
            v = "epsg:" + v
        # ---- Ensure that the colon is present
        if ":" not in v:
            v = "epsg:" + v.split("epsg")[1]
        # ---- Evaluate whether the pre-validator succeeded in finding an acceptable format
        if not re.match(r"^epsg:\d+$", v):
            raise ValueError(
                f"Echopop cannot parse the defined EPSG code ('{v}'). EPSG codes most be formatted "
                f"with strings beginning with 'epsg:' followed by the integer number code (e.g. "
                f"'epsg:4326')."
            )
        # ---- Return the pre-validated entry
        return v


class NASCExports(
    InputModel, arbitrary_types_allowed=True, title="Echoview export processing parameters"
):
    """
    NASC export processing parameters

    Parameters
    ----------
    export_file_directory: str
        Export file directory name.
    nasc_export_directory: str
        Directory name where Echoview export files are located.
    save_file_template: str
        Save file template name.
    save_file_sheetname: str
        File sheetname.
    regions: Dict[str, List[str]]
        Acoustic data region names (list or a single string).
    max_transect_spacing: realposfloat
        Maximum transect spacing (nmi).
    file_columns: List[str]
        File column names included in the final consolidated export *.xlsx file.
    """

    export_file_directory: str
    nasc_export_directory: str
    save_file_template: str
    save_file_sheetname: str
    regions: Dict[str, List[str]]
    max_transect_spacing: realposfloat = Field(ge=0.0, allow_inf_nan=False)
    file_columns: List[str]

    @field_validator("max_transect_spacing", mode="before")
    def validate_realposfloat(cls, v):
        return realposfloat(v)

    @field_validator("save_file_template", mode="after")
    def validate_save_file_template(cls, v):
        # ---- Find all strings contained within curly braces
        template_ids = re.findall(r"{(.*?)}", v)
        # ---- Evaluate valid id's
        if not set(template_ids).issubset(set(["REGION", "YEAR", "GROUP"])):
            # ---- Get the unknown IDs
            unknown_ids = set(template_ids) - set(["REGION", "YEAR", "GROUP"])
            # ---- Raise Error
            raise ValueError(
                f"Haul-to-transect mapping save file template ({v}) contains invalid identifiers "
                f"({list(unknown_ids)}). Valid identifiers within the filename template (bounded "
                f"by curly braces) include: ['YEAR', 'REGION', 'GROUP']."
            )
        # ---- Return values
        return v


class CONFIG_INIT_MODEL(InputModel, arbitrary_types_allowed=True):
    """
    Initialization parameter configuration YAML validator
    """

    stratified_survey_mean_parameters: StratifiedSurveyMeanParameters
    kriging_parameters: KrigingParameters
    bio_hake_age_bin: List[Union[posint, realposfloat]]
    bio_hake_len_bin: List[Union[posint, realposfloat]]
    TS_length_regression_parameters: Dict[str, TSLRegressionParameters]
    geospatial: Geospatial
    nasc_exports: Optional[NASCExports] = Field(default=None)
    haul_to_transect_mapping: Optional[HaulTransectMap] = Field(default=None)
    transect_region_mapping: Optional[TransectRegionMap] = Field(default=None)

    def __init__(self, filename, **kwargs):
        try:
            super().__init__(**kwargs)
        except ValidationError as e:
            # Customize error message
            new_message = str(e).replace(
                self.__class__.__name__,
                f"configured initialization parameters defined in {filename}",
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


class BiologicalFiles(InputModel, title="biological file inputs"):
    """
    Biological data files

    Parameters
    ----------
    length: Union[Dict[str, XLSXFile], XLSXFile]
        An *.xlsx file (or dictionary of files) containing binned length data.
    specimen: Union[Dict[str, XLSXFile], XLSXFile]
        An *.xlsx file (or dictionary of files) containing specimen biodata.
    catch: Union[Dict[str, XLSXFile], XLSXFile]
        An *.xlsx file (or dictionary of files) containing catch/haul biological data.
    haul_to_transect: Union[Dict[str, XLSXFile], XLSXFile]
        An *.xlsx file (or dictionary of files) containing haul-transect mapping data.
    """

    length: Union[Dict[str, XLSXFile], XLSXFile]
    specimen: Union[Dict[str, XLSXFile], XLSXFile]
    catch: Union[Dict[str, XLSXFile], XLSXFile]
    haul_to_transect: Optional[Union[Dict[str, XLSXFile], XLSXFile]] = Field(default=None)


class KrigingFiles(InputModel, title="kriging file inputs"):
    """
    Kriging data files

    Parameters
    ----------
    isobath_200m: XLSXFile
        An *.xlsx file (or dictionary of files) containing 200 m isobath coordinates.
    mesh: XSLXFile
        An *.xlsx file (or dictionary of files) containing kriging mesh node coordinates.
    vario_krig_para: XSLXFile
        An *.xlsx file (or dictionary of files) containing kriging and variogram parameter values.
    """

    vario_krig_para: XLSXFile
    isobath_200m: XLSXFile
    mesh: XLSXFile


class StratificationFiles(InputModel, title="stratification file inputs"):
    """
    Stratification data files

    Parameters
    ----------
    geo_strata: XSLXFile
        An *.xlsx file (or dictionary of files) containing geographically defined strata.
    strata: XLSXFile
        An *.xlsx file (or dictionary of files) containing length-based (e.g. KS) strata
        information.
    """

    strata: XLSXFile
    geo_strata: XLSXFile


class SpeciesDefinition(BaseModel):
    """
    Species definitions
    """

    text_code: Optional[str]
    number_code: Optional[Union[int, float]]


class CONFIG_DATA_MODEL(InputModel):
    """
    Data file configuration YAML validator
    """

    survey_year: int
    biological: BiologicalFiles
    stratification: StratificationFiles
    NASC: Dict[str, XLSXFile]
    species: SpeciesDefinition
    kriging: KrigingFiles
    data_root_dir: Optional[str] = None
    CAN_haul_offset: Optional[int] = None
    ship_id: Optional[Union[int, str, float]] = None
    export_regions: Optional[Dict[str, XLSXFile]] = None

    def __init__(self, filename, **kwargs):
        try:
            super().__init__(**kwargs)
        except ValidationError as e:
            # Customize error message
            new_message = str(e).replace(
                self.__class__.__name__, f"configured data files defined in {filename}"
            )
            raise ValueError(new_message) from e


class VariogramModel(BaseModel, arbitrary_types_allowed=True):
    """
    Base Pydantic model for variogram analysis inputs
    """

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """
        try:
            return cls(**kwargs).model_dump(exclude_none=True)
        except ValidationError as e:
            e.__traceback__ = None
            raise e


class VariogramEmpirical(VariogramModel, arbitrary_types_allowed=True):
    """
    Empirical variogram parameters

    Parameters
    ----------
    azimuth_range: realcircle
        The total azimuth angle range that is allowed for constraining the relative angles between ]
        spatial points, particularly for cases where a high degree of directionality is assumed
    force_lag_zero: bool
        See the `variogram_parameters` argument in
        :fun:`echopop.spatial.variogram.empirical_variogram` for more details on
        `force_lag_zero`
    standardize_coordinates: bool
        When set to `True`, transect coordinates are standardized using reference coordinates

    Returns
    ----------
    VariogramEmpirical: A validated dictionary with the user-defined empirical variogram
    parameter values and default values for any missing parameters/keys.
    """

    azimuth_range: realcircle = Field(default=360.0, ge=0.0, le=360.0, allow_inf_nan=False)
    force_lag_zero: bool = Field(default=True)
    standardize_coordinates: bool = Field(default=True)

    @field_validator("azimuth_range", mode="before")
    def validate_realcircle(cls, v):
        return realcircle(v)

    def __init__(
        self,
        azimuth_range: realcircle = 360.0,
        force_lag_zero: bool = True,
        standardize_coordinates: bool = True,
        **kwargs,
    ):
        """
        Empirical variogram processing parameters
        """

        try:
            super().__init__(
                azimuth_range=azimuth_range,
                force_lag_zero=force_lag_zero,
                standardize_coordinates=standardize_coordinates,
            )
        except ValidationError as e:
            # Drop traceback
            e.__traceback__ = None
            raise e


class VariogramOptimize(VariogramModel, arbitrary_types_allowed=True):
    """
    Variogram optimization (non-linear least squares) parameters

    Parameters
    ----------
    max_fun_evaluations: posint
        The maximum number of evaluations. Defaults to 500.
    cost_fun_tolerance: realposfloat
        Threshold used for determining convergence via incremental changes of the cost function.
        Defaults to 1e-6.
    solution_tolerance: realposfloat
        Threshold used for determining convergence via change of the independent variables.
        Defaults to 1e-8.
    gradient_tolerance: realposfloat
        Threshold used for determining convergence via the gradient norma. Defaults to 1e-8.
    finite_step_size: realposfloat
        The relative step sizes used for approximating the Jacobian via finite differences.
    trust_region_solver: Literal["exact", "base"]
        The method used for solving the trust-region problem by either using the Jacobian
        computed from the first iteration (`"base"`) or via singular value decomposition
        (`"exact"`). Defaults to "exact".
    x_scale: Union[Literal["jacobian"], np.ndarray[realposfloat]]
        When `x_scale="jacobian"`, the characteristic scale is updated across numerical
        iterations via the inverse norms of the Jacobian matrix. Otherwise, a `np.ndarray`
        of the same length as `fit_parameters` can provide a constant scaling factor.
    jacobian_approx: Literal["forward", "central"]
        Indicates whether forward differencing (`"forward"`) or central differencing
        (`"central"`) should be used to approximate the Jacobian matrix.

    Returns
    ----------
    VariogramOptimize: A validated dictionary with the user-defined variogram optimizatization
    parameter values and default values for any missing parameters/keys.
    """

    max_fun_evaluations: posint = Field(default=500, gt=0, allow_inf_nan=False)
    cost_fun_tolerance: realposfloat = Field(default=1e-6, gt=0.0, allow_inf_nan=False)
    gradient_tolerance: realposfloat = Field(default=1e-4, gt=0.0, allow_inf_nan=False)
    solution_tolerance: realposfloat = Field(default=1e-6, gt=0.0, allow_inf_nan=False)
    finite_step_size: realposfloat = Field(default=1e-8, gt=0.0, allow_inf_nan=False)
    trust_region_solver: Literal["base", "exact"] = Field(default="exact")
    x_scale: Union[Literal["jacobian"], np.ndarray[realposfloat]] = Field(default="jacobian")
    jacobian_approx: Literal["forward", "central"] = Field(default="central")

    @field_validator("max_fun_evaluations", mode="before")
    def validate_posint(cls, v):
        return posint(v)

    @field_validator(
        "cost_fun_tolerance",
        "gradient_tolerance",
        "finite_step_size",
        "solution_tolerance",
        mode="before",
    )
    def validate_realposfloat(cls, v):
        return realposfloat(v)

    @field_validator("x_scale", mode="before")
    def validate_xscale(cls, v):
        # Validate `np.ndarray[realposfloat]` case
        if isinstance(v, np.ndarray):
            # ---- Coerce values to a float
            v_float = [float(x) for x in v]
            # ---- Coerce to 'realposfloat', or raise Error
            try:
                v = np.array([realposfloat(x) for x in v_float])
            except ValueError as e:
                e.__traceback__ = None
                raise e
        # Validate `Literal['jacobian']` case
        elif isinstance(v, (int, float, str)) and v != "jacobian":
            raise ValueError(
                "Input should be either the Literal 'jacobian' or a NumPy array of real "
                "positive-only float values."
            )
        # Return 'v'
        return v


class VariogramBase(VariogramModel, arbitrary_types_allowed=True):
    """
    Base variogram model parameters

    Parameters
    ----------
    model: Union[str, List[str]]
        A string or list of model names. A single name represents a single family model. Two
        inputs represent the desired composite model (e.g. the composite J-Bessel and
        exponential model). Defaults to: ['bessel', 'exponential']. Available models and their
        required arguments can be reviewed in the :fun:`echopop.spatial.variogram.variogram`
        function.
    n_lags: posint
        See the `variogram_parameters` argument in
        :fun:`echopop.spatial.variogram.empirical_variogram` for more details on
        `n_lags`.
    sill: Optional[realposfloat]
        See the description of `sill` in
        :fun:`echopop.spatial.variogram.variogram`.
    nugget: Optional[realposfloat]
        See the description of `nugget` in
        :fun:`echopop.spatial.variogram.variogram`.
    correlation_range: Optional[realposfloat]
        See the description of `correlation_range` in
        :fun:`echopop.spatial.variogram.variogram`.
    hole_effect_range: Optional[realposfloat]
        See the description of `hole_effect_range` in
        :fun:`echopop.spatial.variogram.variogram`.
    decay_power: Optional[realposfloat]
        See the description of `decay_power` in
        :fun:`echopop.spatial.variogram.variogram`.
    enhanced_semivariance: Optional[bool]
        See the description of `enhanced_semivariance` in
        :fun:`echopop.spatial.variogram.variogram`.

    Returns
    ----------
    VariogramBase: A validated dictionary with the user-defined variogram parameter values and
    default values for any missing parameters/keys.
    """

    model: Union[str, List[str]] = Field(union_mode="left_to_right")
    n_lags: posint = Field(ge=1, allow_inf_nan=False)
    lag_resolution: Optional[realposfloat] = Field(default=None, gt=0.0, allow_inf_nan=False)
    sill: Optional[realposfloat] = Field(default=None, ge=0.0, allow_inf_nan=False)
    nugget: Optional[realposfloat] = Field(default=None, ge=0.0, allow_inf_nan=False)
    hole_effect_range: Optional[realposfloat] = Field(default=None, ge=0.0, allow_inf_nan=False)
    correlation_range: Optional[realposfloat] = Field(default=None, ge=0.0, allow_inf_nan=False)
    enhance_semivariance: Optional[bool] = Field(default=None)
    decay_power: Optional[realposfloat] = Field(default=None, ge=0.0, allow_inf_nan=False)

    @field_validator("n_lags", mode="before")
    def validate_posint(cls, v):
        return posint(v)

    @field_validator(
        "lag_resolution",
        "sill",
        "nugget",
        "hole_effect_range",
        "correlation_range",
        "decay_power",
        mode="before",
    )
    def validate_realposfloat(cls, v):
        return realposfloat(v)


class InitialValues(VariogramModel, arbitrary_types_allowed=True):
    min: Optional[realposfloat] = Field(default=None)
    value: realposfloat = Field(default=0.0, allow_inf_nan=False)
    max: Optional[posfloat] = Field(default=None)
    vary: bool = Field(default=False)
    """
    Variogram optimization initial values and ranges

    Parameters
    ----------
    min: Optional[realposfloat]
        Minimum value allowed during optimization.
    value: Optional[realposfloat]
        Starting value used for optimization.
    max: Optional[posfloat]
        Maximum value (including infinity) allowed during optimization.
    vary: Optional[bool]
        Boolean value dictating whether a particular parameter will be adjusted
        during optimization ['True'] or held constant ['False']
    """

    @field_validator("min", "value", mode="before")
    def validate_realposfloat(cls, v):
        return realposfloat(v)

    @field_validator("max", mode="before")
    def validate_posfloat(cls, v):
        return posfloat(v)

    @model_validator(mode="after")
    def validate_value_sort(self):

        # Check whether the 'min' and 'max' keys exist
        min = getattr(self, "min", None)
        max = getattr(self, "max", None)
        value = getattr(self, "value", None)
        vary = getattr(self, "vary", None)
        # ---- Group into a dictionary
        value_dict = dict(min=min, value=value, max=max, vary=vary)

        # Evaluate case where 'vary = False'
        if not vary:
            updated_value_dict = InitialValues._DEFAULT_STRUCTURE_EMPTY({"value": value})
            # ---- Update 'min', 'value', 'max'
            self.min = None
            self.max = None
            self.value = updated_value_dict["value"]

        # Evaluate case where 'vary = True' but 'min' and 'max' are not found
        if vary:
            updated_value_dict = InitialValues._DEFAULT_STRUCTURE_OPTIMIZE(
                {k: v for k, v in value_dict.items() if v is not None}
            )
            # ---- Update 'min', 'value', 'max'
            self.min = updated_value_dict["min"]
            self.value = updated_value_dict["value"]
            self.max = updated_value_dict["max"]

            # Ensure valid min-value-max grouping
            # ---- Only apply if parameters are being varied
            if not (self.min <= self.value <= self.max) and self.vary:
                # ---- Raise Error
                raise ValueError(
                    "Optimization minimum, starting, and maximum values  must satisfy the logic: "
                    "`min` <= value` <= `max`."
                )

        # Return values
        return self

    @classmethod
    def _DEFAULT_STRUCTURE_EMPTY(cls, input: Dict[str, Any] = {}):
        return {**dict(vary=False), **input}

    @classmethod
    def _DEFAULT_STRUCTURE_OPTIMIZE(cls, input: Dict[str, Any] = {}):
        return {**dict(min=0.0, value=0.0, max=np.inf, vary=True), **input}


class VariogramInitial(RootModel[InitialValues]):
    root: Dict[str, InitialValues]

    @model_validator(mode="before")
    @classmethod
    def validate_model_params(cls, v):

        # Get valid parameters
        valid_param_names = cls._VALID_PARAMETERS()

        # Compare names of input versus those that are accepted
        if not set(v).issubset(valid_param_names):
            # ---- Get unexpected parameter names
            unexpected = set(v) - set(valid_param_names)
            # ---- Raise Error
            raise ValueError(
                f"Unexpected optimization parameters: {list(unexpected)}. Only the "
                f"following variogram parameters are valid for optimization: "
                f"{valid_param_names}."
            )

        # Return values otherwise
        return v

    @classmethod
    def _VALID_PARAMETERS(cls):
        return ["correlation_range", "decay_power", "hole_effect_range", "nugget", "sill"]

    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """
        try:
            return cls(**kwargs).model_dump(exclude_none=True)
        except ValidationError as e:
            e.__traceback__ = None
            raise e


class MeshCrop(
    BaseModel,
    arbitrary_types_allowed=True,
    title="kriging mesh cropping parameters ('cropping_parameters')",
):
    crop_method: Literal["transect_ends", "convex_hull"] = Field(default="transect_ends")
    num_nearest_transects: posint = Field(gt=0, default=4)
    mesh_buffer_distance: realposfloat = Field(gt=0.0, default=1.25, allow_inf_nan=False)
    latitude_resolution: realposfloat = Field(gt=0.0, default=1.25, allow_inf_nan=False)
    bearing_tolerance: realcircle = Field(gt=0.0, default=15.0, le=180.0, allow_inf_nan=False)

    @field_validator("num_nearest_transects", mode="before")
    def validate_posint(cls, v):
        return posint(v)

    @field_validator("bearing_tolerance", mode="before")
    def validate_realcircle(cls, v):
        return realcircle(v)

    @field_validator("mesh_buffer_distance", "latitude_resolution", mode="before")
    def validate_realposfloat(cls, v):
        return realposfloat(v)

    def __init__(
        self,
        crop_method: Literal["transect_ends", "convex_hull"] = "transect_ends",
        num_nearest_transects: posint = 4,
        mesh_buffer_distance: realposfloat = 1.25,
        latitude_resolution: realposfloat = 1.25,
        bearing_tolerance: realcircle = 15.0,
        **kwargs,
    ):
        """
        Mesh cropping method parameters
        """

        try:
            super().__init__(
                crop_method=crop_method,
                num_nearest_transects=num_nearest_transects,
                mesh_buffer_distance=mesh_buffer_distance,
                latitude_resolution=latitude_resolution,
                bearing_tolerance=bearing_tolerance,
            )
        except ValidationError as e:
            # Drop traceback
            e.__traceback__ = None
            raise e

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method to create a `MeshCrop` instance
        """
        return cls(**kwargs).model_dump(exclude_none=True)


class KrigingParameterInputs(
    BaseModel, arbitrary_types_allowed=True, title="kriging model parameters ('kriging_parameters')"
):
    anisotropy: realposfloat = Field(default=1e-3, allow_inf_nan=False)
    kmin: posint = Field(default=3, ge=3)
    kmax: posint = Field(default=10, ge=3)
    correlation_range: Optional[realposfloat] = Field(default=None, gt=0.0, allow_inf_nan=False)
    search_radius: Optional[realposfloat] = Field(default=None, gt=0.0, allow_inf_nan=False)

    @field_validator("kmin", "kmax", mode="before")
    def validate_posint(cls, v):
        return posint(v)

    @field_validator("anisotropy", "correlation_range", "search_radius", mode="before")
    def validate_realposfloat(cls, v):
        if v is None:
            return v
        else:
            return realposfloat(v)

    @model_validator(mode="before")
    def validate_k_window(cls, values):
        # ---- Get `kmin`
        kmin = values.get("kmin", 3)
        # ---- Get 'kmax'
        kmax = values.get("kmax", 10)
        # ---- Ensure that `kmax >= kmin`
        if kmax < kmin:
            # ---- Raise Error
            raise ValueError(
                f"Defined 'kmax' ({kmax}) must be greater than or equal to 'kmin' ({kmin})."
            )
        # ---- Return values
        return values

    @model_validator(mode="before")
    def validate_spatial_correlation_params(cls, values):
        # ---- Get `correlation_range`
        correlation_range = values.get("correlation_range", None)
        # ---- Get 'search_radius'
        search_radius = values.get("search_radius", None)
        # ---- Ensure that both parameters are not None
        if not correlation_range and not search_radius:
            # ---- Raise Error
            raise ValueError(
                "Both 'correlation_range' and 'search_radius' arguments are missing. At least one "
                "must be defined."
            )
        # ---- Return values
        return values

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method to create a `KrigingParameters` instance
        """

        # Collect errors, if any arise
        try:
            # ---- Test validate
            _ = cls(**kwargs)
            # ---- Edit values if needed
            if kwargs.get("search_radius") is None and kwargs["correlation_range"] is not None:
                kwargs["search_radius"] = kwargs["correlation_range"] * 3
            # ---- Produce the dictionary as an output
            return cls(**kwargs).model_dump(exclude_none=True)
        except ValidationError as e:
            e.__traceback__ = None
            raise e


class KrigingAnalysis(BaseModel, arbitrary_types_allowed=True):
    best_fit_variogram: bool = Field(default=False)
    coordinate_transform: bool = Field(default=True)
    extrapolate: bool = Field(default=False)
    variable: Literal["biomass"] = Field(default="biomass")
    verbose: bool = Field(default=True)

    def __init__(self, **kwargs):
        try:
            super().__init__(**kwargs)
        except ValidationError as e:
            # Drop traceback
            e.__traceback__ = None
            raise e

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method to create a `KrigingAnalysis` instance
        """
        return cls(**kwargs).model_dump(exclude_none=True)
