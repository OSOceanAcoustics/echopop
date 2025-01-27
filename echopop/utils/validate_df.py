import warnings
from typing import Callable, List, Optional, Union

import numpy as np
import pandas as pd
from pandera import DataFrameModel, DataFrameSchema, Field, check, dataframe_check
from pandera.errors import SchemaError, SchemaErrors
from pandera.typing import Series


####################################################################################################
# PANDERA VALIDATORS
# --------------------------------------------------------------------------------------------------
class BaseDataFrame(DataFrameModel):

    class Config:
        strict = "filter"
        metadata = {}

    # Support attributes and methods for coercion/pre-processing
    _DTYPE_COERCION = {
        int: lambda v: v.astype(np.float64, errors="ignore").astype(np.int64, errors="ignore"),
        float: lambda v: v.astype(np.float64, errors="ignore"),
        str: lambda v: v.astype(str, errors="ignore"),
    }

    @classmethod
    def _drop_invalid_rows(
        cls, schema: DataFrameSchema, df: pd.DataFrame, filename: str
    ) -> pd.DataFrame:
        """Drop rows with invalid values"""

        # Create DataFrame copy
        df_copy = df.copy()

        # Get metadata attributes
        schema_meta = schema.get_metadata()[str(cls)]["columns"]

        # Initialize tracking list
        invalid_idx = []

        # Iterate through to find invalid indices
        for param, attr in schema_meta.items():
            # ---- If the metadata attribute is present
            if attr is not None:
                # ---- Only run if `drop_invalid_rows = True`
                if "drop_invalid_rows" in attr and not schema.columns[param].nullable:
                    # ---- Get the violating indices
                    try:
                        if schema.columns[param].regex:
                            # ---- Check for annotation
                            if "annotation" in attr:
                                param = attr["annotation"]
                            # ---- Create mask
                            idx_mask = df_copy.filter(regex=param).isna() | df_copy.filter(
                                regex=param
                            ).isin([np.inf, -np.inf])
                        else:
                            idx_mask = df_copy.filter([param]).isna() | df_copy.filter(
                                [param]
                            ).isin([np.inf, -np.inf])
                        # ---- Parse the indices
                        invalid_idx.extend(df_copy.index[idx_mask.to_numpy().flatten()].to_list())
                    except Exception:
                        pass

        # If there were violating indices, raise a warning to avoid silent data loss
        if len(invalid_idx) > 0:
            warnings.warn(
                f"The following row indices from {filename} contained invalid values and "
                f"were dropped from the validated DataFrame: {invalid_idx}",
                stacklevel=-1,
            )

        # Drop the offending indices and reset the index
        return df_copy.drop(invalid_idx, axis=0).reset_index(drop=True)

    @classmethod
    def _coerce(cls, schema: DataFrameSchema, df: pd.DataFrame) -> pd.DataFrame:
        """Coerce DataFrame datatypes"""

        # Create DataFrame copy
        df_copy = df.copy()

        # Get metadata attributes
        schema_meta = schema.get_metadata()[str(cls)]["columns"]

        # Iterate through to find invalid indices
        for param, attr in schema_meta.items():
            # ---- If the metadata attribute is present
            if attr is not None:
                # ---- Only run if `types` is present (for multi-type columns)
                if "types" in attr and schema.columns[param].coerce:
                    # ---- Check for annotation
                    if "annotation" in attr:
                        param = df_copy.filter(regex=param).columns
                    # ---- Try each ordered datatype
                    for dtype in attr["types"]:
                        try:
                            df_copy[param] = cls._DTYPE_COERCION.get(dtype, None)(df_copy[param])
                            break
                        except Exception:
                            pass

        # Return the coerced DataFrame
        return df_copy

    # Support methods for organizing error codes that are raised
    @classmethod
    def _missing_columns(cls, err: SchemaError) -> List[str]:
        """Extract missing column names from error"""
        return err.failure_cases["failure_case"].tolist()

    @classmethod
    def _inverse_missing_columns(cls, err: SchemaError) -> str:
        """Extract missing regex column names not found in DataFrame"""
        return (
            f"    - '{err.schema.metadata['annotation']}': "
            + f" No valid column name found matching the accepted regex pattern [{err.schema.name}]"
        )

    @classmethod
    def _wrong_dtype(cls, err: SchemaError) -> str:
        """Extract wrong datatype column names from error"""
        return f"    - '{err.schema.name}': {err.__str__().capitalize()}."

    @classmethod
    def _invalid_inequality(cls, err: SchemaError) -> str:
        """Extract invalid inequality column names from error"""
        return f"    - '{err.schema.name}': {err.__str__().capitalize()}."

    @classmethod
    def _mixed_dtypes(cls, err: SchemaError) -> str:
        """Extract mixed datatype column names from error"""
        return f"    - '{err.schema.name}': {err.check.error.capitalize()}."

    @classmethod
    def _nullable(cls, err: SchemaError) -> str:
        """Extract non-nullable column names from error"""
        return (
            f"    - '{err.schema.name}': "
            f"Column is non-nullable and contains invalid values (NaN/Inf)."
        )

    @classmethod
    def _coercion(cls, err: SchemaError) -> str:
        """Extract non-coercible column names from error"""

        # Get the target dtype
        target_dtype = f"'{err.schema._dtype.type.name}'"

        return (
            f"    - '{err.schema.name}': "
            f"Values within column are not coercible to the expected datatype {target_dtype}."
        )

    @classmethod
    def _type_error(cls, err: SchemaError) -> str:
        """Extract column names with a TypeError that has been raised"""
        return f"    - '{err.schema.name}': {err.failure_cases['failure_case'][0]}."

    @classmethod
    def _callable(cls, err: SchemaError, df: pd.DataFrame) -> str:
        """Extract column names from a generating error callable"""
        return f"    - {err.schema.checks[0].error(df)}."

    @classmethod
    def _generic_check(cls, err: SchemaError) -> str:
        """Extract general DATAFRAME and CHECK errors"""
        return f"    - {err.schema.checks[0].error}."

    @classmethod
    def _generic_misc(cls, err: SchemaError) -> str:
        """Extract miscellaneous errors"""
        return f"    - '{err.schema.name}': {err.check.error.capitalize()}."

    @classmethod
    def format_errors(cls, missing_cols: List[str], error_codes: List[str], filename: str) -> str:
        """Format error messages for DataFrame validation"""

        # Format missing columns
        if missing_cols:
            missing_cols_msg = (
                "\n"
                + f"    - The following columns were missing from the DataFrame: {missing_cols}."
            )
        else:
            missing_cols_msg = ""

        # Format the error codes
        if error_codes:
            # ---- Stack
            error_codes_stk = "\n".join(error_codes)
            # ---- Format message
            error_codes_msg = "\n" + f"{error_codes_stk}"
        else:
            error_codes_msg = ""

        # Finalize formatting
        message = (
            f"The following DataFrame validation errors were flagged for {filename}: "
            + f"{missing_cols_msg + error_codes_msg}"
        )

        # Return the error message
        return message

    @classmethod
    def process_errors(
        cls, errors: SchemaError, df: pd.DataFrame, filename: str
    ) -> Union[str, List[str]]:
        """Process errors for DataFrame validation"""

        # Set up lists that will be iteratively populated
        # ---- Generic error codes
        error_codes = []
        # ---- Missing columns (if any; column names will be collated into a single error message)
        missing_cols = []

        # Iterate through the SchemaErrors
        for err in errors.schema_errors:
            # ---- Missing columns
            if err.reason_code.name == "COLUMN_NOT_IN_DATAFRAME":
                missing_cols.extend(cls._missing_columns(err))
            # ---- Invalid column name (i.e. missing column -- generic)
            elif err.reason_code.name == "INVALID_COLUMN_NAME":
                error_codes.extend([cls._inverse_missing_columns(err)])
            # ---- Wrong datatype
            elif err.reason_code.name == "WRONG_DATATYPE":
                error_codes.extend([cls._wrong_dtype(err)])
            # ---- Coercion error
            elif err.reason_code.name == "DATATYPE_COERCION":
                error_codes.extend([cls._coercion(err)])
            # ---- Non-nullable error
            elif err.reason_code.name == "SERIES_CONTAINS_NULLS":
                error_codes.extend([cls._nullable(err)])
            # ---- General attritude-specific and dataframe-generalized checks
            elif err.reason_code.name in ["CHECK_ERROR", "DATAFRAME_CHECK"]:
                # ---- Error-generating function (generally: @dataframe_check)
                if isinstance(err.schema.checks[0].error, Callable):
                    error_codes.extend([cls._callable(err, df)])
                # ---- Inequality checks
                elif err.schema.checks[0].name in [
                    "less_than",
                    "greater_than",
                    "greater_than_or_equal_to",
                    "less_than_or_equal_to",
                ]:
                    # ---- Failure due to a more systematic TypeError issue
                    if "TypeError" in str(err.failure_cases["failure_case"][0]):
                        error_codes.extend([cls._type_error(err)])
                    # ---- General case
                    else:
                        error_codes.extend([cls._invalid_inequality(err)])
                # ---- Mixed datatype check
                elif err.schema.checks[0].name == "mixed_type":
                    error_codes.extend([cls._mixed_dtypes(err)])
                # ---- Generic output from @check and @dataframe_check decorators
                else:
                    error_codes.extend([cls._generic_check(err)])
            # ---- Catch-all miscellaneous errors that do not require specific handling
            else:
                error_codes.extend([cls._generic_misc(err)])

        # Format and raise the SchemaError
        raise SchemaError(cls, df, cls.format_errors(missing_cols, error_codes, filename)) from None

    # Validation methods
    @classmethod
    def pre_validate(cls, df: pd.DataFrame) -> pd.DataFrame:
        """Override in child classes for custom pre-validation"""
        return df.copy()

    @classmethod
    def validate_df(cls, df: pd.DataFrame, filename: str = None) -> pd.DataFrame:
        """Validate DataFrame"""

        # Create copy of DataFrame
        df = df.copy()

        # Exception handling
        # ---- Upon success
        try:
            # ---- Pre-validate the DataFrame (e.g. ambiguous column names)
            df_preproc = cls.pre_validate(df)
            # ---- Remove invalid rows (when nullable)
            df_null = cls._drop_invalid_rows(cls.to_schema(), df_preproc, filename)
            # ---- Coerce multi-type columns
            df_multi = cls._coerce(cls.to_schema(), df_null)
            # ---- Lazy validation of the DataFrame
            df_proc = cls.validate(df_multi, lazy=True)
            # ---- Return the validated DataFrame
            return df_proc
        # ---- Upon failure
        except SchemaErrors as e:
            # ---- Raise the error
            cls.process_errors(e, df, filename)


class LengthBiodata(BaseDataFrame):
    """
    Binned body length DataFrame

    Parameters
    ----------
    haul_num: Union[int, float]
        Haul number.
    length: float
        Body length.
    length_count: int
        Number of animals with the corresponding binned length.
    sex: Union[int, str]
        Animal sex. This can either be represented by an integer value, or a string corresponding
        to male ("male", "m"), female ("female", "f"), or unsexed/unidentified ("unsexed", "u").
    species_id: Union[int, float, str]
        Species identification code.
    """

    haul_num: Series = Field(
        nullable=False, metadata=dict(types=[int, float], drop_invalid_rows=True)
    )
    length: Series[float] = Field(gt=0.0, nullable=False, coerce=True)
    length_count: Series[int] = Field(ge=0, nullable=False, coerce=True)
    sex: Series = Field(nullable=False, metadata=dict(types=[int, str], drop_invalid_rows=True))
    species_id: Series = Field(
        nullable=False, metadata=dict(types=[int, float, str], drop_invalid_rows=True)
    )

    @check(
        "haul_num",
        name="mixed_type",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul_num(cls, haul_num: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in haul_num])

    @check(
        "sex",
        name="mixed_type",
        error=(
            "Column datatype should either be all 'int', or all 'str' contained within "
            "['male', 'm', 'female', 'f', 'unsexed', 'u']"
        ),
    )
    def validate_sex(cls, sex: Series) -> Series[bool]:
        # ---- Check if value is 'int'
        if all(isinstance(x, str) for x in sex):
            return Series([value in ["male", "female", "unsexed", "m", "f", "u"] for value in sex])
        else:
            return Series([isinstance(value, int) for value in sex])

    @check(
        "species_id",
        name="mixed_type",
        error="Column datatype should be either 'int', 'float', or 'str'",
    )
    def validate_species_id(cls, species_id: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float, str)) for value in species_id])


class CatchBiodata(BaseDataFrame):
    """
    Aggregated catch DataFrame

    Parameters
    ----------
    haul_num: Union[int, float]
        Haul number.
    haul_weight: float
        Total weight (kg) of a net haul.
    species_id: Union[int, float, str]
        Species identification code.
    """

    haul_num: Series = Field(
        nullable=False, metadata=dict(types=[int, float], drop_invalid_rows=True)
    )
    haul_weight: Series[float] = Field(
        ge=0.0, coerce=True, nullable=False, metadata=dict(drop_invalid_rows=True)
    )
    species_id: Series = Field(
        nullable=False, metadata=dict(types=[int, float, str], drop_invalid_rows=True)
    )

    @check(
        "haul_num",
        name="mixed_type",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul_num(cls, haul_num: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in haul_num])

    @check(
        "species_id",
        name="mixed_type",
        error="Column datatype should be either 'int', 'float', or 'str'",
    )
    def validate_species_id(cls, species_id: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float, str)) for value in species_id])


class SpecimenBiodata(BaseDataFrame):
    """
    Individual specimen length and weight DataFrame

    Parameters
    ----------
    age: Union[int, float]
        Animal age (in whole or decimal years).
    haul_num: Union[int, float]
        Haul number.
    length: float
        Body length.
    sex: Union[int, str]
        Animal sex. This can either be represented by an integer value, or a string corresponding
        to male ("male", "m"), female ("female", "f"), or unsexed/unidentified ("unsexed", "u").
    species_id: Union[int, float, str]
        Species identification code.
    weight: float
        Individual specimen weight (kg).
    """

    age: Series = Field(ge=0, nullable=True, metadata=dict(types=[int, float]))
    haul_num: Series = Field(
        nullable=False, metadata=dict(types=[int, float], drop_invalid_rows=True)
    )
    length: Series[float] = Field(
        gt=0.0, coerce=True, nullable=False, metadata=dict(drop_invalid_rows=True)
    )
    sex: Series = Field(nullable=False, metadata=dict(types=[int, str], drop_invalid_rows=True))
    species_id: Series = Field(
        nullable=False, metadata=dict(types=[int, float, str], drop_invalid_rows=True)
    )
    weight: Series[float] = Field(gt=0.0, nullable=True, coerce=True)

    @check(
        "age",
        "haul_num",
        name="mixed_type",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul_num(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "sex",
        name="mixed_type",
        error=(
            "Column datatype should either be all 'int', or all 'str' contained within "
            "['male', 'm', 'female', 'f', 'unsexed', 'u']"
        ),
    )
    def validate_sex(cls, sex: Series) -> Series[bool]:
        # ---- Check if value is 'int'
        if all(isinstance(x, str) for x in sex):
            return Series([value in ["male", "female", "unsexed", "m", "f", "u"] for value in sex])
        else:
            return Series([isinstance(value, int) for value in sex])

    @check(
        "species_id",
        name="mixed_type",
        error="Column datatype should be either 'int', 'float', or 'str'",
    )
    def validate_species_id(cls, species_id: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float, str)) for value in species_id])


class HaulTransect(BaseDataFrame):
    """
    Haul-transect map DataFrame

    Parameters
    ----------
    haul_num: Union[int, float]
        Haul number.
    transect_num: Union[int, float]
        Transect number.
    """

    haul_num: Series = Field(nullable=False, metadata=dict(types=[int, float]))
    transect_num: Series = Field(nullable=False, metadata=dict(types=[int, float]))

    @check(
        "haul_num",
        "transect_num",
        name="element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul_num(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])


class KSStrata(BaseDataFrame):
    """
    Length-based stratification DataFrame

    Parameters
    ----------
    fraction*: float
        Column corresponding to the fraction of the total net haul weight that is hake. This column
        must include "fraction" in the name.
    haul*: Union[int, float]
        Haul number. This column must include "haul" in the name.
    stratum*: int
        Length-based stratum index/group. This column must include "stratum" in the name.
    """

    fraction: Series[float] = Field(
        ge=0.0,
        le=1.0,
        nullable=False,
        regex=True,
        coerce=True,
        alias=r".*fraction.*",
        metadata=dict(annotation="fraction", drop_invalid_rows=True),
    )
    haul: Series = Field(
        nullable=False,
        alias=r"(haul|haul_num|haul_start|haul_end)",
        regex=True,
        coerce=True,
        metadata=dict(types=[int, float], drop_invalid_rows=True, annotation="haul"),
    )
    stratum: Series = Field(
        nullable=False,
        regex=True,
        alias=r"stratum(_(?:num|inpfc))?$",
        coerce=True,
        metadata=dict(types=[int, float, str], drop_invalid_rows=True, annotation="stratum"),
    )

    @check(
        "haul",
        name="mixed_type",
        regex=True,
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "stratum",
        name="mixed_type",
        regex=True,
        error="Column datatype should be either 'int', 'float', or 'str'",
    )
    def validate_stratum(cls, stratum: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float, str)) for value in stratum])


class GeoStrata(BaseDataFrame):
    """
    Geographic extents of length- and INPFC-based strata DataFrame

    Parameters
    ----------
    haul*: Union[int, float]
        Haul number. This column must include "haul" in the name.
    northlimit_latitude: float
        Northern limit of the lstratum.
    stratum*: int
        Length-based stratum index/group. This column must include "stratum" in the name.
    """

    haul: Series = Field(
        nullable=False,
        alias=r"(haul|haul_num|haul_start|haul_end)",
        regex=True,
        coerce=True,
        metadata=dict(types=[int, float], drop_invalid_rows=True, annotation="haul"),
    )
    northlimit_latitude: Series[float] = Field(ge=-90.0, le=90.0, coerce=True, nullable=False)
    stratum: Series = Field(
        nullable=False,
        regex=True,
        alias=r"stratum(_(?:num|inpfc))?$",
        coerce=True,
        metadata=dict(types=[int, float, str], drop_invalid_rows=True, annotation="stratum"),
    )

    @check(
        "haul",
        regex=True,
        name="mixed_type",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "stratum",
        regex=True,
        name="mixed_type",
        error="Column datatype should be either 'int', 'float', or 'str'",
    )
    def validate_stratum(cls, stratum: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float, str)) for value in stratum])


class AcousticData(BaseDataFrame):
    """
    Haul-transect map DataFrame

    Parameters
    ----------
    bottom_depth: Optional[float]
        Bottom depth (m).
    haul_num: Union[int, float]
        Haul number. This column must include "haul" in the name.
    latitude: float
        Latitude coordinates.
    layer_height: Optional[float]
        School/aggregation layer height (m).
    layer_mean_depth: Optional[float]
        Mean school/aggregation depth (m).
    longitude: float
        Longitude coordinates.
    nasc: float
        Nautical area scattering coefficient (NASC, m^2 nmi^-2).
    region_id: Optional[Union[int, float, str]]
        Echoview region identification value.
    transect_num: Union[int, float]
        Transect number.
    transect_spacing: float
        Distance (spacing) between transect lines.
    distance_s: float
        Distance from the start of the vessel log to the last ping in the analysis interval. This
        column is automatically generated using the following columns when `distance_s` is absent:
        ['dist_s', 'vl_s', 'vl_start', 'vessel_log_start']
    distance_e: float
        Distance from the start of the vessel log to the first ping in the analysis interval. This
        column is automatically generated using the following columns when `distance_e` is absent:
        ['dist_e', 'vl_e', 'vl_end', 'vessel_log_end']
    """

    bottom_depth: Optional[Series[float]] = Field(nullable=False, coerce=True)
    haul_num: Series = Field(
        nullable=False, metadata=dict(types=[int, float], drop_invalid_rows=True)
    )
    latitude: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, coerce=True)
    layer_height: Optional[Series[float]] = Field(nullable=False, coerce=True)
    layer_mean_depth: Optional[Series[float]] = Field(nullable=False, coerce=True)
    longitude: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, coerce=True)
    nasc: Series[float] = Field(ge=0.0, nullable=False, coerce=True)
    region_id: Optional[Series] = Field(metadata=dict(types=[int, float, str]))
    transect_num: Series = Field(nullable=False, metadata=dict(types=[int, float]))
    transect_spacing: Series[float] = Field(ge=0.0, nullable=False, coerce=True)
    distance_s: Series[float] = Field(ge=0.0, nullable=False, coerce=True)
    distance_e: Series[float] = Field(ge=0.0, nullable=False, coerce=True)

    @classmethod
    def pre_validate(cls, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        distance_s_aliases = ["distance_s", "dist_s", "vl_s", "vl_start", "vessel_log_start"]
        overlap_s = list(set(distance_s_aliases).intersection(df.columns))

        if len(overlap_s) == 1:
            df = df.rename(columns={list(overlap_s)[0]: "distance_s"})

        distance_e_aliases = ["distance_e", "dist_e", "vl_e", "vl_end", "vessel_log_end"]
        overlap_e = list(set(distance_e_aliases).intersection(df.columns))

        if len(overlap_e) == 1:
            df = df.rename(columns={list(overlap_e)[0]: "distance_e"})

        haul_overlap = df.columns[df.columns.str.contains("haul")].tolist()

        if len(haul_overlap) == 1:
            df = df.rename(columns={haul_overlap[0]: "haul_num"})

        return df

    @check(
        "haul_num",
        name="multi_type",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "transect_num",
        name="multi_type",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_transect(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @dataframe_check(
        error=lambda df: (
            f"'distance_s': Multiple overlapping column names associated with the expected column "
            f"`distance_s` were found: "
            f"""{
                list(set(
                    ['dist_s', 'distance_s', 'vl_s', 'vl_start', 'vessel_log_start']
                ).intersection(df.columns))
            }. """
            f"Please either rename or remove overlapping columns to avoid ambiguity"
        ),
        description="distance_s",
    )
    def validate_distance_s_aliases(cls, df: pd.DataFrame) -> bool:
        distance_aliases = ["distance_s", "dist_s", "vl_s", "vl_start", "vessel_log_start"]
        overlap = list(set(distance_aliases).intersection(df.columns))

        return len(overlap) == 1
    
    @dataframe_check(
        error=lambda df: (
            f"'distance_e': Multiple overlapping column names associated with the expected column "
            f"`distance_e` were found: "
            f"""{
                list(set(
                    ['dist_s', 'distance_s', 'vl_s', 'vl_start', 'vessel_log_start']
                ).intersection(df.columns))
            }. """
            f"Please either rename or remove overlapping columns to avoid ambiguity"
        ),
        description="distance_e",
    )
    def validate_distance_e_aliases(cls, df: pd.DataFrame) -> bool:
        distance_aliases = ["distance_e", "dist_e", "vl_e", "vl_end", "vessel_log_end"]
        overlap = list(set(distance_aliases).intersection(df.columns))

        return len(overlap) == 1


class IsobathData(BaseDataFrame):
    """
    Haul-transect map DataFrame

    Parameters
    ----------
    latitude: float
        Latitude coordinates.
    longitude: float
        Longitude coordinates.
    """

    latitude: Series[float] = Field(
        ge=-90.0,
        le=90.0,
        nullable=False,
        alias=r".*latitude.*",
        regex=True,
        coerce=True,
        metadata=dict(annotation="latitude", drop_invalid_rows=True),
    )
    longitude: Series[float] = Field(
        ge=-180.0,
        le=180.0,
        nullable=False,
        alias=r".*longitude.*",
        regex=True,
        coerce=True,
        metadata=dict(annotation="longitude", drop_invalid_rows=True),
    )


class KrigedMesh(BaseDataFrame):
    """
    Haul-transect map DataFrame

    Parameters
    ----------
    fraction*: float
        Fraction of kriging mesh cell that is within the interpolation polygon.
    latitude: float
        Latitude coordinates.
    longitude: float
        Longitude coordinates.
    """

    fraction: Series[float] = Field(
        ge=0.0,
        le=1.0,
        nullable=False,
        alias=r".*fraction.*",
        regex=True,
        coerce=True,
        metadata=dict(annotation="fraction", drop_invalid_rows=True),
    )
    latitude: Series[float] = Field(
        ge=-90.0,
        le=90.0,
        nullable=False,
        alias=r".*latitude.*",
        regex=True,
        coerce=True,
        metadata=dict(annotation="latitude", drop_invalid_rows=True),
    )
    longitude: Series[float] = Field(
        ge=-180.0,
        le=180.0,
        nullable=False,
        alias=r".*longitude.*",
        regex=True,
        coerce=True,
        metadata=dict(annotation="longitude", drop_invalid_rows=True),
    )


class VarioKrigingPara(BaseDataFrame):
    """
    Haul-transect map DataFrame

    Parameters
    ----------
    hole: float
        Length scale or range of the hole effect.
    lscl: float
        The relative length scale, or range at which the correlation between points becomes
        approximately constant.
    nugt: float
        The y-intercept of the variogram representing the short-scale (i.e. smaller than the lag
        resolution) variance.
    powr: float
        The exponent used for variogram models with exponentiated spatial decay terms.
    res: float
        The (scaled) distance between lags.
    sill: float
        The total variance where the change autocorrelation reaches (or nears) 0.0.
    ratio: float
        The directional aspect ratio of anisotropy.
    srad: float
        The adaptive search radius used for kriging.
    kmin: float
        The minimum number of nearest kriging points.
    kmax: float
        The maximum number of nearest kriging points.
    """

    class Config:
        coerce = True

    y_offset: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="dataprep.y_offset")
    corr: Series[float] = Field(ge=0.0, nullable=False, alias="vario.corr")
    hole: Series[float] = Field(ge=0.0, nullable=False, alias="vario.hole")
    lscl: Series[float] = Field(ge=0.0, nullable=False, alias="vario.lscl")
    nugt: Series[float] = Field(ge=0.0, nullable=False, alias="vario.nugt")
    powr: Series[float] = Field(ge=0.0, nullable=False, alias="vario.powr")
    range: Series[float] = Field(ge=0.0, nullable=False, alias="vario.range")
    res: Series[float] = Field(gt=0.0, nullable=False, alias="vario.res")
    sill: Series[float] = Field(ge=0.0, nullable=False, alias="vario.sill")
    ytox_ratio: Series[float] = Field(nullable=False, alias="vario.ytox_ratio")
    ztox_ratio: Series[float] = Field(nullable=False, alias="vario.ztox_ratio")
    blk_nx: Series[int] = Field(gt=0, nullable=False, alias="krig.blk_nx")
    blk_ny: Series[int] = Field(gt=0, nullable=False, alias="krig.blk_ny")
    blk_nz: Series[int] = Field(gt=0, nullable=False, alias="krig.blk_nz")
    dx0: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, alias="krig.dx0")
    dx: Series[float] = Field(nullable=False, alias="krig.dx")
    dy0: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="krig.dy0")
    dy: Series[float] = Field(nullable=False, alias="krig.dy")
    dz: Series[float] = Field(nullable=False, alias="krig.dz")
    elim: Series[float] = Field(nullable=False, alias="krig.elim")
    eps: Series[float] = Field(nullable=False, alias="krig.eps")
    kmax: Series[int] = Field(gt=0, nullable=False, alias="krig.kmax")
    kmin: Series[int] = Field(gt=0, nullable=False, alias="krig.kmin")
    nx: Series[int] = Field(gt=0, nullable=False, alias="krig.nx")
    ny: Series[int] = Field(gt=0, nullable=False, alias="krig.ny")
    nz: Series[int] = Field(gt=0, nullable=False, alias="krig.nz")
    ratio: Series[float] = Field(nullable=False, alias="krig.ratio")
    srad: Series[float] = Field(gt=0.0, nullable=False, alias="krig.srad")
    x_res: Series[float] = Field(nullable=False, alias="krig.x_res")
    xmin: Series[float] = Field(nullable=False, alias="krig.xmin")
    xmax: Series[float] = Field(nullable=False, alias="krig.xmax")
    xmin0: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, alias="krig.xmin0")
    xmax0: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, alias="krig.xmax0")
    y_res: Series[float] = Field(nullable=False, alias="krig.y_res")
    ymin: Series[float] = Field(nullable=False, alias="krig.ymin")
    ymax: Series[float] = Field(nullable=False, alias="krig.ymax")
    ymin0: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="krig.ymin0")
    ymax0: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="krig.ymax0")


####################################################################################################
# ECHOVIEW EXPORT VALIDATION
# --------------------------------------------------------------------------------------------------


class EchoviewCells(BaseDataFrame):
    """
    Echoview ``cells`` ``*.csv`` file

    Parameters
    ----------
    interval: int
        Cell interval number
    layer: int
        Cell layer number
    prc_nasc: float
        Proportioned region-to-cell NASC (m^2 nmi^-2)
    process_id: int
        Unique identification number assigned to each export file
    region_class: str
        Class name
    region_id: int
        Identification number of the export region
    region_name: str
        Region name
    standard_deviation: float
        Standard deviation of sample values (calculated in the linear domain) of the analysis cell
    sv_mean: float
        Mean volumetric backscattering strength of the analysis cell (S_V, dB re. 1 m^-1)

    """

    interval: Series[int] = Field(ge=1, nullable=False)
    layer: Series[int] = Field(ge=1, nullable=False)
    prc_nasc: Series[float] = Field(ge=0.0, nullable=False)
    process_id: Series[int]
    region_class: Series[str]
    region_id: Series[int]
    region_name: Series[str]
    standard_deviation: Series[float] = Field(ge=0.0, nullable=False)
    sv_mean: Series[float] = Field(ge=-999, nullable=True)


class EchoviewIntervals(BaseDataFrame):
    """
    Echoview ``intervals`` ``*.csv`` file

    Parameters
    ----------
    date_s: int
        Date of the first ping in the analysis interval (YYYYMMDD)
    exclude_below_line_depth_mean: float
        Mean bottom depth of the analysis interval (m)
    interval: int
        Cell interval number
    lat_s: float
        Latitude of the first ping in the analysis interval (deicmal degrees)
    lon_s: float
        Longitude of the first ping in the analysis interval (decimal degrees)
    process_id: int
        Unique identification number assigned to each export file
    time_s: str
        Timestamp of the first ping in the analysis interval as a string (HH:MM:OS)
    distance_s: float
        Distance from the start of the vessel log to the last ping in the analysis interval. This
        column is automatically generated using the following columns when `distance_s` is absent:
        ['dist_s', 'vl_s', 'vl_start', 'vessel_log_start']
    distance_e: float
        Distance from the start of the vessel log to the first ping in the analysis interval. This
        column is automatically generated using the following columns when `distance_e` is absent:
        ['dist_e', 'vl_e', 'vl_end', 'vessel_log_end']

    """

    date_s: Series[int]
    distance_e: Series[float] = Field(ge=0.0, nullable=False, coerce=True)
    distance_s: Series[float] = Field(ge=0.0, nullable=False, coerce=True)
    exclude_below_line_depth_mean: Series[float] = Field(coerce=True)
    interval: Series[int] = Field(ge=1, nullable=False)
    lat_s: Series[float] = Field(coerce=True)
    lon_s: Series[float] = Field(coerce=True)
    process_id: Series[int]
    time_s: Series[str]

    @classmethod
    def pre_validate(cls, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        distance_s_aliases = ["distance_s", "dist_s", "vl_s", "vl_start", "vessel_log_start"]
        overlap_s = list(set(distance_s_aliases).intersection(df.columns))

        if len(overlap_s) == 1:
            df = df.rename(columns={list(overlap_s)[0]: "distance_s"})

        distance_e_aliases = ["distance_e", "dist_e", "vl_e", "vl_end", "vessel_log_end"]
        overlap_e = list(set(distance_e_aliases).intersection(df.columns))

        if len(overlap_e) == 1:
            df = df.rename(columns={list(overlap_e)[0]: "distance_e"})

        return df

    @dataframe_check(
        error=lambda df: (
            f"'distance_s': Multiple overlapping column names associated with the expected column "
            f"`distance_s` were found: "
            f"""{
                list(set(
                    ['dist_s', 'distance_s', 'vl_s', 'vl_start', 'vessel_log_start']
                ).intersection(df.columns))
            }. """
            f"Please either rename or remove overlapping columns to avoid ambiguity"
        ),
        description="distance_s",
    )
    def validate_distance_s_aliases(cls, df: pd.DataFrame) -> bool:
        distance_aliases = ["distance_s", "dist_s", "vl_s", "vl_start", "vessel_log_start"]
        overlap = list(set(distance_aliases).intersection(df.columns))

        return len(overlap) == 1

    @dataframe_check(
        error=lambda df: (
            f"'distance_e': Multiple overlapping column names associated with the expected column "
            f"`distance_e` were found: "
            f"""{
                list(set(['dist_e', 'distance_e', 'vl_e', 'vl_end', 'vessel_log_end']).intersection(
                    df.columns
                ))
            }. """
            f"Please either rename or remove overlapping columns to avoid ambiguity"
        ),
        description="distance_e",
    )
    def validate_distance_e_aliases(cls, df: pd.DataFrame) -> bool:
        distance_aliases = ["distance_e", "dist_e", "vl_e", "vl_end", "vessel_log_end"]
        overlap = list(set(distance_aliases).intersection(df.columns))

        return len(overlap) == 1


class EchoviewLayers(BaseDataFrame):
    """
    Echoview ``layers`` ``*.csv`` file

    Parameters
    ----------
    layer: int
        Cell layer number
    layer_depth_max: Union[int, float]
        Maximum depth of the analysis layer (m)
    layer_depth_min: Union[int, float]
        Minimum depth of the analysis layer (m)
    process_id: int
        Unique identification number assigned to each export file

    """

    layer: Series[int] = Field(ge=0, nullable=False)
    layer_depth_max: Series = Field(nullable=False, metadata=dict(types=[int, float]))
    layer_depth_min: Series = Field(nullable=False, metadata=dict(types=[int, float]))
    process_id: Series[int]


####################################################################################################
# CONFIGURATION DATAFRAME MAPPING
# --------------------------------------------------------------------------------------------------
DATASET_DF_MODEL = {
    "biological": {
        "catch": CatchBiodata,
        "haul_to_transect": HaulTransect,
        "length": LengthBiodata,
        "specimen": SpecimenBiodata,
    },
    "stratification": {
        "geo_strata": GeoStrata,
        "inpfc_geo_strata": GeoStrata,
        "strata": KSStrata,
        "inpfc_strata": KSStrata,
    },
    "NASC": {
        "all_ages": AcousticData,
        "no_age1": AcousticData,
    },
    "kriging": {
        "isobath_200m": IsobathData,
        "mesh": KrigedMesh,
        "vario_krig_para": VarioKrigingPara,
    },
}

ECHOVIEW_DF_MODEL = {
    "cells": EchoviewCells,
    "intervals": EchoviewIntervals,
    "layers": EchoviewLayers,
}
