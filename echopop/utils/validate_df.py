import copy
import re
from typing import Optional

import numpy as np
import pandas as pd
from pandera import DataFrameModel, Field, check
from pandera.errors import SchemaError, SchemaErrors
from pandera.typing import Series


####################################################################################################
# UTILITY FUNCTION
# --------------------------------------------------------------------------------------------------
def extract_errors(data, failed_coercion: pd.DataFrame, key="error"):
    """
    Utility function for extract `pandera.SchemaError` and `pandera.SchemaErrors` messages
    """
    errors = []

    if isinstance(data, dict):
        # If the current level is a dictionary, iterate through its keys
        # ---- Check if coercion error occurred
        if "error" in data and "column" in data:
            if data["column"] in failed_coercion.Column.unique():
                # ---- Match the column
                error_msg = failed_coercion.loc[
                    failed_coercion.Column == data["column"], "error"
                ].values[0]
                # ---- Append
                errors.append(f"   -{error_msg}")
            # for k, v in data.items():
            # if "error" in data:
            #     errors.append(f"   -{data["error"].capitalize()}")
            else:
                errors.append(f"   -{data['error'].capitalize()}")
        else:
            for k, v in data.items():
                errors.extend(extract_errors(v, failed_coercion, key))
    elif isinstance(data, list):
        # If the current level is a list, iterate through its elements
        for item in data:
            errors.extend(
                extract_errors(item, failed_coercion, key)
            )  # Recursively search nested lists

    # Drop 'traceback' messages
    return list(set(errors))


####################################################################################################
# PANDERA VALIDATORS
# --------------------------------------------------------------------------------------------------
class BaseDataFrame(DataFrameModel):

    class Config:
        metadata = {}
        strict = False

    _DTYPE_TESTS = {
        int: (
            lambda v: pd.api.types.is_integer_dtype(v)
            or (pd.api.types.is_numeric_dtype(v) and (v % 1 == 0).all())
        ),
        float: lambda v: pd.api.types.is_float_dtype(v),
        str: lambda v: pd.api.types.is_object_dtype(v),
    }

    _DTYPE_COERCION = {
        int: lambda v: v.astype(np.float64, errors="ignore").astype(np.int64, errors="ignore"),
        float: lambda v: v.astype(np.float64, errors="ignore"),
        str: lambda v: v.astype(str, errors="ignore"),
    }

    def __init__(self, *args, **kwargs):
        # Copy class-level __annotations__ to instance
        self.__annotations__ = self.__class__.__annotations__.copy()
        super().__init__(*args, **kwargs)

    @classmethod
    def get_column_types(cls) -> dict:
        # Get the schema of the DataFrameModel
        schema = cls.to_schema()
        column_types = {}

        # Iterate over columns in the schema
        for column_name, column_schema in schema.columns.items():
            # Access the dtype directly
            # ---- If not None:
            if column_schema.dtype:
                column_types[column_name] = column_schema.dtype.type
            # ---- Collect mixed types from metadata
            else:
                column_types[column_name] = cls.get_metadata()[cls.__name__]["columns"][
                    column_name
                ]["types"]

        return column_types

    @classmethod
    def coercion_check(cls, df: pd.DataFrame, col_types: dict) -> pd.DataFrame:

        # Initialize a DataFrame
        errors_coerce = pd.DataFrame(dict(Column=[], error=[]))

        # Coerce the data
        for column_name, dtype in col_types.items():
            # ---- Apply coercion based on column patterns
            for col in df.columns:
                if re.match(column_name, col):
                    # ---- Retrieve the column name
                    # Coerce the column to the appropriate dtype
                    if isinstance(dtype, list):
                        # ---- Initialize the dtype of the validator annotation
                        cls.__annotations__[col] = Series[dtype[0]]
                        # ---- Check for all valid integers
                        for typing in dtype:
                            test = cls._DTYPE_TESTS.get(typing, None)
                            # ---- Adjust typing annotations
                            if test and test(df[col]):
                                cls.__annotations__[col] = Series[typing]
                                # break
                            # ---- Coerce the datatypes
                            try:
                                df[col] = cls._DTYPE_COERCION.get(typing)(df[col])
                                break
                            except Exception as e:
                                e.__traceback__ = None
                                message = (
                                    f"{col.capitalize()} column must be a Series of '{str(dtype)}' "
                                    f"values. Series values could not be automatically coerced."
                                )
                                # errors_coerce.append(e)
                                errors_coerce = pd.concat(
                                    [errors_coerce, pd.DataFrame(dict(Column=col, error=message))]
                                )
                    # ---- If not a List from the metadata attribute
                    else:
                        try:
                            df[col] = df[col].astype(str(dtype))
                        except Exception as e:
                            e.__traceback__ = None
                            message = (
                                f"{col.capitalize()} column must be a Series of '{str(dtype)}' "
                                f"values. Series values could not be automatically coerced."
                            )
                            # errors_coerce.append(e)
                            errors_coerce = pd.concat(
                                [errors_coerce, pd.DataFrame(dict(Column=[col], error=[message]))]
                            )

        # Return the DataFrame
        return errors_coerce

    @classmethod
    def judge(
        cls, df: pd.DataFrame, coercion_failures: pd.DataFrame, filename: str
    ) -> pd.DataFrame:
        # ---- Collect errors
        errors = []
        try:
            return cls.validate(df, lazy=True)
        except SchemaErrors as e:
            e.__traceback__ = None
            errors.append(e)
        # ---- If Errors occur:
        if errors:
            # ---- Join unique errors
            errors_stk = "\n".join(extract_errors(errors[0].message, coercion_failures))
            # ---- Format the error message
            message = (
                f"The following DataFrame validation errors were flagged for {filename}: \n"
                f"{errors_stk}"
            )
            # ---- Raise Error
            raise SchemaError(cls, df, message) from None

    @classmethod
    def validate_df(cls, data: pd.DataFrame, filename: Optional[str] = None) -> pd.DataFrame:
        # ---- Create copy
        df = data.copy()
        # ---- Get the original annotations
        default_annotations = copy.deepcopy(cls.__annotations__)
        # ---- Format the column names
        df.columns = [col.lower() for col in df.columns]
        # ---- Get column types
        column_types = cls.get_column_types()
        # ---- Initialize invalid index
        invalid_idx = {}
        # ---- Initialize column names
        valid_cols = []
        # ---- Find all indices where there are violating NaN/-Inf/Inf values
        for column_name, dtype in column_types.items():
            # ---- Apply coercion based on column patterns
            for col in df.columns:
                # ---- Regular expression matching
                if re.match(column_name, col):
                    # ---- Retrieve the column name
                    col_name = re.match(column_name, col).group(0)
                    # ---- Collect, if valid
                    if cls.to_schema().columns[column_name].regex:
                        valid_cols.append(col)
                    else:
                        valid_cols.append(col_name)
                    # ---- Check if null values are allowed
                    if not cls.to_schema().columns[column_name].nullable:
                        # ---- Get the violating indices
                        invalid_idx.update(
                            {
                                col: df.index[
                                    df[col].isna() | df[col].isin([np.inf, -np.inf])
                                ].to_list()
                            }
                        )
        # ---- Initialize the list
        invalid_lst = []
        # ---- Extend the values
        invalid_lst.extend(
            [
                value
                for key in invalid_idx.keys()
                if invalid_idx[key] is not None and key in valid_cols
                for value in invalid_idx[key]
            ]
        )
        # ---- If the indices are invalid, but can be dropped then drop them
        df.drop(invalid_lst, axis=0, inplace=True)
        # ---- Reset index
        if not df.empty:
            df.reset_index(inplace=True, drop=True)

        # Coercion check
        coercion_failures = cls.coercion_check(df, column_types)

        # Validate
        df_valid = cls.judge(df.filter(valid_cols), coercion_failures, filename)
        # ---- Return to default annotations
        cls.__annotations__ = default_annotations
        # ---- Return the validated DataFrame
        return df_valid


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

    haul_num: Series = Field(nullable=False, metadata=dict(types=[int, float]))
    length: Series[float] = Field(gt=0.0, nullable=False)
    length_count: Series[int] = Field(ge=0, nullable=False)
    sex: Series = Field(nullable=False, metadata=dict(types=[int, str]))
    species_id: Series = Field(nullable=False, metadata=dict(types=[int, float, str]))

    @check(
        "haul_num",
        name="'haul_num' element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul_num(cls, haul_num: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in haul_num])

    @check(
        "sex",
        name="'sex' element-wise datatypes",
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
        name="'species_id' element-wise datatypes",
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

    haul_num: Series = Field(nullable=False, metadata=dict(types=[int, float]))
    haul_weight: Series[float] = Field(ge=0.0, nullable=False)
    species_id: Series = Field(nullable=False, metadata=dict(types=[int, float, str]))

    @check(
        "haul_num",
        name="'haul_num' element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul_num(cls, haul_num: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in haul_num])

    @check(
        "species_id",
        name="'species_id' element-wise datatypes",
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
    haul_num: Series = Field(nullable=False, metadata=dict(types=[int, float]))
    length: Series[float] = Field(gt=0.0, nullable=False)
    sex: Series = Field(nullable=False, metadata=dict(types=[int, str]))
    species_id: Series = Field(nullable=False, metadata=dict(types=[int, float, str]))
    weight: Series[float] = Field(gt=0.0, nullable=True)

    @check(
        "age",
        "haul_num",
        name="element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul_num(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "sex",
        name="'sex' element-wise datatypes",
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
        name="'species_id' element-wise datatypes",
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
        ge=0.0, le=1.0, nullable=False, regex=True, alias=".*fraction.*"
    )
    haul: Series = Field(nullable=False, regex=True, metadata=dict(types=[int, float]))
    stratum: Series = Field(nullable=False, regex=True, metadata=dict(types=[int, float, str]))

    @check(
        "haul",
        name="element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "stratum",
        name="'stratum' element-wise datatypes",
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

    haul: Series = Field(nullable=False, regex=True, metadata=dict(types=[int, float]))
    northlimit_latitude: Series[float] = Field(ge=-90.0, le=90.0, nullable=False)
    stratum: Series = Field(nullable=False, regex=True, metadata=dict(types=[int, float, str]))

    @check(
        "haul",
        name="element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "stratum",
        name="'stratum' element-wise datatypes",
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
    haul*: Union[int, float]
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
    vessel_log_start: float
        Cumulative vessel log distance at the start of each transect interval.
    vessel_log_start: end
        Cumulative vessel log distance at the end of each transect interval.
    """

    bottom_depth: Optional[Series[float]] = Field(nullable=False, coerce=True)
    haul: Series = Field(nullable=False, regex=True, metadata=dict(types=[int, float]))
    latitude: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, regex=True, coerce=True)
    layer_height: Optional[Series[float]] = Field(nullable=False, coerce=True)
    layer_mean_depth: Optional[Series[float]] = Field(nullable=False, coerce=True)
    longitude: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, regex=True, coerce=True)
    nasc: Series[float] = Field(ge=0.0, nullable=False, coerce=True)
    region_id: Optional[Series] = Field(metadata=dict(types=[int, float, str]))
    transect_num: Series = Field(nullable=False, regex=True, metadata=dict(types=[int, float]))
    transect_spacing: Series[float] = Field(ge=0.0, nullable=False, regex=False, coerce=True)
    vessel_log_start: Series[float] = Field(ge=0.0, nullable=False, coerce=True)
    vessel_log_end: Series[float] = Field(ge=0.0, nullable=False, coerce=True)

    @check(
        "haul",
        name="element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_haul(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])

    @check(
        "transect_num",
        name="element-wise datatypes",
        error="Column datatype should be either 'int' or 'float'",
    )
    def validate_transect(cls, v: Series) -> Series[bool]:
        return Series([isinstance(value, (int, float)) for value in v])


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
        ge=-90.0, le=90.0, nullable=False, regex=True, alias=".*latitude.*"
    )
    longitude: Series[float] = Field(
        ge=-180.0, le=180.0, nullable=False, regex=True, alias=".*longitude.*"
    )


class KrigedMesh(IsobathData):
    """
    Haul-transect map DataFrame

    Parameters
    ----------
    fraction*: float
        Fraction of kriging mesh cell that is within the interpolation polygon.
    """

    fraction: Series[float] = Field(
        ge=0.0, le=1.0, nullable=False, regex=True, coerce=True, alias=".*fraction.*"
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
        "inpfc_strata": GeoStrata,
        "strata": KSStrata,
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
