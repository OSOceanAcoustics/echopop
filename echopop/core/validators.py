from typing import Any, Dict, Type, TypeVar

import pandas as pd
import pandera.pandas as pa
from pydantic import BaseModel, ConfigDict, ValidationError

# Type variables
BD = TypeVar("BD", bound="BaseDictionary")
BDF = TypeVar("BDF", bound="BaseDataFrame")


class BaseDictionary(BaseModel):
    """
    Common base Pydantic model providing consistent validation and serialization for file and
    workflow inputs.
    """

    model_config = ConfigDict(
        arbitrary_types_allowed=True, extra="forbid", validate_default=True, populate_by_name=True
    )

    # Validator method
    @classmethod
    def judge(cls: Type[BD], **kwargs) -> BD:
        """
        Safely validate input data and return a model instance.
        """
        try:
            return cls(**kwargs)
        except ValidationError as e:
            # Add schema context to help debugging large validation chains
            raise e

    # Factory method
    @classmethod
    def create(cls: Type[BD], **kwargs) -> Dict[str, Any]:
        """
        Factory creation method that returns a clean dictionary representation of the model.
        """
        return cls.judge(**kwargs).model_dump(exclude_none=True)


class BaseDataFrame(pa.DataFrameModel):
    "Common base Pandera schema model providing consistent validation."

    class Config:
        strict = False
        unique_column_names = True
        coerce = True

    @classmethod
    def pre_validate(cls: Type[BDF], df: pd.DataFrame) -> BDF:
        """
        Hook for child classes to override and perform pre-validation checks or transformations.
        Default: returns input unchanged.
        """
        return df

    @classmethod
    def validate(cls: BDF, df: pd.DataFrame, *args: Any, **kwargs: Any):
        """
        Validate a DataFrame after applying `pre_validate`. Adds an optional uniform error context
        for clearer logs since `pandera` error message formatting can get mangled.
        """
        # Run pre-validation
        try:
            df_checked = cls.pre_validate(df)
            return super().validate(df_checked, *args, **kwargs)
        except Exception as e:
            raise ValueError(f"Validation failed for {cls.__name__}: {e}") from None
