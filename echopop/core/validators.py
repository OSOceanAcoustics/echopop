"""
Core Pydantic and Pandera validators.

This module contains the core Pydantic and Pandera classes that are inherited by all validator
children.
"""

from typing import Any, TypeVar

import pandas as pd
import pandera.pandas as pa
from pydantic import BaseModel, ConfigDict, ValidationError

# Type variables
BD = TypeVar("BD", bound="BaseDictionary")
BDF = TypeVar("BDF", bound="BaseDataFrame")


class BaseDictionary(BaseModel):
    """
    Common base Pydantic model.

    Common base Pydantic model providing consistent validation and serialization for file and
    workflow inputs.
    """

    model_config = ConfigDict(
        arbitrary_types_allowed=True, extra="forbid", validate_default=True, populate_by_name=True
    )

    # Validator method
    @classmethod
    def judge(cls: type[BD], **kwargs) -> BD:
        """Safely validate input data and return a model instance."""
        try:
            return cls(**kwargs)
        except ValidationError as e:
            # Add schema context to help debugging large validation chains
            raise e

    # Factory method
    @classmethod
    def create(cls: type[BD], **kwargs) -> dict[str, Any]:
        """Create a model instance and return a dictionary representation."""
        return cls.judge(**kwargs).model_dump(exclude_none=True)


class BaseDataFrame(pa.DataFrameModel):
    """
    Common base Pandera model.

    Common base Pandera schema model providing consistent validation.
    """

    class Config:
        """Pydantic model configuration."""

        strict = False
        unique_column_names = True
        coerce = True

    @classmethod
    def pre_validate(cls: type[BDF], df: pd.DataFrame) -> BDF:
        """
        Provide a hook for child classes to override for pre-validation.

        This method allows for custom transformations or checks before the main validation logic.
        Returns input unchanged by default.
        """
        return df

    @classmethod
    def validate(cls: BDF, df: pd.DataFrame, *args: Any, **kwargs: Any):
        """
        Validate DataFrame output.

        Validate a DataFrame after applying ``pre_validate``. Adds an optional uniform error context
        for clearer logs since ``pandera`` error message formatting can get mangled.
        """
        # Run pre-validation
        try:
            df_checked = cls.pre_validate(df)
            return super().validate(df_checked, *args, **kwargs)
        except Exception as e:
            raise ValueError(f"Validation failed for {cls.__name__}: {e}") from None
