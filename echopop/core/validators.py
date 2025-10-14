import pandas as pd
import pandera.pandas as pa
from pydantic import BaseModel, ValidationError


class BaseDictionary(BaseModel):
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
            raise e

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """
        return cls.judge(**kwargs).model_dump(exclude_none=True)


class BaseDataFrame(pa.DataFrameModel):
    class Config:
        strict = False
        coerce = True

    @classmethod
    def pre_validate(cls, df):
        """
        Pre-validate DataFrame
        """
        return df

    @classmethod
    def validate(cls, df: pd.DataFrame):
        # Run pre-validation
        df_checked = cls.pre_validate(df)
        # Pass to normal validation
        return super().validate(df_checked)
