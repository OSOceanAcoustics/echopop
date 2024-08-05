import yaml
import re

from functools import reduce
from .sql_methods import SQL
from pathlib import Path
from typing import Union, Tuple, Optional, List

import pandas as pd

import numpy as np

from .live_core import(
    LIVE_FILE_FORMAT_MAP,
    LIVE_INPUT_FILE_CONFIG_MAP
)

def get_unique_identifiers(data_dict: dict,
                           unique_columns: List[str]) -> pd.DataFrame:

    # Gather all dataframes from a dictionary into a list
    df_list = [df for _, df in data_dict.items()]

    # Get unique values of each contrast column across the biological datasets    
    dfs = [pd.DataFrame({col: df[col].unique().tolist()}) for col in unique_columns 
           for df in df_list if not df.empty and isinstance(df, pd.DataFrame)]
    
    # Reduce into a single DataFrame
    if len(unique_columns) > 1:
        return reduce(lambda left, right: pd.merge(left, right, how='cross'), dfs)
    else:
        return reduce(lambda left, right: pd.merge(left, right, how='inner'), dfs)


def query_dataset(db_file: str,
                  data_dict: dict,
                  table_name: str,
                  data_columns: List[str],
                  unique_columns: List[str],
                  constraint: str = None):
    
    # Validate that the desired table exists
    if SQL(db_file, "validate", table_name=table_name):
        # ---- Inspect the SQL table
        inspected_table = SQL(db_file, "inspect", table_name=table_name)
        # ---- Create a list of intersecting column names
        unique_keys = list(set(inspected_table.keys()).intersection(set(unique_columns)))
        # ---- Create list of valid columns
        valid_keys = list(set(inspected_table.keys()).intersection(set(data_columns)))
        # ---- Get unique identifiers
        unique_keys_df = get_unique_identifiers(data_dict, unique_keys)
        # ---- Create conditional string
        conditional_str = (
            " & ".join([f"{col} in {np.unique(unique_keys_df[col])}" 
                        for col in unique_keys_df.columns])
        )
        # ---- Append the additional constraint statement if present
        if constraint is not None:
            conditional_str += f" & {constraint}"
        # ---- SELECT the dataset using the conidtional statement
        data_sql = SQL(db_file, "select", table_name=table_name, columns=valid_keys,
                       condition=conditional_str).filter(data_columns)
    else:
        data_sql = None

    # Return the table DataFrame
    return data_sql
