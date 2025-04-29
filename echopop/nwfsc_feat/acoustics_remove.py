"""
THIS WILL BE REMOVED SOON
KEPT AS REFERENCE FOR NOW
"""
from typing import Union, Literal, Dict
import numpy as np
import pandas as pd


def ts_length_regression(
    length: Union[np.ndarray, float], slope: float, intercept: float
) -> np.ndarray:
    pass


# same as the current aggregate_sigma_bs but with explicit inputs
def get_stratified_sigma_bs(
    df_bio_dict: Dict[pd.DataFrame],  # df_bio_dict from end of load_data.py
    ts_length_regression_dict: dict,  # only 1 species so do not need species code and related checking
    stratum_type: Literal["ks", "inpfc"],  # only do 1 type at each call
) -> pd.DataFrame:
    
    # Organize bio data
    # -- meld the specimen and length dataframes together for downstream calculations
    # -- regroup the length-specific length counts

    # Merge bio and regression dfs
    # -- create DataFrame containing all necessary regression coefficients
    # -- merge with the biological data
    # -- calculate predicted TS from the length values
    # -- convert TS to the linear domain ('sigma_bs')

    # Calculate mean sigma_bs for all hauls and the specified stratum type
    # -- impute sigma_bs values, if necessary, for missing strata
    # -- calculate mean TS for all hauls, KS-strata, and INPFC strata

    df_sigma_bs_haul: pd.DataFrame  # TODO: do we need this for downstream processing?
    df_sigma_bs_stratum: pd.DataFrame

    return df_sigma_bs_stratum
