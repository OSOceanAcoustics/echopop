from typing import Literal
import pandas as pd


# same as the current save_transect_coordinates but with explicit inputs
def get_transect_coordinates(
    df_nasc: pd.DataFrame,
    age_type: Literal["no_age1", "all_ages"],
    stratum_type: Literal["ks", "inpfc"],
) -> pd.DataFrame:
    # NOTE: inpfc stratum renaming should be done in load_data.py::clean_stratification()
    pass
