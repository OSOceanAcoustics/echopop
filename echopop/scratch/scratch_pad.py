import abc
import numpy as np
import pandas as pd
from typing import Union, Dict, List, Optional, Any

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs


# ==============================================================================
# MY IMPLEMENTATIONS OF INVERSION CLASSES
# ==============================================================================
model_parameters = {
    "ts_length_regression": {
        "slope": 20.,
        "intercept": -68.
    },
    "stratify_by": "stratum_ks",
    "strata": df_dict_strata["ks"].stratum_num.unique(),
    "impute_missing_strata": True,
}

strata_options = np.array([1, 2, 3, 4, 5])
sigma_bs_stratum = pd.DataFrame(
    {
        "stratum_num": [1, 2, 3],
        "species_id": np.repeat(94832, 3),
        "sigma_bs_mean": [1.0, 2.0, 3.0],
    }
)


# ==============================================================================
# TRANSECT INTERVAL CORRECTION FUNCTIONS
# ==============================================================================


proportions_dict = dict_df_weight_proportion
group_columns = ["sex"]
group_by = "stratum_ks"
df_nasc = df_nasc_no_age1.copy().set_index([stratify_by])
df_average_weight = df_averaged_weight["all"].copy()
####
# FILTER BY LEVEL
###
# proportions_dict["aged"].xs("male", level="sex")
# proportions_dict["aged"].loc[:, :, 1]

# proportions_dict["aged"].xs("age_bin", axis=)

# value = 1
# target_vals = proportions_dict["aged"].index.get_level_values("age_bin")
# mask = pd.Series([value in interval if hasattr(interval, '__contains__') 
#                     else interval == value 
#                     for interval in target_vals])

# proportions_dict["aged"].iloc[mask]

# level_name = "age_bin"
# level_vals = proportions_dict["aged"].index.get_level_values(level_name)
# mask = [value in val if hasattr(val, '__contains__') else val == value 
#         for val in level_vals]
# proportions_dict["aged"][mask]

###
# ABUNDANCE / BIOMASS CALCULATION
###

df_nasc["abundance"] = np.round(df_nasc["area_interval"] * df_nasc["number_density"])
df_nasc["biomass"] = df_nasc["abundance"] * df_averaged_weight["all"].reindex_like(df_nasc)
df_nasc["biomass_density"] = df_nasc["number_density"] * df_averaged_weight["all"].reindex_like(df_nasc)

