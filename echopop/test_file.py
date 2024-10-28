from echopop import Survey
import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union
import numpy as np
from echopop.graphics import plotting as egp, variogram_interactive as egv
init_config = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config.yml"
file_config = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
survey = Survey(init_config, file_config)
survey.load_survey_data()
survey.load_acoustic_data()
survey.transect_analysis()
survey.fit_variogram()
survey.kriging_analysis(variogram_parameters={"n_lags": 30}, variable="biomass_density")
kind="age_length_distribution"
variable="abundance"
log_base=10
survey.plot(kind="age_length_distribution", variable="biomass", log_base=10)

dataset = survey.analysis["transect"]["biology"]["population"]["tables"]["biomass"]["aged_biomass_df"]
data_dict = survey.analysis["transect"]["biology"]["population"]["tables"]
variable = "biomass"
sex = "all"
dataset_stk = dataset.stack(future_stack=True).sum(axis=1).reset_index(name=variable)
# ---- Convert the length and age columns from intervals into numerics
# -------- Length
dataset_stk["length"] = dataset_stk["length_bin"].apply(lambda x: x.mid).astype(float)
# -------- Age
dataset_stk["age"] = dataset_stk["age_bin"].apply(lambda x: x.mid).astype(float)
# ---- Subset the data based on the sex
dataset_sub = dataset_stk.loc[dataset_stk["sex"] == sex, :]

dataset_sub["biomass"].max()
5e8