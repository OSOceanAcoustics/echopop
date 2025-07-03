import abc
import numpy as np
import pandas as pd
from typing import Union, Dict, List, Optional, Any
from functools import reduce
import pytest

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs
from echopop.nwfsc_feat import utils


# Age-1 weight proportions
nasc_df = df_nasc_no_age1.copy()
stratify_by = ["stratum_ks"]
apportion_by = ["sex"]
number_proportions = dict_df_number_proportion

####

def abundance_to_

group_proportions = utils.create_grouped_table(
    number_proportions,
    group_cols=stratify_by + apportion_by,
    index_cols=apportion_by,
    strat_cols=stratify_by,
    value_col="proportion_overall",
).T

# Index the transect data
nasc_df_idx = nasc_df.set_index(stratify_by)

group_props_expanded = group_proportions.reindex(nasc_df_idx.index)

# Multiply 'number_density' with each column in group_props_expanded
new_cols = nasc_df_idx["number_density"].to_frame().values * group_props_expanded.values

# Assign the new columns to nasc_df_idx with appropriate names
nasc_df_idx[group_props_expanded.columns.map(lambda c: f'number_density_{c}')] = new_cols

nasc_biology_grp["number_density_sex_male"].sum()
nasc_df_idx["number_density_male"].sum()
nasc_biology_grp["number_density_unsexed"].sum()
nasc_df_idx["number_density_unsexed"].sum()
group_proportions.reindex(nasc_df_idx.index)
nasc_df_idx[['stratum_ks']].join(group_proportions, on='stratum_ks')
self.analysis["transect"]["acoustics"]["adult_transect_df"]["biomass"].sum()

df_nasc_no_age1["biomass"].sum()
nasc_biology["interval_area"].mean()

A=df_nasc_no_age1.set_index(["stratum_ks", "transect_num", "longitude", "latitude"])[["nasc", "area_interval", "number_density", "abundance"]]
B=nasc_biology.rename(columns={"stratum_num": "stratum_ks"}).set_index(["stratum_ks", "transect_num", "longitude", "latitude"])[["nasc", "interval_area", "number_density", "abundance"]]

A["abundance"].sum()
np.round(B["abundance"]).sum()

A = df_nasc_no_age1.set_index(["stratum_ks"])
(A["biomass"] * (1-age1_weight_proportions.reindex_like(A))).sum()

df_nasc_no_age1["biomass"].sum()
nasc_biology_grp["biomass"].sum()


A.loc[0.0].loc[29].loc[-124.64068436999999]
B.loc[0.0].loc[29].loc[-124.64068436999999]


D = A["abundance"] - B["abundance"]
np.abs(D).max()


df_nasc_no_age1["abundance"].sum()
nasc_biology["abundance"].sum()
####################################################################################################
from echopop.survey import Survey
from echopop.spatial.transect import correct_transect_intervals
from echopop.biology import age1_metric_proportions

survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data(ingest_exports="echoview")
survey.load_survey_data()
survey.transect_analysis()

self = survey
input_dict, analysis_dict, configuration_dict, settings_dict = self.input, self.analysis["transect"], self.config, self.analysis["settings"]
# Extract the necessary correct strata mean sigma_bs
sigma_bs_strata = analysis_dict["acoustics"]["sigma_bs"]["strata_mean_df"]

# Pull out the length-weight conversion for each stratum
length_weight_strata = analysis_dict["biology"]["weight"]["weight_stratum_df"]

# Get the name of the stratum column
stratum_col = settings_dict["transect"]["stratum_name"]

# Get group-specific columns
age_group_cols = settings_dict["transect"]["age_group_columns"]

# Extract the correct strata dataframe
# ---- Define `strata_df` if KS
if settings_dict["transect"]["stratum"] == "ks":
    strata_df = input_dict["spatial"]["strata_df"].copy()
# Define `inpfc_strata_df` if INPFC
elif settings_dict["transect"]["stratum"] == "inpfc":
    strata_df = input_dict["spatial"]["inpfc_strata_df"].copy()

# Get group-specific column names and create conversion key
name_conversion_key = {age_group_cols["haul_id"]: "haul_num", age_group_cols["nasc_id"]: "nasc"}
# ---- Update if the stratum is not equal to INPFC
if settings_dict["transect"]["stratum"] != "inpfc":
    name_conversion_key.update({age_group_cols["stratum_id"]: stratum_col})

# Rename columns
# ---- Extract NASC data
nasc_data = input_dict["acoustics"]["nasc_df"].copy()
# ---- Change names
nasc_data.rename(columns=name_conversion_key, inplace=True)

# Correct the acoustic survey transect intervals
nasc_interval_df = correct_transect_intervals(nasc_data)

distributions_dict, proportions_dict, TS_L_parameters, settings_dict = (
    input_dict["biology"]["distributions"],
    analysis_dict["biology"]["proportions"],
    configuration_dict["TS_length_regression_parameters"]["pacific_hake"],
    settings_dict
)