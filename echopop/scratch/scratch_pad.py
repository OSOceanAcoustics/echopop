import abc
import numpy as np
import pandas as pd
from typing import Callable, Union, Dict, List, Optional, Any
from functools import reduce
import pytest

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs
from echopop.nwfsc_feat import utils
from typing import Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate

#######

########
from echopop.survey import Survey
from echopop.biology import impute_kriged_values, reallocate_kriged_age1

# survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
#                 survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey = Survey(init_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis(exclude_age1=False)
survey.kriging_analysis()
self = survey
analysis_dict, kriged_mesh, settings_dict = self.analysis, self.results["kriging"]["mesh_results_df"], self.analysis["settings"]["kriging"]

table, settings_dict, variable = kriged_full_table, settings_dict, "biomass_apportioned"

###
apportion.combine_population_tables = combine_population_tables
utils.apply_filters = apply_filters
apportion.redistribute_population_table = redistribute_population_table

###
from echopop.nwfsc_feat import utils

population_table = df_kriged_biomass_table.copy()
exclusion_filter = {"age_bin": 1}
group_by = ["sex"]
redistribute = True
###

# Find any columns that are not in the group_by list
# ---- Get column names
column_names = population_table.columns.names
# ---- Identify extra columns that are not in the group_by list
extra_columns = [col for col in column_names if col not in group_by]
# ---- Stack the population table
stacked_table = population_table.stack(extra_columns, future_stack=True)

# Apply inverse of exclusion filter to get the values being excluded
excluded_grouped_table = utils.apply_filters(stacked_table, 
                                             include_filter=exclusion_filter)

# Replace the excluded values in the full table with 0.
filtered_grouped_table = utils.apply_filters(stacked_table, 
                                             exclude_filter=exclusion_filter,
                                             replace_value=0.)
# filtered_grouped_table1 = filtered_grouped_table.copy()

# Get the sums for each group across the excluded and filtered tables
# ---- Excluded
excluded_grouped_sum = excluded_grouped_table.sum()
# ---- Filtered/included
filtered_grouped_sum = filtered_grouped_table.sum()

# Get the redistributed values that will be added to the filtered table values
adjustment_table = filtered_grouped_table * excluded_grouped_sum / filtered_grouped_sum

# Add the adjustments to the filtered table
filtered_grouped_table += adjustment_table

# Check 
if np.any(filtered_grouped_table.sum() - stacked_table.sum() > 1e-6):
    # ---- If the sums do not match, raise a warning
    check_sums = filtered_grouped_table.sum() - stacked_table.sum() > 1e-6
    # ---- Raise a warning with the indices where the sums do not match
    warnings.warn(
        f"The sums of the table with the redistributed estimates do not match the original table "
        f"filtered table do not match the original table for indices: "
        f"{', '.join(check_sums[check_sums].index.tolist())}"
    )

# Restore the original column structure
redistributed_table = (
    filtered_grouped_table.unstack(extra_columns)
    .reorder_levels(column_names, axis=1)
)



df = stacked_table.copy()
exclude_filter=exclusion_filter
include=False
replace_value=0.
filter_dict = exclude_filter

TEST = filtered_grouped_table - filtered_grouped_table1
TEST["female"].loc[lambda x: x != 0.]

stacked_table.sum() - filtered_grouped_table.sum()

excluded_table
df = stacked_table.copy()
filter_dict = exclusion_filter
include=False
excluded_table

df = population_table.copy()
exclude_filter=exclusion_filter
include_filter: Optional[Dict[str, Any]] = None
# Apply the inverse of exclusion filter

filter_dict = exclude_filter
filter_dict = {"length_bin": [27., 28., 29.]}

replace_value = 0.

kriged_full_table.pivot_table(
    index=["length_bin"],
    columns=["age_bin", "sex"],
    values="biomass_apportioned",
)

standardized_table_sub.loc[5, "male"]
unaged_apportioned_biomass_values.loc[5]

kriged_full_table.set_index(["sex", "length_bin"]).loc["male", 5]
standardized_proportions[name].loc[:, "male"].loc[:, 1]
standardized_table_sub.loc[:, "male"].loc[:, 1]
unaged_apportioned_biomass_values.loc[:, 1]
unaged_apportioned_table.loc[:, 1]
unaged_apportioned_table.loc[6]
(standardized_table[name] - standardized_table_sub).min(axis=1)
standardized_table[name]
standardized_table["unaged"].stack(standardized_table["unaged"].columns.names, future_stack=True)
standardized_table_sub.loc[nonzero_reference_to_table_indices[name][col], col]
