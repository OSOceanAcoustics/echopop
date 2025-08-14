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
from echopop.biology import age1_metric_proportions, impute_kriged_values, reallocate_kriged_age1
from echopop.spatial.transect import correct_transect_intervals
from echopop.analysis import process_transect_data

survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
# survey = Survey(init_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
#                 survey_year_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data(ingest_exports="echoview")
survey.load_survey_data()
survey.transect_analysis()
self = survey
stratum = "ks"
exclude_age1 = True
species_id = 22500
input_dict, analysis_dict, configuration_dict, settings_dict = self.input, self.analysis["transect"], self.config, self.analysis["settings"]


stratify_by = "stratum_ks"
df_average_weight = df_averaged_weight["all"]
ts_length_regression_parameters={"slope": 20., "intercept": -68.}
# number_proportions=dict_df_number_proportion["aged"]
number_proportions = dict_df_number_proportion
length_threshold_min=10.0
weight_proportion_threshold=1e-10
weight_proportions=dict_df_weight_proportion["aged"]
include_filter = {"age_bin": [1]}

###########
df_nasc = df_nasc_no_age1.copy()
number_proportions = dict_df_number_proportion
group_by = ["sex"]
stratify_by = ["stratum_ks"]
exclude_filter = {"sex": "unsexed"}
dataset = df_nasc.copy()

def set_abundance(
    dataset: pd.DataFrame,
    stratify_by: List[str] = [],
    group_by: List[str] = [],
    exclude_filter: Dict[str, str] = {},
    number_proportions: Optional[Dict[str, pd.DataFrame]] = None
):

    # If no grouping, run the simple abundance calculation    
    dataset["abundance"] = dataset["area_interval"] * dataset["number_density"]
    
    # Compute grouped values, if needed
    if number_proportions is not None:      
        # ---- Set the index
        dataset.set_index(stratify_by, inplace=True)
        # ---- Create grouped table from number proportions
        grouped_proportions = utils.create_grouped_table(
            number_proportions,
            group_cols = stratify_by + group_by,
            strat_cols = group_by,
            index_cols = stratify_by,
            value_col = "proportion_overall"
        )
        # ---- Apply exclusion filter, if required
        grouped_proportions_excl = utils.apply_filters(grouped_proportions,
                                                       exclude_filter=exclude_filter)
        # ---- Refine if no grouping
        if len(group_by) == 0:
            grouped_proportions_excl = grouped_proportions_excl["proportion_overall"]
            number_density_vals = dataset["number_density"].values
            abundance_vals = dataset["abundance"].values
        else:
            number_density_vals = dataset["number_density"].values[:, None]
            abundance_vals = dataset["abundance"].values[:, None]     
        # ---- Reindex the table
        grouped_proportions_ridx = grouped_proportions_excl.reindex(dataset.index)
        # ---- Compute number density
        grouped_number_density = (
            number_density_vals * grouped_proportions_ridx
        )
        # ---- Compute abundance
        grouped_abundance = (
            abundance_vals * grouped_proportions_ridx
        )
        # ---- Add the number densities to the dataset
        dataset[
            grouped_number_density.columns.map(lambda c: f"number_density_{c}")
        ] = grouped_number_density.values
        # ---- Add abundances to the dataset
        dataset[
            grouped_abundance.columns.map(lambda c: f"abundance_{c}")
        ] = grouped_abundance.values
        # ---- Reset the index
        dataset.reset_index(inplace=True)

#####
# df_nasc = df_nasc_no_age1.copy()
group_by = ["sex"]
stratify_by = ["stratum_ks"]
df_average_weight = df_averaged_weight["all"].copy()
df_averaged_weight.loc[:, "all"]
def matrix_multiply_grouped_table(
    dataset: pd.DataFrame,
    table: pd.DataFrame,      
    variable: str, 
    output_variable: str,
    group: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    # Create pattern for column filtering
    prefix_pattern = variable + "_"
    
    # Get the overlapping columns
    variable_columns = dataset.filter(like=prefix_pattern).columns

    # Gather the suffixes corresponding to the target group    
    suffixes = variable_columns.str.replace(prefix_pattern, "", regex=False).to_list()

    # Reindex table
    table_ridx = table.reindex(dataset.index)
    
    # Apply an inclusion filter
    target_groups = utils.apply_filters(table_ridx, include_filter={group: suffixes})

    # Get columns that also exist for dataset
    variable_overlap = [prefix_pattern + col for col in target_groups.columns]

    # Reindex the target grouped table
    target_groups_idx = target_groups.reindex(dataset.index)

    # Run the multiplication
    table_matrix = dataset.filter(variable_overlap).to_numpy() * target_groups_idx

    # Set up column names
    column_map = target_groups_idx.columns.map(lambda c: f"{output_variable}_{c}")

    # Add the output variables
    dataset[column_map] = table_matrix.values

    # Calculate the remainder comprising the ungrouped values
    remainder = dataset[variable] - dataset[variable_overlap].sum(axis=1)

    # Calculate the output variable for the ungrouped/excluded values
    remainder_matrix = remainder * table_ridx["all"]

    # Compute the overall output variable
    dataset[output_variable] = dataset[column_map].sum(axis=1) + remainder_matrix

dataset = df_nasc.copy()
table = df_average_weight.copy()
group = table.columns.names[0]
variables = ["number_density", "abundance"]
variable = "number_density"

def set_biomass(
    dataset: pd.DataFrame,
    stratify_by: List[str] = [],
    group_by: List[str] = [],
    df_average_weight: Optional[Union[pd.DataFrame, float]] = None
):

    # Set the index for the dataset
    dataset.set_index(stratify_by, inplace=True)
    
    # Handle stratification and weight alignment
    if isinstance(df_average_weight, (pd.DataFrame, pd.Series)):
        # ---- Ensure weights are properly aligned with the associated dataset
        if (
            hasattr(df_average_weight, "columns") and 
            not set(df_average_weight.index.names).intersection(stratify_by)
        ):
            df_average_weight.set_index(stratify_by, inplace=True)
        elif isinstance(df_average_weight, pd.Series):
            df_average_weight = df_average_weight.to_frame("all")
    else:
        # ---- Create associated Series from a single float
        df_average_weight = pd.DataFrame(
            {
                "all": df_average_weight
            },
            index=dataset.index
        )
        
    # If grouped
    if len(group_by) > 0:
        # ---- Compute the biomass densities across groups
        matrix_multiply_grouped_table(dataset, 
                                      df_average_weight, 
                                      group_by[0], 
                                      "number_density", 
                                      "biomass_density")
        # ---- Compute the biomass densities across groups
        matrix_multiply_grouped_table(dataset, 
                                      df_average_weight, 
                                      group_by[0], 
                                      "abundance", 
                                      "biomass")
    # Ungrouped
    else:
        # ---- Compute biomass densities
        dataset["biomass_density"] = dataset["number_density"] * df_average_weight["all"]
        # ---- Compute biomass
        dataset["biomass"] = dataset["abundance"] * df_average_weight["all"]
        
    # Reset the index
    dataset.reset_index(inplace=True)

######
(
    nasc_biology_df,
    fitted_weight_dict,
    settings_dict,
    population_dict,
    strata_adult_proportions_df
) = (
    nasc_to_biology,
    analysis_dict["biology"]["weight"],
    settings_dict,
    analysis_dict["biology"]["population"],
    strata_adult_proportions
)
######
reapportion_dict = {"nasc": age1_nasc_proportions, 
                    "abundance": age1_number_proportions,
                    "biomass": age1_weight_proportions}
dataset = df_nasc.copy()
#####

def partition_transect_data(
    dataset: pd.DataFrame,
    partition_dict: Dict[str, Union[pd.DataFrame, pd.Series]],
) -> pd.DataFrame:
    
    # Create copy
    dataset = dataset.copy()
    
    # Get the index names
    index_names = list(set().union(*[set(df.index.names) for df in partition_dict.values()]))
    
    # Set the index of the input dataset
    dataset.set_index(index_names, inplace=True)
    
    # NASC, if present
    if "nasc" in partition_dict:
        dataset["nasc"] = dataset["nasc"] * (1 - partition_dict["nasc"].reindex(dataset.index))
        
    # Abundance and number density, if present
    if "abundance" in partition_dict:
        # ---- Get the inverse proportions
        abundance_proportions = (1 - partition_dict["abundance"].reindex(dataset.index))
        # ---- Map the appropriate columns for abundance
        abundance_names = dataset.filter(like="abundance").columns
        # ---- Adjust abundances
        dataset[abundance_names] = (abundance_proportions * dataset[abundance_names].T).T
        # ---- Map the appropriate columns for number density
        number_density_names = dataset.filter(like="number_density").columns
        # ---- Adjust number densities
        dataset[number_density_names] = (abundance_proportions * dataset[number_density_names].T).T
        
    # Biomass and biomass density, if present
    if "biomass" in partition_dict:
        # ---- Get the inverse proportions
        biomass_proportions = (1 - partition_dict["biomass"].reindex(dataset.index))
        # ---- Map the appropriate columns for biomass and biomass density
        biomass_names = dataset.filter(like="biomass").columns
        # ---- Adjust biomass
        dataset[biomass_names] = (biomass_proportions * dataset[biomass_names].T).T    
        
    # Return the repartitioned dataset
    return dataset

################
dataset = df_partitioned_no_age1.copy()
proportions=dict_df_number_proportion["aged"]
stratify_by=["stratum_ks"]
group_by=["sex"]
index=["length_bin"]
columns=["age_bin"]

group_by=["sex", "age_bin", "length_bin"]
variable="biomass"

# Sum over each group
dataset_pvt = dataset.pivot_table(
    index=stratify_by, values=variable, aggfunc="sum", observed=False
)

# Parse the additional columns that are required for grouping
proportions_group_columns = [
    c for c in (list(proportions.index.names) + list(proportions.columns)) 
    if c in group_by + index + columns
]

# Create pivot table
proportions_pvt = utils.create_pivot_table(
    proportions,    
    strat_cols=group_by,
    index_cols=list(set(proportions_group_columns).difference(group_by)) + stratify_by,
    value_col="proportion",
)

# Check if "all" exists for the `group_by` variable
if "all" not in proportions[group_by]:
    # ---- Add 'all'
    proportions_pvt["all"] = proportions_pvt.sum(axis=1)

# Re-pivot based on arguments
proportions_pvt.stack(group_by).unstack(columns + stratify_by)

# Pivot
utils.create_pivot_table(
    proportions,    
    strat_cols=group_by,
    index_cols=list(set(proportions_group_columns).difference(group_by)) + stratify_by,
    value_col="proportion",
)


# Convert to DataFrame(s) to pivot table(s)
proportions_grouped_pvt = {
    k: (
        df
        if utils.is_pivot_table(df)
        else utils.create_pivot_table(
            df,
            index_cols=proportions_group_columns[k],
            strat_cols=stratify_by,
            value_col="proportion",
        )
    )
    for k, df in proportions.items()
}

# Distribute the variable over each table
apportioned_grouped_pvt = {
    k: df.mul(dataset_pvt[variable]).fillna(0.0) for k, df in proportions_grouped_pvt.items()
}

# Re-pivot


TEST = aged_apportioned_biomass.pivot_table(index=["sex", "length_bin"], columns=["age_bin", "stratum_num"], values="biomass_aged", observed=False)
TEST.loc["all"].sum(axis=0).unstack("stratum_num").sum()
TOIT = apportioned_grouped_pvt["aged"].stack("stratum_ks").to_frame("value").reset_index().pivot_table(index=["sex", "length_bin"], 
                                                                                                       columns=["age_bin", "stratum_ks"], 
                                                                                                       values="value", 
                                                                                                       observed=False)

TOIT_ALL = TOIT.loc["male"].sum(axis=0).unstack("stratum_ks").sum() + TOIT.loc["female"].sum(axis=0).unstack("stratum_ks").sum()
TOIT_ALL.sum()

TEST.loc["all"].sum(axis=0).unstack("stratum_num").sum()
proportions_grouped_pvt["aged"]

apportioned_grouped_pvt["aged"].unstack(["age_bin"]).sum().sum() / dataset_pvt.sum()
aged_weight_proportions_all.groupby(["stratum_num"])["weight_proportions"].sum()


nasc_biomass.loc[nasc_biomass.sex == "all", "biomass"].sum()
aged_apportioned_biomass_tbl.sum(axis=1).loc["all"].sum()


aged_apportioned_biomass_tbl.sum(axis=1).loc["male"] - apportioned_grouped_pvt["aged"].unstack(["age_bin"]).sum(axis=1).loc[:, "male"]
# ==================================================================================================
# Get proportions for each stratum specific to age-1
# --------------------------------------------------

# Age-1 NASC proportions
age1_nasc_proportions = get_proportions.get_nasc_proportions_slice(
    number_proportions=dict_df_number_proportion["aged"],
    stratify_by=["stratum_ks"],
    ts_length_regression_parameters={"slope": 20., 
                                     "intercept": -68.},
    include_filter = {"age_bin": [1]}
)

# Age-1 number proportions
age1_number_proportions = get_proportions.get_number_proportions_slice(
    number_proportions=dict_df_number_proportion["aged"],
    stratify_by=["stratum_ks"],
    include_filter = {"age_bin": [1]}
)

# Age-1 weight proportions
age1_weight_proportions = get_proportions.get_weight_proportions_slice(
    weight_proportions=dict_df_weight_proportion["aged"],
    stratify_by=["stratum_ks"],
    include_filter={"age_bin": [1]},
    number_proportions=dict_df_number_proportion,
    length_threshold_min=10.0,
    weight_proportion_threshold=1e-10
)

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
