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

survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data()
survey.load_survey_data()
survey.transect_analysis(exclude_age1=False)
survey.kriging_analysis()

self = survey
analysis_dict = self.analysis
kriged_mesh = self.results["kriging"]["mesh_results_df"]
settings_dict = self.analysis["settings"]["kriging"]
mesh_results_df, transect_dict, settings_dict = kriged_mesh, analysis_dict["transect"], settings_dict 

# Get the stratum column name
stratum_col = settings_dict["stratum_name"]

# Set index for mesh results
mesh_results_df.set_index([stratum_col], inplace=True)

# Get the stratified mean weights
weight_strata = transect_dict["biology"]["weight"]["weight_stratum_df"].copy()
# ---- Sub-select and set index
weight_strata = (
    weight_strata[weight_strata["sex"] == "all"]
    .set_index([stratum_col])
    .reindex(mesh_results_df.index)
)

# Get the average `sigma_bs` per stratum
strata_mean_sigma_bs = (
    (transect_dict["acoustics"]["sigma_bs"]["strata_mean_df"].copy())
    .set_index([stratum_col])
    .reindex(mesh_results_df.index)
)
(transect_dict["acoustics"]["sigma_bs"]["strata_mean_df"].copy()).set_index([stratum_col])

invert_hake.sigma_bs_strata

dict_df_weight_proportion["aged"].unstack("sex").sum(axis=0).unstack("sex")

A = dict_df_weight_proportion["aged"].T.stack("sex", future_stack=True).sum(axis=1).unstack("sex")
AB = aged_props

apportion.standardize_kriged_estimates = standardize_kriged_estimates


df = population_table
DD = dict_kriged_abundance_table["unaged"].sum(axis=1).unstack("sex")
AA = unaged_pivot["abundance_apportioned_unaged"]

AA - DD
AA.sum().sum()
DD.sum().sum()
##########
mesh_data_df = df_kriged_results.copy()
variable = "abundance"
stratify_by = ["stratum_ks"]
group_by = ["sex", "age_bin", "length_bin"]
mesh_proportions_link={"geostratum_ks": "stratum_ks"}
proportions = dict_df_number_proportion["aged"].copy()
# proportions = dict_df_number_proportion.copy()
proportions = dict_df_weight_proportion.copy()
dict_kriged_abundance_table 
dict_kriged_biomass_table

# population_table = dict_kriged_abundance_table.copy()
population_table = dict_kriged_biomass_table.copy()
reference_table = ["aged"]
impute_variable = ["age_bin"]

##########
# Get the table indices
table_indices = {
    k: list(df.index.names)
    for k, df in population_table.items()
}    

# Get the shared index names across the reference tables
reference_index_names = list(
    set.intersection(*(set(v) 
                    for k, v in table_indices.items() 
                    if k in reference_table)) 
)

# Get the indices for each reference table required for summation that must be unstacked
reference_stack_indices = {
    k: list(set(reference_index_names).difference(set(df.index.names)))
    for k, df in population_table.items()
    if k in reference_table
}

# Consolidate the reference tables into a list
reference_list = [
    population_table[k].unstack(reference_stack_indices[k]) for k in reference_table
]

# Sum across all of the references
reference_df = reduce(lambda a, b: a.add(b, fill_value=0), reference_list)

# Get the indices for each target table
target_stack_indices = {
    k: list(set(reference_index_names).difference(set(df.index.names)))
    for k, df in population_table.items()
    if k not in reference_table    
}

# Reorient the columns accordingly 
target_proportions = {k: df.sum(axis=1).unstack(group_by) 
                      for k, df in population_table.items() 
                      if k not in reference_table}

population_table["unaged"].sum(axis=1).unstack(group_by)

# Standardize the proportions
standardized_proportions = {
    k: (
        df *
        reference_df.sum(axis=1).unstack(target_stack_indices[k] + group_by) /
        reference_df.unstack(target_stack_indices[k]).sum(axis=1).unstack(group_by)
    ).fillna(0.)
    for k, df in target_proportions.items()
}

########################
k = "unaged"
df = population_table[k]

df.sum(axis=1).unstack(group_by)
reference_df.sum(axis=1).unstack(target_stack_indices[k] + group_by)

table_input = proportions_standardized["unaged"].copy()
aged_pivot.unstack("sex")["biomass_apportioned"]
aged_pivot["biomass_apportioned"].sum(axis=1).unstack("sex")

aged_age_length_table = aged_pivot["biomass_apportioned"]
unaged_length_table = unaged_pivot["biomass_apportioned_unaged"]
aged_length_totals = aged_length_biomass_totals
unaged_apportioned_table = unaged_apportioned_biomass_values

original_table = target_proportions.copy()
standardized_table = standardized_proportions.copy()

########################
reference_stk = reference_df.sum(axis=1).unstack(impute_variable).sum(axis=1).unstack(group_by)
reference_group_stk = reference_df.sum(axis=1).unstack(group_by + impute_variable)

# Apply similar transformation to the standardized proportions
table_trans = {
    k: df.stack(group_by, future_stack=True).sum(axis=1).unstack(group_by)
    for k, df in standardized_table.items()
}

# Get the mask for all 0.0's and non-zeros
ref_zero_mask = reference_stk == 0.
ref_nonzero_mask = reference_stk != 0.

# Create translation for row numbers
interval_to_numeric = {interval: i for i, interval in enumerate(reference_stk.index)}

# Gather the indices for each column
ref_zero_indices = {col: reference_stk.index[ref_zero_mask[col]].tolist()
                    for col in reference_stk.columns}
ref_nonzero_indices = {col: reference_stk.index[ref_nonzero_mask[col]].tolist() 
                       for col in reference_stk.columns}

# Convert to row numbers
ref_zero_rows = {col: np.array([interval_to_numeric[ival] for ival in intervals])
                 for col, intervals in ref_zero_indices.items()}
ref_nonzero_rows = {col: np.array([interval_to_numeric[ival] for ival in intervals])
                 for col, intervals in ref_nonzero_indices.items()}

# Apply the non-zero reference indices to the original target tables
table_nonzeros_mask = {
    k: {
        col: df[col].loc[ref_zero_indices[col]] != 0.
        for col in df.columns
    }
    for k, df in original_table.items()
}

# Get the actual indices of the masked values
nonzero_reference_to_table_indices = {
    k: {
        col: np.array(ref_zero_indices[col])[table_nonzeros_mask[k][col]]
        for col in table_nonzeros_mask[k]
    }
    for k in table_nonzeros_mask
}

# Convert to numeric indices for compatibility with `iloc`
nonzero_reference_to_table_rows = {
    k: {
        col: np.array([interval_to_numeric[ival] for ival in intervals])
        for col, intervals in nonzero_reference_to_table_indices[k].items()
    }
    for k in nonzero_reference_to_table_indices
}

# Iterate through the tables to impute
name = "unaged"
table = table_nonzeros_mask[name].copy()

for name, table in table_nonzeros_mask.items():
    # ---- Get the column name indices
    column_indices = target_stack_indices[name] + group_by
    # ---- Restack the standardized dataset    
    standardized_table_sub = standardized_table[name].copy()
    standardized_table_sub.columns = standardized_table_sub.columns.reorder_levels(list(reversed(column_indices)))
    # ---- Get the nearest-neighbor rows and recompute the indices
    imputed_rows = {
        col: arr[
            np.argmin(
                np.abs(ref_zero_rows[col][table[col]][:, np.newaxis] - ref_nonzero_rows[col]), 
            axis=1)
        ]
        for col, arr in ref_nonzero_rows.items()
    }
    # ---- Impute to replace these values
    imputed_values = {
        col: (
            original_table[name][col].loc[nonzero_reference_to_table_indices[name][col]].to_numpy() * 
            reference_group_stk.iloc[imputed_rows[col]][col].T / 
            reference_stk.iloc[imputed_rows[col]][col]
        ).T
        for col in table.keys()
    }
    # ---- Update the standardized values
    for col in table.keys():
        standardized_table_sub.iloc[
            nonzero_reference_to_table_rows[name][col],
            standardized_table_sub.columns.get_loc(col)
        ] = (
            imputed_values[col]
        )
    # ---- Update the original dictionary
    standardized_table[name] = standardized_table_sub
    
standardized_table_sub.loc[5, "male"]
unaged_apportioned_biomass_values.loc[5]

kriged_full_table.set_index(["sex", "length_bin"]).loc["male", 5]
standardized_table[name].loc[:, "male"].loc[:, 1]
standardized_table_sub.loc[:, "male"].loc[:, 1]
unaged_apportioned_biomass_values.loc[:, 1]
unaged_apportioned_table.loc[:, 1]
unaged_apportioned_table.loc[6]
(standardized_table[name] - standardized_table_sub).min(axis=1)
standardized_table[name]
standardized_table["unaged"].stack(standardized_table["unaged"].columns.names, future_stack=True)
standardized_table_sub.loc[nonzero_reference_to_table_indices[name][col], col]
unaged_values_pvt.iloc[
    female_nonzero_unaged_idx, unaged_values_pvt.columns.get_loc("female")
] = 0
unaged_apportioned_table
unaged_apportioned_biomass_values - unaged_apportioned_table


standardized_table_sub.loc[nonzero_reference_to_table_indices[name][col]][col] = (
    original_table[name][col].loc[nonzero_reference_to_table_indices[name][col]].to_numpy() * 
    reference_group_stk.iloc[imputed_rows[col]][col].T / 
    reference_stk.iloc[imputed_rows[col]][col]
).T.to_numpy()
nonzero_reference_to_table_indices[name][col]
standardized_table[name].loc[6]
standardized_table_sub.loc[6]
df = standardized_table[name]
df.columns.reorder_levels(["sex", "age_bin"])
type(df)

df = pd.DataFrame({
        'A':list('abcdef'),
         'B':[4,5,4,5,5,4],
         'C':[7,8,9,4,2,3],
         'D':[1,3,5,7,1,0],
         'E':[5,3,6,9,2,4],
         'F':list('aaabbb')
}).set_index(['A','B','C','D'])


df = df.reorder_levels(['B','C','A','D']).sort_index()

standardized_table["unaged"].loc[:, [1, 2, 3]].loc["female"]

TEWST = standardized_table["unaged"].columns.names
standardized_table["unaged"].xs(tuple(TEWST), level=["(0.5, 1.5]", "male"], axis=1)
standardized_table["unaged"].xs(col, level=group_by[0], axis=1)
standardized_table["unaged"].xs((), level="age_bin", axis=1)
result = {
    k: {
        key: df.xs(key=key, level=group_by, axis=1)
        for key in {tuple(col[i] for i in level_pos) for col in df.columns}
    }
    for k, df in table_nonzeros_mask.items()
}
standardized_table["unaged"].columns.get_loc("male")
standardized_table["unaged"].iloc[ref_nonzero_rows["male"]]

table[col] = table[col].apply(lambda x: True)
arr = ref_nonzero_rows[col]
level = columns.names.index(level_name)
mask = columns.get_level_values(level) == value
return columns[mask]

interval_to_numeric[np.array(ref_zero_indices[col])[table_nonzeros_mask[k][col]]]
TEST_C[np.argmin(np.abs(TEST_A[TEST_B][:, np.newaxis] - TEST_C), axis=1)]
TEST_A = ref_zero_rows["male"]
TEST_B = table_nonzeros_mask["unaged"]["male"]
TEST_C = ref_nonzero_rows["male"]

test = nonzero_reference_to_table_indices["unaged"]["male"]

male_zero_aged == TEST_A
male_nonzero_unaged == TEST_B

male_zero_aged[male_nonzero_unaged]
TEST_A[TEST_B]

if len(imputation_indices) > 0:
    # ---- Find the nearest-neighbor non-zero index
    


table_nonzeros_mask.keys()
table_nonzeros_mask.values()

k = "unaged"
col = "male"
np.array(ref_zero_indices[col])[table_nonzeros_mask[k][col]]

k = "unaged"
col = "female"
np.array(ref_zero_indices[col])[table_nonzeros_mask[k][col]]

mask = (df.loc[pd.Index(ref_zero_indices[col]), "male"] != 0.)
df.loc[pd.Index(ref_zero_indices[col])[mask]]
df

ref_zero_indices["male"].
ref_zero_indices["male"].to_numpy()[table_nonzeros_mask["unaged"]["male"].to_numpy()]

np.array(unaged_length_table.index.values)[male_nonzero_unaged_idx]

k = "unaged"
col = "male"
df = original_table[k]


unaged_length_table["male"].iloc[male_zero_aged].shape


dict_kriged_biomass_table["unaged"].sum(axis=1).unstack("sex")


ref_zero_intervals = np.array([interval.mid for interval in ref_zero_indices[col]])
ref_nonzero_intervals = np.array([interval.mid for interval in ref_nonzero_indices[col]])

zero_mid := np.array([interval.mid for interval in zero_intervals]),
nonzero_mid := np.array([interval.mid for interval in nonzero_intervals])

####
def nearest_neighbor_imputation(
    apportioned_tables,  # Dictionary of DataFrames to impute
    zero_indices_dict,  # Dictionary of zero indices: {'dataset': {'column': [indices]}}
    length_tables,  # Dictionary of length tables: {'dataset': DataFrame}
    age_length_tables,  # Dictionary of age-length tables: {'dataset': DataFrame}
    length_totals,  # Dictionary of length totals: {'dataset': DataFrame or Series}
):
    # Stack all datasets into a single dictionary comprehension
    return {
        dataset: (
            lambda df, zero_dict, length_table, age_length_table, length_total: (
                # Vectorized imputation for all columns
                [
                    (
                        # Get zero and nonzero indices for this column
                        zero_intervals := zero_dict[col],
                        nonzero_intervals := df.index[df[col] != 0.0].to_numpy(),
                        # Skip if no nonzero indices
                        None if len(nonzero_intervals) == 0 or len(zero_intervals) == 0 else (
                            # Convert intervals to midpoints
                            zero_mid := np.array([interval.mid for interval in zero_intervals]),
                            nonzero_mid := np.array([interval.mid for interval in nonzero_intervals]),
                            # Find nearest nonzero index for each zero index
                            imputed_idx := np.argmin(np.abs(zero_mid[:, None] - nonzero_mid), axis=1),
                            imputed_nonzero_intervals := nonzero_intervals[imputed_idx],
                            # Get row positions in df for zero and imputed indices
                            zero_row_idx := [df.index.get_loc(idx) for idx in zero_intervals],
                            imputed_row_idx := [df.index.get_loc(idx) for idx in imputed_nonzero_intervals],
                            col_idx := df.columns.get_loc(col),
                            # Perform the imputation calculation
                            df.values[zero_row_idx, col_idx] := (
                                length_table[col].iloc[zero_row_idx].to_numpy()
                                * age_length_table.loc[col].iloc[imputed_row_idx].T
                                / length_total[col].iloc[imputed_row_idx]
                            ).T
                        )
                    )
                    for col in zero_dict
                ],
                # Return the updated DataFrame
                df
            )[1]
        )(
            apportioned_tables[dataset],
            zero_indices_dict[dataset],
            length_tables[dataset],
            age_length_tables[dataset],
            length_totals[dataset],
        )
        for dataset in zero_indices_dict
    }
col = "male"
df = table_trans["unaged"]

table_input["unaged"].index[table_input["unaged"] == 0.]

cal_vals = (aged_pivot["biomass_apportioned"], 
            unaged_pivot["biomass_apportioned_unaged"],
            aged_length_biomass_totals,
            unaged_apportioned_biomass_values,
            settings_dict,
            "biomass")

(aged_age_length_table,
unaged_length_table,
aged_length_totals,
unaged_apportioned_table,
settings_dict,
variable) = cal_vals




analysis_dict, kriged_mesh, settings_dict = self.analysis, self.results["kriging"]["mesh_results_df"], self.analysis["settings"]["kriging"]


# Sum the kriged weights for each stratum
# ---- Extract stratum column name
stratum_col = settings_dict["stratum_name"]
# ---- Extract the biological variable (independent of area)
# biology_col = settings_dict["variable"].replace("_density", "")
# ---- Sum abundance for each stratum
summed_abundance = kriged_mesh.groupby([stratum_col], observed=False)["abundance"].sum()
# ---- Sum biomass for each stratum
summed_biomass = kriged_mesh.groupby([stratum_col], observed=False)["biomass"].sum()

# Extract the weight proportions from the analysis object
proportions_dict = analysis_dict["transect"]["biology"]["proportions"]

# Prepare the number/abundance proportions
# ---- Aged
aged_abundance_proportions = proportions_dict["number"]["aged_length_proportions_df"].copy()
# ---- Unaged
unaged_abundance_proportions = proportions_dict["number"]["unaged_length_proportions_df"].copy()

# Prepare the weight/biomass proportions
# ---- Aged
aged_biomass_proportions = proportions_dict["weight"]["aged_weight_proportions_df"].copy()
# ---- Unaged
unaged_biomass_proportions = proportions_dict["weight"]["unaged_weight_proportions_df"].copy()
# ---- Aged-unaged sexed weight proportions
unaged_sex_biomass_proportions = proportions_dict["weight"][
    "aged_unaged_sex_weight_proportions_df"
].copy()[[stratum_col, "sex", "weight_proportion_overall_unaged"]]

# Apportion abundances
# ---- Pivot unaged
unaged_abundance_proportions_pvt = unaged_abundance_proportions.pivot_table(
    columns=["sex", "length_bin"],
    index=[stratum_col],
    values="proportion_number_overall_unaged",
    observed=False,
)
# ---- Apportion the abundances
unaged_apportioned_abundance = (
    ((summed_abundance * unaged_abundance_proportions_pvt.transpose()).fillna(0.0))
    .stack()
    .reset_index(name="abundance_apportioned_unaged")
)
# ---- Set index for latter merging
unaged_apportioned_abundance.set_index([stratum_col, "sex", "length_bin"], inplace=True)
# ---- Pivot aged
aged_abundance_proportions_pvt = aged_abundance_proportions.pivot_table(
    columns=["sex", "age_bin", "length_bin"],
    index=[stratum_col],
    values="proportion_number_overall_aged",
    observed=False,
)
# ---- Apportion the abundances
aged_apportioned_abundance = (
    ((summed_abundance * aged_abundance_proportions_pvt.transpose()).fillna(0.0))
    .stack()
    .reset_index(name="abundance_apportioned")
)
# ---- Set index for latter merging
aged_apportioned_abundance.set_index(
    [stratum_col, "sex", "age_bin", "length_bin"], inplace=True
)

# Compute the apportioned unaged kriged biological values per stratum
# ---- Merge the unaged proportions
unaged_sexed_apportioned = unaged_biomass_proportions.merge(unaged_sex_biomass_proportions)
# ---- Set index to stratum, sex, length_bin columns
unaged_sexed_apportioned.set_index([stratum_col, "sex", "length_bin"], inplace=True)
# ---- Merge
unaged_sexed_apportioned["abundance_apportioned_unaged"] = unaged_apportioned_abundance
# ---- Reset the index
unaged_sexed_apportioned.reset_index(["sex", "length_bin"], inplace=True)
# ---- Set the index based on `summed_biomass`
summed_biomass_indexed = summed_biomass.reindex(unaged_sexed_apportioned.index)
# ---- Append the stratum-aggregated biomass values
unaged_sexed_apportioned["biomass_apportioned_unaged"] = (
    unaged_sexed_apportioned["weight_proportion"]
    * unaged_sexed_apportioned["weight_proportion_overall_unaged"]
    * summed_biomass_indexed
)

# Distribute biological values over the overall proportions (i.e. relative to aged and unaged
# fish) for aged fish
# ---- Set index to stratum column
aged_biomass_proportions.set_index([stratum_col, "sex", "age_bin", "length_bin"], inplace=True)
# ---- Compute the distributed abundance values
aged_biomass_proportions["abundance_apportioned"] = aged_apportioned_abundance
# ---- Reset the index
aged_biomass_proportions.reset_index(["sex", "age_bin", "length_bin"], inplace=True)
# ---- Compute the distributed biomass values
aged_biomass_proportions["biomass_apportioned"] = (
    aged_biomass_proportions["weight_proportion_overall"] * summed_biomass
).fillna(0.0)

# Distribute the aged biological distributions over unaged length distributions to estimate
# aged distributions
# ---- Pivot aged data
aged_pivot = aged_biomass_proportions.reset_index().pivot_table(
    index=["sex", "length_bin"],
    columns=["age_bin"],
    values=["abundance_apportioned", "biomass_apportioned"],
    aggfunc="sum",
    observed=False,
)
# ---- Calculate the total biomass values for each sex per length bin
aged_length_biomass_totals = aged_pivot["biomass_apportioned"].sum(axis=1).unstack("sex")
# ---- Pivot unaged data
unaged_pivot = unaged_sexed_apportioned.reset_index().pivot_table(
    index=["length_bin"],
    columns=["sex"],
    values=["abundance_apportioned_unaged", "biomass_apportioned_unaged"],
    aggfunc="sum",
    observed=False,
)
# ---- Calculate the new unaged biomass values distributed over age
unaged_apportioned_biomass_values = (
    unaged_pivot["biomass_apportioned_unaged"]
    * aged_pivot.unstack("sex")["biomass_apportioned"]
    / aged_length_biomass_totals
).fillna(0)
