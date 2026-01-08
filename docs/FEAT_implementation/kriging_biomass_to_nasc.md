(kriged-biomass-to-nasc)=
# Apportion kriged biomass, abundance, and NASC

## Estimate kriged abundance and $\textit{NASC}$

Kriged biomass values can be converted to abundance and subsequently $\textit{NASC}$ using the `mesh_to_nasc` function from the `nwfsc_feat.apportionment` module. This is a function that can be configured in a variety of ways and has the following arguments:

- `mesh_data`: A `pandas.DataFrame` containing the gridded {ref}`kriged biomass estimates<kriging-imp>`.
- `biodata`: A dictionary of `xarray.DataArray`s that contain the [calculated weight proportions](../example_notebooks/weight_proportions.ipynb).
- `mesh_biodata_link`: A dictionary that maps column names from `mesh_data` to those in `biodata`. In this case, we can use it to match the haul-based strata (e.g., `stratum_ks`) with the latitude-based ones in the gridded data (e.g., `geostratum_ks`).
- `stratum_weights`: An `xarray.DataArray` that contains the {ref}`average weights per stratum <stratum-weights>` required for converting biomass to abundance. However, the overall means are used here as opposed to the sex-specific estimates.
- `stratum_sigma_bs`: A `pandas.DataFrame` containing the [mean $\sigma_\text{bs}$ for each stratum](sigma-bs-strata).

```python
from echopop.workflows.nwfsc_feat import apportionment

apportionment.mesh_biomass_to_nasc(
    mesh_data=unextrapolated_results,
    biodata=dict_ds_weight_proportion,
    group_columns=["sex", "stratum_ks"],
    mesh_biodata_link={"geostratum_ks": "stratum_ks"},
    stratum_weights=da_averaged_weight.sel(sex="all"),
    stratum_sigma_bs=invert_hake.sigma_bs_strata,
)
```

The `mesh_biomass_to_nasc` function modifies the DataFrame in-place to add the appropriate columns for sex-apportioned abundance and subsequently $\textit{NASC}$. 

## Distributing kriged abundance and biomass over age and length

Kriged abundance and biomass are both distributed over age ($\alpha$) and length ($\ell$) distributions using the `apportionment.distribute_population_estimates` function, which has the following arguments:

- `data`: `pandas.DataFrame` containing the population estimates.
- `proportions`: Either an `xarray.Dataset` or dictionary containing the number or weight proportions (depending on whether abundance or biomass is being distributed).
- `variable`: Name of the column in `data` containing the values to distribute (e.g., `"abundance"`, `"biomass"`).
- `group_columns`: List of column names that define any biological groups for the distribution (e.g., `["sex", "age_bin", "length_bin"]`).
- `data_proportions_link`: Dictionary that links colum names in `data` to those in `proportions`. For instance, the dictionary `{"stratum_A": "stratum_B"}` links `"stratum_A"` in `data` with `stratum_B` in `proportions`. 

```python
from echopop.workflows.nwfsc_feat import apportionment

# Abundance 
dict_ds_kriged_abundance_table = apportionment.distribute_population_estimates(
    data = df_kriged_results,
    proportions = dict_ds_number_proportion,
    variable = "abundance",
    group_columns = ["sex", "age_bin", "length_bin", "stratum_ks"],
    data_proportions_link={"geostratum_ks": "stratum_ks"},
)

# Biomass
dict_ds_kriged_biomass_table = apportionment.distribute_population_estimates(
    data=df_kriged_results,
    proportions = dict_da_weight_proportion,
    variable = "biomass",
    group_columns = ["sex", "age_bin", "length_bin", "stratum_ks"],
    data_proportions_link={"geostratum_ks": "stratum_ks"},
)
```

## Redistributing unaged estimates over age

The next step involves taking the unaged abundance and biomass distributions and redistributing them over age. This is done with the `apportionment.distribute_unaged_from_aged` function, which has the arguments:

- `population_table`: Population estimates to be redistributed.
- `reference_table`: Reference population table used for the redistribution.
- `collapse_dims`: List of dimension names to collapse (sum over) during standardization (e.g., `["stratum"]`). These dimensions are summed over in the initial step, and the resulting standardized table will not have these dimensions.
- `impute`: When `True`, perform a nearest-neighbor imputation for missing joint $\alpha$-$\ell$ values within each `group_columns` variable.
- `impute_variable`: List of variables used for imputation and required when `impute = True`. This typically refers to the dimension being imputed (e.g., `["age_bin"]`).

```python

# Abundance
dict_ds_kriged_abundance_table["redistributed_unaged"] = apportionment.distribute_unaged_from_aged(
    population_table = dict_ds_kriged_abundance_table["unaged"],
    reference_table = dict_ds_kriged_abundance_table["aged"],
    collapse_dims = ["stratum_ks"],
    impute = False,
)

# Biomass
dict_ds_kriged_biomass_table["redistributed_unaged"] = apportionment.distribute_unaged_from_aged(
    population_table = dict_ds_kriged_biomass_table["unaged"],
    reference_table = dict_ds_kriged_biomass_table["aged"],
    collapse_dims = ["stratum_ks"],
    impute=True,
    impute_variable=["age_bin"],
)
```

With these redistributed values, we can now appropriately sum the aged and unaged distributions since they share the same dimensions: `"sex"`, `"age_bin"`, and `"length_bin"`. This consolidation is done using the `apportionment.sum_population_tables` function, which has the following arguments:

- `population_tables`: The dictionary of `xarray.DataArray`s containing population estimate tables to be combined. Keys are the table names (e.g., `"aged"` and `"redistributed_unaged"`).

```python
# Abundance
da_kriged_abundance_table = apportionment.sum_population_tables(
    population_tables={
        "aged": dict_ds_kriged_abundance_table["aged"],
        "unaged": dict_ds_kriged_abundance_table["redistributed_unaged"]
    },
)

# Biomass
da_kriged_biomass_table = apportionment.sum_population_tables(
    population_tables={
        "aged": dict_ds_kriged_biomass_table["aged"],
        "unaged": dict_ds_kriged_biomass_table["redistributed_unaged"]
    },
)
```

## Re-allocating estimates for removed groups

In some cases, we may exclude certain $\ell$, $\alpha$, or other contrasts from earlier in the workflow. An expected example of this would be removing age-1 fish from the entire workflow. Despite these groups being removed from the workflow in general, there is some "leakage" that occurs in the distributions due to how the biological data are processed. This can be accounted by using the `reallocate_excluded_estimates` function, which has the arguments:

- `population_table`: The consolidated population table as an `xarray.DataArray` with biological groups corresponding coordinates (i.e., `"sex"`, `"age_bin"`, `"length_bin"`). 
- `exclusion_filter`: Dictionary specifying which table to exclude and redistribute. Keys are column and index names, values are the categories to exclude. For example, `exclusion_filter = {"age_bin": [1]}` would exclude age-1 fish.
- `group_columns`: List of column names that define the grouping variables for redistribution. For example, `group_columns = ["sex"]` would redistribute the age-1 estimates within each sex.

So if we excluded age-1 fish earlier in the analysis, then we can remove these via:

```python
# Abundance
ds_kriged_abundance_table_noage1 = apportionment.reallocate_excluded_estimates(
    population_table=da_kriged_abundance_table,
    exclusion_filter={"age_bin": [1]},
    group_columns=["sex"],
)

# Biomass
ds_kriged_biomass_table_noage1 = apportionment.reallocate_excluded_estimates(
    population_table=da_kriged_biomass_table,
    exclusion_filter={"age_bin": [1]},
    group_columns=["sex"],
)
```