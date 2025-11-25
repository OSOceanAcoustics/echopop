(kriged-biomass-to-nasc)=
# Apportioning kriged biomass, abundance, and NASC

## Estimating kriged abundance and $\textit{NASC}$

Kriged biomass values can be converted to abundance and subsequently $\textit{NASC}$ using the `mesh_to_nasc` function from the `nwfsc_feat.apportionment` module. This is a function that can be configured in a variety of ways and has the following arguments:

- `mesh_data_df`: A `pandas.DataFrame` containing the gridded {ref}`kriged biomass estimates<kriging-imp>`.
- `biodata`: A `pandas.DataFrame` or dictionary that contains the [calculated weight proportions](../example_notebooks/weight_proportions.ipynb).
- `mesh_biodata_link`: A dictionary that maps column names from `mesh_data_df` to those in `biodata`. In this case, we can use it to match the haul-based strata (e.g., `stratum_ks`) with the latitude-based ones in the gridded data (e.g., `geostratum_ks`).
- `stratum_weights_df`: A `pandas.DataFrame` that contains the {ref}`average weights per stratum <stratum-weights>` required for converting biomass to abundance. However, the overall means are used here as opposed to the sex-specific estimates.
- `stratum_sigma_bs_df`: A `pandas.DataFrame` containing the [mean $\sigma_\text{bs}$ for each stratum](sigma-bs-strata).

```python
from echopop.workflows.nwfsc_feat import apportionment

apportionment.mesh_biomass_to_nasc(
    mesh_data_df=df_kriged_results,
    biodata=dict_df_weight_proportion,
    group_by=["sex"],
    mesh_biodata_link={"geostratum_ks": "stratum_ks"},
    stratum_weights_df=df_averaged_weight["all"],
    stratum_sigma_bs_df=invert_hake.sigma_bs_strata,
)
```
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>geostratum_ks</th>
      <th>et_id</th>
      <th>area (km^2)</th>
      <th>latitude</th>
      <th>longitude</th>
      <th>fraction</th>
      <th>cell cut by polygon?</th>
      <th>geostratum_inpfc</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>kriged_variance</th>
      <th>sample_variance</th>
      <th>cell_cv</th>
      <th>biomass</th>
      <th>biomass_female</th>
      <th>biomass_male</th>
      <th>abundance_female</th>
      <th>abundance_male</th>
      <th>abundance</th>
      <th>nasc</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>7</td>
      <td>85670</td>
      <td>21.4369</td>
      <td>49.099727</td>
      <td>-126.024144</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.125747</td>
      <td>0.200900</td>
      <td>...</td>
      <td>0.500427</td>
      <td>NaN</td>
      <td>0.034444</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>7</td>
      <td>85427</td>
      <td>21.4369</td>
      <td>49.057959</td>
      <td>-126.024127</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.119084</td>
      <td>0.198853</td>
      <td>...</td>
      <td>0.131102</td>
      <td>NaN</td>
      <td>0.017630</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>7</td>
      <td>85184</td>
      <td>21.4369</td>
      <td>49.016196</td>
      <td>-126.024110</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.109859</td>
      <td>0.196807</td>
      <td>...</td>
      <td>0.414271</td>
      <td>NaN</td>
      <td>0.031339</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>7</td>
      <td>84941</td>
      <td>21.4369</td>
      <td>48.974438</td>
      <td>-126.024093</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.099092</td>
      <td>0.194760</td>
      <td>...</td>
      <td>0.646779</td>
      <td>NaN</td>
      <td>0.039158</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>7</td>
      <td>84698</td>
      <td>21.4369</td>
      <td>48.932686</td>
      <td>-126.024076</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.088307</td>
      <td>0.192714</td>
      <td>...</td>
      <td>0.732735</td>
      <td>1.373507</td>
      <td>0.041679</td>
      <td>284647.842003</td>
      <td>198307.321715</td>
      <td>86340.520288</td>
      <td>258852.039798</td>
      <td>112700.931062</td>
      <td>371552.970860</td>
      <td>1791.545291</td>
    </tr>
    <tr>
      <th>5</th>
      <td>8</td>
      <td>84455</td>
      <td>21.4369</td>
      <td>48.890939</td>
      <td>-126.024060</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.077506</td>
      <td>0.190669</td>
      <td>...</td>
      <td>0.757142</td>
      <td>NaN</td>
      <td>0.042368</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>6</th>
      <td>8</td>
      <td>84212</td>
      <td>21.4369</td>
      <td>48.849198</td>
      <td>-126.024043</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.066959</td>
      <td>0.188623</td>
      <td>...</td>
      <td>0.718963</td>
      <td>1.812516</td>
      <td>0.041286</td>
      <td>337598.042729</td>
      <td>296729.780259</td>
      <td>40868.262470</td>
      <td>256175.630352</td>
      <td>35282.784527</td>
      <td>291458.414879</td>
      <td>1760.043038</td>
    </tr>
    <tr>
      <th>7</th>
      <td>8</td>
      <td>83969</td>
      <td>21.4369</td>
      <td>48.807461</td>
      <td>-126.024026</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.056558</td>
      <td>0.186578</td>
      <td>...</td>
      <td>0.809815</td>
      <td>NaN</td>
      <td>0.043817</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>8</th>
      <td>8</td>
      <td>83726</td>
      <td>21.4369</td>
      <td>48.765730</td>
      <td>-126.024009</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.046140</td>
      <td>0.184533</td>
      <td>...</td>
      <td>0.325678</td>
      <td>6.237652</td>
      <td>0.027787</td>
      <td>1618.000216</td>
      <td>1422.131611</td>
      <td>195.868604</td>
      <td>1227.768449</td>
      <td>169.099182</td>
      <td>1396.867631</td>
      <td>8.435327</td>
    </tr>
    <tr>
      <th>9</th>
      <td>8</td>
      <td>83483</td>
      <td>21.4369</td>
      <td>48.724004</td>
      <td>-126.023992</td>
      <td>1.0</td>
      <td>N</td>
      <td>6</td>
      <td>0.035707</td>
      <td>0.182488</td>
      <td>...</td>
      <td>0.085441</td>
      <td>0.242226</td>
      <td>0.014232</td>
      <td>46618.385601</td>
      <td>40974.951169</td>
      <td>5643.434433</td>
      <td>35374.892049</td>
      <td>4872.144521</td>
      <td>40247.036570</td>
      <td>243.041590</td>
    </tr>
  </tbody>
</table>
<p>10 rows × 22 columns</p>
</div>

The `mesh_biomass_to_nasc` function modifies the DataFrame in-place to add the appropriate columns for sex-apportioned abundance and subsequently $\textit{NASC}$. 

## Distributing kriged abundance and biomass over age and length

Kriged abundance and biomass are both distributed over age ($\alpha$) and length ($\ell$) distributions using the `apportionment.distribute_population_estimates` function, which has the following arguments:

- `data`: `pandas.DataFrame` containing the population estimates.
- `proportions`: Either a `pandas.DataFrame` or dictionary containing the number or weight proportions (depending on whether abundance or biomass is being distributed).
- `variable`: Name of the column in `data` containing the values to distribute (e.g., `"abundance"`, `"biomass"`).
- `group_by`: List of column names that define any biological groups for the distribution (e.g., `["sex", "age_bin", "length_bin"]`).
- `stratify_by`: List of columns used for stratification (e.g., `["stratum_name"]`).
- `data_proportions_link`: Dictionary that links colum names in `data` to those in `proportions`. For instance, the dictionary `{"stratum_A": "stratum_B"}` links `"stratum_A"` in `data` with `stratum_B` in `proportions`. 

```python
from echopop.workflows.nwfsc_feat import apportionment

# Abundance 
dict_kriged_abundance_table = apportionment.distribute_population_estimates(
    data=df_kriged_results,
    proportions=dict_df_number_proportion,
    variable="abundance",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
    data_proportions_link={"geostratum_ks": "stratum_ks"},
)

# Biomass
dict_kriged_biomass_table = apportionment.distribute_population_estimates(
    data=df_kriged_results,
    proportions=dict_df_weight_proportion,
    variable="biomass",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
    data_proportions_link={"geostratum_ks": "stratum_ks"},
)
```

## Redistributing unaged estimates over age

The next step involves taking the unaged abundance and biomass distributions and redistributing them over age. This is done with the `apportionment.distribute_unaged_from_aged` function, which has the arguments:

- `population_table`: Population estimates to be redistributed.
- `reference_table`: Reference population table used for the redistribution.
- `group_by`: List of columns that define the grouping variables for this rediistribution (e.g., `["sex"]`)
- `impute`: When `True`, perform a nearest-neighbor imputation for missing joint $\alpha$-$\ell$ values within each `group_by` variable.
- `impute_variable`: List of variables used for imputation and required when `impute = True`. This typically refers to the dimension being imputed (e.g., `["age_bin"]`).

```python

# Abundance
dict_kriged_abundance_table["redistributed_unaged"] = apportionment.distribute_unaged_from_aged(
    population_table=dict_kriged_abundance_table["unaged"],
    reference_table=dict_kriged_abundance_table["aged"],
    group_by=["sex"],
    impute=False,
)

# Biomass
dict_kriged_biomass_table["redistributed_unaged"] = apportionment.distribute_unaged_from_aged(
    population_table=dict_kriged_biomass_table["unaged"],
    reference_table=dict_kriged_biomass_table["aged"],
    group_by=["sex"],
    impute=True,
    impute_variable=["age_bin"],
)
```

With these redistributed values, we can now appropriately sum the aged and unaged distributions since they share the same dimensions: `"sex"`, `"age_bin"`, and `"length_bin"`. This consolidation is done using the `apportionment.sum_population_tables` function, which has the following arguments:

- `population_table`: The dictionary of population estimate tables to combined. Keys are the table names (e.g., `"aged"` and `"redistributed_unaged"`) and values are the `pandas.DataFrame`s with the associated population estimates (i.e., abundance and biomass).
- `table_names`: List of table names (keys) from `population_table` to include in the summation. Tables not in this list will be ignored. 
- `table_index`: List of column names to use as the index in the final combined table (e.g., `["length_bin"]`). This means that the rows of the output table will correspond to length $\ell$.
- `table_columns`: List of column names to use as the column indices in the final combined table (e.g., `["age_bin", "sex"]`). This means that the columns of the output table will correspond to age $\alpha$ and sex $s$.

```python
# Abundance
df_kriged_abundance_table = apportionment.sum_population_tables(
    population_table=dict_kriged_abundance_table,
    table_names=["aged", "redistributed_unaged"],
    table_index=["length_bin"],
    table_columns=["age_bin", "sex"],
)

# Biomass
df_kriged_biomass_table = apportionment.sum_population_tables(
    population_table=dict_kriged_biomass_table,
    table_names=["aged", "redistributed_unaged"],
    table_index=["length_bin"],
    table_columns=["age_bin", "sex"],
)
```

So the final table output for `df_kriged_abundance_table` would look like:

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead tr th {
        text-align: left;
    }

    .dataframe thead tr:last-of-type th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th>age_bin</th>
      <th colspan="2" halign="left">(0.5, 1.5]</th>
      <th colspan="2" halign="left">(1.5, 2.5]</th>
      <th colspan="2" halign="left">(2.5, 3.5]</th>
      <th colspan="2" halign="left">(3.5, 4.5]</th>
      <th colspan="2" halign="left">(4.5, 5.5]</th>
      <th>...</th>
      <th colspan="2" halign="left">(17.5, 18.5]</th>
      <th colspan="2" halign="left">(18.5, 19.5]</th>
      <th colspan="2" halign="left">(19.5, 20.5]</th>
      <th colspan="2" halign="left">(20.5, 21.5]</th>
      <th colspan="2" halign="left">(21.5, 22.5]</th>
    </tr>
    <tr>
      <th>sex</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
      <th>...</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
      <th>female</th>
      <th>male</th>
    </tr>
    <tr>
      <th>length_bin</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>(1.0, 3.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(3.0, 5.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(5.0, 7.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(7.0, 9.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(9.0, 11.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(11.0, 13.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(13.0, 15.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(15.0, 17.0]</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(17.0, 19.0]</th>
      <td>0.000000e+00</td>
      <td>3.201968e+06</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>(19.0, 21.0]</th>
      <td>7.294352e+07</td>
      <td>8.713349e+07</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>10 rows × 44 columns</p>
</div>

## Re-allocating estimates for removed groups

In some cases, we may exclude certain $\ell$, $\alpha$, or other contrasts from earlier in the workflow. An expected example of this would be removing age-1 fish from the entire workflow. Despite these groups being removed from the workflow in general, there is some "leakage" that occurs in the distributions due to how the biological data are processed. This can be accounted by using the `reallocate_excluded_estimates` function, which has the arguments:

- `population_table`: The consolidated population table with biological groups as the index and columns (i.e., `"sex"`, `"age_bin"`, `"length_bin"`). 
- `exclusion_filter`: Dictionary specifying which table to exclude and redistribute. Keys are column and index names, values are the categories to exclude. For example, `exclusion_filter = {"age_bin": [1]}` would exclude age-1 fish.
- `group_by`: List of column names that define the grouping variables for redistribution. For example, `group_by = ["sex"]` would redistribute the age-1 estimates within each sex.

So if we excluded age-1 fish earlier in the analysis, then we can remove these via:

```python
# Abundance
df_kriged_abundance_table_noage1 = apportionment.reallocate_excluded_estimates(
    population_table=df_kriged_abundance_table,
    exclusion_filter={"age_bin": [1]},
    group_by=["sex"],
)

# Biomass
df_kriged_biomass_table_noage1 = apportionment.reallocate_excluded_estimates(
    population_table=df_kriged_biomass_table,
    exclusion_filter={"age_bin": [1]},
    group_by=["sex"],
)
```

So we can contrast `df_kriged_abundance_table_noage1` with `df_kriged_abundance_table` by seeing how the age-1 estimates were re-distributed across $\alpha > 1$. So first we can check the difference in total biomass for each $\alpha$ (only the first 10 rows are displayed below):

```python
df_kriged_abundance_table_noage1.sum() - df_kriged_abundance_table.sum()

age_bin     sex   
(0.5, 1.5]  female   -3.210495e+08
            male     -3.774863e+08
(1.5, 2.5]  female    4.976895e+07
            male      5.085311e+07
(2.5, 3.5]  female    9.509865e+07
            male      1.111761e+08
(3.5, 4.5]  female    3.913087e+06
            male      6.893927e+06
(4.5, 5.5]  female    8.905839e+07
            male      1.112596e+08
dtype: float64
```

But if we sum the entire tables, we can see that the overall summed biomass did not significantly change.

```python
df_kriged_abundance_table_noage1.sum().sum() - df_kriged_abundance_table.sum().sum()

np.float64(4.76837158203125e-07)
```


