(data-structure)=
# Core data structures

Input data, intermediate variables and calculations, and results are stored within several attributes associated with the `Survey`-class object. The attributes and stored data are detailed in the sections below. See [the example notebook](echopop-example-workflow) for information on how these attributes are used in the general `echopop` workflow. 

## Class attributes

There are currently **five** unique class attributes that store different types of information and data. 

* `Survey.meta`: Necessary information such as the date the object was created and general data workflow/provenance would be collected.
* `Survey.config`: Background configuration settings. 
* `Survey.input`: Imported acoustic, biological, kriging, and stratification [data detailed elsewhere](input-files).
* `Survey.analysis`: Intermediate data products, variables, and other calculations that may be both relevant to users (e.g. number proportions) or are required for computing the final results. 
* `Survey.results`: Overall results for each analysis.

Each of these attributes are organized using nested dictionaries that may contain an uneven number of levels. For instance, calculated number and weight proportions computed via the transect analysis are nested via the following structure:
* `Survey`-class object
  * `analysis`: parent node containing intermediate calculations for all analyses in a nested dictionary structure
    * `transect`: intermediate calculations computed via the `Survey.transect_analysis(...)` method
      * `biology`: biological calculations 
        * `proportions`: proportion calculations
          * `number`: number proportions
            * Various number proportion `DataFrame`'s
          * `weight`: weight proportions
            * Various weight proportion `DataFrame`'s
which can be accessed using the following Python code:
```python
# Accessing number proportion dataframes
Survey.analysis['transect']['biology']['proportions']['number']
# Accessing weight proportion dataframes
Survey.analysis['transect']['biology']['proportions']['weight']
```

### Variable naming convention

In the `Survey` class data structure, the following suffixes are used to denote the data type of each entry:
- `_df`: A standard `DataFrame` indexed by row number with at least two columns.
- `_tbl`: A pivot table, which is a `DataFrame` whose rows and columns are both indexed by different variables besides the row and column numbers/names.
- `_config`: A dictionary containing specific parameter arguments for a particular analysis or model.
- ` `: No suffix indicates either a dictionary with additional nested levels, or a scalar number/string (e.g. a single survey-wide mean estimate)

These suffix conventions can help aid in navigating the nested dictionaries via `.keys()` within each `Survey`-class attribute to see what is included at each level. For instance, users can run the following command:
```python 
# If the Survey-class object is named 'survey':
survey.input['biology'].keys()
```
which yields the following output:
```python
dict_keys(['catch_df', 'distributions', 'length_df', 'haul_to_transect_df', 'specimen_df'])
```

## `Survey`-class data structure

This data structure is meant to provide an intentional approach toward organizing the data to limit any ambiguities concerning the identity, provenance, or interpretation of the stored data and results contained within the `Survey`-class object. When the `Survey`-class object is ran using all available analyses, the overall data structure of each attribute (and their stored variables, `DataFrames`'s, etc.) can be mapped out as indicated below.

```{contents}
:local:
:depth: 4
```

### **Metadata** (`.meta`)
This is a currently underdeveloped attribute that will be adapted for more thourough use in later `echopop` versions. It comprises two keys: 
- `date`: Datetime the `Survey`-class was generated.
- `provenance`: An iteratively populated dictionary that tracks what analyses/processing steps have been performed on the object.

### **Configuration** (`.config`)
Configuration settings defined in the two `.yml` files are imported, reformatted, and stored in this attribute. This attribute comprises 11 keys: 
- `biological`: Biological data file and sheet names.
- `CAN_haul_offset`: A constant that is applied to Canadian haul numbers to differentiate them from US trawls.
- `data_root_dir`: The root directory containing all data files.
- `geospatial`: Geospatial variables (i.e. EPSG string).
- `kriging`: Kriging data file and sheet names.
- `NASC`: Georeferenced acoustic trawl data file and sheet names.
- `species`: Species that will be analyzed that includes both the text and equivalent number codes.
- `stratification`: Stratification data file and sheet names.
- `stratified_survey_mean_parameters`: Algorithm parameters used for the stratified analysis that incorporates the Jolly and Hampton (1990) algorithm.
- `survey_year`: The year that the survey was conducted.
- `TS_length_regression_parameters`: Regression coefficient terms for species-specific TS-length relationships.

### **Input data** (`.input`)
Acoustic, biological, spatial, and statistical data defined in the configuration `.yml` file undergo initial processing after being read in and are stored within a single attribute. This attribute comprises four keys: `'acoustics'`, `'biology'`, `'spatial'`, and `'statistics'`. Each of these keys contain other sub-keys that are sometimes nested. 

#### Acoustic data (`['acoustics']`)
- `nasc_df`: Georeferenced acoustic trawl data for all-aged and age-2+ fish.

#### Biological data (`['biology']`)
- `catch_df`: Total unaged weights for each haul.
- `distributions`: Age and length distributions.
  - `age_bins_df`: Age distribution/histogram bins.
  - `length_bins_df`: Length distribution/histogram bins.
- `haul_to_transect_df`: Haul-to-transect key that links haul numbers to their respective transects.
- `length_df`: Unaged (Station 1) length measurements.
- `specimen_df`: Aged (Station 2) length and weight measurements.

#### Spatial data (`['spatial']`)
- `strata_df`: The length-based 'KS' stratum definitions and fraction of hake for each haul.
- `geo_strata_df`: The latitudinal extents/ranges of each KS stratum throughout the survey region.
- `inpfc_strata_df`: The INPFC stratum definitions and their respective latitudinal limits

#### Statistical data (`['statistics']`)
- `kriging`: Data and parameters required for the kriging analysis.
  - `mesh_df`: The survey kriging mesh used for interpolation.
  - `isobath_200m_df`: The 200 m isobath coordinates.
  - `model_config`: A dictionary comprising all required parameter values for the kriging analysis.
- `variogram`: Data and parameters required for the variogram analysis.
  - `model_config`: Dictionary comprising all required arguments for the variogram analysis.

### **Analysis** (`.analysis`)
This attribute is to store all intermediate data products, variables, and other calculations necessary to compute the final resutls of each analysis. It comprises four keys: `transect`, `settings`, `stratified`, and `kriging`. Each of these keys contain other sub-keys that are sometimes nested. 

#### Analysis settings (`['settings']`)
There are currently four sub-keys that population the `settings` key: `kriging`, `stratified`, `transect`, and `variogram`. All of these sub-keys include user-defined arguments (or the default values) used for each analysis. This provides a recipe approach that helps orchestrate how functions contained within each analysis operate and further tracks the details of each analysis for replicability. 

#### Transect variables (`['transect']`)
- `acoustics`: Intermediate acoustic calculations.
    - `sigma_bs`: Aggregate $\sigma_{bs}$ calculations.
        - `haul_mean_df`: The mean $\sigma_{bs}$ for each haul.
        - `strata_mean_df`: The mean $\sigma_{bs}$ for each stratum (KS or INPFC).
    - `adult_transect_df`: Georeferenced population estimates (e.g. abundance, biomass) distributed across the original survey transect longitude and latitude coordinates.
- `biology`: Intermediate biological calculations.
  - `distributions`:
    - `binned_aged_counts_df`: The number of aged fish distributed across length- and age-bins.
    - `binned_aged_counts_filtered_df`: The number of aged fish distributed across length- and age-bins when explicitly excluding unsexed fish or those missing weight estimates.
    - `binned_uanged_counts_df`: The number of unaged fish distributed across length-bins.
    - `weight`: The summed fish weights distributed across length- and age-bins (when applicable).
      - `aged_length_weight_tbl`: Summed aged fish weights distributed across age- and length-bins.
      - `unaged_length_weight_tbl`: Summed unaged fish weights distributed across length-bins.
  - `population`: Population estimates.
    - `tables`: Total population estimates distributed across age- and length-bins.
      - `abundance`: Abundance estimates.
        - `aged_abundance_df`: The summed aged fish abundance distributed across age- and length-bins.
        - `unaged_abundance_df`: The summed unaged fish abundance distributed across length-bins.
        - `unaged_age1_abundance_df`: The summed unaged fish abundance accounting for age-1 proportions when age-1 fish are excluded from parts of the analysis.
      - `biomass`: Biomass estimates.
        - `aged_biomass_df`: The total biomass estimated from the transect data distributed across age- and length-bins.
  - `proportions`: Number and weight proportions.
    - `number`: Number (count) proportions of age- and length-bins.
      - `age_proportions_df`:
      - `sex_proportions_df`: 
      - `aged_length_proportions_df`: 
      - `unaged_length_proportions_df`:
    - `weight`: Weight proportions of age- and length-bins.
      - `aged_weight_proportions_df`: Weight proportions of aged fish across age- and length-bins.
      - `unaged_weight_proportions_df`: Weight proportions of unaged fish across length-bins.
      - `aged_unaged_sex_weight_proportions_df`: The relative weight proportions of males and females belonging to both aged and unaged fish groups. 
      - `aged_unaged_weight_proportions_df`: Weight proportions of unaged fish across length-bins distributed across age-bins based on the relative age distributions estimated from the specifically aged fish.
  - `weight`: Weight biometrics.
    - `length_weight_regression`: Length-weight regression parameters and fitted values.
      - `parameters_df`: Regression coefficients.
      - `weight_fitted_df`: Fitted weights for each length-bin.
    - `weight_stratum_df`: Mean weight for male, female, and all fish specific to each KS or INPFC stratum.
- `coordinates`: The stored coordinates including the original latitude, longitude, and stratum (KS and INPFC) for each interval along every transect so any spatial transformations can be readily converted back to their original forms.

#### Stratified variables (`['stratified']`)
This key can comprise two sub-keys depending on whether `kriging_analysis(...)` and/or `transect_analysis(...)` are run: `kriging` and `transect`, respectivvely. Otherwise, the contents of both of these sub-keys contain an identically named `DataFrame` called `stratified_replicates_df` that represents the various statistics computed from the Jolly and Hampton (1990) algorithm over all bootstrapped iterations.

#### Kriging variables (`['kriging']`)
- `mesh_df`: A copy of the kriging mesh located in `.analysis['statistics']['kriging']['mesh_df']` that can be updated (e.g. cropping, standardizing) depending on user arguments in `kriging_analysis(...)`. 
- `transect_df`: A reduced copy of the transect-specific population results (`.analysis['acoustics']['adult_transect_df']`) with the target variable that will be kriged that can be updated (e.g. coordinate transformation).

### **Final results** (`.results`)
The results of each analysis are stored in this attribute. This currently can include three keys: `transect`, `stratified`, and `kriging`.

#### Transect results (`['transect']`)
- `biomass_summary_df`: A DataFrame that estimates the total biomass of age-1, age-2+, and all fish partitioned by sex (`male`, `female`, `all`, `unsexed`) and haul species composition (`mixed`).

#### Stratified results (`['stratified']`)
Similar to `.analysis['stratified']`, two sub-keys can be located here: `transect` and `kriging`. These represent nested sub-keys with stratified statistics computed using the transect- and kriged-specific data, respectively. Each of these sub-keys include the following entries:

- `variable`: The population estimate used for the stratified analysis (e.g. biomass).
- `ci_percentile`: The confidence interval percentile.
- `num_transects`: The number of (virtual) transects represented in the analysis.
- `total_area`: The total area coverage of the (virtual) transects.
- `mean`: Various statistics computed for the mean `variable`.
- `variance`: Various statistics computed for the variance of `variable`.
- `cv`: Various statistics computed for the coefficient of variation of `variable`.
- `total`: Various statistics computed for the total `variable`.

#### Kriging results (`['kriging']`)
- `variable`: The population estimate used for the kriging analysis (e.g. biomass density).
- `survey_mean`: The overall mean `variable` estimate per unit area (e.g. biomass density).
- `survey_estimate`: The total `variable` estimate integrated over area (e.g. biomass density converted to biomass).
- `survey_cv`: The mean kriged `variable` coefficient of variation.
- `mesh_results_df`: A DataFrame that includes the georeferenced kriged mean, variance, and other metrics for the entire mesh grid.
- `tables`: Apportioned kriging results integrated over area (i.e. abundance, biomass).
  - `aged_tbl`: Total kriged values apportioend over age- and length-bins for each sex among aged fish.
  - `unaged_tbl`: Total kriged values apportioned over length-bins for each sex among unaged fish.
  - `overall_apportionment_df`: Total kriged values apportioned over age- and length-bins inclusive of both aged and unaged fish.