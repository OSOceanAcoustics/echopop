# Classes and data structures

## Classes

Python EchoPro is designed using an object-oriented approach where data are stored in core classes and computations are called as methods (functions) on those classes.

- Core class
  - survey.Survey
  - All actions and computations are orchestrated from this class
- computation
  - computation.transect_results.ComputeTransectVariables
  - computation.semivariogram.SemiVariogram
  - computation.kriging.Kriging
  - computation.kriging_variables.ComputeKrigingVariables
  - computation.bootstrapping.Bootstrapping
- data loaders
  - data_loader.kriging_mesh.KrigingMesh
  - data_loader.biological_data.LoadBioData
  - data_loader.stratification_data.LoadStrataData
- report generation
  - reports.Reports

### Other notes

- `bio_calc`: From `computation.transect_results.ComputeTransectVariables`; and also from Kriging?
  - Categorize the long list of df's and gdf's available on `bio_calc`; eg, `kriging`, `kriging_bin`, `transect`, `weight_fraction`, etc. These can be categorized either by prefix or suffix
  - Only the 6 "results" tables are gdf's: transect vs kriging, and sex (all[empty]/male/female). eg, `transect_results_male_gdf`
  - `bin_ds` is an xarray dataset
  - Most of the other objects are dataframes (`_df`), but there are a few others that don't have a clarifying suffix: `strata_sig_b`, `all_strata`, `missing_strata`, `percentage_transects_selected`, `sel_tran_strata_choice` (a dict), `stratum_choices` (a dict)


## Data structures

- df's acessible from the `survey`
    - innputs
    - others
- `bin_ds` 
- results gdf's