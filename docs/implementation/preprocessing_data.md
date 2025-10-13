(preprocessing-data)=
# Loading and preprocessing data

Acoustic, biological, and other data loaded into `Echopop` are preprocessed in various ways to enable full functionality with the software. This includes changes to dataset names, column names within file-specific dataframes, and binning data (e.g. age, latitude, length). 

## Configuring files
### `survey_year_{YEAR}_config.yml`
The `survey_year_{YEAR}_config.yml` configuration file defines general metadata and the filepaths for the datasets that will be loaded into `Echopop`.

#### Parameters
* `survey_year`: metadata entry detailing the survey year.
* `species`: metadata entries detailing target species.
    * `text_code`: target species -- string (e.g. "blue_whale").
    * `number_code`: target species -- numeric (e.g. 123456789)
* `CAN_haul_offset`: an offset value applied to Canadian haul numbers [**specific to the FEAT hake survey**].
* `data_root_dir`: root data directory path where all data are stored.
* `report_path`: an optional directory path where report tables for all results are stored [**specific to the FEAT hake survey**].

#### Input data files
* `biological`: all biodata.
  * `catch`/`length`/`specimen`: aggregated haul data, binned length data, and specimen-specific length and weight data.
    * `{REGION}`: a region code (e.g. `US`, `CAN`) that separates similar biodata collected from different geographical regions or sources (e.g. ships).
      * `filename`: name of the `*.xlsx` file.
      * `sheetname`: sheet containing required data. *Note*: `sheetname` is typically written as a string (e.g. `Sheet1`), but it can also be expressed as a list (e.g. `[Sheet1, Sheet2]`) when data from multiple sheets within the same file are required.
* `stratification`: stratification defininitions.
  * `strata`: length-based stratification definitions.
      * `filename`: *same as above*.
      * `sheetname`: *same as above*.
  * `geo_strata`: latitude-based (i.e. INPFC) stratificaiton definitions.
      * `filename`: *same as above*.
      * `sheetname`: *same as above*.
* `NASC`: vertically integrated acoustic backscatter data.
  * `{GROUP}`: a group-specific code for differentiating backscatter inputs (e.g. `age-1`, `all_ages`).
      * `filename`: *same as above*.
      * `sheetname`: *same as above*.
* `kriging`: kriging and variogram parameters.
  * `mesh`: mesh grid for which kriged population estimates are projected.
      * `filename`: *same as above*.
      * `sheetname`: *same as above*.
  * `isobath_200m`: latitude and longitude coordinates of the 200 m isobath [**specific to the FEAT hake survey**].
      * `filename`: *same as above*.
      * `sheetname`: *same as above*.
  * `vario_krig_para`: kriging and variogram model parameters.
      * `filename`: *same as above*.
      * `sheetname`: *same as above*.
  
### `initialization_config.yml`
The `initialization_config.yml` configuration file defines various parameters that dictate how biological data are binned (i.e. age and length), how Echoview exports are processed and consolidated into single files, TS-length regression parameters, and additional geospatial processing inputs.

#### Parameters
* `stratified_survey_mean_parameters`: parameters required for applying the Jolly and Hampton (1990) analysis.
  * `strata_transect_proportion`: percent of transects resampled within each stratum.
  * `num_replicates`: number of bootstrap/resampling replicates are conducted.
  * `mesh_transects_per_latitude`: the number of virtual/synthetic transects that are generated for running this analysis on kriged mesh values.
* `bio_hake_len_bin`: a list of values representing the binned length distribution that has the following format: `[minimum length, maximum length, number of bins]` [**specific to the FEAT hake survey**].
* `bio_hake_age_bin`: a list of values representing the binned age distribution that has the following format: `[minimum age, maximum age, number of bins]` [**specific to the FEAT hake survey**].
* `nasc_exports`:
  * `export_file_directory`: directory path containing Echoview exports. *Note*: this assumes that the directory is contained within `data_root_dir`.
  * `nasc_export_directory`: directory path where processed exports will be saved. *Note*: this assumes that the directory is contained within `data_root_dir`.
  * `save_file_template`: the filename format of the processed export files that can contain the following tags: `{REGION}`, `{YEAR}`, and/or `{GROUP}`.
  * `save_file_sheetname`: sheetname of the processed export file.
  * `regions`: acoustic data region names.
    * `{GROUP}`: a single string (e.g. `age-1 fish`) or a list of strings that represent how backscatter measurements are integrated for each `{GROUP}` region.
  * `max_transect_spacing`: the maximum spacing between transects (nmi).
  * `file_columns`: a list file columns that will be included in the processed export file.
  * `transect_region_mapping`: definitions for parsing acoustic region names to parse critical metadata.
    * `pattern`: a regular expression pattern that can contain the following tags: `{REGION_CLASS}`, `{HAUL_NUM}`, and `{COUNTRY}`.
    * `parts`: each tag in `pattern` can be broken down by various regular expressions to extract different information, such as different `{REGION_CLASS}` (e.g. `Hake`, `Age-1 Hake`, `Non-hake`). Each component within `parts` must contain a `pattern` and `label` entry. The tags contained within `parts` **must** match those defined in the overall `pattern` (e.g. `"{REGION_CLASS}{HAUL_NUM}{COUNTRY}"`). For instance, the three entries for this could be formatted via:
       ```yaml
        REGION_CLASS:
        - pattern: ^[hH](?![a-zA-Z]|1a)
          label: Hake
        - pattern: ^[hH]1[aA][mM]
          label: Age-1 Hake Mix
        HAUL_NUM:
        - pattern: '[0-9]+'
          label: None
        COUNTRY:
        - pattern: ^[cC]
          label: CAN
        - pattern: ^[uU]
          label: US
        ```
  * `TS_length_regression_parameters`: TS-length regression parameters.
    * `{SPECIES}`: this should match the `species` parameter defined in the `survey_year_{YEAR}_config.yml` configuration file.
      * `number_code`: target species -- numeric (e.g. 123456789)
      * `TS_L_slope`: regression slope.
      * `TS_L_intercept`: regression *y*-intercept.
      * `length_units`: units for length used for computing the regression.
  * `geospatial`: geospatial settings
    * `init`: EPSG integer code for geodetic paramaterization (e.g. `epsg:4326`).
  * `kriging_parameters`: additional kriging parameters.
    * `A0`: base area of each mesh grid cell.
    * `longitude_reference`: reference longitude for kriging mesh coordinate standardization.
    * `longitude_offset`: longitude offset applied for kriging mesh coordinate standardization.
    * `latitude_offset`: latitude offset applied for kriging mesh coordinate standardization.

## Loading data into `Echopop`
```{figure} /images/Echopop_input_schematic.png
:width: 600px
:name: schematic

A color-coded schematic that provides a visual overview of how data are loaded and preprocessed. Click on schematic to zoom in.
```

**<span style="color:#00A200">Data are loaded</span>** into `Echopop` via `Survey.load_survey_data()` and `Survey.load_acoustic_data()` for the `[Biological, Kriging, Stratification]` and `[NASC]` datasets, respectively. The names of these datasets are mutated and stored slightly differently within the `Survey.input` attribute:

- `Biological` ➡️ `Survey.input["biology"]`
  - `catch` ➡️ `Survey.input["biology"]["catch_df"]`
  - `length` ➡️ `Survey.input["biology"]["length_df"]`
  - `specimen` ➡️ `Survey.input["biology"]["specimen_df"]`
- `Configuration` ➡️ `Survey.input["biology"]`
  - `bio_hake_len_bin`/`bio_hake_age_bin` ➡️ `Survey.input["biology"]["distributions]`
- `Kriging` ➡️ `Survey.input["statistics"]` 
  - `mesh` ➡️ `Survey.input["statistics"]["kriging"]["mesh_df"]`
  - `isobath_200m` ➡️ `Survey.input["statistics"]["kriging"]["isobath_200m_df"]`
  - `vario_krig_para` ➡️ `Survey.input["statistics"]["kriging"]["vario_krig_para"]`
- `NASC` ➡️ `Survey.input["acoustics"]`
  - `{GROUP}` (all) ➡️ `Survey.input["acoustics"]["nasc_df"]`
- `Stratification` ➡️ `Survey.input["spatial"]`
  - `strata` ➡️: `Survey.input["spatial"]["strata_df"]`
  - `strata` ➡️: `Survey.input["spatial"]["inpfc_strata_df"]`
  - `geo_strata` ➡️ `Survey.input["spatial"]["geo_strata_df"]`
  - `geo_strata` ➡️ `Survey.input["spatial"]["inpfc_geo_strata_df"]`

Echoview exports can be **<span style="color:#6666FF">alternatively processed and loaded</span>** into `Echopop` by incorporating the `nasc_exports` parameters within `initialization_config.yml`. These files can also processed outside of the same `Echopop` workflow whereby the processed exports can then be saved and used to parameterize the `NASC` dataset definiations within the `survey_year_{YEAR}_config.yml` configuration file.

Once the kriging and variogram parameters are loaded, the specific model parameterizations are **<span style="color:#BBB982">reformatted and split</span>** into separate dictionaries:

`Survey.input["statistics"]["kriging"]["vario_krig_para"]` ➡️ `["kriging"]["model_config"]`
<span style="display:inline-block; width: 51ch;"></span> ↘️ `["variogram"]["model_config"]`

Once split, `Survey.input["statistics"]["kriging"]["vario_krig_para"]` is removed from `Survey.input["statistics"]["kriging"]`.

## Preprocessing data
Several **<span style="color:#FF8000">preprocessing steps</span>** are conducted to ensure that the data comprise all of the necessary information required for computing the transect results via `Survey.transect_results()`.

1. When `Survey.input["spatial"]` is fully validated and loaded, the associated dataframes are used to **<span style="color:#FF9494">stratify data</span>** within both `Survey.input["acoustics"]` and `Survey.input["biology"]`. Since the stratification data files are loaded separately from the acoustic data, acoustic transect data are **only** stratified once both datasets are loaded. Acoustic and biological data that contain the column `haul_num` are stratified using both the length- and INPFC-based strata.
2. When `Survey.biology` is fully validated and loaded, the species information provided in the configuration file is used to identify all valid hauls where that species was captured across all of the biological datasets (i.e. `["catch_df"]`, `["length_df"]`, and `["specimen_df"]`). This is then used to **<span style="color:#319FFF">zero out NASC values</span>** `Survey.input["acoustics"]["nasc_df"]` that do not correspond to hauls associated with the target species.


