(input-files)=
# Input files

Input files used in an Echopop run, grouped by data type. The tables below describe the data columns required by Echopop; other columns will be ignored. All input files are in Excel format. File paths, names, and Excel tab names are specified in the survey year configuration file (e.g., `survey_year_2019_config.yml`).

Biological data are always separated into US vs Canada files. All other data files combine US and Canadian data.

:::{admonition} Dataset structures
:class: note
*See page for [configuration dataset file organization](implementation/preprocessing_data) for more details.*
:::

```{contents}
:local:
:depth: 3
```

:::{admonition} `EchoPro` column names
:class: note
There may be some inconsistencies in the columns used by files in `EchoPro` for previous years. These column names are based on the files associated with the 2019 survey.
:::

## Biological data

### Length

`Echopop` column | `EchoPro` column | Data type | Units | Description
--- | --- | --- | ----- | --- 
haul_num | Haul | integer<br>float |   | Haul number. <br> Rows with missing values are removed
species_id | Species_Code | integer<br>float<br>string |   | Species identification code (ID)<br>Rows with missing values are removed
sex | Sex | integer | | Sex of the animal <br> Male: `1`/`"m"`/`"male"`, Female: `2`/`"f"`/`"female"`, Unsexed: `3`/`"u"`/`"unsexed"` <br>Missing values are replaced with `"unsexed"`
length | Length | float | cm <br> (0.0, ∞) | Animal fork length <br> Missing values are replaced with `NaN`
length_count | Frequency | integer | count <br> [0, ∞)  | Number of animals with the corresponding binned fork length <br> Missing values are replaced with `0`

### Specimen

`Echopop` column | `EchoPro` column | Data type | Units &nbsp; | Description
--- | --- | --- | --- | --- 
haul_num | Haul | integer | | Haul number. <br> Rows with missing values are removed
species_id | Species_Code | integer |    | Species identification code (ID) <br> Rows with missing values are removed
sex | Sex | integer | | Sex of the animal <br> Male: `1`/`"m"`/`"male"`, Female: `2`/`"f"`/`"female"`, Unsexed: `3`/`"u"`/`"unsexed"` <br> Missing values are replaced with `"unsexed"`
length | Length | float | cm <br> (0.0, ∞) | Animal fork length <br> Missing values are replaced with `NaN`
weight | Weight | float | kg <br> (0.0, ∞) | Specimen weight <br> Missing values are replaced with `NaN`
age | Age | float <br> integer | years <br> [0.0, ∞)  | Age of the animal <br> Missing values are replaced with `NaN`

### Catch

`Echopop` column | `EchoPro` column | Data type | Units |  Description
--- | --- | --- | --- | --- 
haul_num | Haul | integer | | Haul number <br> Rows with missing values are removed
species_id | Species_Code | integer | | Species identification code (ID) <br> Rows with missing values are removed
haul_weight | Weight_In_Haul | float | kg <br> [0.0, ∞) | Haul weight <br> Rows with missing values are removed

## Stratification

Strata may be based on age-length (`KS`, Kolmogorov-Smirnov test) or regional (`INPFC`, International North Pacific Fisheries Commission) stratifications. Each file contains two tabs, one for each strata type.

### Strata

File that relates the stratification to the haul.


`Echopop` column | `EchoPro` column | Data type | Units | Description
--- | --- | --- | --- | --- 
stratum_num | Strata index | integer | | Index/grouping representing the stratum identifier based on <br> either length (KS) or latitude (INPFC) <br> Rows with missing values are removed
haul_num | Haul | integer | | Haul number <br> Rows with missing values are removed
fraction_hake | wt | float | proportion <br> [0.0-1.0] | Fraction of the haul weight that is hake <br> Missing values are replaced with `0.0`

### Geostrata

File that defines the geographic definition of strata.


`Echopop` column | `EchoPro` column | Data type | Units | Description
--- | --- | --- | --- | --- 
stratum_num | Strata index | integer | | Index/grouping representing the stratum identifier based on <br> either length (KS) or latitude (INPFC) <br> Rows with missing values are removed
northlimit_latitude | Latitude (upper limit) | float | decimal degrees <br> [-90.0, 90.0] | Northern limit of stratum <br> Rows with missing values are removed 


## NASC

### No age-1 fish

NASC (Nautical Area Scattering Coefficient) values that do not include age-1 values. Values are defined along transects at cells with an approximately 0.5 nmi spacing, 

`Echopop` column | `EchoPro` column | Data type | Units | Description
--- | --- | --- | --- |  --- 
transect_num | Transect | integer <br> float | | Transect number <br> Rows with missing values are removed
distance_s | VL start | float | nmi <br> [0.0, ∞) |  Vessel log cumulative distance at start of transect interval <br> Missing values are replaced with `NaN`
distance_e | VL end | float | nmi <br> [0.0, ∞) | Vessel log cumulative distance at end of transect interval <br> Missing values are replaced with `NaN`
latitude | Latitude | float | decimal degrees <br> [-90.0, 90.0] | Transect interval center latitude <br> Missing values are replaced with `NaN`
longitude | Longitude | float | decimal degrees <br> [-180.0, 180.0] | Transect interval center longitude <br> Missing values are replaced with `NaN`
transect_spacing | Spacing | float | nmi <br> [0.0, ∞) | Distance (spacing) between transects <br> Missing values are replaced with the user-defined <br> value for `max_transect_spacing` in the <br>`initialization_config.yml` configuration file
NASC | NASC | float | m<sup>2</sup> nmi<sup>-2</sup> <br> [0.0, ∞) | Nautical area scattering strength ($\textit{NASC}$) <br> Missing values are replaced with `0.0`
haul_num | Haul | integer <br> float | | Assigned haul number <br> A value of `0` is used for transect intervals where no haul <br> was present or used <br> Rows with missing values are removed


### All aged fish

NASC values that include all ages. The file structure is the same as for the files excluding age-1 fish.

## Kriging

### Mesh

The "Mesh" file containing the centroids of the Kriging grid cells. The default grid size is 2.5 nmi by 2.5 nmi.

`Echopop` column | `EchoPro` column | Data type | Units |  Description
--- | --- | --- | --- | --- 
centroid_latitude | Latitude of centroid | float | decimal degrees <br> [-90.0, 90.0] | Mesh cell centroid latitude <br> Rows with missing values are removed
centroid_longitude | Longitude of centroid | float | decimal degrees <br> [-180.0, 180.0] | Mesh cell centroid longitude <br> Rows with missing values are removed
fraction_cell_in_polygon | Cell portion | float | proportion <br> [0.0-1.0] | Fraction of mesh cell that is within the interpolation polygon that delineates the mesh <br> Missing values are replaced with `0.0`

### Smoothed shelf-break contour

Smoothed isobath contour used to transform the mesh points. A set of point locations delineating the 200 meter bathymetric contour that represents the shelf break.

`Echopop` column | `EchoPro` column | Data type | Units | Description
--- | --- | --- | --- | --- 
latitude | Latitude | float | decimal degrees <br> [-90.0, 90.0] | Isobath latitude <br> Rows with missing values are removed
longitude | Longitude | float | decimal degrees <br> [-180.0, 180.0] | Isobath longitude <br> Rows with missing values are removed

### Kriging and variogram parameters

This file comprises two columns: 1) the parameter names and 2) their associated values. 

`Echopop` parameter | `EchoPro` parameter | Data type | Valid range | Description
--- | --- | --- | --- | ---
correlation_range | vario.lscl | float | [0.0, ∞) | The relative length scale, or range at which the correlation between points becomes approximately constant
sill | vario.sill | float | [0.0, ∞) | The total variance where the change autocorrelation reaches (or nears) 0.0
nugget | vario.nugt | float | [0.0, ∞) | The $y$-intercept of the variogram representing the short-scale (i.e. smaller than the lag resolution) variance
decay_power | vario.pwr | float | [0.0, ∞) | The exponent used for variogram models with exponentiated spatial decay terms
hole_effect_range | vario.hole | float | [0.0, ∞) | Length scale or range of the hole effect
lag_resolution | vario.res | float | (0.0, ∞) | The (scaled) distance between lags
anisotropy | vario.ratio | float | (0.0, ∞) | The directional aspect ratio of anisotropy
search_radius | krig.srad | float | (0.0, ∞) | The adaptive search radius used for kriging
kmin | kmin | integer | (0, ∞) | The minimum number of nearest kriging points
kmax | kmax | integer | (kmin, ∞) | The maximum number of nearest kriging points
