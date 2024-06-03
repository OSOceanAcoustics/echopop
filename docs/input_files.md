(input-files)=
# Input files

Input files used in an Echopop run, grouped by data type. The tables below describe the data columns required by Echopop; other columns will be ignored. All input files are in Excel format. File paths, names, and Excel tab names are specified in the survey year configuration file (e.g., `survey_year_2019_config.yml`).

Biological data are always separated into US vs Canada files. All other data files combine US and Canadian data.

To minimize duplication in the data file description tables below, additional definitions and information for some variables found in multiple files is provided here, especiallly for column names:

- `haul_num`: Haul number. Identifies the haul the collected data come from. A haul is usally described as a collection of trawls for a certain section of the survey.
- `transect_num`: Transect number.
- `species_id`: Species identification code (ID). Identifies what species is associated with the collected data. Pacific hake is 22500.
- `N/P`: Empty value Not Permitted.
- `nmi`: Nautical miles.
- `Old name`: Column name used previously in the Matlab EchoPro program

```{contents}
:local:
:depth: 3
```


## Biological (trawl) data

**Current base directory** used with the sample files: `Biological`

Data files from the US and Canada are found in subdirectories `US` and `CAN`, respectively.

### Length

**Current sample file (US data)** relative to base directory: `US/2019_biodata_length.xls`, sheet `biodata_length`

Column name | Old name | Data type | Units | Empty value   | Description
--- | --- | --- | --- | --- | --- 
haul_num | Haul | integer | | N/P | Haul number
species_id | Species_Code | integer | | N/P | Species identification code (ID)
sex | Sex | integer | | N/P | Sex of the animal. 1=Male, 2=Female, 3=Unknown/Not determined ("unsexed")
length | Length | float | cm | | Length of the animal
length_count | Frequency | float | | empty (blank) | Number of animals in the haul, of a particular species, and of a certain sex and length. For example, we have 5 Hake from haul 1 that are males with length 20cm

### Specimen

**Current sample file (US data)** relative to base directory: `US/2019_biodata_specimen_AGES.xls`, sheet `biodata_specimen`

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
haul_num | Haul | integer | | N/P | Haul number
species_id | Species_Code | integer | | N/P | Species identification code (ID)
sex | Sex | integer | | N/P | Sex of the animal. 1=Male, 2=Female, 3=Unknown/Not determined ("unsexed")
length | Length | float | cm | | Length of the animal
weight | Weight | float | kg | empty (blank) | Weight of the animal
age | Age | float | years | empty (blank) | Age of the animal

### Catch

**Current sample file (US data)** relative to base directory: `US/2019_biodata_catch.xls`, sheet `biodata_catch`

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
haul_num | Haul | integer | | N/P | Haul number
species_id | Species_Code | integer | | N/P | Species identification code (ID)
haul_weight | Weight_In_Haul | float | kg | N/P | Haul weight

### Haul vs transect

File containing the mapping between hauls and transects. This is a new file that replaces the sole information that was being used from the gear file. Note that rows with empty `transect_num` must be omitted.

**Current sample file (US data)** relative to base directory: `US/haul_to_transect_mapping_2019.xls`, single sheet
Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
haul_num | Haul | integer | | N/P | Haul number
transect_num | Transect | integer | | N/P | Transect number


## Stratification

**Current base directory** used with the sample files: `Stratification`

Strata may be based on age-length (`KS`, Kolmogorov-Smirnov test) or regional (`INPFC`, International North Pacific Fisheries Commission) stratifications. Each file contains two tabs, one for each strata type.

### Strata

File that relates the stratification to the haul.

**Current sample file (US data)** relative to base directory: `US_CAN strata 2019_final.xlsx`, sheets `Base KS` and `INPC`

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
stratum_num | Cluster name / INPFC | integer | | N/P | Stratum number for KS or INPC strata (`Base KS` or `INPC` tab, respectively). For `Base KS`, 0 = Low sample size. The Old names listed are for the `Base KS` and `INPFC` tabs, respectively
haul_num | Haul | integer | | N/P | Haul number
fraction_hake | wt | float | 0-1 | N/P | Fraction of the haul weight that is hake

### Geo-strata

File that defines the geographic definition of strata.

**Current sample file (US data)** relative to base directory: `Stratification_geographic_Lat_2019_final.xlsx`, sheets `stratification1` and `INPC`

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
stratum_num | Strata index | integer | | N/P | Stratum number for KS or INPC strata (`stratification1` or `INPC` tab, respectively)
northlimit_latitude | Latitude (upper limit) | float | decimal degrees | N/P | Northern limit of stratum


## NASC

**Current base directory** used with the sample files: `Exports`

### No Age 1

NASC (Nautical Area Scattering Coefficient) values that do not include age1 values. Values are defined along transects at cells with an approximately 0.5 nmi spacing, 

**Current sample file (US data)** relative to base directory: `US_CAN_detailsa_2019_table2y+_ALL_final - updated.xlsx`, single sheet

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
transect_num | Transect | integer | | N/P | Transect number
vessel_log_start | VL start | float | nmi | N/P | Vessel log cumulative distance at start of transect cell
vessel_log_end | VL end | float | nmi | N/P | Vessel log cumulative distance at end of transect cell
latitude | Latitude | float | decimal degrees | N/P | Transect cell center latitude
longitude | Longitude | float | decimal degrees | N/P | Transect cell center longitude
stratum_num | Stratum | integer | | N/P | Base KS stratum number
transect_spacing | Spacing | float | nmi | N/P | Distance (spacing) between transects
NASC | NASC | float | m<sup>2</sup> nmi<sup>-2</sup> | N/P | Nautical Area Scattering Coefficient
haul_num | Assigned haul | integer | | N/P | Haul number. A value of 0 is used for transect cells where a haul was not present or used.

The following columns are currently not used in core computations. They are used in reports and in some plots (plots not implemented yet). The column names are the original names and have not been "sanitized".

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
Region ID | Region ID | int | | | 
Bottom depth | Bottom depth | float | meters? | | 
Layer mean depth | Layer mean depth | float | meters? | | 
Layer height | Layer height | float | meters? | | 

### All ages

NASC values that include all ages. The file structure is the same as for the "No Age 1" file.

**Current sample file (US data)** relative to base directory: `US_CAN_detailsa_2019_table1y+_ALL_final - updated.xlsx`, single sheet


## Kriging

**Current base directory** used with the sample files: `Kriging_files`

### Mesh

The "Mesh" file containing the centroids of the Kriging grid cells. Grid size is 2.5 nmi by 2.5 nmi.

**Current sample file (US data)** relative to base directory: `Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx`, sheet `krigedgrid2_5nm_forChu`

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
centroid_latitude | Latitude of centroid | float | decimal degrees | N/P | Cell centroid latitude
centroid_longitude | Longitude of centroid | float | decimal degrees | N/P | Cell centroid longitude
fraction_cell_in_polygon | Cell portion | float | 0-1 | N/P | Fraction of mesh cell that is within the interpolation polygon that delineates the mesh

### Smoothed shelf-break contour

Smoothed isobath contour used to transform the mesh points. A set of point locations delineating the 200 meter bathymetric contour that represents the shelf break.

**Current sample file (US data)** relative to base directory: `Kriging_grid_files/transformation_isobath_coordinates.xlsx`, sheet `Smoothing_EasyKrig`

Column name | Old name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- | --- 
latitude | Latitude | float | decimal degrees | N/P | Point latitude
longitude | Longitude | float | decimal degrees | N/P | Point longitude
