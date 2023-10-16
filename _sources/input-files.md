# Input files


## Biological (trawl) data

### Length

**Current sample file (US data):** `Biological/US/2019_biodata_length.xls`, sheet `biodata_length`

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
haul_num | integer | | | Describes what haul the collected data comes from. A haul is usally described as a collection of trawls for a certain section of the survey.
species_id | string | | | Identifies what species is associated with the collected data.
sex | string | | | The sex of the animal. This is put into three categories: male, female, and  unsexed.
length | float | cm | | The length of the animal.
length_count | integer | | | The number of animals in the haul, of a particular species, and of a certain sex and length. For example, we have 5 Hake from haul 1 that are males with length 20cm.

### Specimen

**Current sample file (US data):** `Biological/US/2019_biodata_specimen_AGES.xls`, sheet `biodata_specimen`

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
haul_num | integer | | | Describes what haul the collected data comes from. A haul is usally described as a collection of trawls for a certain section of the survey.
species_id | string | | | Identifies what species is associated with the collected data.
sex | string | | | The sex of the animal. This is put into three categories: male, female, and unsexed.
length | float | cm | | The length of the animal. 
weight | float | kg | | The weight of the animal.
age | float | years | | The age of the animal. I believe that this is usually based on the length.

### Catch

**Current sample file (US data):** `Biological/US/2019_biodata_catch.xls`, sheet `biodata_catch`

**Note:** We are currently assessing which columns are actually being used. The original column names were not modified previously, but will be revisited.

### Haul vs transect

File containing the mapping between hauls and transects. **This is a new file that replaces the sole information that was being used from the gear file.**

**Current sample file (US data):** `Biological/US/haul_to_transect_mapping_2019.xls`, single sheet
Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
haul_num | integer | | | Describes what haul the collected data comes from. A haul is usally described as a collection of trawls for a certain section of the survey.
transect_num | integer | | | Transect number (From Net Configuration Form?)


## Stratification

### Strata

The "Strata" file relates the stratification to the haul.

**Current sample file:** `Stratification/US_CAN strata 2019_final.xlsx`, sheet `Base KS`

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
year | integer | | |
stratum_num | integer | | |
*strata name* | string | | |
haul_num | integer | | |
fraction_hake | float | 0-1 | |

### Geo-strata

The "Geo-strata" file defines the geographic definition of strata.

**Current sample file:** `Stratification/Stratification_geographic_Lat_2019_final.xlsx`, `stratification1` sheet

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
stratum_num | integer | | |
Latitude (upper limit) | float | decimal degrees | |
*haul start* | integer | | |
*haul end* | integer | | |
*Strata Name* | string | | |


## NASC

### No Age 1

NASC values that do not include age1 values.

**Current sample file:**  `Exports/US_CAN_detailsa_2019_table2y+_ALL_final - updated.xlsx`, single sheet

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
transect_num | integer | | |
vessel_log_start | float | | |
vessel_log_end | float | | |
latitude | float | decimal degrees | |
longitude | float | decimal degrees | |
stratum_num | integer | | |
transect_spacing | float | nmi? | |
NASC | float | | |

### All ages

NASC values that include all ages.

**Current sample file:**  `Exports/US_CAN_detailsa_2019_table1y+_ALL_final - updated.xlsx`, single sheet

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
transect_num | integer | | |
vessel_log_start | float | nmi? | |
vessel_log_end | float | | nmi? |
latitude | float | decimal degrees | |
longitude | float | decimal degrees | |
stratum_num | integer | | |
transect_spacing | float | | |
NASC | float | | |


## Kriging

### Mesh

The "Mesh" file contains the centroids of the Kriging grid.

**Current sample file:** `Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx`, sheet `krigedgrid2_5nm_forChu`

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
centroid_latitude | float | decimal degrees | |
centroid_longitude | float | decimal degrees | |
fraction_cell_in_polygon | float | 0-1? | |


### Smoothed contours

Smoothed isobath (shelf break, 200 meters?) contour used to transform the mesh points.

**Current sample file:** `Kriging_files/Kriging_grid_files/transformation_isobath_coordinates.xlsx`, sheet `Smoothing_EasyKrig`

Column name | Data type | Units | Empty value | Description
--- | --- | --- | --- | --- 
latitude | float | decimal degrees | |
longitude | float | decimal degrees | |
