#################
biodata_filepath = Path("C:/Users/Brandyn/Documents/GitHub/Data/Biological/1995-2023_biodata_redo.xlsx")
biodata_sheet_map = {
    "catch": "biodata_catch", 
    "length": "biodata_length",
    "specimen": "biodata_specimen",
}
species_code = 22500
subset_dict = {
    "ships": {
        160: {
            "survey": 201906
        },
        584: {
            "survey": 2019097,
            "haul_offset": 200
        }
    },
    "species_code": [22500]
}

FEAT_TO_ECHOPOP_BIODATA_COLUMNS = {
    "frequency": "length_count",
    "haul": "haul_num",
    "weight_in_haul": "haul_weight",
}

column_name_map = FEAT_TO_ECHOPOP_BIODATA_COLUMNS
sheet_name = sheet_map["catch"]
biodata_label_map = {
    "sex": {
        1: "male",
        2: "female",
        3: "unsexed"
    }
}

ROOT_PATH = Path("C:/Users/Brandyn/Documents/GitHub/Data")

biodata_filepath = ROOT_PATH / "Biological/1995-2023_biodata_redo.xlsx"

dict_df_bio = load_biological_data(biodata_filepath, biodata_sheet_map, FEAT_TO_ECHOPOP_BIODATA_COLUMNS, subset_dict, biodata_label_map)

ROOT_PATH = Path("C:/Users/Brandyn/Documents/GitHub/Data")
FEAT_TO_ECHOPOP_STRATA_COLUMNS = {
    "fraction_hake": "nasc_proportion",
    "haul": "haul_num",
    "stratum": "stratum_num",
}

strata_filepath = ROOT_PATH / "Stratification/US_CAN strata 2019_final.xlsx"
strata_sheet_map = {
    "inpfc": "INPFC",
    "ks": "Base KS",
}
column_name_map = FEAT_TO_ECHOPOP_STRATA_COLUMNS

dict_df_strata = load_strata(strata_filepath, strata_sheet_map, column_name_map)

FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS = {
    "latitude (upper limit)": "northlimit_latitude",
    "stratum": "stratum_num",
}
geostrata_filepath = ROOT_PATH / "Stratification/Stratification_geographic_Lat_2019_final.xlsx"
geostrata_sheet_map = {
    "inpfc": "INPFC",
    "ks": "stratification1",
}
column_name_map = FEAT_TO_ECHOPOP_STRATA_COLUMNS

dict_df_geostrata = load_geostrata(geostrata_filepath, geostrata_sheet_map, FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS)

data = dict_df_bio
strata_df = dict_df_strata["ks"].copy()

join_strata_by_haul(nasc_all_ages_df, dict_df_strata["inpfc"])
join_strata_by_haul(nasc_all_ages_df, dict_df_strata["ks"])
join_strata_by_haul(dict_df_bio, dict_df_strata["inpfc"])
join_strata_by_haul(dict_df_bio, dict_df_strata["ks"])

data = nasc_all_ages_df.copy()
geostrata_df = dict_df_geostrata["ks"].copy()
FEAT_TO_ECHOPOP_MESH_COLUMNS = {
    "centroid_latitude": "latitude",
    "centroid_longitude": "longitude",
    "fraction_cell_in_polygon": "fraction",
}

column_name_map = FEAT_TO_ECHOPOP_MESH_COLUMNS
mesh_filepath = ROOT_PATH / "Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx"
mesh_sheet_name = "krigedgrid2_5nm_forChu"

mesh_df = load_mesh_data(mesh_filepath, mesh_sheet_name, FEAT_TO_ECHOPOP_MESH_COLUMNS)

join_geostrata_by_latitude(mesh_df, dict_df_geostrata["inpfc"])

geostatistic_params_filepath = ROOT_PATH / "Kriging_files/default_vario_krig_settings_2019_US_CAN.xlsx"
geostatistic_params_sheet_name = "Sheet1"
FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS = {
    "hole": "hole_effect_range",
    "lscl": "correlation_range",
    "nugt": "nugget",
    "powr": "decay_power",
    "ratio": "anisotropy",
    "res": "lag_resolution",
    "srad": "search_radius",
}
column_name_map = FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS

kriging_params_dict, variogram_params_dict = load_kriging_variogram_params(
    geostatistic_params_filepath,
    geostatistic_params_sheet_name,
    FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS
)

####################################################################################################








