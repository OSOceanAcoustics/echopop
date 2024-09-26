"""
Core API/components for class structures, class methods, and utility/helper methods
"""

from datetime import datetime

import numpy as np
import pandas as pd

# Name configuration dictionary
NAME_CONFIG = {
    "Age": "age",
    "Cell portion": "fraction_cell_in_polygon",
    "Frequency": "length_count",
    "haul": "haul_num",
    "haul end": "haul_end",
    "haul start": "haul_start",
    "Haul": "haul_num",
    "Latitude": "latitude",
    "Latitude (upper limit)": "northlimit_latitude",
    "Latitude of centroid": "centroid_latitude",
    "Length": "length",
    "Longitude": "longitude",
    "Longitude of centroid": "centroid_longitude",
    "strata": "stratum_num",
    "Sex": "sex",
    "Ship": "ship_id",
    "Spacing": "transect_spacing",
    "Species_Code": "species_id",
    "Species_Name": "species_name",
    "Strata Index": "stratum_num",
    "Stratum": "stratum_num",
    "Transect": "transect_num",
    "VL start": "vessel_log_start",
    "VL end": "vessel_log_end",
    "wt": "fraction_hake",
    "Weight": "weight",
    "Weight_In-Haul": "haul_weight",
}

ECHOVIEW_EXPORT_MAP = {
    "Region_ID": dict(name="region_id", type=int),
    "Region_name": dict(name="region_name", type=str),
    "Region_class": dict(name="region_class", type=str),
    "Interval": dict(name="interval", type=int),
    "Layer": dict(name="layer", type=int),
    "Sv_mean": dict(name="sv_mean", type=float),
    "PRC_NASC": dict(name="nasc", type=float),
    "transect_num": dict(name="transect_num", type=int),
    "VL_start": dict(name="vessel_log_start", type=float),
    "VL_end": dict(name="vessel_log_end", type=float),
    "Lat_S": dict(name="latitude", type=float),
    "Lat_M": dict(name="latitude", type=float),
    "Lat_E": dict(name="latitude", type=float),
    "Lon_S": dict(name="longitude", type=float),
    "Lon_M": dict(name="longitude", type=float),
    "Lon_E": dict(name="longitude", type=float),
    "Layer_depth_min": dict(name="layer_depth_min", type=float),
    "Layer_depth_max": dict(name="layer_depth_max", type=float),
    "Exclude_below_line_depth_mean": dict(name="max_depth", type=float),
}

REGION_EXPORT_MAP = {
    "no_age1": {
        "Tranect": dict(name="transect_num", type=int),
        "Region ID": dict(name="region_id", type=int),
        "Trawl #": dict(name="haul_num", type=int),
        "Region Name": dict(name="region_name", type=str),
        "Region Calss": dict(name="region_class", type=str),
    },
    "all_ages": {
        "Tranect": dict(name="transect_num", type=int),
        "Region ID": dict(name="region_id", type=int),
        "Trawl #": dict(name="haul_num", type=int),
        "Region Name": dict(name="region_name", type=str),
        "Region Calss": dict(name="region_class", type=str),
    },
}

BIODATA_HAUL_MAP = {
    "Haul": dict(name="haul_num", type=int),
    "Transect": dict(name="transect_num", type=int),
}

# ``LAYER_NAME_MAP`` is a hard-coded dictionary that aids in re-mapping the data tree
# structure of various attributes contained within a generated ``Survey`` class object.
# structure of various attributes contained within a generated ``Survey`` class object.
# Each dictionary comprises the same dataset names contained within the configuration file.
# Nested within each of these dictionaries include three variables:
# Nested within each of these dictionaries include three variables:
#     name (string): the name of the attribute appended to a ``Survey`` class object.
#     data (list[string]): the name of the expected data layer, which is formatted as a dataframe.
#     superlayer (list[string]): layer names that are prepended onto the nested dictionary structure
#     of the data attribute.
#     superlayer (list[string]): layer names that are prepended onto the nested dictionary structure
#     of the data attribute.
# While this reference library is only required in the sense that it is hard-coded into the
# the current implementation concerning loading data into a ``Survey`` class object, it may be
# the current implementation concerning loading data into a ``Survey`` class object, it may be
# helpful for organizing data with more parsimonious and representative names. It may also be
# helpful for future extensions enabling users to set additional data attributes that aren't
# included in the ``Survey`` class initialization and parameterization (when the object is
# generated).
DATA_STRUCTURE = {
    "meta": {
        "provenance": {"imported_datasets": set(), "date": f"{datetime.now():%Y-%m-%d %H:%M:%S%z}"},
    },
    "input": {
        "acoustics": {
            "nasc_df": pd.DataFrame(),
        },
        "biology": {
            "catch_df": pd.DataFrame(),
            "distributions": {
                "age_bins_df": pd.DataFrame(),
                "length_bins_df": pd.DataFrame(),
            },
            "length_df": pd.DataFrame(),
            "haul_to_transect_df": pd.DataFrame(),
            "specimen_df": pd.DataFrame(),
        },
        "spatial": {
            "strata_df": pd.DataFrame(),
            "geo_strata_df": pd.DataFrame(),
            "inpfc_strata_df": pd.DataFrame(),
        },
        "statistics": {
            "kriging": {
                "mesh_df": pd.DataFrame(),
                "isobath_200m_df": pd.DataFrame(),
                "vario_krig_para_df": pd.DataFrame(),
                "model_config": dict(),
            },
            "variogram": {
                "model_config": dict(),
            },
        },
    },
    "analysis": {
        "transect": {
            "acoustics": {"sigma_bs": dict()},
            "biology": {
                "counts": dict(),
                "distributions": dict(),
                "population": {"tables": dict()},
                "proportions": {"number": dict(), "weight": dict()},
                "weight": {"length_weight_regression": dict()},
            },
            "coordinates": pd.DataFrame(),
        },
        "settings": {
            "kriging": dict(),
            "stratified": dict(),
            "transect": dict(),
            "variogram": dict(),
        },
        "stratified": dict(),
    },
    "results": {"transect": dict(), "stratified": dict(), "kriging": dict(), "variogram": dict()},
}

LAYER_NAME_MAP = {
    "biological": {
        "name": "biology",
        "data": ["length", "specimen", "catch", "haul_to_transect"],
        "data_label": [
            "length:unaged lengths",
            "specimen:aged lengths",
            "catch:unaged catch weights",
            "haul_to_transect:haul-to-transect key",
        ],
        "superlayer": [],
        "data_tree": {
            "catch_df": pd.DataFrame(),
            "distributions": {
                "length_bins_arr": np.array([]),
                "age_bins_arr": np.array([]),
            },
            "length_df": pd.DataFrame(),
            "haul_to_transect_df": pd.DataFrame(),
            "specimen_df": pd.DataFrame(),
        },
    },
    "stratification": {
        "name": "spatial",
        "data": ["strata", "geo_strata", "inpfc_strata"],
        "data_label": [
            "strata:KS strata",
            "geo_strata:Georeferenced KS strata",
            "inpfc_strata:INPFC strata",
        ],
        "superlayer": [],
        "data_tree": {
            "strata_df": pd.DataFrame(),
            "geo_strata_df": pd.DataFrame(),
            "inpfc_strata_df": pd.DataFrame(),
        },
    },
    "NASC": {
        "name": "nasc",
        "data": ["no_age1", "all_ages"],
        "data_label": ["no_age1:Survey NASC", "all_ages:no_age1:Survey NASC"],
        "superlayer": ["acoustics"],
        "data_tree": {
            "nasc": {
                "nasc_df": pd.DataFrame(),
            },
        },
    },
    "kriging": {
        "name": "kriging",
        "data": ["mesh", "isobath_200m", "vario_krig_para"],
        "data_label": ["mesh:Kriging mesh", "isobath_200m:200m isobath"],
        "superlayer": ["statistics"],
        "data_tree": {
            "kriging": {
                "mesh_df": pd.DataFrame(),
                "isobath_200m_df": pd.DataFrame(),
                "vario_krig_para_df": pd.DataFrame(),
                "model_config": dict(),
            },
            "variogram": {
                "model_config": dict(),
            },
        },
    },
}
