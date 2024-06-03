"""
Core API/components for class structures, class methods, and utility/helper methods
"""

from datetime import datetime

import numpy as np
import pandas as pd

# `Survey` object data structure
CONFIG_MAP = {
    "biological": {
        # BIOLOGICAL DATASET -- LENGTH
        "length": {
            "haul_num": int,
            "species_id": int,
            "sex": np.uint8,
            "length": np.float64,
            "length_count": int,
        },
        # BIOLOGICAL DATASET -- CATCH
        "catch": {
            "haul_num": int,
            "species_id": int,
            "haul_weight": np.float64,
        },
        # BIOLOGICAL DATASET -- SPECIMEN
        "specimen": {
            "haul_num": int,
            "species_id": int,
            "sex": np.int8,
            "length": np.float64,
            "weight": np.float64,
            "age": np.float64,
        },
        # BIOLOGICAL DATASET -- HAUL TO TRANSECT
        "haul_to_transect": {
            "haul_num": int,
            "transect_num": int,
        },
    },
    "stratification": {
        # STRATIFICATION DATASET -- STRATA
        "strata": {
            "stratum_num": int,
            "haul_num": int,
            "fraction_hake": np.float64,
        },
        # STRATIFICATION DATASET -- GEOSTRATA
        "geo_strata": {
            "stratum_num": int,
            "northlimit_latitude": np.float64,
            "haul start": int,
            "haul end": int,
        },
        "inpfc_strata": {
            "stratum_num": int,
            "northlimit_latitude": np.float64,
            "haul start": int,
            "haul end": int,
        },
    },
    "NASC": {
        # ACOUSTIC DATASET -- NASC EXCLUDING AGE-1
        "no_age1": {
            "transect_num": int,
            "vessel_log_start": np.float64,
            "vessel_log_end": np.float64,
            "latitude": np.float64,
            "longitude": np.float64,
            "stratum_num": int,
            "transect_spacing": np.float64,
            "NASC": np.float64,
            "haul_num": int,
        },
        # ACOUSTIC DATASET -- NASC EXCLUDING AGE-1
        "all_ages": {
            "transect_num": int,
            "vessel_log_start": np.float64,
            "vessel_log_end": np.float64,
            "latitude": np.float64,
            "longitude": np.float64,
            "stratum_num": int,
            "transect_spacing": np.float64,
            "NASC": np.float64,
            "haul_num": int,
        },
    },
    "kriging": {
        # STATISTICAL DATASET -- KRIGING MESH
        "mesh": {
            "centroid_latitude": np.float64,
            "centroid_longitude": np.float64,
            "fraction_cell_in_polygon": np.float64,
        },
        # STATISTICAL DATASET -- 200 M ISOBATH
        "isobath_200m": {
            "latitude": np.float64,
            "longitude": np.float64,
        },
        # STATISTICAL DATASET -- KRIGING VARIOGRAM PARAMETERIZATION
        "vario_krig_para": {
            "dataprep.y_offset": np.float64,
            "vario.range": np.float64,
            "vario.res": np.float64,
            "vario.vario": np.float64,
            "vario.corr": np.float64,
            "vario.nugt": np.float64,
            "vario.sill": np.float64,
            "vario.lscl": np.float64,
            "vario.powr": np.float64,
            "vario.hole": np.float64,
            "vario.model": np.float64,
            "vario.ytox_ratio": np.float64,
            "vario.ztox_ratio": np.float64,
            "vario.dim": np.float64,
            "krig.x_res": np.float64,
            "krig.y_res": np.float64,
            "krig.proc_opt": np.float64,
            "krig.load_para": np.float64,
            "krig.vario_para": np.float64,
            "krig.krig_para": np.float64,
            "krig.both_para": np.float64,
            "krig.load_griddata_file": np.float64,
            "krig.xmin0": np.float64,
            "krig.xmax0": np.float64,
            "krig.dx0": np.float64,
            "krig.ymin0": np.float64,
            "krig.ymax0": np.float64,
            "krig.dy0": np.float64,
            "krig.nx": np.float64,
            "krig.ny": np.float64,
            "krig.nz": np.float64,
            "krig.dx": np.float64,
            "krig.dy": np.float64,
            "krig.dz": np.float64,
            "krig.model": np.float64,
            "krig.scheme": np.float64,
            "krig.blk_nx": np.float64,
            "krig.blk_ny": np.float64,
            "krig.blk_nz": np.float64,
            "krig.srad": np.float64,
            "krig.kmin": int,
            "krig.kmax": int,
            "krig.elim": np.float64,
            "krig.eps": np.float64,
            "krig.ratio": np.float64,
            "krig.xmin": np.float64,
            "krig.xmax": np.float64,
            "krig.ymin": np.float64,
            "krig.ymax": np.float64,
        },
    },
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
    "meta": {"provenance": dict(), "date": f"{datetime.now():%Y-%m-%d %H:%M:%S%z}"},
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
    "results": {"transect": dict(), "stratified": dict(), "kriging": dict()},
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
