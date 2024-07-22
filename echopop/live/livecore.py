from datetime import datetime

import pandas as pd

LIVE_DATA_STRUCTURE = {
    "meta": {
        "provenance": dict(),
        "date": list(),
    },
    "database": {
        "acoustics": None,
        "biology": None,
        "population": None,
    },
    "input": {
        "acoustics": {
            "nasc_df": pd.DataFrame(),
        },
        "biology": {
            "catch_df": pd.DataFrame(),
            "distributions": {
                "length_bins_df": pd.DataFrame(),
            },
            "length_df": pd.DataFrame(),
            "specimen_df": pd.DataFrame(),
        },
    },
    "results": {
        "acoustics": dict(),
        "biology": dict(),
        "stratified": dict(),        
    },
}

# TODO: Update structure with additional information (as needed)
# TODO: Documentation
LIVE_INPUT_FILE_CONFIG_MAP = {
    "acoustics": {
        "xarray_coordinates": {
            "distance": float,
            "depth": float,
        },
        "xarray_variables": {
            "NASC": float,
            "frequency_nominal": float, 
            "latitude": float,
            "longitude": float,
            "ping_time": "datetime64[ns]",
        }
    },
    "biology": {
        "catch": {
            "dtypes": {
                "partition": str,
                "species_code": int,
                "sample_weight_kg": float,
                "catch_perc": float,
            },
            "names": {
                "partition": "trawl_partition",
                "species_code": "species_id",
                "sample_weight_kg": "haul_weight",
                "catch_perc": "catch_percentage",
            }
        },
        "length": {
            "dtypes": {
                "sex": str,
                "rounded_length": int,
                "frequency": int,
            },
            "names": {
                "sex": "sex",
                "rounded_length": "length",
                "frequency": "length_count",
            }
        },
        "specimen": {
            "dtypes": {
                "rounded_length": int,
                "organism_weight": float,
                "sex": str,
            },
            "names": {
                "sex": "sex",
                "rounded_length": "length",
                "organism_weight": "weight"
            },
        },
    },
}

LIVE_FILE_FORMAT_MAP = {
    "DATE:YYYYMM": {
        "name": "date",
        "dtype": "datetime[ns]",
        "expression": r"(?P<DATE>\d{6})",
    },
    "DATE:YYYYMMDD": {
        "name": "date",
        "dtype": "datetime[ns]",
        "expression": r"(?P<DATE>\d{8})",
    },
    "HAUL": {
        "name": "haul_num",
        "dtype": int,
        "expression": r"(?P<HAUL>\d+)",
    },
    "SPECIES_CODE": {
        "name": "species_id",
        "dtype": int,
        "expression": r"(?P<SPECIES_CODE>\d+)"
    },
    "FILE_ID": {
        "name": "file_id",
        "dtype": str,
        "expression": r"(?P<FILE_ID>.+)"
    },
}
